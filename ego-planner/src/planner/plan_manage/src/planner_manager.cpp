// #include <fstream>
#include <plan_manage/planner_manager.h>
#include <thread>

namespace ego_planner
{

  // SECTION interfaces for setup and query

  EGOPlannerManager::EGOPlannerManager() {}

  EGOPlannerManager::~EGOPlannerManager() { std::cout << "des manager" << std::endl; }

  void EGOPlannerManager::initPlanModules(ros::NodeHandle &nh, PlanningVisualization::Ptr vis)
  //初始化规划参数
  {
    /* read algorithm parameters */

    nh.param("manager/max_vel", pp_.max_vel_, -1.0);
    nh.param("manager/max_acc", pp_.max_acc_, -1.0);
    nh.param("manager/max_jerk", pp_.max_jerk_, -1.0);
    nh.param("manager/feasibility_tolerance", pp_.feasibility_tolerance_, 0.0);
    nh.param("manager/control_points_distance", pp_.ctrl_pt_dist, -1.0);
    nh.param("manager/planning_horizon", pp_.planning_horizen_, 5.0);//规划的视野范围

    local_data_.traj_id_ = 0;//局部轨迹段index
    grid_map_.reset(new GridMap);//生成地图
    grid_map_->initMap(nh);//初始化

    bspline_optimizer_rebound_.reset(new BsplineOptimizer);
    bspline_optimizer_rebound_->setParam(nh);
    bspline_optimizer_rebound_->setEnvironment(grid_map_);
    bspline_optimizer_rebound_->a_star_.reset(new AStar);
    bspline_optimizer_rebound_->a_star_->initGridMap(grid_map_, Eigen::Vector3i(100, 100, 100));//(100, 100, 100)局部地图的大小

    visualization_ = vis;
  }

  // !SECTION

  // SECTION rebond replanning

  bool EGOPlannerManager::reboundReplan(Eigen::Vector3d start_pt, Eigen::Vector3d start_vel,
                                        Eigen::Vector3d start_acc, Eigen::Vector3d local_target_pt,
                                        Eigen::Vector3d local_target_vel, bool flag_polyInit, bool flag_randomPolyTraj)
  {

    static int count = 0;
    std::cout << endl
              << "[rebo replan]: -------------------------------------" << count++ << std::endl;
    cout.precision(3);//输出浮点精度为3
    cout << "start: " << start_pt.transpose() << ", " << start_vel.transpose() << "\ngoal:" << local_target_pt.transpose() << ", " << local_target_vel.transpose()
         << endl;

    if ((start_pt - local_target_pt).norm() < 0.2)//起点和局部规划的终点距离太短
    {
      cout << "Close to goal" << endl;
      continous_failures_count_++;//错误的次数+1
      return false;
    }

    ros::Time t_start = ros::Time::now();//开始时间
    ros::Duration t_init, t_opt, t_refine;

    /*** STEP 1: INIT初始化 ***/
    //什么意思？？？？？？？
    double ts = (start_pt - local_target_pt).norm() > 0.1 ? pp_.ctrl_pt_dist / pp_.max_vel_ * 1.2 : pp_.ctrl_pt_dist / pp_.max_vel_ * 5; // pp_.ctrl_pt_dist / pp_.max_vel_ is too tense, and will surely exceed the acc/vel limits
    
    vector<Eigen::Vector3d> point_set, start_end_derivatives;//存放轨迹点以及开始结束点的速度加速度
    static bool flag_first_call = true, flag_force_polynomial = false;
    bool flag_regenerate = false;
    do
    {
      point_set.clear();//轨迹点清空
      start_end_derivatives.clear();//开始结束点的速度加速度清空
      flag_regenerate = false;//
     //初始的轨迹来自一个minimumsnap轨迹
      if (flag_first_call || flag_polyInit || flag_force_polynomial /*|| ( start_pt - local_target_pt ).norm() < 1.0*/) // Initial path generated from a min-snap traj by order.
      //如果是第一次调用或者要求使用多项式优化
      {
        flag_first_call = false;
        flag_force_polynomial = false;

        PolynomialTraj gl_traj;//多项式轨迹

        double dist = (start_pt - local_target_pt).norm();//起点和局部规划目标点的距离
        //这段轨迹时间的确定，依据梯形原则，参考ppt362
        double time = pow(pp_.max_vel_, 2) / pp_.max_acc_ > dist ? sqrt(dist / pp_.max_acc_) : (dist - pow(pp_.max_vel_, 2) / pp_.max_acc_) / pp_.max_vel_ + 2 * pp_.max_vel_ / pp_.max_acc_;

        if (!flag_randomPolyTraj)//如果不加入随机点
        {
          gl_traj = PolynomialTraj::one_segment_traj_gen(start_pt, start_vel, start_acc, local_target_pt, local_target_vel, Eigen::Vector3d::Zero(), time);
        }//生成从起点到局部轨迹目标点的minimumsnap轨迹

        else//加入随机点后再minimumsnap
        {
          //水平方向向量：和（0 0 1）叉乘
          Eigen::Vector3d horizen_dir = ((start_pt - local_target_pt).cross(Eigen::Vector3d(0, 0, 1))).normalized();
          //垂直方向向量
          Eigen::Vector3d vertical_dir = ((start_pt - local_target_pt).cross(horizen_dir)).normalized();
          //随机插值点 为什么这样定义？？？？？？？？？？？？？？？？？
          Eigen::Vector3d random_inserted_pt = (start_pt + local_target_pt) / 2 +
                                               (((double)rand()) / RAND_MAX - 0.5) * (start_pt - local_target_pt).norm() * horizen_dir * 0.8 * (-0.978 / (continous_failures_count_ + 0.989) + 0.989) +
                                               (((double)rand()) / RAND_MAX - 0.5) * (start_pt - local_target_pt).norm() * vertical_dir * 0.4 * (-0.978 / (continous_failures_count_ + 0.989) + 0.989);

          Eigen::MatrixXd pos(3, 3);//位置矩阵 分别存放起点 局部规划目标点和随机插值点
          pos.col(0) = start_pt;
          pos.col(1) = random_inserted_pt;
          pos.col(2) = local_target_pt;
          Eigen::VectorXd t(2);
          t(0) = t(1) = time / 2;//均匀分配两段时间
          gl_traj = PolynomialTraj::minSnapTraj(pos, start_vel, local_target_vel, start_acc, Eigen::Vector3d::Zero(), t);
        }//获得三个点的minimumsnap轨迹

        double t;
        bool flag_too_far;
        ts *= 1.5; //  在接下来的时间里，ts将除以1.5   ts will be divided by 1.5 in the next
        do
        {
          ts /= 1.5;
          point_set.clear();//清空位置点
          flag_too_far = false;
          Eigen::Vector3d last_pt = gl_traj.evaluate(0);//初始化上一个点的位置为起点
          for (t = 0; t < time; t += ts)//通过循环时间循环轨迹 存放位置点向量中
          {
            Eigen::Vector3d pt = gl_traj.evaluate(t);//当前的位置
            if ((last_pt - pt).norm() > pp_.ctrl_pt_dist * 1.5)//如果两个点之间的距离>1.5*相邻B样条控制点之间的距离
            {
              flag_too_far = true;//距离太远标志为true
              break;//跳出循环
            }
            last_pt = pt;//更新位置点
            point_set.push_back(pt);//存放位置点到向量中
          }
        } while (flag_too_far || point_set.size() < 7); //确保初始路径有足够的点 To make sure the initial path has enough points.
       
        t -= ts;
        start_end_derivatives.push_back(gl_traj.evaluateVel(0));
        start_end_derivatives.push_back(local_target_vel);
        start_end_derivatives.push_back(gl_traj.evaluateAcc(0));
        start_end_derivatives.push_back(gl_traj.evaluateAcc(t));//开始结束点的速度加速度
      }
      else //  由先前轨迹生成的初始路径  Initial path generated from previous trajectory.
      {

        double t;
        double t_cur = (ros::Time::now() - local_data_.start_time_).toSec();//当前的时间

        vector<double> pseudo_arc_length;//伪弧长
        vector<Eigen::Vector3d> segment_point;//这段轨迹中的点
        pseudo_arc_length.push_back(0.0);
        for (t = t_cur; t < local_data_.duration_ + 1e-3; t += ts)//从当前时间到局部规划终止时间
        {
          segment_point.push_back(local_data_.position_traj_.evaluateDeBoorT(t));//把b样条的轨迹点位置加入向量中
          if (t > t_cur)
          {
            pseudo_arc_length.push_back((segment_point.back() - segment_point[segment_point.size() - 2]).norm() + pseudo_arc_length.back());//加这一段曲线的长度
          }
        }
        t -= ts;

        
        double poly_time = (local_data_.position_traj_.evaluateDeBoorT(t) - local_target_pt).norm() / pp_.max_vel_ * 2;//结束点到局部目标local_target_pt 以max_vel_/2 的速度行驶的时间
        if (poly_time > ts)
        {
          PolynomialTraj gl_traj = PolynomialTraj::one_segment_traj_gen(local_data_.position_traj_.evaluateDeBoorT(t),
                                                                        local_data_.velocity_traj_.evaluateDeBoorT(t),
                                                                        local_data_.acceleration_traj_.evaluateDeBoorT(t),
                                                                        local_target_pt, local_target_vel, Eigen::Vector3d::Zero(), poly_time);

          for (t = ts; t < poly_time; t += ts)
          {
            if (!pseudo_arc_length.empty())
            {
              segment_point.push_back(gl_traj.evaluate(t));
              pseudo_arc_length.push_back((segment_point.back() - segment_point[segment_point.size() - 2]).norm() + pseudo_arc_length.back());
            }
            else
            {
              ROS_ERROR("pseudo_arc_length is empty, return!");
              continous_failures_count_++;
              return false;
            }
          }
        }//如果该时间大于ts ，则用PolynomialTraj::one_segment_traj_gen() 规划从local_data_ 结束点到局部目标local_target_pt 的轨迹gl_traj
        //接着在gl_traj上等时间间隔ts的距离累加赋值到pseudo_arc_length 并将点保存到segment_point 。

        double sample_length = 0;
        double cps_dist = pp_.ctrl_pt_dist * 1.5; // cps_dist will be divided by 1.5 in the next
        size_t id = 0;
        do
        {
          cps_dist /= 1.5;
          point_set.clear();
          sample_length = 0;
          id = 0;
          while ((id <= pseudo_arc_length.size() - 2) && sample_length <= pseudo_arc_length.back())//循环segment_point 中每两个相邻点
          {
            if (sample_length >= pseudo_arc_length[id] && sample_length < pseudo_arc_length[id + 1])
            {
              point_set.push_back((sample_length - pseudo_arc_length[id]) / (pseudo_arc_length[id + 1] - pseudo_arc_length[id]) * segment_point[id + 1] +
                                  (pseudo_arc_length[id + 1] - sample_length) / (pseudo_arc_length[id + 1] - pseudo_arc_length[id]) * segment_point[id]);
                                  //在两点中均匀插值
              sample_length += cps_dist;
            }
            else
              id++;
          }
          point_set.push_back(local_target_pt);
        } while (point_set.size() < 7); // If the start point is very close to end point, this will help
        //令cps_dist =ctrl_pt_dist*1.5接着对于segment_point 中每两个相邻点，如果距离大于cps_dist ，则在这两点间加入几个点，使得相邻两点距离小于等于cps_dist 。
        //如果最终路径点的个数小于7，则减小cps_dist 重新加点，直到路径点的个数大于等于7。

        start_end_derivatives.push_back(local_data_.velocity_traj_.evaluateDeBoorT(t_cur));
        start_end_derivatives.push_back(local_target_vel);
        start_end_derivatives.push_back(local_data_.acceleration_traj_.evaluateDeBoorT(t_cur));
        start_end_derivatives.push_back(Eigen::Vector3d::Zero());//将local_data_现在的速度加速度和局部终点的速度加速度压入start_end_derivatives

       // 如果路径点的个数大于planning_horizen_/pp_.ctrl_pt_dist*3 则重新使用规则1 进行初始化。
        if (point_set.size() > pp_.planning_horizen_ / pp_.ctrl_pt_dist * 3) // point_set >56 初始路径通常太长The initial path is unnormally too long!
        {
          flag_force_polynomial = true;
          flag_regenerate = true;//
        }
      }
    } while (flag_regenerate);

    Eigen::MatrixXd ctrl_pts;
    UniformBspline::parameterizeToBspline(ts, point_set, start_end_derivatives, ctrl_pts);//从路径点得到控制点ctrl_pts

    vector<vector<Eigen::Vector3d>> a_star_pathes;//A*规划路径
    a_star_pathes = bspline_optimizer_rebound_->initControlPoints(ctrl_pts, true);//得到pv对 障碍物段和搜索的无碰撞路径
    //在这一步骤中，已经将当前的控制点存入了cps_pt中，即在类中已经进行了初始化
    t_init = ros::Time::now() - t_start;

    static int vis_id = 0;
    visualization_->displayInitPathList(point_set, 0.2, 0);//可视化路径点
    visualization_->displayAStarList(a_star_pathes, vis_id);//可视化搜索到的A*无碰撞路径

    t_start = ros::Time::now();

    /*** STEP 2: OPTIMIZE 优化***/
    bool flag_step_1_success = bspline_optimizer_rebound_->BsplineOptimizeTrajRebound(ctrl_pts, ts);//规划是否成功(rebound所得到的无碰撞路径)
    cout << "first_optimize_step_success=" << flag_step_1_success << endl;
    if (!flag_step_1_success)//如果没成功
    {
      // visualization_->displayOptimalList( ctrl_pts, vis_id );
      continous_failures_count_++;
      return false;
    }
    //visualization_->displayOptimalList( ctrl_pts, vis_id );

    t_opt = ros::Time::now() - t_start;
    t_start = ros::Time::now();

    /*** STEP 3: REFINE(RE-ALLOCATE TIME) IF NECESSARY  如果有必要 进行时间重置***/
    UniformBspline pos = UniformBspline(ctrl_pts, 3, ts);//由ts和当前控制点产生的均匀b样条曲线
    pos.setPhysicalLimits(pp_.max_vel_, pp_.max_acc_, pp_.feasibility_tolerance_);//设置物理约束

    double ratio;
    bool flag_step_2_success = true;
    if (!pos.checkFeasibility(ratio, false))//如果动力学不可行
    {
      cout << "Need to reallocate time." << endl;

      Eigen::MatrixXd optimal_control_points;
      flag_step_2_success = refineTrajAlgo(pos, start_end_derivatives, ratio, ts, optimal_control_points);//进行refine后是否成功
      if (flag_step_2_success)
        pos = UniformBspline(optimal_control_points, 3, ts);//如果成功了，更新控制点
    }

    if (!flag_step_2_success)//如果没成功
    {
      printf("\033[34mThis refined trajectory hits obstacles. It doesn't matter if appeares occasionally. But if continously appearing, Increase parameter \"lambda_fitness\".\n\033[0m");
      continous_failures_count_++;
      return false;
    }

    t_refine = ros::Time::now() - t_start;

    // save planned results
    updateTrajInfo(pos, ros::Time::now());//更新local_data_

    cout << "total time:\033[42m" << (t_init + t_opt + t_refine).toSec() << "\033[0m,optimize:" << (t_init + t_opt).toSec() << ",refine:" << t_refine.toSec() << endl;

    // success. YoY
    continous_failures_count_ = 0;
    return true;
  }

//给出停止位置，控制点为6个这个停止点，规划一条b样条轨迹并且存储到local.data中（相同的控制点生成的轨迹起始就是一个点，不管何时都在这个位置。
  bool EGOPlannerManager::EmergencyStop(Eigen::Vector3d stop_pos)//规划出一条在stop_pos处停止的三阶均匀B样条轨迹并赋值给local_data_
  {
    Eigen::MatrixXd control_points(3, 6);
    for (int i = 0; i < 6; i++)
    {
      control_points.col(i) = stop_pos;
    }

    updateTrajInfo(UniformBspline(control_points, 3, 1.0), ros::Time::now());

    return true;
  }

  bool EGOPlannerManager::planGlobalTrajWaypoints(const Eigen::Vector3d &start_pos, const Eigen::Vector3d &start_vel, const Eigen::Vector3d &start_acc,
                                                  const std::vector<Eigen::Vector3d> &waypoints, const Eigen::Vector3d &end_vel, const Eigen::Vector3d &end_acc)
  {

    // 生成全局参考轨迹 generate global reference trajectory

    vector<Eigen::Vector3d> points;//点集合
    points.push_back(start_pos);//先压入起点

    for (size_t wp_i = 0; wp_i < waypoints.size(); wp_i++)//循环路径点
    {
      points.push_back(waypoints[wp_i]);//加到代点集合中
    }

    double total_len = 0;//初始化路径总长度
    total_len += (start_pos - waypoints[0]).norm();
    for (size_t i = 0; i < waypoints.size() - 1; i++)
    {
      total_len += (waypoints[i + 1] - waypoints[i]).norm();
    }//计算路径总长度

    //  如果太远，插入中间点  insert intermediate points if too far
    vector<Eigen::Vector3d> inter_points;
    double dist_thresh = max(total_len / 8, 4.0);

    for (size_t i = 0; i < points.size() - 1; ++i)//循环所有点
    {
      inter_points.push_back(points.at(i));
      double dist = (points.at(i + 1) - points.at(i)).norm();

      if (dist > dist_thresh)//如果路标点间距太大
      {
        int id_num = floor(dist / dist_thresh) + 1;

        for (int j = 1; j < id_num; ++j)
        {
          Eigen::Vector3d inter_pt =
              points.at(i) * (1.0 - double(j) / id_num) + points.at(i + 1) * double(j) / id_num;
          inter_points.push_back(inter_pt);//等间隔插入几个点
        }
      }
    }

    inter_points.push_back(points.back());//压入最后一个点

    // for ( int i=0; i<inter_points.size(); i++ )
    // {
    //   cout << inter_points[i].transpose() << endl;
    // }

    // write position matrix
    int pt_num = inter_points.size();
    Eigen::MatrixXd pos(3, pt_num);
    for (int i = 0; i < pt_num; ++i)
      pos.col(i) = inter_points[i];//写成矩阵形式

    Eigen::Vector3d zero(0, 0, 0);
    Eigen::VectorXd time(pt_num - 1);//每一段的时间，段数是点数-1
    for (int i = 0; i < pt_num - 1; ++i)
    {
      time(i) = (pos.col(i + 1) - pos.col(i)).norm() / (pp_.max_vel_);
    }//以最大速度飞行

    time(0) *= 2.0;
    time(time.rows() - 1) *= 2.0;//始末段分配的时间为2倍原来的时间

    PolynomialTraj gl_traj;//多项式轨迹
    if (pos.cols() >= 3)
      gl_traj = PolynomialTraj::minSnapTraj(pos, start_vel, end_vel, start_acc, end_acc, time);
    //如果路标点大于等于3，则用PolynomialTraj::minSnapTraj minimum snap方法计算出连续轨迹
    else if (pos.cols() == 2)
      gl_traj = PolynomialTraj::one_segment_traj_gen(start_pos, start_vel, start_acc, pos.col(1), end_vel, end_acc, time(0));
    //否则使用PolynomialTraj::one_segment_traj_gen多项式拟合计算出连续轨迹
    else
      return false;

    auto time_now = ros::Time::now();
    global_data_.setGlobalTraj(gl_traj, time_now);

    return true;
  }

  bool EGOPlannerManager::planGlobalTraj(const Eigen::Vector3d &start_pos, const Eigen::Vector3d &start_vel, const Eigen::Vector3d &start_acc,
                                         const Eigen::Vector3d &end_pos, const Eigen::Vector3d &end_vel, const Eigen::Vector3d &end_acc)
  {

    // 和上边类似只不过这里只有始末两个点generate global reference trajectory

    vector<Eigen::Vector3d> points;
    points.push_back(start_pos);
    points.push_back(end_pos);

    // insert intermediate points if too far
    vector<Eigen::Vector3d> inter_points;
    const double dist_thresh = 4.0;

    for (size_t i = 0; i < points.size() - 1; ++i)
    {
      inter_points.push_back(points.at(i));
      double dist = (points.at(i + 1) - points.at(i)).norm();

      if (dist > dist_thresh)
      {
        int id_num = floor(dist / dist_thresh) + 1;

        for (int j = 1; j < id_num; ++j)
        {
          Eigen::Vector3d inter_pt =
              points.at(i) * (1.0 - double(j) / id_num) + points.at(i + 1) * double(j) / id_num;
          inter_points.push_back(inter_pt);
        }
      }
    }

    inter_points.push_back(points.back());

    // write position matrix
    int pt_num = inter_points.size();
    Eigen::MatrixXd pos(3, pt_num);
    for (int i = 0; i < pt_num; ++i)
      pos.col(i) = inter_points[i];

    Eigen::Vector3d zero(0, 0, 0);
    Eigen::VectorXd time(pt_num - 1);
    for (int i = 0; i < pt_num - 1; ++i)
    {
      time(i) = (pos.col(i + 1) - pos.col(i)).norm() / (pp_.max_vel_);
    }

    time(0) *= 2.0;
    time(time.rows() - 1) *= 2.0;

    PolynomialTraj gl_traj;
    if (pos.cols() >= 3)
      gl_traj = PolynomialTraj::minSnapTraj(pos, start_vel, end_vel, start_acc, end_acc, time);
    else if (pos.cols() == 2)
      gl_traj = PolynomialTraj::one_segment_traj_gen(start_pos, start_vel, start_acc, end_pos, end_vel, end_acc, time(0));
    else
      return false;

    auto time_now = ros::Time::now();
    global_data_.setGlobalTraj(gl_traj, time_now);

    return true;
  }

  bool EGOPlannerManager::refineTrajAlgo(UniformBspline &traj, vector<Eigen::Vector3d> &start_end_derivative, double ratio, double &ts, Eigen::MatrixXd &optimal_control_points)
  //refine的过程
  {
    double t_inc;

    Eigen::MatrixXd ctrl_pts; //控制点 = traj.getControlPoint()

    // std::cout << "ratio: " << ratio << std::endl;
    reparamBspline(traj, start_end_derivative, ratio, ctrl_pts, ts, t_inc);//放大时间后的b样条曲线控制点在ctrl_pts

    traj = UniformBspline(ctrl_pts, 3, ts);//在新的控制点上再次生成b样条曲线

    double t_step = traj.getTimeSum() / (ctrl_pts.cols() - 3);
    bspline_optimizer_rebound_->ref_pts_.clear();//清空ref向量
    for (double t = 0; t < traj.getTimeSum() + 1e-4; t += t_step)//循环轨迹
      bspline_optimizer_rebound_->ref_pts_.push_back(traj.evaluateDeBoorT(t));//Φs上的点，即初始的轨迹。迭代过程中Φf会逐渐变化，但是它不变，是基准。要衡量优化后的点和它之间的拟合度

    bool success = bspline_optimizer_rebound_->BsplineOptimizeTrajRefine(ctrl_pts, ts, optimal_control_points);//这一步中将当前的控制点ctrl_pts赋予到了cps_.pt中

    return success;//返回是否优化成功
  }

  void EGOPlannerManager::updateTrajInfo(const UniformBspline &position_traj, const ros::Time time_now)
  //更新参数
  {
    local_data_.start_time_ = time_now;
    local_data_.position_traj_ = position_traj;
    local_data_.velocity_traj_ = local_data_.position_traj_.getDerivative();
    local_data_.acceleration_traj_ = local_data_.velocity_traj_.getDerivative();
    local_data_.start_pos_ = local_data_.position_traj_.evaluateDeBoorT(0.0);
    local_data_.duration_ = local_data_.position_traj_.getTimeSum();
    local_data_.traj_id_ += 1;
  }

  void EGOPlannerManager::reparamBspline(UniformBspline &bspline, vector<Eigen::Vector3d> &start_end_derivative, double ratio,
                                         Eigen::MatrixXd &ctrl_pts, double &dt, double &time_inc)
   //在原来的b样条轨迹的基础上拉长时间，然后遍历新的时间下的曲线寻找路径点。再由当前的路径点重新生成控制点
  {
    double time_origin = bspline.getTimeSum();//初始的b样条曲线时间
    int seg_num = bspline.getControlPoint().cols() - 3;//
    // double length = bspline.getLength(0.1);
    // int seg_num = ceil(length / pp_.ctrl_pt_dist);

    bspline.lengthenTime(ratio);//延长时间
    double duration = bspline.getTimeSum();//延长后的b样条曲线时间
    dt = duration / double(seg_num);//等间隔取N个点，N为bspline原来的控制点数-3 即seg_num
    time_inc = duration - time_origin;//轨迹增加的时长

    vector<Eigen::Vector3d> point_set;//路径点
    for (double time = 0.0; time <= duration + 1e-4; time += dt)
    {
      point_set.push_back(bspline.evaluateDeBoorT(time));//加入路径点
    }
    UniformBspline::parameterizeToBspline(dt, point_set, start_end_derivative, ctrl_pts);//由路径点转为控制点
  }

} // namespace ego_planner
