
#include <plan_manage/ego_replan_fsm.h>

namespace ego_planner
{

  void EGOReplanFSM::init(ros::NodeHandle &nh) //初始化，给定ros节点
  {
    current_wp_ = 0;
    exec_state_ = FSM_EXEC_STATE::INIT;
    have_target_ = false;
    have_odom_ = false;

    /*  fsm param  */
    nh.param("fsm/flight_type", target_type_, -1);
    nh.param("fsm/thresh_replan", replan_thresh_, -1.0);
    nh.param("fsm/thresh_no_replan", no_replan_thresh_, -1.0);
    nh.param("fsm/planning_horizon", planning_horizen_, -1.0);
    nh.param("fsm/planning_horizen_time", planning_horizen_time_, -1.0);
    nh.param("fsm/emergency_time_", emergency_time_, 1.0);

    nh.param("fsm/waypoint_num", waypoint_num_, -1); //固定的路径点个数
    for (int i = 0; i < waypoint_num_; i++)          //循环，将这些固定的路径点加入矩阵中
    {
      nh.param("fsm/waypoint" + to_string(i) + "_x", waypoints_[i][0], -1.0);
      nh.param("fsm/waypoint" + to_string(i) + "_y", waypoints_[i][1], -1.0);
      nh.param("fsm/waypoint" + to_string(i) + "_z", waypoints_[i][2], -1.0);
    }

    /* initialize main modules */
    //p.reset(q) //q为智能指针要指向的新对象
    //会令智能指针p中存放指针q，即p指向q的空间，而且会释放原来的空间。（默认是delete）
    visualization_.reset(new PlanningVisualization(nh));   //可视化类 
    planner_manager_.reset(new EGOPlannerManager);         // planner类
    planner_manager_->initPlanModules(nh, visualization_); //初始化planner类相关参数

    /* callback */
    exec_timer_ = nh.createTimer(ros::Duration(0.01), &EGOReplanFSM::execFSMCallback, this);          //定时器状态机切换 周期0.01秒
    safety_timer_ = nh.createTimer(ros::Duration(0.05), &EGOReplanFSM::checkCollisionCallback, this); //安全检查 周期0.05秒

    odom_sub_ = nh.subscribe("/odom_world", 1, &EGOReplanFSM::odometryCallback, this); //保存无人机当前里程计信息，包括位置、速度和姿态

    bspline_pub_ = nh.advertise<ego_planner::Bspline>("/planning/bspline", 10);
    data_disp_pub_ = nh.advertise<ego_planner::DataDisp>("/planning/data_display", 100);

    if (target_type_ == TARGET_TYPE::MANUAL_TARGET)
      waypoint_sub_ = nh.subscribe("/waypoint_generator/waypoints", 1, &EGOReplanFSM::waypointCallback, this);
    else if (target_type_ == TARGET_TYPE::PRESET_TARGET)
    {
      ros::Duration(1.0).sleep();
      while (ros::ok() && !have_odom_)
        ros::spinOnce();
      planGlobalTrajbyGivenWps();
    }
    else
      cout << "Wrong target_type_ value! target_type_=" << target_type_ << endl;
  }

  void EGOReplanFSM::planGlobalTrajbyGivenWps() //给出路径点后规划全局轨迹
  {
    std::vector<Eigen::Vector3d> wps(waypoint_num_); //所有路径点
    for (int i = 0; i < waypoint_num_; i++)
    {
      wps[i](0) = waypoints_[i][0];
      wps[i](1) = waypoints_[i][1];
      wps[i](2) = waypoints_[i][2];

      end_pt_ = wps.back(); //终点
    }
    bool success = planner_manager_->planGlobalTrajWaypoints(odom_pos_, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), wps, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
    //从里程计给出的当前位置作为起点 ，起点速度加速度为0 终点速度加速度为0规划一条通过路径点的minimumsnap轨迹,存储在全局轨迹global_data_中

    for (size_t i = 0; i < (size_t)waypoint_num_; i++) //循环所有路径点将其可视化
    {
      visualization_->displayGoalPoint(wps[i], Eigen::Vector4d(0, 0.5, 0.5, 1), 0.3, i);
      ros::Duration(0.001).sleep();
    }

    if (success) //如果成功规划出了minimumsnap轨迹
    {

      /*** display ***/
      constexpr double step_size_t = 0.1;
      int i_end = floor(planner_manager_->global_data_.global_duration_ / step_size_t);
      std::vector<Eigen::Vector3d> gloabl_traj(i_end);
      for (int i = 0; i < i_end; i++)
      {
        gloabl_traj[i] = planner_manager_->global_data_.global_traj_.evaluate(i * step_size_t);
      } //以step_size_t 为时间间隔 得到全局轨迹点

      end_vel_.setZero(); //终点速度设置为0
      have_target_ = true;
      have_new_target_ = true;

      /*** FSM ***/
      // if (exec_state_ == WAIT_TARGET)
      changeFSMExecState(GEN_NEW_TRAJ, "TRIG");
      // else if (exec_state_ == EXEC_TRAJ)
      //   changeFSMExecState(REPLAN_TRAJ, "TRIG");

      // visualization_->displayGoalPoint(end_pt_, Eigen::Vector4d(1, 0, 0, 1), 0.3, 0);
      ros::Duration(0.001).sleep();
      visualization_->displayGlobalPathList(gloabl_traj, 0.1, 0);
      ros::Duration(0.001).sleep();
    }
    else
    {
      ROS_ERROR("Unable to generate global trajectory!");
    }
  }

  void EGOReplanFSM::waypointCallback(const nav_msgs::PathConstPtr &msg) // const nav_msgs::PathConstPtr &msg 代表里程计的一个class指针，存储里程计的信息
  {
    if (msg->poses[0].pose.position.z < -0.1)
      return;

    cout << "Triggered!" << endl;
    trigger_ = true;
    init_pt_ = odom_pos_; //现在里程计的位置作为起点位置

    bool success = false;
    end_pt_ << msg->poses[0].pose.position.x, msg->poses[0].pose.position.y, 1.0; //终点位置
    success = planner_manager_->planGlobalTraj(odom_pos_, odom_vel_, Eigen::Vector3d::Zero(), end_pt_, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
    //规划一条从起点到终点的全局minimumsnap轨迹（注意起点加速度和终点速度加速度为0）,存储在全局轨迹global_data_中

    visualization_->displayGoalPoint(end_pt_, Eigen::Vector4d(0, 0.5, 0.5, 1), 0.3, 0); //可视化

    if (success)
    {

      /*** display ***/
      constexpr double step_size_t = 0.1;
      int i_end = floor(planner_manager_->global_data_.global_duration_ / step_size_t);
      vector<Eigen::Vector3d> gloabl_traj(i_end);
      for (int i = 0; i < i_end; i++)
      {
        gloabl_traj[i] = planner_manager_->global_data_.global_traj_.evaluate(i * step_size_t);
      } //以step_size_t 为时间间隔 得到全局轨迹点

      end_vel_.setZero();
      have_target_ = true;
      have_new_target_ = true;

      /*** FSM ***/
      if (exec_state_ == WAIT_TARGET) //如果现在的状态是等待目标点，就将状态变为正在生成轨迹
        changeFSMExecState(GEN_NEW_TRAJ, "TRIG");
      else if (exec_state_ == EXEC_TRAJ) //如果现在的状态是正在执行轨迹，就将状态变为重规划轨迹
        changeFSMExecState(REPLAN_TRAJ, "TRIG");

      // visualization_->displayGoalPoint(end_pt_, Eigen::Vector4d(1, 0, 0, 1), 0.3, 0);
      visualization_->displayGlobalPathList(gloabl_traj, 0.1, 0);
    }
    else
    {
      ROS_ERROR("Unable to generate global trajectory!");
    }
  }

  void EGOReplanFSM::odometryCallback(const nav_msgs::OdometryConstPtr &msg)
  { //由里程计得到当前的位置速度和方向
    //方向是什么意思？？？？？？？？？？？？？？？？？
    odom_pos_(0) = msg->pose.pose.position.x;
    odom_pos_(1) = msg->pose.pose.position.y;
    odom_pos_(2) = msg->pose.pose.position.z;

    odom_vel_(0) = msg->twist.twist.linear.x;
    odom_vel_(1) = msg->twist.twist.linear.y;
    odom_vel_(2) = msg->twist.twist.linear.z;

    // odom_acc_ = estimateAcc( msg );

    odom_orient_.w() = msg->pose.pose.orientation.w;
    odom_orient_.x() = msg->pose.pose.orientation.x;
    odom_orient_.y() = msg->pose.pose.orientation.y;
    odom_orient_.z() = msg->pose.pose.orientation.z;

    have_odom_ = true;
  }

  void EGOReplanFSM::changeFSMExecState(FSM_EXEC_STATE new_state, string pos_call)
  //改变状态机的状态 切换到new_state模式，若新状态与当前状态一致，则continously_called_times_++ ，否则continously_called_times_置1。
  {

    if (new_state == exec_state_)
      continously_called_times_++;
    else
      continously_called_times_ = 1;

    static string state_str[7] = {"INIT", "WAIT_TARGET", "GEN_NEW_TRAJ", "REPLAN_TRAJ", "EXEC_TRAJ", "EMERGENCY_STOP"};
    int pre_s = int(exec_state_);
    exec_state_ = new_state;
    cout << "[" + pos_call + "]: from " + state_str[pre_s] + " to " + state_str[int(new_state)] << endl;
  }

  std::pair<int, EGOReplanFSM::FSM_EXEC_STATE> EGOReplanFSM::timesOfConsecutiveStateCalls()
  {
    return std::pair<int, FSM_EXEC_STATE>(continously_called_times_, exec_state_);
  } //返回当前状态被连续切换的次数

  void EGOReplanFSM::printFSMExecState()
  {
    static string state_str[7] = {"INIT", "WAIT_TARGET", "GEN_NEW_TRAJ", "REPLAN_TRAJ", "EXEC_TRAJ", "EMERGENCY_STOP"};

    cout << "[FSM]: state: " + state_str[int(exec_state_)] << endl;
  } //打印出当前的状态

  void EGOReplanFSM::execFSMCallback(const ros::TimerEvent &e)
  //进入状态机时会类似于关中断，防止在中断过程中响应其他中断。
  //每秒会判断是否有里程计信息和目标位置。
  {

    static int fsm_num = 0;
    fsm_num++;
    if (fsm_num == 100)
    {
      printFSMExecState();
      if (!have_odom_) //没获得当前状态
        cout << "no odom." << endl;
      if (!trigger_) //没获得目标点
        cout << "wait for goal." << endl;
      fsm_num = 0;
    }

    switch (exec_state_) //对当前状态进行判断
    {
    case INIT: //如果是INIT
    {
      if (!have_odom_) //没获得当前状态直接返回
      {
        return;
      }
      if (!trigger_) //没获得目标点直接返回
      {
        return;
      }
      changeFSMExecState(WAIT_TARGET, "FSM"); //改变当前状态为WAIT_TARGET
      break;
    }

    case WAIT_TARGET: //如果是WAIT_TARGET
    {
      if (!have_target_) //全局轨迹没生成直接返回
        return;
      else
      {
        changeFSMExecState(GEN_NEW_TRAJ, "FSM"); //改变当前状态为GEN_NEW_TRAJ
      }
      break;
    }

    case GEN_NEW_TRAJ: //如果是GEN_NEW_TRAJ
    {
      start_pt_ = odom_pos_; //起点为当前的位置
      start_vel_ = odom_vel_;
      start_acc_.setZero(); //起点位置速度加速度

      // Eigen::Vector3d rot_x = odom_orient_.toRotationMatrix().block(0, 0, 3, 1);
      // start_yaw_(0)         = atan2(rot_x(1), rot_x(0));
      // start_yaw_(1) = start_yaw_(2) = 0.0;

      bool flag_random_poly_init;
      if (timesOfConsecutiveStateCalls().first == 1) //如果状态连续切换的次数为1
        flag_random_poly_init = false;
      else
        flag_random_poly_init = true;

      bool success = callReboundReplan(true, flag_random_poly_init); //进行优化
      if (success)                                                   //如果成功了
      {

        changeFSMExecState(EXEC_TRAJ, "FSM"); // 状态变为EXEC_TRAJ
        flag_escape_emergency_ = true;
      }
      else
      {
        changeFSMExecState(GEN_NEW_TRAJ, "FSM"); //否则状态还是GEN_NEW_TRAJ
      }
      break;
    }

    case REPLAN_TRAJ: //如果REPLAN_TRAJ
    {

      if (planFromCurrentTraj()) //如果根据当前的轨迹进行重优化成功了
      {
        changeFSMExecState(EXEC_TRAJ, "FSM"); //改变状态为EXEC_TRAJ
      }
      else
      {
        changeFSMExecState(REPLAN_TRAJ, "FSM"); //否则状态还是REPLAN_TRAJ
      }

      break;
    }

    case EXEC_TRAJ: //如果EXEC_TRAJ
    {
      /* determine if need to replan 看是否需要重规划*/
      LocalTrajData *info = &planner_manager_->local_data_;//局部规划的轨迹
      ros::Time time_now = ros::Time::now();
      double t_cur = (time_now - info->start_time_).toSec(); //当前时间
      t_cur = min(info->duration_, t_cur);//比较当前时间和局部轨迹的结束时间的大小，取其中小的

      Eigen::Vector3d pos = info->position_traj_.evaluateDeBoorT(t_cur); //局部轨迹b样条曲线位置

      /* && (end_pt_ - pos).norm() < 0.5 */
      if (t_cur > info->duration_ - 1e-2)
      {
        have_target_ = false;

        changeFSMExecState(WAIT_TARGET, "FSM"); //超时了，即本次的局部规划已经到了终点,状态改为WAIT_TARGET
        return;
      }

      //如果到起点的距离太小（刚出发）或者到全局终点的距离太小就不需要重规划了,继续执行轨迹
      else if ((end_pt_ - pos).norm() < no_replan_thresh_) //如果全局终点到当前点的距离<设定的不重规划的距离
      {
        // cout << "near end" << endl;
        return;
      }
      else if ((info->start_pos_ - pos).norm() < replan_thresh_) //局部规划起点到当前点的距离<设定的重规划的距离
      {
        // cout << "near start" << endl;
        return;
      }
      else
      {
        changeFSMExecState(REPLAN_TRAJ, "FSM"); //否则改变状态为REPLAN_TRAJ
      }
      break;
    }

    case EMERGENCY_STOP: //如果EMERGENCY_STOP
    {

      if (flag_escape_emergency_) // 避免重复调用 Avoiding repeated calls
      {
        callEmergencyStop(odom_pos_); //在当前位置紧急停止
      }
      else
      {
        if (odom_vel_.norm() < 0.1)                //当前速度很小
          changeFSMExecState(GEN_NEW_TRAJ, "FSM"); //改变状态为GEN_NEW_TRAJ
      }

      flag_escape_emergency_ = false; //紧急停止为false
      break;
    }
    }

    data_disp_.header.stamp = ros::Time::now();
    data_disp_pub_.publish(data_disp_);
  }

  bool EGOReplanFSM::planFromCurrentTraj()
  {

    LocalTrajData *info = &planner_manager_->local_data_;
    ros::Time time_now = ros::Time::now();
    double t_cur = (time_now - info->start_time_).toSec(); //当前时间

    // cout << "info->velocity_traj_=" << info->velocity_traj_.get_control_points() << endl;

    start_pt_ = info->position_traj_.evaluateDeBoorT(t_cur);
    start_vel_ = info->velocity_traj_.evaluateDeBoorT(t_cur);
    start_acc_ = info->acceleration_traj_.evaluateDeBoorT(t_cur); //获得起点位置速度加速度

    bool success = callReboundReplan(false, false); // rebound优化

    if (!success) //没成功
    {
      success = callReboundReplan(true, false); // rebound优化
      // changeFSMExecState(EXEC_TRAJ, "FSM");
      if (!success) //没成功
      {
        success = callReboundReplan(true, true); // rebound优化
        if (!success)                            //没成功
        {
          return false; //返回没成功
        }
      }
    }

    return true;
  }

  void EGOReplanFSM::checkCollisionCallback(const ros::TimerEvent &e)
  //
  {
    LocalTrajData *info = &planner_manager_->local_data_; //获得局部轨迹
    auto map = planner_manager_->grid_map_;               //获得局部地图

    if (exec_state_ == WAIT_TARGET || info->start_time_.toSec() < 1e-5)
      return;

    /* ---------- check trajectory 检查轨迹---------- */
    constexpr double time_step = 0.01;
    double t_cur = (ros::Time::now() - info->start_time_).toSec();
    double t_2_3 = info->duration_ * 2 / 3;
    for (double t = t_cur; t < info->duration_; t += time_step) //循环轨迹
    {
      if (t_cur < t_2_3 && t >= t_2_3) // If t_cur < t_2_3, only the first 2/3 partition of the trajectory is considered valid and will get checked.
        break;

      if (map->getInflateOccupancy(info->position_traj_.evaluateDeBoorT(t))) //如果发生碰撞
      {
        if (planFromCurrentTraj()) // 进行重规划并且成功了 Make a chance
        {
          changeFSMExecState(EXEC_TRAJ, "SAFETY"); //改变当前状态为EXEC_TRAJ
          return;
        }
        else
        {
          if (t - t_cur < emergency_time_) // 如果从当前的位置到这个发现的障碍物所经历的时间<emergency_time_就停止（说明眼前就有个障碍物，需要停止）0.8s of emergency time
          {
            ROS_WARN("Suddenly discovered obstacles. emergency stop! time=%f", t - t_cur); //发现了紧急障碍，停止
            changeFSMExecState(EMERGENCY_STOP, "SAFETY");
          }
          else
          {
            // ROS_WARN("current traj in collision, replan.");
            changeFSMExecState(REPLAN_TRAJ, "SAFETY"); //否则重规划
          }
          return;
        }
        break;
      }
    }
  }

  bool EGOReplanFSM::callReboundReplan(bool flag_use_poly_init, bool flag_randomPolyTraj)
  {

    getLocalTarget(); //获得局部规划终点

    bool plan_success =
        planner_manager_->reboundReplan(start_pt_, start_vel_, start_acc_, local_target_pt_, local_target_vel_, (have_new_target_ || flag_use_poly_init), flag_randomPolyTraj);
    // rebound规划的结果是否成功,并更新local_data中的数据
    have_new_target_ = false;

    cout << "final_plan_success=" << plan_success << endl;

    if (plan_success) //如果成功
    {

      auto info = &planner_manager_->local_data_; //局部规划轨迹的数据

      /* publish traj */
      ego_planner::Bspline bspline; // b样条曲线消息，定义在plan_manage中的msg文件夹中
      bspline.order = 3;
      bspline.start_time = info->start_time_;
      bspline.traj_id = info->traj_id_;

      Eigen::MatrixXd pos_pts = info->position_traj_.getControlPoint(); //获得轨迹控制点
      bspline.pos_pts.reserve(pos_pts.cols());
      for (int i = 0; i < pos_pts.cols(); ++i) //循环控制点
      {
        geometry_msgs::Point pt;
        pt.x = pos_pts(0, i);
        pt.y = pos_pts(1, i);
        pt.z = pos_pts(2, i);
        bspline.pos_pts.push_back(pt);
      }

      Eigen::VectorXd knots = info->position_traj_.getKnot(); //获得时间knot
      bspline.knots.reserve(knots.rows());
      for (int i = 0; i < knots.rows(); ++i)
      {
        bspline.knots.push_back(knots(i));
      }

      bspline_pub_.publish(bspline);

      visualization_->displayOptimalList(info->position_traj_.get_control_points(), 0); //可视化
    }

    return plan_success;
  }

  bool EGOReplanFSM::callEmergencyStop(Eigen::Vector3d stop_pos)
  //输入当前的位置
  {

    planner_manager_->EmergencyStop(stop_pos); //生成紧急停止的b样条轨迹（其实就是一个点）

    auto info = &planner_manager_->local_data_; // local_data_的指针

    /* publish traj */
    ego_planner::Bspline bspline;
    bspline.order = 3;
    bspline.start_time = info->start_time_;
    bspline.traj_id = info->traj_id_;

    Eigen::MatrixXd pos_pts = info->position_traj_.getControlPoint(); //得到控制点
    bspline.pos_pts.reserve(pos_pts.cols());
    for (int i = 0; i < pos_pts.cols(); ++i)
    {
      geometry_msgs::Point pt;
      pt.x = pos_pts(0, i);
      pt.y = pos_pts(1, i);
      pt.z = pos_pts(2, i);
      bspline.pos_pts.push_back(pt);
    } //控制点

    Eigen::VectorXd knots = info->position_traj_.getKnot();
    bspline.knots.reserve(knots.rows());
    for (int i = 0; i < knots.rows(); ++i)
    {
      bspline.knots.push_back(knots(i));
    }

    bspline_pub_.publish(bspline);

    return true;
  }

  void EGOReplanFSM::getLocalTarget() //获得局部规划的target点
  {
    double t;

    double t_step = planning_horizen_ / 20 / planner_manager_->pp_.max_vel_;

    double dist_min = 9999, dist_min_t = 0.0;
    for (t = planner_manager_->global_data_.last_progress_time_; t < planner_manager_->global_data_.global_duration_; t += t_step)
    //循环全局轨迹上的点
    {
      Eigen::Vector3d pos_t = planner_manager_->global_data_.getPosition(t); //获得在全局轨迹上的位置
      double dist = (pos_t - start_pt_).norm();                              //计算到起点的距离

      if (t < planner_manager_->global_data_.last_progress_time_ + 1e-5 && dist > planning_horizen_)
      {
        // todo
        ROS_ERROR("last_progress_time_ ERROR !!!!!!!!!");
        ROS_ERROR("last_progress_time_ ERROR !!!!!!!!!");
        ROS_ERROR("last_progress_time_ ERROR !!!!!!!!!");
        ROS_ERROR("last_progress_time_ ERROR !!!!!!!!!");
        ROS_ERROR("last_progress_time_ ERROR !!!!!!!!!");
        return;
      }
      if (dist < dist_min) //更新最小距离
      {
        dist_min = dist;
        dist_min_t = t;
      }
      if (dist >= planning_horizen_) //如果超出了规划范围
      {
        local_target_pt_ = pos_t;                                        //当前点作为局部规划终点
        planner_manager_->global_data_.last_progress_time_ = dist_min_t; //last_progress_time_设置为全局轨迹中距离当前点最近点所对应的时间
        break;
      }
    }
    if (t > planner_manager_->global_data_.global_duration_) //如果超出了时间范围 Last global point
    {
      local_target_pt_ = end_pt_; //全局规划终点作为局部规划的终点
    }

    //如果全局终点到局部轨迹的目标点的距离小于无人机以做大速度和加速度可飞行的距离，则局部轨迹的目标点速度设为0；否则局部轨迹的目标点速度为全局轨迹在该点的速度
    if ((end_pt_ - local_target_pt_).norm() < (planner_manager_->pp_.max_vel_ * planner_manager_->pp_.max_vel_) / (2 * planner_manager_->pp_.max_acc_))
    {
      // local_target_vel_ = (end_pt_ - init_pt_).normalized() * planner_manager_->pp_.max_vel_ * (( end_pt_ - local_target_pt_ ).norm() / ((planner_manager_->pp_.max_vel_*planner_manager_->pp_.max_vel_)/(2*planner_manager_->pp_.max_acc_)));
      // cout << "A" << endl;
      local_target_vel_ = Eigen::Vector3d::Zero();
    }
    else
    {
      local_target_vel_ = planner_manager_->global_data_.getVelocity(t);
      // cout << "AA" << endl;
    }
  }

} // namespace ego_planner
