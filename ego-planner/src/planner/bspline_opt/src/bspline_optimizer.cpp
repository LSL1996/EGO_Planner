#include "bspline_opt/bspline_optimizer.h"
#include "bspline_opt/gradient_descent_optimizer.h"
// using namespace std;

namespace ego_planner
{

  void BsplineOptimizer::setParam(ros::NodeHandle &nh)//设置参数
  {
    //四个优化项前的系数lamada
    nh.param("optimization/lambda_smooth", lambda1_, -1.0);
    nh.param("optimization/lambda_collision", lambda2_, -1.0);
    nh.param("optimization/lambda_feasibility", lambda3_, -1.0);
    nh.param("optimization/lambda_fitness", lambda4_, -1.0);

    //安全距离 最大加速度和速度
    nh.param("optimization/dist0", dist0_, -1.0);
    nh.param("optimization/max_vel", max_vel_, -1.0);
    nh.param("optimization/max_acc", max_acc_, -1.0);
//b样条曲线的次数
    nh.param("optimization/order", order_, 3);
  }

  void BsplineOptimizer::setEnvironment(const GridMap::Ptr &env)//
  {
    this->grid_map_ = env;//加载地图
  }

  void BsplineOptimizer::setControlPoints(const Eigen::MatrixXd &points)
  {
    cps_.points = points;//给控制点类中的点矩阵赋值
  }

  void BsplineOptimizer::setBsplineInterval(const double &ts) { bspline_interval_ = ts; }
//设置b样条曲线的knot span（时间间隔）

  /* This function is very similar to check_collision_and_rebound(). 
   * It was written separately, just because I did it once and it has been running stably since March 2020.
   * But I will merge then someday.*/

  //初始化控制点
  std::vector<std::vector<Eigen::Vector3d>> BsplineOptimizer::initControlPoints(Eigen::MatrixXd &init_points, bool flag_first_init /*= true*/)
  {

    if (flag_first_init)//如果是第一次初始化
    {
      cps_.clearance = dist0_;//安全距离
      cps_.resize(init_points.cols());//resize
      cps_.points = init_points;//初始化控制点矩阵
    }

    /*** Segment the initial trajectory according to obstacles 根据障碍物分割初始轨迹***/
    constexpr int ENOUGH_INTERVAL = 2;//值不会改变并且在编译过程就能得到计算结果的表达式https://blog.csdn.net/janeqi1987/article/details/103542802
    //ENOUGH_INTERVAL代表空白点的个数。为了保证障碍物轨迹不是很离散。
    //规定障碍物轨迹两边必须空出3无碰撞（free）个点。
    //即一段障碍物轨迹如果要开始，障碍物轨迹起点前必须有三个free。
    //如果要结束，必须连续3个free出现
    //（但是注意，如果在第一段或者最后一段的话可以不遵循此规律。
    //最开始一旦有障碍物出现就开始，最后一旦有free点出现就结束。否则可能会出现无法开始或结束的现象

    //每一次迭代的精度  代表这两个控制点的连线直线离散化成很多的点
    double step_size = grid_map_->getResolution() / ((init_points.col(0) - init_points.rightCols(1)).norm() / (init_points.cols() - 1)) / 2;
    //grid_map_->getResolution()得到地图比例尺 目的是从实际地图转换到grid地图
    //((init_points.col(0) - init_points.rightCols(1)).norm()第一列和最后一列差再取模 代表初始点和最后一个点的欧氏距离
    //(init_points.cols() - 1)列数-1
    //((init_points.col(0) - init_points.rightCols(1)).norm() / (init_points.cols() - 1))/2代表起点到终点的距离/段数/2，即每一段的长度的一半
    //

    
    vector<std::pair<int, int>> segment_ids;//存储障碍物段的起点和终点index
    int same_occ_state_times = ENOUGH_INTERVAL + 1;//
    bool occ, last_occ = false;//当前点和上一个点是否为障碍物的标志
    bool flag_got_start = false, flag_got_end = false, flag_got_end_maybe = false;
    //是否获得起点index标志
    // 是否可能获得终点index标志（因为在障碍点后出现空白点不一定就要结束，只是可能，还要看是不是满足轨迹后连续3个free点出现的准则） 
    //是否获得终点index的标志（真的结束）

    //获得终止控制点的角标,仅仅检查前三分之二的点（整体的定义是不考虑开头和结束的order_个点）
    int i_end = (int)init_points.cols() - order_ - ((int)init_points.cols() - 2 * order_) / 3; // only check closed 2/3 points.
    //整体结束的角标-（控制点的个数-2* b样条曲线阶数即整体控制点的个数）/3
    
    
    for (int i = order_; i <= i_end; ++i)//循环需要考虑的控制点  i_start就是order_
    {
      for (double a = 1.0; a >= 0.0; a -= step_size)//循环迭代,将两个控制点的连线直线离散化成很多的点
      {
        occ = grid_map_->getInflateOccupancy(a * init_points.col(i - 1) + (1 - a) * init_points.col(i));
        //看在栅格地图中所对应的是否是障碍物 occ=1是障碍物 occ=0不是障碍物
        //(a * init_points.col(i - 1) + (1 - a) * init_points.col(i))代表这两个控制点的连线上的点，a是负责将直线连线离散化成很多的点
        // cout << setprecision(5);
        // cout << (a * init_points.col(i-1) + (1-a) * init_points.col(i)).transpose() << " occ1=" << occ << endl;

        if (occ && !last_occ)//如果这个点是障碍物，并且上一个点不是障碍物
        {
          if (same_occ_state_times > ENOUGH_INTERVAL || i == order_)//前面出现了3个空白点或者在第一段（第一段内出现障碍点那么直接把第一个点作为障碍物起点）
          {
            in_id = i - 1;//得到起始点的角标id
            flag_got_start = true;//找到起点的标志符设置为true
          }
          same_occ_state_times = 0;//free点的个数设置为0
          flag_got_end_maybe = false; // terminate in advance 提前终止标志符为false
        }
        else if (!occ && last_occ)//如果这个点不是障碍物，并且上一个点是障碍物
        {
          out_id = i;//得到结束点的角标id
          flag_got_end_maybe = true;//提前终止标志符为true
          same_occ_state_times = 0;//free点的个数设置为0
        }
        else//当都不是障碍物或者都是障碍物的时候，+1
        {
          ++same_occ_state_times;//障碍点或者free点+1
          //注意：其实这个变量也代表连续的障碍点个数，但是由于不管碰到开始还是结束点，它都会归0
          //没有用连续障碍物的个数去判断任何事情，只用到了free点的个数
        }

        if (flag_got_end_maybe && (same_occ_state_times > ENOUGH_INTERVAL || (i == (int)init_points.cols() - order_)))
        //如果提前终止标志符为true 并且后面出现了3个空白点或者在最后一段
        {
          flag_got_end_maybe = false;//提前终止标志符为false
          flag_got_end = true;//终止标志符为true
        }

        last_occ = occ;//当前点的occ作为上一个occ 迭代

        if (flag_got_start && flag_got_end)//如果找到起点并且找到终点的标志符均为true
        {
          flag_got_start = false;
          flag_got_end = false;
          segment_ids.push_back(std::pair<int, int>(in_id, out_id));//把每一段障碍轨迹的起点id和终点id放入其中
        }
      }
    }

    /*** a star search  A*算法搜索***/
    vector<vector<Eigen::Vector3d>> a_star_pathes;//定义一个向量，每一个向量中放着一段A*算法求得的路径（每段路径也是一个vector，里面存者很多点的三维坐标）
    for (size_t i = 0; i < segment_ids.size(); ++i)//循环每一段障碍物路径
    {
      //cout << "in=" << in.transpose() << " out=" << out.transpose() << endl;
      Eigen::Vector3d in(init_points.col(segment_ids[i].first)), out(init_points.col(segment_ids[i].second));
      //起点in为障碍物段的起点，终点out为障碍物段的终点
      if (a_star_->AstarSearch(/*(in-out).norm()/10+0.05*/ 0.1, in, out))//如果A*算法寻找出来了一条路径
      {
        a_star_pathes.push_back(a_star_->getPath());//把这条路径加入a_star_pathes中
      }
      else
      {
        ROS_ERROR("a star error, force return!");
        return a_star_pathes;
      }
    }

//虽然写了这段代码，但是好像并没用，因为MINIMUM_PERCENT设置成了0?????????
    /*** calculate bounds 计算边界    目的是什么？？？***/
    int id_low_bound, id_up_bound;//定义边界上下限
    vector<std::pair<int, int>> bounds(segment_ids.size());//定义和碰撞段数相同的边界向量
    for (size_t i = 0; i < segment_ids.size(); i++)//循环每一段碰撞段
    {

      if (i == 0) // first segment第一段
      {
        id_low_bound = order_;//下限为第一个点的角标
        if (segment_ids.size() > 1)//如果障碍物轨迹段数>1
        {
          id_up_bound = (int)(((segment_ids[0].second + segment_ids[1].first) - 1.0f) / 2); // id_up_bound : -1.0f fix()
        }//上限为（第一段障碍轨迹终点+第二条障碍轨迹起点-1）/2
        else//如果障碍物轨迹段数<=1
        {
          id_up_bound = init_points.cols() - order_ - 1;//上限为最后一个考虑的控制点
        }
      }
      else if (i == segment_ids.size() - 1) // 最后一段last segment, i != 0 here
      {
        id_low_bound = (int)(((segment_ids[i].first + segment_ids[i - 1].second) + 1.0f) / 2); // id_low_bound : +1.0f ceil()
        id_up_bound = init_points.cols() - order_ - 1;
      }
      else
      {
        id_low_bound = (int)(((segment_ids[i].first + segment_ids[i - 1].second) + 1.0f) / 2); // id_low_bound : +1.0f ceil()
        id_up_bound = (int)(((segment_ids[i].second + segment_ids[i + 1].first) - 1.0f) / 2);  // id_up_bound : -1.0f fix()
      }

      bounds[i] = std::pair<int, int>(id_low_bound, id_up_bound);
    }//存储每一段的上下限

    // cout << "+++++++++" << endl;
    // for ( int j=0; j<bounds.size(); ++j )
    // {
    //   cout << bounds[j].first << "  " << bounds[j].second << endl;
    // }

    /*** Adjust segment length 调整每一段的长度***/
    vector<std::pair<int, int>> final_segment_ids(segment_ids.size());//最终障碍物段的起点和终点index
    constexpr double MINIMUM_PERCENT = 0.0; // 每一段都保证有足够的点来产生足够的推力 Each segment is guaranteed to have sufficient points to generate sufficient thrust
    int minimum_points = round(init_points.cols() * MINIMUM_PERCENT), num_points;
    for (size_t i = 0; i < segment_ids.size(); i++)//循环每一段障碍物路径
    {
      /*** Adjust segment length ***/
      num_points = segment_ids[i].second - segment_ids[i].first + 1;//每一段障碍物轨迹中控制点的个数
      //cout << "i = " << i << " first = " << segment_ids[i].first << " second = " << segment_ids[i].second << endl;
      if (num_points < minimum_points)//比最少的控制点数还要少
      {
        double add_points_each_side = (int)(((minimum_points - num_points) + 1.0f) / 2);//在障碍物一端增加的控制点数目

        final_segment_ids[i].first = segment_ids[i].first - add_points_each_side >= bounds[i].first ? segment_ids[i].first - add_points_each_side : bounds[i].first;

        final_segment_ids[i].second = segment_ids[i].second + add_points_each_side <= bounds[i].second ? segment_ids[i].second + add_points_each_side : bounds[i].second;
      }
      else//否则不用调整
      {
        final_segment_ids[i].first = segment_ids[i].first;
        final_segment_ids[i].second = segment_ids[i].second;
      }

      //cout << "final:" << "i = " << i << " first = " << final_segment_ids[i].first << " second = " << final_segment_ids[i].second << endl;
    }

    /*** Assign data to each segment 为每个段分配数据***/
    for (size_t i = 0; i < segment_ids.size(); i++)//循环每一段障碍物路径
    {
      // step 1
      for (int j = final_segment_ids[i].first; j <= final_segment_ids[i].second; ++j)//循环这段碰撞轨迹上的控制点
        cps_.flag_temp[j] = false;//这些控制点的flag标志为false

      // step 2
      int got_intersection_id = -1;

      //当这段控制点数>=3时
      for (int j = segment_ids[i].first + 1; j < segment_ids[i].second; ++j)//循环这段碰撞轨迹上除了开头结尾其他的控制点
      {
        Eigen::Vector3d ctrl_pts_law(cps_.points.col(j + 1) - cps_.points.col(j - 1)), intersection_point;
        //ctrl_pts_law:第j个控制点的切向量  intersection_point:第j个控制点的法向量和A*规划的轨迹交点
        int Astar_id = a_star_pathes[i].size() / 2, last_Astar_id; // Let "Astar_id = id_of_the_most_far_away_Astar_point" will be better, but it needs more computation
        //Astar_id:这段碰撞轨迹用A*算法寻路出来的节点个数/2
        double val = (a_star_pathes[i][Astar_id] - cps_.points.col(j)).dot(ctrl_pts_law), last_val = val;
        //val=A*算法轨迹点到第j个控制点的向量点乘切向量
        while (Astar_id >= 0 && Astar_id < (int)a_star_pathes[i].size())
        {
          last_Astar_id = Astar_id;

          if (val >= 0)
            --Astar_id;
          else
            ++Astar_id;

          val = (a_star_pathes[i][Astar_id] - cps_.points.col(j)).dot(ctrl_pts_law);

          if (val * last_val <= 0 && (abs(val) > 0 || abs(last_val) > 0)) // val = last_val = 0.0 is not allowed
          {
            intersection_point =
                a_star_pathes[i][Astar_id] +
                ((a_star_pathes[i][Astar_id] - a_star_pathes[i][last_Astar_id]) *
                 (ctrl_pts_law.dot(cps_.points.col(j) - a_star_pathes[i][Astar_id]) / ctrl_pts_law.dot(a_star_pathes[i][Astar_id] - a_star_pathes[i][last_Astar_id])) // = t
                );//得到第j个控制点的法向量和A*规划的轨迹交点的坐标,含义证明在pad上

            //cout << "i=" << i << " j=" << j << " Astar_id=" << Astar_id << " last_Astar_id=" << last_Astar_id << " intersection_point = " << intersection_point.transpose() << endl;

            got_intersection_id = j;//代表找到了和第j个控制点的法向量和A*规划的轨迹交点的坐标
            break;
          }
        }

        if (got_intersection_id >= 0)//如果找到了第j个控制点的法向量和A*规划的轨迹交点
        {
          cps_.flag_temp[j] = true;//这个控制点的flag标志改为true
          double length = (intersection_point - cps_.points.col(j)).norm();
          //法向量的长度
          if (length > 1e-5)//长度>1e-5
          {
            for (double a = length; a >= 0.0; a -= grid_map_->getResolution())//当前长度一点一点减去比例尺,目的是离散化从j控制点到垂直交点这个向量
            {//开始是在外面A*轨迹上，a越小越靠近j控制点
              occ = grid_map_->getInflateOccupancy((a / length) * intersection_point + (1 - a / length) * cps_.points.col(j));
              //这个法向量上的点是否为障碍物

              if (occ || a < grid_map_->getResolution())//如果是障碍物，或者长度已经到了最小（下一次就循环结束了）
              {
                if (occ)//如果是障碍物
                  a += grid_map_->getResolution();//往外找
                cps_.base_point[j].push_back((a / length) * intersection_point + (1 - a / length) * cps_.points.col(j));//第j个控制点所对应方向向量上的碰撞起始点
                cps_.direction[j].push_back((intersection_point - cps_.points.col(j)).normalized());//第j个控制点所对应的方向向量
                break;
              }
            }
          }
        }
      }

      /* Corner case: the segment length is too short. Here the control points may outside the A* path, leading to opposite gradient direction. So I have to take special care of it */
      if (segment_ids[i].second - segment_ids[i].first == 1)
      //当这段控制点数=2时
      {
        Eigen::Vector3d ctrl_pts_law(cps_.points.col(segment_ids[i].second) - cps_.points.col(segment_ids[i].first)), intersection_point;
        //切向量的新定义
        Eigen::Vector3d middle_point = (cps_.points.col(segment_ids[i].second) + cps_.points.col(segment_ids[i].first)) / 2;
        //这段轨迹中点坐标
        int Astar_id = a_star_pathes[i].size() / 2, last_Astar_id; // Let "Astar_id = id_of_the_most_far_away_Astar_point" will be better, but it needs more computation
        double val = (a_star_pathes[i][Astar_id] - middle_point).dot(ctrl_pts_law), last_val = val;
        while (Astar_id >= 0 && Astar_id < (int)a_star_pathes[i].size())
        {
          last_Astar_id = Astar_id;

          if (val >= 0)
            --Astar_id;
          else
            ++Astar_id;

          val = (a_star_pathes[i][Astar_id] - middle_point).dot(ctrl_pts_law);

          if (val * last_val <= 0 && (abs(val) > 0 || abs(last_val) > 0)) // val = last_val = 0.0 is not allowed
          {
            intersection_point =
                a_star_pathes[i][Astar_id] +
                ((a_star_pathes[i][Astar_id] - a_star_pathes[i][last_Astar_id]) *
                 (ctrl_pts_law.dot(middle_point - a_star_pathes[i][Astar_id]) / ctrl_pts_law.dot(a_star_pathes[i][Astar_id] - a_star_pathes[i][last_Astar_id])) // = t
                );
                //这块的计算同上
            if ((intersection_point - middle_point).norm() > 0.01) //如果法向量的长度> 1cm.
            {
              cps_.flag_temp[segment_ids[i].first] = true;//这个控制点的flag标志改为true
              cps_.base_point[segment_ids[i].first].push_back(cps_.points.col(segment_ids[i].first));//该控制点所对应方向向量上的碰撞起始点
              cps_.direction[segment_ids[i].first].push_back((intersection_point - middle_point).normalized());//方向向量

              got_intersection_id = segment_ids[i].first;//代表找到了这段只包含两个控制点的碰撞轨迹中起点控制点所对应的的法向量和A*规划的轨迹交点的坐标
            }
            break;
          }
        }
      }

      //step 3//证明理解在pad上。
      //如果障碍物控制点对之间有控制点没有碰撞起始点和方向向量，则将最后一次得到的牵引点和方向作为碰撞起始点和方向向量，压入对应的cps_.base_point和cps_.direction。
      if (got_intersection_id >= 0)//找到了这个点的法向量
      {
        for (int j = got_intersection_id + 1; j <= final_segment_ids[i].second; ++j)//从这个控制点的下一个开始循环，直到这段碰撞轨迹的最后一个控制点
          if (!cps_.flag_temp[j])//如果后面的控制点没找到对应法向量
          {
            //就把上一个控制点对应的碰撞点和方向向量放进去
            cps_.base_point[j].push_back(cps_.base_point[j - 1].back());
            cps_.direction[j].push_back(cps_.direction[j - 1].back());
          }

        for (int j = got_intersection_id - 1; j >= final_segment_ids[i].first; --j)//从这个控制点的上一个开始循环，直到这段碰撞轨迹的第一个控制点
          if (!cps_.flag_temp[j])//如果前面的控制点没找到对应法向量
          {
            //就把下一个控制点对应的碰撞点和方向向量放进去
            cps_.base_point[j].push_back(cps_.base_point[j + 1].back());
            cps_.direction[j].push_back(cps_.direction[j + 1].back());
          }
      }
      else
      {
        // Just ignore, it does not matter ^_^.
        // ROS_ERROR("Failed to generate direction! segment_id=%d", i);
      }
    }//只有碰撞段连续的那些点才会存储多个

    return a_star_pathes;//返回A*搜索的路径
  }

  int BsplineOptimizer::earlyExit(void *func_data, const double *x, const double *g, const double fx, const double xnorm, const double gnorm, const double step, int n, int k, int ls)
  {
    BsplineOptimizer *opt = reinterpret_cast<BsplineOptimizer *>(func_data);
    // cout << "k=" << k << endl;
    // cout << "opt->flag_continue_to_optimize_=" << opt->flag_continue_to_optimize_ << endl;
    return (opt->force_stop_type_ == STOP_FOR_ERROR || opt->force_stop_type_ == STOP_FOR_REBOUND);
  //返回是否要提前停止
  }

  double BsplineOptimizer::costFunctionRebound(void *func_data, const double *x, double *grad, const int n)
  {
    BsplineOptimizer *opt = reinterpret_cast<BsplineOptimizer *>(func_data);//使*opt引用func_data
    //http://c.biancheng.net/view/410.html
    double cost;//
    opt->combineCostRebound(x, grad, cost, n);

    opt->iter_num_ += 1;
    return cost;//得到rebound过程的cost
  }

  double BsplineOptimizer::costFunctionRefine(void *func_data, const double *x, double *grad, const int n)
  {
    BsplineOptimizer *opt = reinterpret_cast<BsplineOptimizer *>(func_data);

    double cost;
    opt->combineCostRefine(x, grad, cost, n);

    opt->iter_num_ += 1;
    return cost;//得到refine过程的cost
  }

  void BsplineOptimizer::calcDistanceCostRebound(const Eigen::MatrixXd &q, double &cost,
                                                 Eigen::MatrixXd &gradient, int iter_num, double smoothness_cost)
  {
    cost = 0.0;//初始化cost
    int end_idx = q.cols() - order_;//最后一个要考虑的控制点的角标
    double demarcation = cps_.clearance;//安全距离
    double a = 3 * demarcation, b = -3 * pow(demarcation, 2), c = pow(demarcation, 3);
    //定义距离函数时所需要的系数，参考论文公式5

    force_stop_type_ = DONT_STOP;//强制停止类型为不停止
    if (iter_num > 3 && smoothness_cost / (cps_.size - 2 * order_) < 0.1) // 0.1 is an experimental value that indicates the trajectory is smooth enough.
    //如果迭代次数>3并且光滑项已经满足了需求
    {
      check_collision_and_rebound();//检查是否有新的碰撞物
    }

    /*** calculate distance cost and gradient ***/
    for (auto i = order_; i < end_idx; ++i)//循环控制点
    {
      for (size_t j = 0; j < cps_.direction[i].size(); ++j)//循环其中所存储的方向向量集
      {
        double dist = (cps_.points.col(i) - cps_.base_point[i][j]).dot(cps_.direction[i][j]);//公式5中的dij，注意这里是从base_point指向控制点Q
        double dist_err = cps_.clearance - dist;//公式5中的cij(实际距离-安全距离)
        Eigen::Vector3d dist_grad = cps_.direction[i][j];//即公式7中的vij

        //梯度的定义，参考公式7
        if (dist_err < 0)
        {
          /* do nothing */
        }
        else if (dist_err < demarcation)//<安全距离
        {
          cost += pow(dist_err, 3);
          gradient.col(i) += -3.0 * dist_err * dist_err * dist_grad;
        }
        else
        {
          cost += a * dist_err * dist_err + b * dist_err + c;
          gradient.col(i) += -(2.0 * a * dist_err + b) * dist_grad;
        }
      }
    }
  }

  void BsplineOptimizer::calcFitnessCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient)//计算光滑项
  //给出两条控制点线，第一条是时间分配后的(Φf)，第二个是refine（优化）前的控制点线(Φs)
  {

    cost = 0.0;//初始化cost值

    int end_idx = q.cols() - order_;//最后要考虑的控制点

    // def: f = |x*v|^2/a^2 + |x×v|^2/b^2
    double a2 = 25, b2 = 1;//系数
    for (auto i = order_ - 1; i < end_idx + 1; ++i)//循环控制点
    {
      //Φf:(q.col(i - 1) + 4 * q.col(i) + q.col(i + 1)) 新的轨迹 / 6.0  Φs:ref_pts_[i - 1]旧的轨迹
      Eigen::Vector3d x = (q.col(i - 1) + 4 * q.col(i) + q.col(i + 1)) / 6.0 - ref_pts_[i - 1];
      //Φ̇s:(ref_pts_[i] - ref_pts_[i - 2])  Φ̇s其实就是这个点切向量
      Eigen::Vector3d v = (ref_pts_[i] - ref_pts_[i - 2]).normalized();

      double xdotv = x.dot(v);//点乘
      Eigen::Vector3d xcrossv = x.cross(v);//叉乘

      double f = pow((xdotv), 2) / a2 + pow(xcrossv.norm(), 2) / b2;//cost值 公式18
      cost += f;

      //梯度值，具体为什么这么算不是很清楚！！！！！！
      Eigen::Matrix3d m;
      m << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
      Eigen::Vector3d df_dx = 2 * xdotv / a2 * v + 2 / b2 * m * xcrossv;

      gradient.col(i - 1) += df_dx / 6;
      gradient.col(i) += 4 * df_dx / 6;
      gradient.col(i + 1) += df_dx / 6;
    }
  }

  void BsplineOptimizer::calcSmoothnessCost(const Eigen::MatrixXd &q, double &cost,
                                            Eigen::MatrixXd &gradient, bool falg_use_jerk /* = true*/)//光滑项的计算
  {

    cost = 0.0;//初始化cost的值

    if (falg_use_jerk)//如果使用minimumjerk(3阶)
    {
      Eigen::Vector3d jerk, temp_j;

      for (int i = 0; i < q.cols() - 3; i++)//
      {
        /* evaluate jerk */
        jerk = q.col(i + 3) - 3 * q.col(i + 2) + 3 * q.col(i + 1) - q.col(i);//由Qi点得到Ji点ts是常数，所以其实可以不用管
        cost += jerk.squaredNorm();//cost+这个Ji点的二范数
        temp_j = 2.0 * jerk;//因为是取的二范数模的平方，所以求导后有有一个2
        /* jerk gradient 对控制点求导 */
        // jerk = q.col(i + 3) - 3 * q.col(i + 2) + 3 * q.col(i + 1) - q.col(i)根据这个式子获得导数
        //具体证明见pad
        gradient.col(i + 0) += -temp_j;
        gradient.col(i + 1) += 3.0 * temp_j;
        gradient.col(i + 2) += -3.0 * temp_j;
        gradient.col(i + 3) += temp_j;
      }
    }
    else//如果使用minimumacc(2阶)   和3阶原理一样
    {
      Eigen::Vector3d acc, temp_acc;

      for (int i = 0; i < q.cols() - 2; i++)
      {
        /* evaluate acc */
        acc = q.col(i + 2) - 2 * q.col(i + 1) + q.col(i);
        cost += acc.squaredNorm();
        temp_acc = 2.0 * acc;
        /* acc gradient */
        gradient.col(i + 0) += temp_acc;
        gradient.col(i + 1) += -2.0 * temp_acc;
        gradient.col(i + 2) += temp_acc;
      }
    }
  }

  void BsplineOptimizer::calcFeasibilityCost(const Eigen::MatrixXd &q, double &cost,
                                             Eigen::MatrixXd &gradient)
  {

    #define SECOND_DERIVATIVE_CONTINOUS 
//条件编译，如果宏定义了这个变量，就执行以下编译，否则不动https://blog.csdn.net/wordwarwordwar/article/details/84932183
#ifdef SECOND_DERIVATIVE_CONTINOUS

    cost = 0.0;
    double demarcation = 1.0; // 界线1m/s, 1m/s/s
    double ar = 3 * demarcation, br = -3 * pow(demarcation, 2), cr = pow(demarcation, 3);
    double al = ar, bl = -br, cl = cr;

    /* abbreviation */
    double ts, ts_inv2, ts_inv3;
    ts = bspline_interval_;
    ts_inv2 = 1 / ts / ts;
    ts_inv3 = 1 / ts / ts / ts;

    /* velocity feasibility */
    for (int i = 0; i < q.cols() - 1; i++)
    {
      Eigen::Vector3d vi = (q.col(i + 1) - q.col(i)) / ts;

      for (int j = 0; j < 3; j++)
      {
        if (vi(j) > max_vel_ + demarcation)
        {
          double diff = vi(j) - max_vel_;
          cost += (ar * diff * diff + br * diff + cr) * ts_inv3; // multiply ts_inv3 to make vel and acc has similar magnitude

          double grad = (2.0 * ar * diff + br) / ts * ts_inv3;
          gradient(j, i + 0) += -grad;
          gradient(j, i + 1) += grad;
        }
        else if (vi(j) > max_vel_)
        {
          double diff = vi(j) - max_vel_;
          cost += pow(diff, 3) * ts_inv3;
          ;

          double grad = 3 * diff * diff / ts * ts_inv3;
          ;
          gradient(j, i + 0) += -grad;
          gradient(j, i + 1) += grad;
        }
        else if (vi(j) < -(max_vel_ + demarcation))
        {
          double diff = vi(j) + max_vel_;
          cost += (al * diff * diff + bl * diff + cl) * ts_inv3;

          double grad = (2.0 * al * diff + bl) / ts * ts_inv3;
          gradient(j, i + 0) += -grad;
          gradient(j, i + 1) += grad;
        }
        else if (vi(j) < -max_vel_)
        {
          double diff = vi(j) + max_vel_;
          cost += -pow(diff, 3) * ts_inv3;

          double grad = -3 * diff * diff / ts * ts_inv3;
          gradient(j, i + 0) += -grad;
          gradient(j, i + 1) += grad;
        }
        else
        {
          /* nothing happened */
        }
      }
    }

    /* acceleration feasibility */
    for (int i = 0; i < q.cols() - 2; i++)
    {
      Eigen::Vector3d ai = (q.col(i + 2) - 2 * q.col(i + 1) + q.col(i)) * ts_inv2;

      for (int j = 0; j < 3; j++)
      {
        if (ai(j) > max_acc_ + demarcation)
        {
          double diff = ai(j) - max_acc_;
          cost += ar * diff * diff + br * diff + cr;

          double grad = (2.0 * ar * diff + br) * ts_inv2;
          gradient(j, i + 0) += grad;
          gradient(j, i + 1) += -2 * grad;
          gradient(j, i + 2) += grad;
        }
        else if (ai(j) > max_acc_)
        {
          double diff = ai(j) - max_acc_;
          cost += pow(diff, 3);

          double grad = 3 * diff * diff * ts_inv2;
          gradient(j, i + 0) += grad;
          gradient(j, i + 1) += -2 * grad;
          gradient(j, i + 2) += grad;
        }
        else if (ai(j) < -(max_acc_ + demarcation))
        {
          double diff = ai(j) + max_acc_;
          cost += al * diff * diff + bl * diff + cl;

          double grad = (2.0 * al * diff + bl) * ts_inv2;
          gradient(j, i + 0) += grad;
          gradient(j, i + 1) += -2 * grad;
          gradient(j, i + 2) += grad;
        }
        else if (ai(j) < -max_acc_)
        {
          double diff = ai(j) + max_acc_;
          cost += -pow(diff, 3);

          double grad = -3 * diff * diff * ts_inv2;
          gradient(j, i + 0) += grad;
          gradient(j, i + 1) += -2 * grad;
          gradient(j, i + 2) += grad;
        }
        else
        {
          /* nothing happened */
        }
      }
    }

#else//没有宏定义这个变量SECOND_DERIVATIVE_CONTINOUS

    cost = 0.0;//初始化cost=0
    /* abbreviation */
    double ts, /*vm2, am2, */ ts_inv2;
    // vm2 = max_vel_ * max_vel_;
    // am2 = max_acc_ * max_acc_;

    ts = bspline_interval_;//b样条曲线时间间隔
    ts_inv2 = 1 / ts / ts;//  1/时间间隔的平方

    /* velocity feasibility 速度可行性*/
    for (int i = 0; i < q.cols() - 1; i++)
    {
      Eigen::Vector3d vi = (q.col(i + 1) - q.col(i)) / ts;//速度点Vi

      //cout << "temp_v * vi=" ;
      for (int j = 0; j < 3; j++)//循环三个纬度xyz
      {
        if (vi(j) > max_vel_)//如果超过了最大速度限制
        {
          // cout << "fuck VEL" << endl;
          // cout << vi(j) << endl;
          cost += pow(vi(j) - max_vel_, 2) * ts_inv2; //cost增加，为当前速度和最大速度差的平方 ts_inv2为了使速度和加速度有相同的阶次
          //梯度
          gradient(j, i + 0) += -2 * (vi(j) - max_vel_) / ts * ts_inv2;
          gradient(j, i + 1) += 2 * (vi(j) - max_vel_) / ts * ts_inv2;
          // multiply ts_inv2 to make vel and acc has similar magnitude
          //速度和加速度都是控制点除了t的四次方，就只有和Q相关了。为了防止时间会对梯度产生影响，所以统一次数
        }
        else if (vi(j) < -max_vel_)//同理，小于了最小速度
        {
          cost += pow(vi(j) + max_vel_, 2) * ts_inv2;
          gradient(j, i + 0) += -2 * (vi(j) + max_vel_) / ts * ts_inv2;
          gradient(j, i + 1) += 2 * (vi(j) + max_vel_) / ts * ts_inv2;
        }
        else
        {
          /* code */
        }
      }
    }

    /* acceleration feasibility加速度可行性 */
    for (int i = 0; i < q.cols() - 2; i++)
    {
      Eigen::Vector3d ai = (q.col(i + 2) - 2 * q.col(i + 1) + q.col(i)) * ts_inv2;//加速度点Ai

      //cout << "temp_a * ai=" ;
      for (int j = 0; j < 3; j++)//循环三个纬度xyz
      {
        if (ai(j) > max_acc_)//如果超过了最大加速度限制
        {
          // cout << "fuck ACC" << endl;
          // cout << ai(j) << endl;
          cost += pow(ai(j) - max_acc_, 2);//cost增加，为当前加速度和最大加速度差的平方
          //梯度
          gradient(j, i + 0) += 2 * (ai(j) - max_acc_) * ts_inv2;
          gradient(j, i + 1) += -4 * (ai(j) - max_acc_) * ts_inv2;
          gradient(j, i + 2) += 2 * (ai(j) - max_acc_) * ts_inv2;
        }
        else if (ai(j) < -max_acc_)//同理，小于了最小加速度
        {
          cost += pow(ai(j) + max_acc_, 2);

          gradient(j, i + 0) += 2 * (ai(j) + max_acc_) * ts_inv2;
          gradient(j, i + 1) += -4 * (ai(j) + max_acc_) * ts_inv2;
          gradient(j, i + 2) += 2 * (ai(j) + max_acc_) * ts_inv2;
        }
        else
        {
          /* code */
        }
      }
      //cout << endl;
    }

#endif
  }

  bool BsplineOptimizer::check_collision_and_rebound(void)//
  {

    int end_idx = cps_.size - order_;//最后一个考虑的控制点的角标

    /*** Check and segment the initial trajectory according to obstacles ***/
    int in_id, out_id;//起始id和结束id
    vector<std::pair<int, int>> segment_ids;//存储每一段起始和结束控制点的index
    bool flag_new_obs_valid = false;//出现新的障碍物标志为false
    int i_end = end_idx - (end_idx - order_) / 3;//只考虑前三分之二的控制点
    for (int i = order_ - 1; i <= i_end; ++i)//循环控制点
    {

      bool occ = grid_map_->getInflateOccupancy(cps_.points.col(i));//控制点是否在障碍物内

      /*** check if the new collision will be valid   检查新的障碍物是否是有效的 ***/
      if (occ)//如果控制点在障碍物内
      {
        for (size_t k = 0; k < cps_.direction[i].size(); ++k)//循环这个控制点的deriction[i]所存放的所有方向向量
        {
          cout.precision(2);//https://blog.csdn.net/huangchijun11/article/details/72934222
          //输出的时候设定输出值以新的浮点数精度值显示，即小数点后保留2位

          if ((cps_.points.col(i) - cps_.base_point[i][k]).dot(cps_.direction[i][k]) < 1 * grid_map_->getResolution()) // 当前点位于所有碰撞点之外   current point is outside all the collision_points.
          {
            occ = false; // Not really takes effect, just for better hunman understanding.
            break;
          }
        }
      }//什么意思？？??????

      if (occ)//如果是有效的碰撞
      {
        flag_new_obs_valid = true;//新的障碍物是有效的标志

        int j;
        for (j = i - 1; j >= 0; --j)//循环这个有效的碰撞控制点前的控制点
        {
          occ = grid_map_->getInflateOccupancy(cps_.points.col(j));//检查其是不是在障碍内
          if (!occ)//如果不是障碍
          {
            in_id = j;//起始id
            break;
          }
        }
        if (j < 0) //没有找到 fail to get the obs free point
        {
          ROS_ERROR("ERROR! the drone is in obstacle. This should not happen.");
          in_id = 0;
        }

        for (j = i + 1; j < cps_.size; ++j)//循环这个有效的碰撞控制点后所有控制点
        {
          occ = grid_map_->getInflateOccupancy(cps_.points.col(j));//检查其是不是在障碍内

          if (!occ)//如果不是障碍
          {
            out_id = j;//终止id
            break;
          }
        }
        if (j >= cps_.size) // 没有找到 fail to get the obs free point
        {
          ROS_WARN("WARN! terminal point of the current trajectory is in obstacle, skip this planning.");

          force_stop_type_ = STOP_FOR_ERROR;//因为发生错误而停止
          return false;
        }

        i = j + 1;

        segment_ids.push_back(std::pair<int, int>(in_id, out_id));//碰撞段的起始和终止控制点id
      }
    }

    if (flag_new_obs_valid)//如果新的障碍物是有效的
    {
      vector<vector<Eigen::Vector3d>> a_star_pathes;//A*路径集
      for (size_t i = 0; i < segment_ids.size(); ++i)//循环碰撞段
      {
        /*** a star search ***/
        Eigen::Vector3d in(cps_.points.col(segment_ids[i].first)), out(cps_.points.col(segment_ids[i].second));
        if (a_star_->AstarSearch(/*(in-out).norm()/10+0.05*/ 0.1, in, out))//A*寻找到了路径
        {
          a_star_pathes.push_back(a_star_->getPath());
        }//把路径加进去
        else
        {
          ROS_ERROR("a star error");
          segment_ids.erase(segment_ids.begin() + i);//删除这个碰撞段
          i--;
        }
      }

      /*** Assign parameters to each segment ***/
      for (size_t i = 0; i < segment_ids.size(); ++i)//循环碰撞段
      {
        // step 1
        for (int j = segment_ids[i].first; j <= segment_ids[i].second; ++j)//循环这一段上的所有控制点
          cps_.flag_temp[j] = false;//所有控制点的标志为false

        // step 2
        int got_intersection_id = -1;
        for (int j = segment_ids[i].first + 1; j < segment_ids[i].second; ++j)//循环这一段上的所有控制点
        {
          Eigen::Vector3d ctrl_pts_law(cps_.points.col(j + 1) - cps_.points.col(j - 1)), intersection_point;
          int Astar_id = a_star_pathes[i].size() / 2, last_Astar_id; // Let "Astar_id = id_of_the_most_far_away_Astar_point" will be better, but it needs more computation
          double val = (a_star_pathes[i][Astar_id] - cps_.points.col(j)).dot(ctrl_pts_law), last_val = val;
          while (Astar_id >= 0 && Astar_id < (int)a_star_pathes[i].size())
          {
            last_Astar_id = Astar_id;

            if (val >= 0)
              --Astar_id;
            else
              ++Astar_id;

            val = (a_star_pathes[i][Astar_id] - cps_.points.col(j)).dot(ctrl_pts_law);

            // cout << val << endl;

            if (val * last_val <= 0 && (abs(val) > 0 || abs(last_val) > 0)) // val = last_val = 0.0 is not allowed
            {
              intersection_point =
                  a_star_pathes[i][Astar_id] +
                  ((a_star_pathes[i][Astar_id] - a_star_pathes[i][last_Astar_id]) *
                   (ctrl_pts_law.dot(cps_.points.col(j) - a_star_pathes[i][Astar_id]) / ctrl_pts_law.dot(a_star_pathes[i][Astar_id] - a_star_pathes[i][last_Astar_id])) // = t
                  );

              got_intersection_id = j;
              break;
            }
          }//得到第j个控制点的法向量和A*规划的轨迹交点的坐标

          if (got_intersection_id >= 0)//如果找到了法向量和A*规划的轨迹交点
          {
            cps_.flag_temp[j] = true;//标记为true
            double length = (intersection_point - cps_.points.col(j)).norm();
            if (length > 1e-5)
            {
              for (double a = length; a >= 0.0; a -= grid_map_->getResolution())
              {
                bool occ = grid_map_->getInflateOccupancy((a / length) * intersection_point + (1 - a / length) * cps_.points.col(j));

                if (occ || a < grid_map_->getResolution())
                {
                  if (occ)
                    a += grid_map_->getResolution();
                  cps_.base_point[j].push_back((a / length) * intersection_point + (1 - a / length) * cps_.points.col(j));//第j个控制点所对应方向向量上的碰撞起始点
                  cps_.direction[j].push_back((intersection_point - cps_.points.col(j)).normalized());//第j个控制点所对应的方向向量
                  break;
                }
              }
            }
            else
            {
              got_intersection_id = -1;//长度太短（A*算法路径和障碍物太近了）
            }
          }
        }

        //step 3
        if (got_intersection_id >= 0)//这一段中最后一个找到法向量和A*规划的轨迹交点的控制点
        {
          for (int j = got_intersection_id + 1; j <= segment_ids[i].second; ++j)//遍历此控制点后的点
            if (!cps_.flag_temp[j])
            {
              cps_.base_point[j].push_back(cps_.base_point[j - 1].back());
              cps_.direction[j].push_back(cps_.direction[j - 1].back());
            }//如果没找到，就用它前一个存入其中当做其方向向量

          for (int j = got_intersection_id - 1; j >= segment_ids[i].first; --j)
            if (!cps_.flag_temp[j])
            {
              cps_.base_point[j].push_back(cps_.base_point[j + 1].back());
              cps_.direction[j].push_back(cps_.direction[j + 1].back());
            }
        }//每个点存储方向向量集
        else
          ROS_WARN("Failed to generate direction. It doesn't matter.");
      }

      force_stop_type_ = STOP_FOR_REBOUND;
      return true;//如果新的障碍物有效，返回true
    }

    return false;//新的障碍物无效，返回false
  }

  bool BsplineOptimizer::BsplineOptimizeTrajRebound(Eigen::MatrixXd &optimal_points, double ts)
  {
    setBsplineInterval(ts);//设置b样条曲线的时间间隔

    bool flag_success = rebound_optimize();

    optimal_points = cps_.points;//曲线的控制点

    return flag_success;//返回是否规划成功
  }

  bool BsplineOptimizer::BsplineOptimizeTrajRefine(const Eigen::MatrixXd &init_points, const double ts, Eigen::MatrixXd &optimal_points)
  {

    setControlPoints(init_points);//初始化控制点 
    setBsplineInterval(ts);//得到时间间隔

    bool flag_success = refine_optimize();//得到当前的轨迹是否安全

    optimal_points = cps_.points;

    return flag_success;
  }

  bool BsplineOptimizer::rebound_optimize()
  {
    iter_num_ = 0;//迭代次数
    int start_id = order_;//开始考虑的控制点index（从order_开始考虑）
    int end_id = this->cps_.size - order_;//最后考虑的控制点（总个数-order_)
    //说白了就是不考虑开头和结尾的order_个点，只考虑中间

    variable_num_ = 3 * (end_id - start_id);//变量总数，考虑的控制点数*3（三维）
    double final_cost;//最终的cost值

    ros::Time t0 = ros::Time::now(), t1, t2;
    int restart_nums = 0, rebound_times = 0;
    ;
    bool flag_force_return, flag_occ, success;
    new_lambda2_ = lambda2_;
    constexpr int MAX_RESART_NUMS_SET = 3;
    do
    {
      /* ---------- prepare ---------- */
      min_cost_ = std::numeric_limits<double>::max();//int类型所能表达的最大值。其实就是初始化cost，给一个巨大的数值
     // https://blog.csdn.net/fengbingchun/article/details/77922558/
      
      iter_num_ = 0;//迭代次数
      flag_force_return = false;
      flag_occ = false;
      success = false;

      double q[variable_num_];//建立一个和所需要优化变量数目相同的数组q，用来存储优化变量的值
      memcpy(q, cps_.points.data() + 3 * start_id, variable_num_ * sizeof(q[0]));//函数使用参考https://blog.csdn.net/qq_26747797/article/details/82794333
      //其实就是把需要考虑的控制点的3纬坐标（即要优化的变量）拷贝到q向量中作为优化变量的初始值

      lbfgs::lbfgs_parameter_t lbfgs_params;//lbfgs梯度优化的参数
      lbfgs::lbfgs_load_default_parameters(&lbfgs_params);//初始化梯度优化参数
      lbfgs_params.mem_size = 16;
      lbfgs_params.max_iterations = 200;
      lbfgs_params.g_epsilon = 0.01;//参数的一些改动

      /* ---------- optimize ----------优化过程 */
      t1 = ros::Time::now();//开始时间
      int result = lbfgs::lbfgs_optimize(variable_num_, q, &final_cost, BsplineOptimizer::costFunctionRebound, NULL, BsplineOptimizer::earlyExit, this, &lbfgs_params);
      t2 = ros::Time::now();//结束时间
      double time_ms = (t2 - t1).toSec() * 1000;//优化时间
      double total_time_ms = (t2 - t0).toSec() * 1000;//整体时间

      /* ---------- success temporary, check collision again 暂时成功，再次检查碰撞  ---------- */
      if (result == lbfgs::LBFGS_CONVERGENCE ||
          result == lbfgs::LBFGSERR_MAXIMUMITERATION ||
          result == lbfgs::LBFGS_ALREADY_MINIMIZED ||
          result == lbfgs::LBFGS_STOP)//优化成功
      {
        //ROS_WARN("Solver error in planning!, return = %s", lbfgs::lbfgs_strerror(result));
        flag_force_return = false;

        UniformBspline traj = UniformBspline(cps_.points, 3, bspline_interval_);//得到均匀b样条曲线的knot span
        double tm, tmp;
        traj.getTimeSpan(tm, tmp);//得到开始和结束的时间
        //tm和tmp取决于u_向量，而定义的时候已经将时间定义在了[0, duration]上，即有效时间段
        double t_step = (tmp - tm) / ((traj.evaluateDeBoorT(tmp) - traj.evaluateDeBoorT(tm)).norm() / grid_map_->getResolution());
        //计算时间step   后头一项代表开始时和结束时的距离 以此为基准离散化时间 
        //把遍历轨迹变为遍历时间

        for (double t = tm; t < tmp * 2 / 3; t += t_step) //只检查前2/3  Only check the closest 2/3 partition of the whole trajectory.
        {
          flag_occ = grid_map_->getInflateOccupancy(traj.evaluateDeBoorT(t));//t时间轨迹是否碰到障碍物
          if (flag_occ)//如果碰到了
          {
            //cout << "hit_obs, t=" << t << " P=" << traj.evaluateDeBoorT(t).transpose() << endl;

            if (t <= bspline_interval_) // First 3 control points in obstacles!
            //如果还在第一个有效时间段内，即[up_,up_+1]发生了碰撞。说明前四个控制点会碰撞（前四个控制点的基础函数决定着这个时间段内的点）
            {
              cout << cps_.points.col(1).transpose() << "\n"
                   << cps_.points.col(2).transpose() << "\n"
                   << cps_.points.col(3).transpose() << "\n"
                   << cps_.points.col(4).transpose() << endl;//
              ROS_WARN("First 3 control points in obstacles! return false, t=%f", t);
              return false;
            }

            break;//如果碰到了障碍就跳出循环
          }
        }

        if (!flag_occ)//全程没有碰撞
        {
          printf("\033[32miter(+1)=%d,time(ms)=%5.3f,total_t(ms)=%5.3f,cost=%5.3f\n\033[0m", iter_num_, time_ms, total_time_ms, final_cost);
          success = true;//成功规划成功规划标志
        }
        else // restart 发生了碰撞，重新开始规划
        {
          restart_nums++;//重新规划次数+1
          initControlPoints(cps_.points, false);
          //根据现在的控制点重新去规划，找到碰撞段得到无碰撞A*路径
          new_lambda2_ *= 2;//增大碰撞项系数

          printf("\033[32miter(+1)=%d,time(ms)=%5.3f,keep optimizing\n\033[0m", iter_num_, time_ms);
        }
      }
      else if (result == lbfgs::LBFGSERR_CANCELED)//最小化过程已被取消
      {
        flag_force_return = true;//强制返回标志为true
        rebound_times++;//弹回次数+1
        cout << "iter=" << iter_num_ << ",time(ms)=" << time_ms << ",rebound." << endl;
      }
      else
      {
        ROS_WARN("Solver error. Return = %d, %s. Skip this planning.", result, lbfgs::lbfgs_strerror(result));
        // while (ros::ok());
      }

    } while ((flag_occ && restart_nums < MAX_RESART_NUMS_SET) ||
             (flag_force_return && force_stop_type_ == STOP_FOR_REBOUND && rebound_times <= 20));
             //当(发生碰撞并且重规划次数<最大重规划次数)||(强制返回并且停止类型固定并且弹回次数<=20) 循环规划
    return success;
  }

  bool BsplineOptimizer::refine_optimize()
  {
    iter_num_ = 0;//迭代次数
    int start_id = order_;//开始考虑的控制点index（从order_开始考虑）
    int end_id = this->cps_.points.cols() - order_;//最后考虑的控制点（总个数-order_)
    variable_num_ = 3 * (end_id - start_id);//变量总数

    double q[variable_num_];///建立一个和所需要优化变量数目相同的数组q，用来存储优化变量的值
    double final_cost;//cost数值

    memcpy(q, cps_.points.data() + 3 * start_id, variable_num_ * sizeof(q[0]));//初始变量值拷贝

    double origin_lambda4 = lambda4_;//拟合项的初始系数lambda4
    bool flag_safe = true;//安全标志为true
    int iter_count = 0;
    do
    {
      lbfgs::lbfgs_parameter_t lbfgs_params;
      lbfgs::lbfgs_load_default_parameters(&lbfgs_params);//初始化lbfgs参数
      lbfgs_params.mem_size = 16;
      lbfgs_params.max_iterations = 200;
      lbfgs_params.g_epsilon = 0.001;//参数的一些修改

      int result = lbfgs::lbfgs_optimize(variable_num_, q, &final_cost, BsplineOptimizer::costFunctionRefine, NULL, NULL, this, &lbfgs_params);//优化后的结果，不同数字代表优化的不同情况
      if (result == lbfgs::LBFGS_CONVERGENCE ||
          result == lbfgs::LBFGSERR_MAXIMUMITERATION ||
          result == lbfgs::LBFGS_ALREADY_MINIMIZED ||
          result == lbfgs::LBFGS_STOP)//如果优化成功，跳过
      {
        //pass
      }
      else
      {
        ROS_ERROR("Solver error in refining!, return = %d, %s", result, lbfgs::lbfgs_strerror(result));
      }//报错

      UniformBspline traj = UniformBspline(cps_.points, 3, bspline_interval_);//得到均匀b样条曲线的knot span
      double tm, tmp;//开始终止时间
      traj.getTimeSpan(tm, tmp);//时间间隔
      double t_step = (tmp - tm) / ((traj.evaluateDeBoorT(tmp) - traj.evaluateDeBoorT(tm)).norm() / grid_map_->getResolution()); // Step size is defined as the maximum size that can passes throgth every gird.
      for (double t = tm; t < tmp * 2 / 3; t += t_step)//遍历前2/3轨迹
      {
        if (grid_map_->getInflateOccupancy(traj.evaluateDeBoorT(t)))//如果碰到了障碍
        {
          // cout << "Refined traj hit_obs, t=" << t << " P=" << traj.evaluateDeBoorT(t).transpose() << endl;

          Eigen::MatrixXd ref_pts(ref_pts_.size(), 3);//新建立一个矩阵，存储所有ref_pts_，每一行是一个点，三列分别为三个坐标
          for (size_t i = 0; i < ref_pts_.size(); i++)
          {
            ref_pts.row(i) = ref_pts_[i].transpose();//每行放一个点
          }

          flag_safe = false;//是不安全的
          break;//跳出前2/3轨迹这一循环
        }
      }

      if (!flag_safe)//如果轨迹不安全
        lambda4_ *= 2;//重新改变lamada4的数值

      iter_count++;//循环次数+1
    } while (!flag_safe && iter_count <= 0);//进行一次就跳出（因为iter_count>0了）

    lambda4_ = origin_lambda4;//重新再赋予lamada4值

    //cout << "iter_num_=" << iter_num_ << endl;

    return flag_safe;//返回是否安全
  }

  void BsplineOptimizer::combineCostRebound(const double *x, double *grad, double &f_combine, const int n)
  {

    memcpy(cps_.points.data() + 3 * order_, x, n * sizeof(x[0]));
    //将指针x所指示的位置存储的数据，拷贝到cps_.points中

    /* ---------- evaluate cost and gradient ---------- */
    double f_smoothness, f_distance, f_feasibility;

    Eigen::MatrixXd g_smoothness = Eigen::MatrixXd::Zero(3, cps_.size);
    Eigen::MatrixXd g_distance = Eigen::MatrixXd::Zero(3, cps_.size);
    Eigen::MatrixXd g_feasibility = Eigen::MatrixXd::Zero(3, cps_.size);

    calcSmoothnessCost(cps_.points, f_smoothness, g_smoothness);//计算光滑项和梯度 
    calcDistanceCostRebound(cps_.points, f_distance, g_distance, iter_num_, f_smoothness);//计算碰撞项和梯度
    calcFeasibilityCost(cps_.points, f_feasibility, g_feasibility);//计算可行性项和梯度

    f_combine = lambda1_ * f_smoothness + new_lambda2_ * f_distance + lambda3_ * f_feasibility;//整体的cost
    //printf("origin %f %f %f %f\n", f_smoothness, f_distance, f_feasibility, f_combine);

    Eigen::MatrixXd grad_3D = lambda1_ * g_smoothness + new_lambda2_ * g_distance + lambda3_ * g_feasibility;//整体的梯度矩阵，注意碰撞项的系数变成了new_lamada2
    memcpy(grad, grad_3D.data() + 3 * order_, n * sizeof(grad[0]));
  }

  void BsplineOptimizer::combineCostRefine(const double *x, double *grad, double &f_combine, const int n)
  {

    memcpy(cps_.points.data() + 3 * order_, x, n * sizeof(x[0]));

    /* ---------- evaluate cost and gradient ---------- */
    double f_smoothness, f_fitness, f_feasibility;

    Eigen::MatrixXd g_smoothness = Eigen::MatrixXd::Zero(3, cps_.points.cols());
    Eigen::MatrixXd g_fitness = Eigen::MatrixXd::Zero(3, cps_.points.cols());
    Eigen::MatrixXd g_feasibility = Eigen::MatrixXd::Zero(3, cps_.points.cols());

    //time_satrt = ros::Time::now();

    calcSmoothnessCost(cps_.points, f_smoothness, g_smoothness);//计算光滑项和梯度
    calcFitnessCost(cps_.points, f_fitness, g_fitness);// 计算拟合项和梯度 
    calcFeasibilityCost(cps_.points, f_feasibility, g_feasibility);//计算可行性项和梯度

    /* ---------- convert to solver format...---------- */
    f_combine = lambda1_ * f_smoothness + lambda4_ * f_fitness + lambda3_ * f_feasibility;//整体的cost
    // printf("origin %f %f %f %f\n", f_smoothness, f_fitness, f_feasibility, f_combine);

    Eigen::MatrixXd grad_3D = lambda1_ * g_smoothness + lambda4_ * g_fitness + lambda3_ * g_feasibility;//整体的梯度矩阵
    memcpy(grad, grad_3D.data() + 3 * order_, n * sizeof(grad[0]));
  }

} // namespace ego_planner