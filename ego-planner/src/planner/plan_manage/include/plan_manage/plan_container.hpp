#ifndef _PLAN_CONTAINER_H_
#define _PLAN_CONTAINER_H_

#include <Eigen/Eigen>
#include <vector>
#include <ros/ros.h>

#include <bspline_opt/uniform_bspline.h>
#include <traj_utils/polynomial_traj.h>

using std::vector;

namespace ego_planner
{

  class GlobalTrajData//全局轨迹
  {
  private:
  public:
    PolynomialTraj global_traj_;//全局:多项式（minimumsnap）轨迹
    vector<UniformBspline> local_traj_;//局部:b样条曲线轨迹，这个其实没用到

    double global_duration_;//全局时间
    ros::Time global_start_time_;//全局规划的开始时间
    double local_start_time_, local_end_time_;//局部规划的开始时间和结束时间
    double time_increase_;
    double last_time_inc_;
    double last_progress_time_;

    GlobalTrajData(/* args */) {}

    ~GlobalTrajData() {}

    bool localTrajReachTarget() { return fabs(local_end_time_ - global_duration_) < 0.1; }
    //看局部规划的终点时间和全局规划的时间差是否<0.1

    void setGlobalTraj(const PolynomialTraj &traj, const ros::Time &time)
    //设置全局轨迹
    {
      global_traj_ = traj;//全局轨迹
      global_traj_.init();//初始化全局轨迹段数和总时间
      global_duration_ = global_traj_.getTimeSum();//全局规划的总时间
      global_start_time_ = time;//全局规划开始时间
      //这块其实没啥用，没有用到
      local_traj_.clear();//局部轨迹清空
      local_start_time_ = -1;
      local_end_time_ = -1;//局部规划开始结束时间
      time_increase_ = 0.0;
      last_time_inc_ = 0.0;
      last_progress_time_ = 0.0;
    }

    void setLocalTraj(UniformBspline traj, double local_ts, double local_te, double time_inc)
   //设置局部轨迹
    {
      local_traj_.resize(3);//大小设置为3，分别存放位置速度加速度
      local_traj_[0] = traj;//位置轨迹 
      local_traj_[1] = local_traj_[0].getDerivative();//速度轨迹
      local_traj_[2] = local_traj_[1].getDerivative();//加速度轨迹

      local_start_time_ = local_ts;
      local_end_time_ = local_te;//局部开始结束时间
      global_duration_ += time_inc;
      time_increase_ += time_inc;
      last_time_inc_ = time_inc;
    }

    Eigen::Vector3d getPosition(double t)
    //获得位置 
    //这几个变量什么意思？？？？？？？？？time_inc time_increase_
    {
      //时间在全局规划的范围内，获得全局轨迹的位置
      if (t >= -1e-3 && t <= local_start_time_)
      {
        return global_traj_.evaluate(t - time_increase_ + last_time_inc_);//
      }
      else if (t >= local_end_time_ && t <= global_duration_ + 1e-3)
      {
        return global_traj_.evaluate(t - time_increase_);
      }
      //时间在局部规划的范围内，获得局部轨迹的位置 
      else
      {
        double tm, tmp;
        local_traj_[0].getTimeSpan(tm, tmp);//局部轨迹的起始和终止时间
        return local_traj_[0].evaluateDeBoor(tm + t - local_start_time_);
      }
    }
    //速度加速度同理
    Eigen::Vector3d getVelocity(double t)
    {
      if (t >= -1e-3 && t <= local_start_time_)
      {
        return global_traj_.evaluateVel(t);
      }
      else if (t >= local_end_time_ && t <= global_duration_ + 1e-3)
      {
        return global_traj_.evaluateVel(t - time_increase_);
      }
      else
      {
        double tm, tmp;
        local_traj_[0].getTimeSpan(tm, tmp);
        return local_traj_[1].evaluateDeBoor(tm + t - local_start_time_);
      }
    }

    Eigen::Vector3d getAcceleration(double t)
    {
      if (t >= -1e-3 && t <= local_start_time_)
      {
        return global_traj_.evaluateAcc(t);
      }
      else if (t >= local_end_time_ && t <= global_duration_ + 1e-3)
      {
        return global_traj_.evaluateAcc(t - time_increase_);
      }
      else
      {
        double tm, tmp;
        local_traj_[0].getTimeSpan(tm, tmp);
        return local_traj_[2].evaluateDeBoor(tm + t - local_start_time_);
      }
    }

    // 获取球体内局部轨迹的Bspline参数化数据   get Bspline paramterization data of a local trajectory within a sphere
    // start_t: 轨迹的开始时间 start time of the trajectory
    // dist_pt: 离散点之间的距离distance between the discretized points
    void getTrajByRadius(const double &start_t, const double &des_radius, const double &dist_pt,
                         vector<Eigen::Vector3d> &point_set, vector<Eigen::Vector3d> &start_end_derivative,
                         double &dt, double &seg_duration)
    //在一定的半径范围内获取局部轨迹
    {
      double seg_length = 0.0; //  截断段的长度  length of the truncated segment
      double seg_time = 0.0;   //  截断段的时间  duration of the truncated segment
      double radius = 0.0;     // 到这段局部轨迹第一个点的距离distance to the first point of the segment

      double delt = 0.2;
      Eigen::Vector3d first_pt = getPosition(start_t); // 起始点的位置 first point of the segment
      Eigen::Vector3d prev_pt = first_pt;              // 上一个点 previous point
      Eigen::Vector3d cur_pt;                          // 现在这个点 current point

      //  继续前进，直到轨迹超过半径或全局时间  go forward until the traj exceed radius or global time

      while (radius < des_radius && seg_time < global_duration_ - start_t - 1e-3)
      //当前点到起始点的距离<设置的半径 并且时间在全局规划的时间内
      {
        seg_time += delt;//时间前进delta
        seg_time = min(seg_time, global_duration_ - start_t);//局部规划的时间（注意不能超出全局时间所以取了min）

        cur_pt = getPosition(start_t + seg_time);//现在的位置
        seg_length += (cur_pt - prev_pt).norm();//长度
        prev_pt = cur_pt;//更新点
        radius = (cur_pt - first_pt).norm();//当前点到起始点的距离
      }

      //  通过所需的点密度获得参数dt  get parameterization dt by desired density of points
      int seg_num = floor(seg_length / dist_pt);//这一段的点数目：总长度/离散点的距离

      // get outputs

      seg_duration = seg_time; // 局部轨迹的时间 duration of the truncated segment
      dt = seg_time / seg_num; //两点之间的时间长度dt=局部轨迹的时间/点个数 time difference between two points

      for (double tp = 0.0; tp <= seg_time + 1e-4; tp += dt)
      {
        cur_pt = getPosition(start_t + tp);
        point_set.push_back(cur_pt);
      }//细化离散时间得到轨迹上的点

      start_end_derivative.push_back(getVelocity(start_t));
      start_end_derivative.push_back(getVelocity(start_t + seg_time));
      start_end_derivative.push_back(getAcceleration(start_t));
      start_end_derivative.push_back(getAcceleration(start_t + seg_time));//得到开始结束点的速度加速度
    }

    //  获取固定持续时间局部轨迹的Bspline参数化数据  get Bspline paramterization data of a fixed duration local trajectory
    // start_t:局部轨迹开始时间  start time of the trajectory 
    // duration: 局部轨迹持续时间  time length of the segment
    // seg_num: 把本段轨迹离散化成seg_num段  discretized the segment into *seg_num* parts
    void getTrajByDuration(double start_t, double duration, int seg_num,
                           vector<Eigen::Vector3d> &point_set,
                           vector<Eigen::Vector3d> &start_end_derivative, double &dt)
    {
      dt = duration / seg_num;//时间/段数
      Eigen::Vector3d cur_pt;
      for (double tp = 0.0; tp <= duration + 1e-4; tp += dt)
      {
        cur_pt = getPosition(start_t + tp);
        point_set.push_back(cur_pt);
      }////离散时间得到轨迹上的点

      start_end_derivative.push_back(getVelocity(start_t));
      start_end_derivative.push_back(getVelocity(start_t + duration));
      start_end_derivative.push_back(getAcceleration(start_t));
      start_end_derivative.push_back(getAcceleration(start_t + duration));
    }//得到开始结束点的速度加速度
  };

  struct PlanParameters//规划参数
  {
    /* planning algorithm parameters */
    double max_vel_, max_acc_, max_jerk_; // physical limits
    double ctrl_pt_dist;                  //  相邻B样条控制点之间的距离  distance between adjacient B-spline control points 0.4
    double feasibility_tolerance_;        // permitted ratio of vel/acc exceeding limits
    double planning_horizen_;//7.5
    /* processing time */
    double time_search_ = 0.0;
    double time_optimize_ = 0.0;
    double time_adjust_ = 0.0;
  };

  struct LocalTrajData//局部轨迹
  {
    /* info of generated traj */

    int traj_id_;
    double duration_;
    double global_time_offset; // This is because when the local traj finished and is going to switch back to the global traj, the global traj time is no longer matches the world time.
    //这是因为当本地轨迹完成并将切换回全局轨迹时，全局轨迹时间不再与世界时间匹配。 
    ros::Time start_time_;
    Eigen::Vector3d start_pos_;
    UniformBspline position_traj_, velocity_traj_, acceleration_traj_;
  };

} // namespace ego_planner

#endif