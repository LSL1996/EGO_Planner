#include <ros/ros.h>
#include <visualization_msgs/Marker.h>

#include <plan_manage/ego_replan_fsm.h>

using namespace ego_planner;
//规划部分的主函数
int main(int argc, char **argv)
{
  ros::init(argc, argv, "ego_planner_node");//ros节点初始化
  ros::NodeHandle nh("~");//创建节点句柄

  EGOReplanFSM rebo_replan;//定义一个fsm类

  rebo_replan.init(nh);//初始化节点

  ros::Duration(1.0).sleep();//关于time和duration的讲解https://blog.csdn.net/datase/article/details/80156409
  ros::spin();
  //关于回调函数和spin的知识参考https://blog.csdn.net/qq_33898609/article/details/105935613

  return 0;
}
