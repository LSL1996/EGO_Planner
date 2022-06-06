#include <ros/ros.h>
#include <visualization_msgs/Marker.h>

#include <plan_manage/ego_replan_fsm.h>

using namespace ego_planner;
//�滮���ֵ�������
int main(int argc, char **argv)
{
  ros::init(argc, argv, "ego_planner_node");//ros�ڵ��ʼ��
  ros::NodeHandle nh("~");//�����ڵ���

  EGOReplanFSM rebo_replan;//����һ��fsm��

  rebo_replan.init(nh);//��ʼ���ڵ�

  ros::Duration(1.0).sleep();//����time��duration�Ľ���https://blog.csdn.net/datase/article/details/80156409
  ros::spin();
  //���ڻص�������spin��֪ʶ�ο�https://blog.csdn.net/qq_33898609/article/details/105935613

  return 0;
}
