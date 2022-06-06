#ifndef _GRID_MAP_H
#define _GRID_MAP_H

#include <Eigen/Eigen>
#include <Eigen/StdVector>
#include <cv_bridge/cv_bridge.h>
#include <geometry_msgs/PoseStamped.h>
#include <iostream>
#include <random>
#include <nav_msgs/Odometry.h>
#include <queue>
#include <ros/ros.h>
#include <tuple>
#include <visualization_msgs/Marker.h>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>

#include <message_filters/subscriber.h>
#include <message_filters/sync_policies/approximate_time.h>
#include <message_filters/sync_policies/exact_time.h>
#include <message_filters/time_synchronizer.h>

#include <plan_env/raycast.h>

#define logit(x) (log((x) / (1 - (x))))//logit 函数 自变量从(0,1)连续单调地映射到正负无穷

using namespace std;

// voxel hashing
template <typename T>
struct matrix_hash : std::unary_function<T, size_t> {
  std::size_t operator()(T const& matrix) const {
    size_t seed = 0;
    for (size_t i = 0; i < matrix.size(); ++i) {
      auto elem = *(matrix.data() + i);
      seed ^= std::hash<typename T::Scalar>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

// constant parameters

struct MappingParameters {//地图参数结构体

  /* map properties */
  Eigen::Vector3d map_origin_, map_size_;
  Eigen::Vector3d map_min_boundary_, map_max_boundary_;  // map range in pos 实际位置的最大范围
  Eigen::Vector3i map_voxel_num_; // map range in index 角标grid地图中的最大范围
  Eigen::Vector3d local_update_range_;
  double resolution_, resolution_inv_;
  double obstacles_inflation_;
  string frame_id_;
  int pose_type_;

  /* camera parameters 摄像机参数  */
  double cx_, cy_, fx_, fy_;
  //fx fy为两个方向的焦距，至于为什么会有两个不同的，参考https://blog.csdn.net/xujie126/article/details/118026380?spm=1001.2101.3001.6650.1&utm_medium=distribute.pc_relevant.none-task-blog-2%7Edefault%7ECTRLIST%7ERate-1.pc_relevant_default&depth_1-utm_source=distribute.pc_relevant.none-task-blog-2%7Edefault%7ECTRLIST%7ERate-1.pc_relevant_default&utm_relevant_index=

  /* depth image projection filtering 深度图像投影滤波  */
  double depth_filter_maxdist_, depth_filter_mindist_, depth_filter_tolerance_;
  int depth_filter_margin_;
  bool use_depth_filter_;
  double k_depth_scaling_factor_;
  int skip_pixel_;

  /* raycasting */
  double p_hit_, p_miss_, p_min_, p_max_, p_occ_;  // 占据概率 occupancy probability
  double prob_hit_log_, prob_miss_log_, clamp_min_log_, clamp_max_log_,
      min_occupancy_log_;                   // 概率取logic函数 logit of occupancy probability
  double min_ray_length_, max_ray_length_;  //激光探测的最短和最长距离  range of doing raycasting

  /* local map update and clear */
  int local_map_margin_;//局部地图的大小

  /* visualization and computation time display 可视化和计算时间显示  */
  double visualization_truncate_height_, virtual_ceil_height_, ground_height_;
  bool show_occ_time_;

  /* active mapping */
  double unknown_flag_;
};

// intermediate mapping data for fusion

struct MappingData {
  // main map data, occupancy of each voxel and Euclidean distance 
  //主地图数据,每个体素是否被占据（是障碍物）和欧几里德距离 

  std::vector<double> occupancy_buffer_;//一维向量（把整个地图的三纬信息存储为一维，同A*data）是否被占据（是障碍物）
  std::vector<char> occupancy_buffer_inflate_;//是否膨胀（根据无人机的大小对障碍物进行膨胀，因为最终要把无人机当成质点）

  // camera position and pose data相机的位置和姿态信息

  Eigen::Vector3d camera_pos_, last_camera_pos_;//相机这次的位置和上次的位置
  Eigen::Quaterniond camera_q_, last_camera_q_;//相机坐标系到世界坐标系的四元数 当前相机的和上一次的

  // depth image data 深度图像信息

  cv::Mat depth_image_, last_depth_image_;//上一次和这一次的深度图像
  int image_cnt_;

  Eigen::Matrix4d cam2body_;//四维矩阵，从相机到机体的转换矩阵

  // flags of map state 地图的状态标志

  bool occ_need_update_, local_updated_;
  bool has_first_depth_;
  bool has_odom_, has_cloud_;

  // depth image projected point cloud 
  //点云图提供的深度图像信息

  vector<Eigen::Vector3d> proj_points_;
  int proj_points_cnt;

  // flag buffers for speeding up raycasting用于加速光线投射的标志存储区

  vector<short> count_hit_, count_hit_and_miss_;
  vector<char> flag_traverse_, flag_rayend_;
  char raycast_num_;
  queue<Eigen::Vector3i> cache_voxel_;

  // range of updating grid 更新网格的范围

  Eigen::Vector3i local_bound_min_, local_bound_max_;//局部地图的最大位置和最小位置

  // computation time 计算时间

  double fuse_time_, max_fuse_time_;
  int update_num_;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW//大致意思是为了指针对齐，但是还是不太懂
};

class GridMap {//定义Gridmap类
public:
  GridMap() {}//啥意思？？？？
  ~GridMap() {}

  enum { POSE_STAMPED = 1, ODOMETRY = 2, INVALID_IDX = -10000 };//枚举变量,https://www.cnblogs.com/ForFreeDom/archive/2012/03/22/2412055.html

  // occupancy map management  占用地图管理   
  //许多函数的声明
  void resetBuffer();//初始化
  void resetBuffer(Eigen::Vector3d min, Eigen::Vector3d max);//初始化

  inline void posToIndex(const Eigen::Vector3d& pos, Eigen::Vector3i& id);//实际位置到角标index
  inline void indexToPos(const Eigen::Vector3i& id, Eigen::Vector3d& pos);//角标index到实际位置
  inline int toAddress(const Eigen::Vector3i& id);//从格点的角标得到其在实际存储空间（一维）中的位置
  inline int toAddress(int& x, int& y, int& z);////从三个常数得到其在实际存储空间（一维）中的位置(和上面函数实现的内容一样)
  inline bool isInMap(const Eigen::Vector3d& pos);//输入实际位置看是否在地图中
  inline bool isInMap(const Eigen::Vector3i& idx);//输入角标index位置看是否在地图中

  inline void setOccupancy(Eigen::Vector3d pos, double occ = 1);//函数中有不懂的量 设置是否为障碍，用0和1填充occupancy_buffer_向量
  inline void setOccupied(Eigen::Vector3d pos);//函数中有不懂的量 设置障碍,在occupancy_buffer_inflate_向量中赋值为1 
  inline int getOccupancy(Eigen::Vector3d pos);//函数中有不懂的量 输入位置，返回1 0 或-1（不在地图）
  inline int getOccupancy(Eigen::Vector3i id);//函数中有不懂的量 输入角标index，返回1 0 或-1（不在地图）
  inline int getInflateOccupancy(Eigen::Vector3d pos);//函数中有不懂的量 输入位置,返回occupancy_buffer_inflate_中相对应元素的值

  inline void boundIndex(Eigen::Vector3i& id);//输入index，把它限制在grid的最大和最小index之间
  inline bool isUnknown(const Eigen::Vector3i& id);//函数中有不懂的量 判断格点是否已知 输入index,返回teue or false
  inline bool isUnknown(const Eigen::Vector3d& pos);//输入位置 判断格点是否已知 转换为index后使用上一个函数,返回teue or false
  inline bool isKnownFree(const Eigen::Vector3i& id);//函数中有不懂的量 判断格点是否已知并且不是障碍物 输入index,返回teue or false
  inline bool isKnownOccupied(const Eigen::Vector3i& id);//函数中有不懂的量 判断格点是否是障碍物 输入index,返回teue or false

  void initMap(ros::NodeHandle& nh);// 里面关于视觉的函数不是很明白  输入ros节点初始化地图 

  void publishMap();
  void publishMapInflate(bool all_info = false);

  void publishUnknown();
  void publishDepth();

  bool hasDepthObservation();
  bool odomValid();
  void getRegion(Eigen::Vector3d& ori, Eigen::Vector3d& size);
  inline double getResolution();
  Eigen::Vector3d getOrigin();
  int getVoxelNum();

  typedef std::shared_ptr<GridMap> Ptr;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  MappingParameters mp_;
  MappingData md_;

  // get depth image and camera pose
  void depthPoseCallback(const sensor_msgs::ImageConstPtr& img,
                         const geometry_msgs::PoseStampedConstPtr& pose);
  void depthOdomCallback(const sensor_msgs::ImageConstPtr& img, const nav_msgs::OdometryConstPtr& odom);
  void cloudCallback(const sensor_msgs::PointCloud2ConstPtr& img);
  void odomCallback(const nav_msgs::OdometryConstPtr& odom);

  // update occupancy by raycasting
  void updateOccupancyCallback(const ros::TimerEvent& /*event*/);
  void visCallback(const ros::TimerEvent& /*event*/);

  // main update process
  void projectDepthImage();
  void raycastProcess();
  void clearAndInflateLocalMap();

  inline void inflatePoint(const Eigen::Vector3i& pt, int step, vector<Eigen::Vector3i>& pts);
  int setCacheOccupancy(Eigen::Vector3d pos, int occ);
  Eigen::Vector3d closetPointInMap(const Eigen::Vector3d& pt, const Eigen::Vector3d& camera_pt);

  // typedef message_filters::sync_policies::ExactTime<sensor_msgs::Image,
  // nav_msgs::Odometry> SyncPolicyImageOdom; typedef
  // message_filters::sync_policies::ExactTime<sensor_msgs::Image,
  // geometry_msgs::PoseStamped> SyncPolicyImagePose;
  typedef message_filters::sync_policies::ApproximateTime<sensor_msgs::Image, nav_msgs::Odometry>
      SyncPolicyImageOdom;
  typedef message_filters::sync_policies::ApproximateTime<sensor_msgs::Image, geometry_msgs::PoseStamped>
      SyncPolicyImagePose;
  typedef shared_ptr<message_filters::Synchronizer<SyncPolicyImagePose>> SynchronizerImagePose;
  typedef shared_ptr<message_filters::Synchronizer<SyncPolicyImageOdom>> SynchronizerImageOdom;

  ros::NodeHandle node_;
  shared_ptr<message_filters::Subscriber<sensor_msgs::Image>> depth_sub_;
  shared_ptr<message_filters::Subscriber<geometry_msgs::PoseStamped>> pose_sub_;
  shared_ptr<message_filters::Subscriber<nav_msgs::Odometry>> odom_sub_;
  SynchronizerImagePose sync_image_pose_;
  SynchronizerImageOdom sync_image_odom_;

  ros::Subscriber indep_cloud_sub_, indep_odom_sub_;
  ros::Publisher map_pub_, map_inf_pub_;
  ros::Publisher unknown_pub_;
  ros::Timer occ_timer_, vis_timer_;

  //
  uniform_real_distribution<double> rand_noise_;
  normal_distribution<double> rand_noise2_;
  default_random_engine eng_;
};

/* ============================== definition of inline function
 * ============================== */

inline int GridMap::toAddress(const Eigen::Vector3i& id) {//得到在内存中的位置，这个和A*算法中data数据一样，和内存存储方式有关.类似data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z]
  return id(0) * mp_.map_voxel_num_(1) * mp_.map_voxel_num_(2) + id(1) * mp_.map_voxel_num_(2) + id(2);
}

inline int GridMap::toAddress(int& x, int& y, int& z) {
  return x * mp_.map_voxel_num_(1) * mp_.map_voxel_num_(2) + y * mp_.map_voxel_num_(2) + z;
}

inline void GridMap::boundIndex(Eigen::Vector3i& id) {//当前index和最大index中取一个小的
//（不可以超过最大），然后再和0取一个最大的（不可以比0还小）。其实就是把index范围限制在了最大最小之间
  Eigen::Vector3i id1;
  id1(0) = max(min(id(0), mp_.map_voxel_num_(0) - 1), 0);//
  id1(1) = max(min(id(1), mp_.map_voxel_num_(1) - 1), 0);
  id1(2) = max(min(id(2), mp_.map_voxel_num_(2) - 1), 0);
  id = id1;
}

inline bool GridMap::isUnknown(const Eigen::Vector3i& id) {//输入index
  Eigen::Vector3i id1 = id;
  boundIndex(id1);//得到被grid大小限制的index
  return md_.occupancy_buffer_[toAddress(id1)] < mp_.clamp_min_log_ - 1e-3;//<成立返回true，否则返回false
}

inline bool GridMap::isUnknown(const Eigen::Vector3d& pos) {
  Eigen::Vector3i idc;
  posToIndex(pos, idc);
  return isUnknown(idc);
}

inline bool GridMap::isKnownFree(const Eigen::Vector3i& id) {//输入index
  Eigen::Vector3i id1 = id;
  boundIndex(id1);//限制index后的index
  int adr = toAddress(id1);//转换为一维向量中的位置

  // return md_.occupancy_buffer_[adr] >= mp_.clamp_min_log_ &&
  //     md_.occupancy_buffer_[adr] < mp_.min_occupancy_log_;
  return md_.occupancy_buffer_[adr] >= mp_.clamp_min_log_ && md_.occupancy_buffer_inflate_[adr] == 0;//前一个判断代表格点已知，后一个代表不是障碍物
}

inline bool GridMap::isKnownOccupied(const Eigen::Vector3i& id) {//输入index
  Eigen::Vector3i id1 = id;
  boundIndex(id1);
  int adr = toAddress(id1);

  return md_.occupancy_buffer_inflate_[adr] == 1;//是障碍物
}

inline void GridMap::setOccupied(Eigen::Vector3d pos) {
  if (!isInMap(pos)) return;//不在地图里，错误

  Eigen::Vector3i id;
  posToIndex(pos, id);//转化为角标index

  md_.occupancy_buffer_inflate_[id(0) * mp_.map_voxel_num_(1) * mp_.map_voxel_num_(2) +
                                id(1) * mp_.map_voxel_num_(2) + id(2)] = 1;//内存中这个格点的位置设置为1
}

inline void GridMap::setOccupancy(Eigen::Vector3d pos, double occ) {//输入实际位置和是否是障碍的标志（occ）
  if (occ != 1 && occ != 0) {//不是0也不是1，是错误的
    cout << "occ value error!" << endl;
    return;
  }

  if (!isInMap(pos)) return;//不在地图里，错误

  Eigen::Vector3i id;
  posToIndex(pos, id);//实际位置转化为角标index

  md_.occupancy_buffer_[toAddress(id)] = occ;//内存中这个格点的位置设置为occ（1：是障碍；0：不是障碍）
}

inline int GridMap::getOccupancy(Eigen::Vector3d pos) {
  if (!isInMap(pos)) return -1;//不在地图里，return -1

  Eigen::Vector3i id;
  posToIndex(pos, id);//实际位置转化为角标index

  return md_.occupancy_buffer_[toAddress(id)] > mp_.min_occupancy_log_ ? 1 : 0;//如果occupancy_buffer_[toAddress(id)] > min_occupancy_log_ return 1;else return 0
}

inline int GridMap::getInflateOccupancy(Eigen::Vector3d pos) {
  if (!isInMap(pos)) return -1;//不在地图里，return -1

  Eigen::Vector3i id;
  posToIndex(pos, id);//实际位置转化为角标index

  return int(md_.occupancy_buffer_inflate_[toAddress(id)]);
}

inline int GridMap::getOccupancy(Eigen::Vector3i id) {
  if (id(0) < 0 || id(0) >= mp_.map_voxel_num_(0) || id(1) < 0 || id(1) >= mp_.map_voxel_num_(1) ||
      id(2) < 0 || id(2) >= mp_.map_voxel_num_(2))
    return -1;//不在地图里，return -1

  return md_.occupancy_buffer_[toAddress(id)] > mp_.min_occupancy_log_ ? 1 : 0;//如果occupancy_buffer_[toAddress(id)] > min_occupancy_log_ return 1;else return 0
}

inline bool GridMap::isInMap(const Eigen::Vector3d& pos) {//输入实际位置看是否在地图中
  if (pos(0) < mp_.map_min_boundary_(0) + 1e-4 || pos(1) < mp_.map_min_boundary_(1) + 1e-4 ||
      pos(2) < mp_.map_min_boundary_(2) + 1e-4) {
    // cout << "less than min range!" << endl;
    return false;//比最小的小
  }
  if (pos(0) > mp_.map_max_boundary_(0) - 1e-4 || pos(1) > mp_.map_max_boundary_(1) - 1e-4 ||
      pos(2) > mp_.map_max_boundary_(2) - 1e-4) {
    return false;//比最大的大
  }
  return true;
}

inline bool GridMap::isInMap(const Eigen::Vector3i& idx) {//输入角标index看是否在地图中
  if (idx(0) < 0 || idx(1) < 0 || idx(2) < 0) {
    return false;
  }
  if (idx(0) > mp_.map_voxel_num_(0) - 1 || idx(1) > mp_.map_voxel_num_(1) - 1 ||
      idx(2) > mp_.map_voxel_num_(2) - 1) {
    return false;
  }
  return true;
}

inline void GridMap::posToIndex(const Eigen::Vector3d& pos, Eigen::Vector3i& id) {//从实际位置转换到格点坐标
  for (int i = 0; i < 3; ++i) id(i) = floor((pos(i) - mp_.map_origin_(i)) * mp_.resolution_inv_);//实际位置减去地图的最小位置向下取整，再乘比例尺
}

inline void GridMap::indexToPos(const Eigen::Vector3i& id, Eigen::Vector3d& pos) {//从格点坐标转换到实际位置
  for (int i = 0; i < 3; ++i) pos(i) = (id(i) + 0.5) * mp_.resolution_ + mp_.map_origin_(i);
}

inline void GridMap::inflatePoint(const Eigen::Vector3i& pt, int step, vector<Eigen::Vector3i>& pts) {//点的膨胀
  int num = 0;

  /* ---------- + shape inflate ---------- */
  // for (int x = -step; x <= step; ++x)
  // {
  //   if (x == 0)
  //     continue;
  //   pts[num++] = Eigen::Vector3i(pt(0) + x, pt(1), pt(2));
  // }
  // for (int y = -step; y <= step; ++y)
  // {
  //   if (y == 0)
  //     continue;
  //   pts[num++] = Eigen::Vector3i(pt(0), pt(1) + y, pt(2));
  // }
  // for (int z = -1; z <= 1; ++z)
  // {
  //   pts[num++] = Eigen::Vector3i(pt(0), pt(1), pt(2) + z);
  // }

  /* ---------- all inflate ---------- */
  for (int x = -step; x <= step; ++x)
    for (int y = -step; y <= step; ++y)
      for (int z = -step; z <= step; ++z) {
        pts[num++] = Eigen::Vector3i(pt(0) + x, pt(1) + y, pt(2) + z);//获得这个点膨胀后的点集合，以此点为中心，step为半径向两边扩展
      }
}

inline double GridMap::getResolution() { return mp_.resolution_; }

#endif
