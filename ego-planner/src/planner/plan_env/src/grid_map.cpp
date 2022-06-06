#include "plan_env/grid_map.h"

// #define current_img_ md_.depth_image_[image_cnt_ & 1]
// #define last_img_ md_.depth_image_[!(image_cnt_ & 1)]

void GridMap::initMap(ros::NodeHandle &nh)//输入节点nh
{
  node_ = nh;//

  /* get parameter 从文件中获得node中的参数数据*/
  double x_size, y_size, z_size;
  //从前到后依次为launch中的关键字 实际变量名 默认值（获取参数失败时会自动设置默认参数值）https://blog.csdn.net/sinat_16643223/article/details/112441563
  node_.param("grid_map/resolution", mp_.resolution_, -1.0);//以下数据在plan_manage/launch/advanced_param.xml72-79
  node_.param("grid_map/map_size_x", x_size, -1.0);//地图的大小 40 40 3
  node_.param("grid_map/map_size_y", y_size, -1.0);
  node_.param("grid_map/map_size_z", z_size, -1.0);
  node_.param("grid_map/local_update_range_x", mp_.local_update_range_(0), -1.0);//地图的局部更新范围 5.5 5.5 4.5
  node_.param("grid_map/local_update_range_y", mp_.local_update_range_(1), -1.0);
  node_.param("grid_map/local_update_range_z", mp_.local_update_range_(2), -1.0);
  node_.param("grid_map/obstacles_inflation", mp_.obstacles_inflation_, -1.0);//障碍物膨胀系数  0.099
//以下数据在plan_manage/launch/advanced_param.xml83-86
// camera parameter 相机的参数
  node_.param("grid_map/fx", mp_.fx_, -1.0);
  node_.param("grid_map/fy", mp_.fy_, -1.0);//相机焦距 
  node_.param("grid_map/cx", mp_.cx_, -1.0);
  node_.param("grid_map/cy", mp_.cy_, -1.0);//像素平面中心点移动到左上角的偏移量

//这些参数什么意思？？？？？？？和视觉有关？？是否使用深度滤波器
  node_.param("grid_map/use_depth_filter", mp_.use_depth_filter_, true);//是否使用深度滤波器 true
  node_.param("grid_map/depth_filter_tolerance", mp_.depth_filter_tolerance_, -1.0);//0.15
  node_.param("grid_map/depth_filter_maxdist", mp_.depth_filter_maxdist_, -1.0);//5
  node_.param("grid_map/depth_filter_mindist", mp_.depth_filter_mindist_, -1.0);//0.2
  node_.param("grid_map/depth_filter_margin", mp_.depth_filter_margin_, -1);//1
  node_.param("grid_map/k_depth_scaling_factor", mp_.k_depth_scaling_factor_, -1.0);//1000
  node_.param("grid_map/skip_pixel", mp_.skip_pixel_, -1);// 一个像素的边长（这里认为一个像素是一个正方形，所以只需要边长一个值）
  
  //以下数据在plan_manage/launch/advanced_param.xml96-102
   // <!-- local fusion -->这些概率和ray-cast代表什么意思？???概率栅格地图的变量
  node_.param("grid_map/p_hit", mp_.p_hit_, 0.70);//障碍物概率 0.65
  node_.param("grid_map/p_miss", mp_.p_miss_, 0.35);//碰撞概率 0.35
  node_.param("grid_map/p_min", mp_.p_min_, 0.12);//不碰撞概率0.12
  node_.param("grid_map/p_max", mp_.p_max_, 0.97);//最大概率0.90
  node_.param("grid_map/p_occ", mp_.p_occ_, 0.80);//障碍物概率 0.8
  //raycast射线长度
  node_.param("grid_map/min_ray_length", mp_.min_ray_length_, -0.1);//0.1
  node_.param("grid_map/max_ray_length", mp_.max_ray_length_, -0.1);//4.5

 //以下数据在plan_manage/launch/advanced_param.xml104-105 什么意思？???
  node_.param("grid_map/visualization_truncate_height", mp_.visualization_truncate_height_, 999.0);//截断高度 2.4
  node_.param("grid_map/virtual_ceil_height", mp_.virtual_ceil_height_, -0.1);//屋顶高度（主要在局部地图建立的时候限制此高度以上为障碍物，即不可以飞行超过此高度）2.5
 
  //以下数据在plan_manage/launch/advanced_param.xml106-107 什么意思？???
  node_.param("grid_map/show_occ_time", mp_.show_occ_time_, false);
  node_.param("grid_map/pose_type", mp_.pose_type_, 1);

  node_.param("grid_map/frame_id", mp_.frame_id_, string("world"));
  node_.param("grid_map/local_map_margin", mp_.local_map_margin_, 1);//局部地图的边缘 30
  node_.param("grid_map/ground_height", mp_.ground_height_, 1.0);//地面高度 -0.01

  mp_.resolution_inv_ = 1 / mp_.resolution_;//比例尺求倒数 is 10
  mp_.map_origin_ = Eigen::Vector3d(-x_size / 2.0, -y_size / 2.0, mp_.ground_height_);//地图起始位置 -20 -20 -0.01
  mp_.map_size_ = Eigen::Vector3d(x_size, y_size, z_size);//地图大小 40 40 3

//对这些概率数值取logit函数值log((x) / (1 - (x))) 把自变量从(0,1)连续单调地映射到正负无穷https://blog.csdn.net/songyunli1111/article/details/82285156
//概率栅格地图建立 
  mp_.prob_hit_log_ = logit(mp_.p_hit_);
  mp_.prob_miss_log_ = logit(mp_.p_miss_);
  mp_.clamp_min_log_ = logit(mp_.p_min_);
  mp_.clamp_max_log_ = logit(mp_.p_max_);
  mp_.min_occupancy_log_ = logit(mp_.p_occ_);
  mp_.unknown_flag_ = 0.01;
//输出
  cout << "hit: " << mp_.prob_hit_log_ << endl;
  cout << "miss: " << mp_.prob_miss_log_ << endl;
  cout << "min log: " << mp_.clamp_min_log_ << endl;
  cout << "max: " << mp_.clamp_max_log_ << endl;
  cout << "thresh log: " << mp_.min_occupancy_log_ << endl;

  for (int i = 0; i < 3; ++i)
    mp_.map_voxel_num_(i) = ceil(mp_.map_size_(i) / mp_.resolution_);//栅格地图的最大角标

  mp_.map_min_boundary_ = mp_.map_origin_;//实际地图的最小边界 -20 -20 -0.01
  mp_.map_max_boundary_ = mp_.map_origin_ + mp_.map_size_;//实际地图的最大边界20 20 2.99

  // initialize data buffers

  int buffer_size = mp_.map_voxel_num_(0) * mp_.map_voxel_num_(1) * mp_.map_voxel_num_(2);//初始化grid变量存储空间的大小即x*y*z的大小

  md_.occupancy_buffer_ = vector<double>(buffer_size, mp_.clamp_min_log_ - mp_.unknown_flag_);//概率栅格地图建立时的向量，当激光探测后，根据其中元素的大小和判断是否为障碍物的阈值来判断是否为障碍物 初始化数值为mp_.clamp_min_log_ - mp_.unknown_flag_即所有栅格状态是未知的
  md_.occupancy_buffer_inflate_ = vector<char>(buffer_size, 0);// 最终确认是不是障碍物的向量 即内部元素为0 或1 初始化数值为0 初始认为多所有点都是free的

  md_.count_hit_and_miss_ = vector<short>(buffer_size, 0);
  md_.count_hit_ = vector<short>(buffer_size, 0);
  md_.flag_rayend_ = vector<char>(buffer_size, -1);
  md_.flag_traverse_ = vector<char>(buffer_size, -1);//道理相同

  md_.raycast_num_ = 0;//激光雷达探测次数

  md_.proj_points_.resize(640 * 480 / mp_.skip_pixel_ / mp_.skip_pixel_);//初始化存储一个深度图像元素的矩阵，640和480为一张图片的长和宽 skip_pixel_为一个像素的边长
  //resize设置大小，参考https://blog.csdn.net/jiayizhenzhenyijia/article/details/97898577?spm=1001.2101.3001.6661.1&utm_medium=distribute.pc_relevant_t0.none-task-blog-2%7Edefault%7ECTRLIST%7ERate-1.pc_relevant_paycolumn_v3&depth_1-utm_source=distribute.pc_relevant_t0.none-task-blog-2%7Edefault%7ECTRLIST%7ERate-1.pc_relevant_paycolumn_v3&utm_relevant_index=1
  md_.proj_points_cnt = 0;
  md_.cam2body_ << 0.0, 0.0, 1.0, 0.0,
      -1.0, 0.0, 0.0, 0.0,
      0.0, -1.0, 0.0, -0.02,
      0.0, 0.0, 0.0, 1.0;//初始化从相机到飞机机体的旋转矩阵

  /* init callback *///这块具体是在干什么？？？
//给定不同的pose_type_ 会有着不同的结果，但是这是什么意思？
  depth_sub_.reset(new message_filters::Subscriber<sensor_msgs::Image>(node_, "/grid_map/depth", 50));//订阅深度图像信息

  if (mp_.pose_type_ == POSE_STAMPED)//此状态订阅相机在是世界坐标系的坐标
  {
    pose_sub_.reset(
        new message_filters::Subscriber<geometry_msgs::PoseStamped>(node_, "/grid_map/pose", 25));

    sync_image_pose_.reset(new message_filters::Synchronizer<SyncPolicyImagePose>(
        SyncPolicyImagePose(100), *depth_sub_, *pose_sub_));
    sync_image_pose_->registerCallback(boost::bind(&GridMap::depthPoseCallback, this, _1, _2));
  }
  else if (mp_.pose_type_ == ODOMETRY)//此状态订阅机体在是世界坐标系的坐标
  {
    odom_sub_.reset(new message_filters::Subscriber<nav_msgs::Odometry>(node_, "/grid_map/odom", 100));

    sync_image_odom_.reset(new message_filters::Synchronizer<SyncPolicyImageOdom>(
        SyncPolicyImageOdom(100), *depth_sub_, *odom_sub_));
    sync_image_odom_->registerCallback(boost::bind(&GridMap::depthOdomCallback, this, _1, _2));
  }

  // use odometry and point cloud使用里程计和点云地图这块具体是在干什么？？？
  //subscribe参考https://blog.csdn.net/weixin_38258767/article/details/104243642
  indep_cloud_sub_ =
      node_.subscribe<sensor_msgs::PointCloud2>("/grid_map/cloud", 10, &GridMap::cloudCallback, this);
  indep_odom_sub_ =
      node_.subscribe<nav_msgs::Odometry>("/grid_map/odom", 10, &GridMap::odomCallback, this);

  occ_timer_ = node_.createTimer(ros::Duration(0.05), &GridMap::updateOccupancyCallback, this);
  vis_timer_ = node_.createTimer(ros::Duration(0.05), &GridMap::visCallback, this);

  map_pub_ = node_.advertise<sensor_msgs::PointCloud2>("/grid_map/occupancy", 10);
  map_inf_pub_ = node_.advertise<sensor_msgs::PointCloud2>("/grid_map/occupancy_inflate", 10);

  unknown_pub_ = node_.advertise<sensor_msgs::PointCloud2>("/grid_map/unknown", 10);

  md_.occ_need_update_ = false;
  md_.local_updated_ = false;
  md_.has_first_depth_ = false;
  md_.has_odom_ = false;
  md_.has_cloud_ = false;
  md_.image_cnt_ = 0;

  md_.fuse_time_ = 0.0;
  md_.update_num_ = 0;
  md_.max_fuse_time_ = 0.0;

  // rand_noise_ = uniform_real_distribution<double>(-0.2, 0.2);
  // rand_noise2_ = normal_distribution<double>(0, 0.2);
  // random_device rd;
  // eng_ = default_random_engine(rd());
}

void GridMap::resetBuffer()
{
  Eigen::Vector3d min_pos = mp_.map_min_boundary_;//实际地图最小位置
  Eigen::Vector3d max_pos = mp_.map_max_boundary_;//实际地图最大位置

  resetBuffer(min_pos, max_pos);

  md_.local_bound_min_ = Eigen::Vector3i::Zero();//栅格地图的最小边界角标
  md_.local_bound_max_ = mp_.map_voxel_num_ - Eigen::Vector3i::Ones();//栅格地图的最大边界角标
}

void GridMap::resetBuffer(Eigen::Vector3d min_pos, Eigen::Vector3d max_pos)//给定最大位置和最小位置.初始化地图的相关信息
{

  Eigen::Vector3i min_id, max_id;//定义栅格地图最大最小角标index
  posToIndex(min_pos, min_id);//得到栅格地图最小角标index
  posToIndex(max_pos, max_id);//得到栅格地图最大角标index

  boundIndex(min_id);//限制角标大小，得到真正的最小角标index
  boundIndex(max_id);//限制角标大小，得到真正的最大角标index

  /* reset occ and dist buffer *///初始化地图障碍和距离
  for (int x = min_id(0); x <= max_id(0); ++x)//循环遍历从最小到最大的角标index
    for (int y = min_id(1); y <= max_id(1); ++y)
      for (int z = min_id(2); z <= max_id(2); ++z)
      {
        md_.occupancy_buffer_inflate_[toAddress(x, y, z)] = 0;//初始化所有边界内的栅格状态均为free
      }
}

int GridMap::setCacheOccupancy(Eigen::Vector3d pos, int occ)//输入一个实际位置和是否发生碰撞的标志（occ）
{
  if (occ != 1 && occ != 0)
    return INVALID_IDX;//返回无效角标

  Eigen::Vector3i id;
  posToIndex(pos, id);//从实际位置转化到角标index
  int idx_ctns = toAddress(id);//得到在存储空间中的位置

  md_.count_hit_and_miss_[idx_ctns] += 1;//探测到一次就+1


  if (md_.count_hit_and_miss_[idx_ctns] == 1)//如果是第一次探索到这个voxel
  {
    md_.cache_voxel_.push(id);//把当前voxel的角标index向量加到该queue队列中
  }

  if (occ == 1)//如果occ==1（碰撞了）
    md_.count_hit_[idx_ctns] += 1;//碰撞就+1

  return idx_ctns;
}

void GridMap::projectDepthImage()//通过深度图像得到，将图像上的点转化为世界坐标系，存储在md_.proj_points_
{
  // md_.proj_points_.clear();
  md_.proj_points_cnt = 0;//初始化proj_points_cnt为0

  uint16_t *row_ptr;//定义一个unit16_t类型的指针
  // int cols = current_img_.cols, rows = current_img_.rows;
  int cols = md_.depth_image_.cols;//深度图像列数
  int rows = md_.depth_image_.rows;//深度图像行数

  double depth;

  Eigen::Matrix3d camera_r = md_.camera_q_.toRotationMatrix();//从四元数变成相机坐标系到世界坐标系的旋转矩阵
//
  // cout << "rotate: " << md_.camera_q_.toRotationMatrix() << endl;
  // std::cout << "pos in proj: " << md_.camera_pos_ << std::endl;

  if (!mp_.use_depth_filter_)//判断是否使用深度过滤器，如果不使用
  {
    for (int v = 0; v < rows; v++)//循环像素坐标的行
    {
      row_ptr = md_.depth_image_.ptr<uint16_t>(v);//得到这一行的指针

      for (int u = 0; u < cols; u++)//循环像素坐标的列
      {

        Eigen::Vector3d proj_pt;
        depth = (*row_ptr++) / mp_.k_depth_scaling_factor_;//得到相机坐标系下的深度
        //
        proj_pt(0) = (u - mp_.cx_) * depth / mp_.fx_;//u - mp_.cx_这是从像素坐标系变到了图像坐标系 后面是从图像坐标系到相机坐标系 参考pdf第六页
        proj_pt(1) = (v - mp_.cy_) * depth / mp_.fy_;
        proj_pt(2) = depth;//从像素坐标系到相机坐标系

        proj_pt = camera_r * proj_pt + md_.camera_pos_;//从相机坐标系到世界坐标系
        //

        if (u == 320 && v == 240)//如果到了最后一个像素
          std::cout << "depth: " << depth << std::endl;
        md_.proj_points_[md_.proj_points_cnt++] = proj_pt;//把相机拍下来的点转换成世界坐标系下的坐标并存储到该vector中
      }
    }
  }
  /* use depth filter */
  else//如果使用深度过滤器这块的理论知识不是很理解？？？？？？？？？？？？？？？？？？？？？？？？？？
  {

    if (!md_.has_first_depth_)//如果这个值为false就把它设置成true
      md_.has_first_depth_ = true;
    else//如果本来就是true
    {
      Eigen::Vector3d pt_cur, pt_world, pt_reproj;//定义三个向量

      Eigen::Matrix3d last_camera_r_inv;//逆四元数：上一个相机到世界的四元数旋转的逆四元数
      last_camera_r_inv = md_.last_camera_q_.inverse();//从四元数到旋转矩阵
      const double inv_factor = 1.0 / mp_.k_depth_scaling_factor_;//1/1000

      for (int v = mp_.depth_filter_margin_; v < rows - mp_.depth_filter_margin_; v += mp_.skip_pixel_)
      //循环深度图像的行，从即depth_filter_margin_1开始，到 rows - mp_.depth_filter_margin_即320-1=319，每次增加 skip_pixel_即2
      {
        row_ptr = md_.depth_image_.ptr<uint16_t>(v) + mp_.depth_filter_margin_;//这一行的指针

        for (int u = mp_.depth_filter_margin_; u < cols - mp_.depth_filter_margin_;
             u += mp_.skip_pixel_)//同理，循环深度图像的列
        {

          depth = (*row_ptr) * inv_factor;//得到深度
          row_ptr = row_ptr + mp_.skip_pixel_;//下一列元素的指针

          // filter depth
          // depth += rand_noise_(eng_);
          // if (depth > 0.01) depth += rand_noise2_(eng_);

          if (*row_ptr == 0)//深度为0
          {
            depth = mp_.max_ray_length_ + 0.1;//深度值改为4.5+0.1=4.6
          }
          else if (depth < mp_.depth_filter_mindist_)//如果深度比最小的0.2还小
          {
            continue;//跳过
          }
          else if (depth > mp_.depth_filter_maxdist_)//如果深度比最大的5还大
          {
            depth = mp_.max_ray_length_ + 0.1;//深度值改为4.5+0.1=4.6
          }

          // project to world frame
          pt_cur(0) = (u - mp_.cx_) * depth / mp_.fx_;
          pt_cur(1) = (v - mp_.cy_) * depth / mp_.fy_;
          pt_cur(2) = depth;//从像素坐标系到相机坐标系

          pt_world = camera_r * pt_cur + md_.camera_pos_;//得到pt_world
          // if (!isInMap(pt_world)) {
          //   pt_world = closetPointInMap(pt_world, md_.camera_pos_);
          // }

          md_.proj_points_[md_.proj_points_cnt++] = pt_world;//从相机坐标系到世界坐标系

          // check consistency with last image, disabled... 检查与上一张图片的一致性 
          //false后是不是没用？？？？
          if (false)
          {
            pt_reproj = last_camera_r_inv * (pt_world - md_.last_camera_pos_);
            double uu = pt_reproj.x() * mp_.fx_ / pt_reproj.z() + mp_.cx_;
            double vv = pt_reproj.y() * mp_.fy_ / pt_reproj.z() + mp_.cy_;

            if (uu >= 0 && uu < cols && vv >= 0 && vv < rows)
            {
              if (fabs(md_.last_depth_image_.at<uint16_t>((int)vv, (int)uu) * inv_factor -
                       pt_reproj.z()) < mp_.depth_filter_tolerance_)
              {
                md_.proj_points_[md_.proj_points_cnt++] = pt_world;
              }
            }
            else
            {
              md_.proj_points_[md_.proj_points_cnt++] = pt_world;
            }
          }
        }
      }
    }
  }

  /* maintain camera pose for consistency check ,迭代相机的相关数据*/

  md_.last_camera_pos_ = md_.camera_pos_;
  md_.last_camera_q_ = md_.camera_q_;
  md_.last_depth_image_ = md_.depth_image_;
}

void GridMap::raycastProcess()//对深度图像拍摄到的所有点进行raycast，判断是不是障碍物 更新栅格地图的概率
{
  // if (md_.proj_points_.size() == 0)
  if (md_.proj_points_cnt == 0)
    return;

  ros::Time t1, t2;

  md_.raycast_num_ += 1;

  int vox_idx;
  double length;

  // bounding box of updated region 更新区域的边界 
  double min_x = mp_.map_max_boundary_(0);
  double min_y = mp_.map_max_boundary_(1);
  double min_z = mp_.map_max_boundary_(2);

  double max_x = mp_.map_min_boundary_(0);
  double max_y = mp_.map_min_boundary_(1);
  double max_z = mp_.map_min_boundary_(2);

  RayCaster raycaster;
  Eigen::Vector3d half = Eigen::Vector3d(0.5, 0.5, 0.5);//voxel的中点
  Eigen::Vector3d ray_pt, pt_w;

  for (int i = 0; i < md_.proj_points_cnt; ++i)//循环深度图像中得到的每个点
  {
    pt_w = md_.proj_points_[i];//这个图像中点的世界坐标

    // set flag for projected point

    if (!isInMap(pt_w))//如果图像中点不在地图
    {
      pt_w = closetPointInMap(pt_w, md_.camera_pos_);//图像中点拉回到地图边缘

      length = (pt_w - md_.camera_pos_).norm();//图像中点到相机的距离
      if (length > mp_.max_ray_length_)//超出了最大探测长度
      {
        pt_w = (pt_w - md_.camera_pos_) / length * mp_.max_ray_length_ + md_.camera_pos_;//把这个点拉到raycast能探测到的最大位置
      }
      vox_idx = setCacheOccupancy(pt_w, 0);//得到其所在voxel的index在存储buffer中的角标
    }
    else//图像中点本来就在地图中
    {
      length = (pt_w - md_.camera_pos_).norm();

      if (length > mp_.max_ray_length_)//超过了最大范围
      {
        pt_w = (pt_w - md_.camera_pos_) / length * mp_.max_ray_length_ + md_.camera_pos_;
        vox_idx = setCacheOccupancy(pt_w, 0);
      }
      else//本来就在地图中且没超过raycast的探测范围
      {
        vox_idx = setCacheOccupancy(pt_w, 1);//此点设置为障碍物，md_.count_hit_[vox_idx] += 1
      }
    }

    max_x = max(max_x, pt_w(0));
    max_y = max(max_y, pt_w(1));
    max_z = max(max_z, pt_w(2));//第一次比较是图像中点和实际地图最小边界较大的一个
    //后面便是比较两个图像中的点大小，然后取其中较大的一个。最终得到图像中所有点中xyz的最大坐标

    min_x = min(min_x, pt_w(0));
    min_y = min(min_y, pt_w(1));
    min_z = min(min_z, pt_w(2));//第一次比较是图像中点和实际地图最大边界较小的一个
    //后面便是比较两个图像中的点大小，然后取其中较小的一个。最终得到图像中所有点中xyz的最小坐标



    // raycasting between camera center and point 相机中心和点之间进行raycast

    if (vox_idx != INVALID_IDX)//不是无效位置
    {
      if (md_.flag_rayend_[vox_idx] == md_.raycast_num_)
      //本次激光雷达的探测次数==这个点上一次被激光雷达探测的次数（即在这次雷达探测中，这个voxel被反复探测了）
      {
        continue;
      }
      else
      {
        md_.flag_rayend_[vox_idx] = md_.raycast_num_;
      }
      //将本次激光探测的次数设置为这个点上一次被激光雷达探测的次数
    }

    raycaster.setInput(pt_w / mp_.resolution_, md_.camera_pos_ / mp_.resolution_);//设置raycast的起点和终点

    while (raycaster.step(ray_pt))//当还没有到本次raycast的终点时
    {
      Eigen::Vector3d tmp = (ray_pt + half) * mp_.resolution_;//现在点的实际位置  voxel角标+0.5，即体素的中心再乘上比例尺 得到实际的位置
      length = (tmp - md_.camera_pos_).norm();//现在点到相机中心的距离

      // if (length < mp_.min_ray_length_) break;

      vox_idx = setCacheOccupancy(tmp, 0);//得到其所在voxel的index在存储buffer中的角标

      if (vox_idx != INVALID_IDX)
      {
        if (md_.flag_traverse_[vox_idx] == md_.raycast_num_)
        {
          break;
        }
        else
        {
          md_.flag_traverse_[vox_idx] = md_.raycast_num_;//设置激光雷达的探测次数
        }
      }
    }
  }

  min_x = min(min_x, md_.camera_pos_(0));
  min_y = min(min_y, md_.camera_pos_(1));
  min_z = min(min_z, md_.camera_pos_(2));//图像中所有点再和相机比较一下，取最小的

  max_x = max(max_x, md_.camera_pos_(0));
  max_y = max(max_y, md_.camera_pos_(1));
  max_z = max(max_z, md_.camera_pos_(2));
  max_z = max(max_z, mp_.ground_height_);//图像中所有点再和相机比较一下，取最大的

  posToIndex(Eigen::Vector3d(max_x, max_y, max_z), md_.local_bound_max_);//得到本次视野所探测到的最大voxel的index
  posToIndex(Eigen::Vector3d(min_x, min_y, min_z), md_.local_bound_min_);//得到本次视野所探测到的最小voxel的index
  boundIndex(md_.local_bound_min_);
  boundIndex(md_.local_bound_max_);

  md_.local_updated_ = true;

  // update occupancy cached in queue
  Eigen::Vector3d local_range_min = md_.camera_pos_ - mp_.local_update_range_;
  Eigen::Vector3d local_range_max = md_.camera_pos_ + mp_.local_update_range_;//
  //局部地图更新的范围 是5.5 5.5 4.5，用当前的相机位置+上为局部更新的最大位置，-则为最小位置
  
  Eigen::Vector3i min_id, max_id;
  posToIndex(local_range_min, min_id);
  posToIndex(local_range_max, max_id);
  boundIndex(min_id);
  boundIndex(max_id);//转换成index

  // std::cout << "cache all: " << md_.cache_voxel_.size() << std::endl;

  while (!md_.cache_voxel_.empty())//所有被探测到的voxel的队列不为空
  {

    Eigen::Vector3i idx = md_.cache_voxel_.front();//取出第一个voxel
    int idx_ctns = toAddress(idx);
    md_.cache_voxel_.pop();//删除

    double log_odds_update =
        md_.count_hit_[idx_ctns] >= md_.count_hit_and_miss_[idx_ctns] - md_.count_hit_[idx_ctns] ? mp_.prob_hit_log_ : mp_.prob_miss_log_;
        //判断为碰撞的次数是否>检测到这个voxel的次数-碰撞的次数
    
    md_.count_hit_[idx_ctns] = md_.count_hit_and_miss_[idx_ctns] = 0;//清0

    if (log_odds_update >= 0 && md_.occupancy_buffer_[idx_ctns] >= mp_.clamp_max_log_)
    {
      continue;
    }
    else if (log_odds_update <= 0 && md_.occupancy_buffer_[idx_ctns] <= mp_.clamp_min_log_)
    {
      md_.occupancy_buffer_[idx_ctns] = mp_.clamp_min_log_;
      continue;
    }

    bool in_local = idx(0) >= min_id(0) && idx(0) <= max_id(0) && idx(1) >= min_id(1) &&
                    idx(1) <= max_id(1) && idx(2) >= min_id(2) && idx(2) <= max_id(2);//清0
    if (!in_local)//不在,概率设置为最低
    {
      md_.occupancy_buffer_[idx_ctns] = mp_.clamp_min_log_;
    }
    //否则加上log_odds_update
    md_.occupancy_buffer_[idx_ctns] =
        std::min(std::max(md_.occupancy_buffer_[idx_ctns] + log_odds_update, mp_.clamp_min_log_),
                 mp_.clamp_max_log_);//
  }
}

Eigen::Vector3d GridMap::closetPointInMap(const Eigen::Vector3d &pt, const Eigen::Vector3d &camera_pt)
//把地图范围外的点pt拉到地图边缘。
{
  Eigen::Vector3d diff = pt - camera_pt;
  Eigen::Vector3d max_tc = mp_.map_max_boundary_ - camera_pt;
  Eigen::Vector3d min_tc = mp_.map_min_boundary_ - camera_pt;

  double min_t = 1000000;

  for (int i = 0; i < 3; ++i)
  {
    if (fabs(diff[i]) > 0)
    {

      double t1 = max_tc[i] / diff[i];
      if (t1 > 0 && t1 < min_t)
        min_t = t1;

      double t2 = min_tc[i] / diff[i];
      if (t2 > 0 && t2 < min_t)
        min_t = t2;
    }
  }

  return camera_pt + (min_t - 1e-3) * diff;
}//把地图范围外的点pt拉到地图边缘

void GridMap::clearAndInflateLocalMap()//清空并膨胀局部地图
{
  /*clear outside local*/
  const int vec_margin = 5;
  // Eigen::Vector3i min_vec_margin = min_vec - Eigen::Vector3i(vec_margin,
  // vec_margin, vec_margin); Eigen::Vector3i max_vec_margin = max_vec +
  // Eigen::Vector3i(vec_margin, vec_margin, vec_margin);
  Eigen::Vector3i min_cut = md_.local_bound_min_ -
                            Eigen::Vector3i(mp_.local_map_margin_, mp_.local_map_margin_, mp_.local_map_margin_);
                            //本次raycast探测到的最大voxel三个维度分别-30个voxel
  Eigen::Vector3i max_cut = md_.local_bound_max_ +
                            Eigen::Vector3i(mp_.local_map_margin_, mp_.local_map_margin_, mp_.local_map_margin_);
                            //本次raycast探测到的最大voxel三个维度分别+30个voxel
  boundIndex(min_cut);
  boundIndex(max_cut);

  Eigen::Vector3i min_cut_m = min_cut - Eigen::Vector3i(vec_margin, vec_margin, vec_margin);//-5个voxel
  Eigen::Vector3i max_cut_m = max_cut + Eigen::Vector3i(vec_margin, vec_margin, vec_margin);//+5个voxel
  boundIndex(min_cut_m);//min_cut_m里头是在 min_cut里头有一圈
  boundIndex(max_cut_m);//max_cut_m是在 max_cut外头有一圈 

  // clear data outside the local range 

  for (int x = min_cut_m(0); x <= max_cut_m(0); ++x)
    for (int y = min_cut_m(1); y <= max_cut_m(1); ++y)
    {

      for (int z = min_cut_m(2); z < min_cut(2); ++z)
      {
        int idx = toAddress(x, y, z);
        md_.occupancy_buffer_[idx] = mp_.clamp_min_log_ - mp_.unknown_flag_;
      }

      for (int z = max_cut(2) + 1; z <= max_cut_m(2); ++z)
      {
        int idx = toAddress(x, y, z);
        md_.occupancy_buffer_[idx] = mp_.clamp_min_log_ - mp_.unknown_flag_;
      }
    }

  for (int z = min_cut_m(2); z <= max_cut_m(2); ++z)
    for (int x = min_cut_m(0); x <= max_cut_m(0); ++x)
    {

      for (int y = min_cut_m(1); y < min_cut(1); ++y)
      {
        int idx = toAddress(x, y, z);
        md_.occupancy_buffer_[idx] = mp_.clamp_min_log_ - mp_.unknown_flag_;
      }

      for (int y = max_cut(1) + 1; y <= max_cut_m(1); ++y)
      {
        int idx = toAddress(x, y, z);
        md_.occupancy_buffer_[idx] = mp_.clamp_min_log_ - mp_.unknown_flag_;
      }
    }

  for (int y = min_cut_m(1); y <= max_cut_m(1); ++y)
    for (int z = min_cut_m(2); z <= max_cut_m(2); ++z)
    {

      for (int x = min_cut_m(0); x < min_cut(0); ++x)
      {
        int idx = toAddress(x, y, z);
        md_.occupancy_buffer_[idx] = mp_.clamp_min_log_ - mp_.unknown_flag_;
      }

      for (int x = max_cut(0) + 1; x <= max_cut_m(0); ++x)
      {
        int idx = toAddress(x, y, z);
        md_.occupancy_buffer_[idx] = mp_.clamp_min_log_ - mp_.unknown_flag_;
      }
    }
//cut 和cut_m之间的那一圈清零 
//local_bound是最小的立方体，cut是包围它的一个小立方体，cut_m是一个大的立方体。本次将cut和cut_m两个立方体中间的部分occupancy_buffer_清零


  // inflate occupied voxels to compensate robot size膨胀障碍物（为了补偿无人机的尺寸）

  int inf_step = ceil(mp_.obstacles_inflation_ / mp_.resolution_);//膨胀的比例
  // int inf_step_z = 1;
  vector<Eigen::Vector3i> inf_pts(pow(2 * inf_step + 1, 3));
  // inf_pts.resize(4 * inf_step + 3);
  Eigen::Vector3i inf_pt;

  // clear outdated data
  for (int x = md_.local_bound_min_(0); x <= md_.local_bound_max_(0); ++x)
    for (int y = md_.local_bound_min_(1); y <= md_.local_bound_max_(1); ++y)
      for (int z = md_.local_bound_min_(2); z <= md_.local_bound_max_(2); ++z)
      {
        md_.occupancy_buffer_inflate_[toAddress(x, y, z)] = 0;//local_bound这个立方体内occupancy_buffer_inflate_清零
      }

  // inflate obstacles
  for (int x = md_.local_bound_min_(0); x <= md_.local_bound_max_(0); ++x)
    for (int y = md_.local_bound_min_(1); y <= md_.local_bound_max_(1); ++y)
      for (int z = md_.local_bound_min_(2); z <= md_.local_bound_max_(2); ++z)//循环local_bound立方体
      {

        if (md_.occupancy_buffer_[toAddress(x, y, z)] > mp_.min_occupancy_log_)//如果这个点的概率>是障碍物的最小概率
        {
          inflatePoint(Eigen::Vector3i(x, y, z), inf_step, inf_pts);//得到膨胀后的点集inf_pts

          for (int k = 0; k < (int)inf_pts.size(); ++k)//循环膨胀点集
          {
            inf_pt = inf_pts[k];
            int idx_inf = toAddress(inf_pt);
            if (idx_inf < 0 ||
                idx_inf >= mp_.map_voxel_num_(0) * mp_.map_voxel_num_(1) * mp_.map_voxel_num_(2))
            {
              continue;
            }
            md_.occupancy_buffer_inflate_[idx_inf] = 1;//把这个点设置为障碍物
          }
        }
      }

  // add virtual ceiling to limit flight height 加入虚拟屋顶限制飞行高度
  if (mp_.virtual_ceil_height_ > -0.5)
  {
    int ceil_id = floor((mp_.virtual_ceil_height_ - mp_.map_origin_(2)) * mp_.resolution_inv_);//限制的屋顶高度
    for (int x = md_.local_bound_min_(0); x <= md_.local_bound_max_(0); ++x)
      for (int y = md_.local_bound_min_(1); y <= md_.local_bound_max_(1); ++y)
      {
        md_.occupancy_buffer_inflate_[toAddress(x, y, ceil_id)] = 1;//这个高度全设置为障碍物（不可以超过这个高度）
      }
  }
}

void GridMap::visCallback(const ros::TimerEvent & /*event*/)
{

  publishMap();
  publishMapInflate(true);
}

void GridMap::updateOccupancyCallback(const ros::TimerEvent & /*event*/)
{
  if (!md_.occ_need_update_)//不需要更新，直接返回
    return;

  /* update occupancy */
  // ros::Time t1, t2, t3, t4;
  // t1 = ros::Time::now();

  projectDepthImage();//得到深度图像
  // t2 = ros::Time::now();
  raycastProcess();
  // t3 = ros::Time::now();

  if (md_.local_updated_)
    clearAndInflateLocalMap();//更新局部地图

  // t4 = ros::Time::now();

  // cout << setprecision(7);
  // cout << "t2=" << (t2-t1).toSec() << " t3=" << (t3-t2).toSec() << " t4=" << (t4-t3).toSec() << endl;;

  // md_.fuse_time_ += (t2 - t1).toSec();
  // md_.max_fuse_time_ = max(md_.max_fuse_time_, (t2 - t1).toSec());

  // if (mp_.show_occ_time_)
  //   ROS_WARN("Fusion: cur t = %lf, avg t = %lf, max t = %lf", (t2 - t1).toSec(),
  //            md_.fuse_time_ / md_.update_num_, md_.max_fuse_time_);

  md_.occ_need_update_ = false;
  md_.local_updated_ = false;
}

void GridMap::depthPoseCallback(const sensor_msgs::ImageConstPtr &img,
                                const geometry_msgs::PoseStampedConstPtr &pose)
{
  /* get depth image */
  cv_bridge::CvImagePtr cv_ptr;
  cv_ptr = cv_bridge::toCvCopy(img, img->encoding);

  if (img->encoding == sensor_msgs::image_encodings::TYPE_32FC1)
  {
    (cv_ptr->image).convertTo(cv_ptr->image, CV_16UC1, mp_.k_depth_scaling_factor_);
  }
  cv_ptr->image.copyTo(md_.depth_image_);

  // std::cout << "depth: " << md_.depth_image_.cols << ", " << md_.depth_image_.rows << std::endl;

  /* get pose */
  md_.camera_pos_(0) = pose->pose.position.x;
  md_.camera_pos_(1) = pose->pose.position.y;
  md_.camera_pos_(2) = pose->pose.position.z;//得到相机的位置
  md_.camera_q_ = Eigen::Quaterniond(pose->pose.orientation.w, pose->pose.orientation.x,
                                     pose->pose.orientation.y, pose->pose.orientation.z);//相机到世界坐标系的四元数
  if (isInMap(md_.camera_pos_))
  {
    md_.has_odom_ = true;//得到了相机的位置
    md_.update_num_ += 1;//地图更新次数+1
    md_.occ_need_update_ = true;
  }//如果相机在地图范围内，则md_.occ_need_update_=true，表示需要更新栅格地图
  else
  {
    md_.occ_need_update_ = false;//相机不在地图范围内则不需要更新地图
  }
}
void GridMap::odomCallback(const nav_msgs::OdometryConstPtr &odom)
{
  if (md_.has_first_depth_)
    return;

  md_.camera_pos_(0) = odom->pose.pose.position.x;
  md_.camera_pos_(1) = odom->pose.pose.position.y;
  md_.camera_pos_(2) = odom->pose.pose.position.z;//得到相机的位置
//如果还没有收到第一份深度信息，则把相机位置保存在md_.camera_pos_ 。
  md_.has_odom_ = true;
}

void GridMap::cloudCallback(const sensor_msgs::PointCloud2ConstPtr &img)
{

  pcl::PointCloud<pcl::PointXYZ> latest_cloud;
  pcl::fromROSMsg(*img, latest_cloud);

  md_.has_cloud_ = true;//得到点云数据

  if (!md_.has_odom_)
  {
    std::cout << "no odom!" << std::endl;
    return;
  }

  if (latest_cloud.points.size() == 0)
    return;

  if (isnan(md_.camera_pos_(0)) || isnan(md_.camera_pos_(1)) || isnan(md_.camera_pos_(2)))
    return;

  this->resetBuffer(md_.camera_pos_ - mp_.local_update_range_,
                    md_.camera_pos_ + mp_.local_update_range_);//局部地图数据的初始化，以相机位置为中心，local_update_range_为更新范围

  pcl::PointXYZ pt;
  Eigen::Vector3d p3d, p3d_inf;

  int inf_step = ceil(mp_.obstacles_inflation_ / mp_.resolution_);//膨胀系数
  int inf_step_z = 1;

  double max_x, max_y, max_z, min_x, min_y, min_z;

  min_x = mp_.map_max_boundary_(0);
  min_y = mp_.map_max_boundary_(1);
  min_z = mp_.map_max_boundary_(2);

  max_x = mp_.map_min_boundary_(0);
  max_y = mp_.map_min_boundary_(1);
  max_z = mp_.map_min_boundary_(2);

  for (size_t i = 0; i < latest_cloud.points.size(); ++i)//循化点云数据中的点
  {
    pt = latest_cloud.points[i];
    p3d(0) = pt.x, p3d(1) = pt.y, p3d(2) = pt.z;//三个坐标

    /* point inside update range */
    Eigen::Vector3d devi = p3d - md_.camera_pos_;//点云中的点到相机的距离
    Eigen::Vector3i inf_pt;

    if (fabs(devi(0)) < mp_.local_update_range_(0) && fabs(devi(1)) < mp_.local_update_range_(1) &&
        fabs(devi(2)) < mp_.local_update_range_(2))//这个点在局部地图范围内
    {

      /* inflate the point膨胀这个点 */
      for (int x = -inf_step; x <= inf_step; ++x)
        for (int y = -inf_step; y <= inf_step; ++y)
          for (int z = -inf_step_z; z <= inf_step_z; ++z)
          {

            p3d_inf(0) = pt.x + x * mp_.resolution_;
            p3d_inf(1) = pt.y + y * mp_.resolution_;
            p3d_inf(2) = pt.z + z * mp_.resolution_;

            max_x = max(max_x, p3d_inf(0));
            max_y = max(max_y, p3d_inf(1));
            max_z = max(max_z, p3d_inf(2));

            min_x = min(min_x, p3d_inf(0));
            min_y = min(min_y, p3d_inf(1));
            min_z = min(min_z, p3d_inf(2));

            posToIndex(p3d_inf, inf_pt);

            if (!isInMap(inf_pt))
              continue;

            int idx_inf = toAddress(inf_pt);

            md_.occupancy_buffer_inflate_[idx_inf] = 1;//设置为障碍物
          }
    }
  }

  min_x = min(min_x, md_.camera_pos_(0));
  min_y = min(min_y, md_.camera_pos_(1));
  min_z = min(min_z, md_.camera_pos_(2));

  max_x = max(max_x, md_.camera_pos_(0));
  max_y = max(max_y, md_.camera_pos_(1));
  max_z = max(max_z, md_.camera_pos_(2));

  max_z = max(max_z, mp_.ground_height_);//所有点云数据和相机位置的最大最小位置

  posToIndex(Eigen::Vector3d(max_x, max_y, max_z), md_.local_bound_max_);
  posToIndex(Eigen::Vector3d(min_x, min_y, min_z), md_.local_bound_min_);//得到本次点云数据的最大最小voxel

  boundIndex(md_.local_bound_min_);
  boundIndex(md_.local_bound_max_);
}

void GridMap::publishMap()//occupancy_buffer_中比最小的障碍物概率大的存在点云中
{

  if (map_pub_.getNumSubscribers() <= 0)
    return;

  pcl::PointXYZ pt;
  pcl::PointCloud<pcl::PointXYZ> cloud;

  Eigen::Vector3i min_cut = md_.local_bound_min_;
  Eigen::Vector3i max_cut = md_.local_bound_max_;

  int lmm = mp_.local_map_margin_ / 2;
  min_cut -= Eigen::Vector3i(lmm, lmm, lmm);
  max_cut += Eigen::Vector3i(lmm, lmm, lmm);

  boundIndex(min_cut);
  boundIndex(max_cut);

  for (int x = min_cut(0); x <= max_cut(0); ++x)
    for (int y = min_cut(1); y <= max_cut(1); ++y)
      for (int z = min_cut(2); z <= max_cut(2); ++z)
      {
        if (md_.occupancy_buffer_[toAddress(x, y, z)] < mp_.min_occupancy_log_)
          continue;

        Eigen::Vector3d pos;
        indexToPos(Eigen::Vector3i(x, y, z), pos);
        if (pos(2) > mp_.visualization_truncate_height_)
          continue;
        pt.x = pos(0);
        pt.y = pos(1);
        pt.z = pos(2);
        cloud.push_back(pt);
      }

  cloud.width = cloud.points.size();
  cloud.height = 1;
  cloud.is_dense = true;
  cloud.header.frame_id = mp_.frame_id_;
  sensor_msgs::PointCloud2 cloud_msg;

  pcl::toROSMsg(cloud, cloud_msg);
  map_pub_.publish(cloud_msg);
}

void GridMap::publishMapInflate(bool all_info)//发布膨胀栅格地图 occupancy_buffer_inflate_中为1的点存在点云中
{

  if (map_inf_pub_.getNumSubscribers() <= 0)
    return;

  pcl::PointXYZ pt;
  pcl::PointCloud<pcl::PointXYZ> cloud;

  Eigen::Vector3i min_cut = md_.local_bound_min_;
  Eigen::Vector3i max_cut = md_.local_bound_max_;

  if (all_info)
  {
    int lmm = mp_.local_map_margin_;
    min_cut -= Eigen::Vector3i(lmm, lmm, lmm);
    max_cut += Eigen::Vector3i(lmm, lmm, lmm);
  }

  boundIndex(min_cut);
  boundIndex(max_cut);

  for (int x = min_cut(0); x <= max_cut(0); ++x)
    for (int y = min_cut(1); y <= max_cut(1); ++y)
      for (int z = min_cut(2); z <= max_cut(2); ++z)
      {
        if (md_.occupancy_buffer_inflate_[toAddress(x, y, z)] == 0)
          continue;

        Eigen::Vector3d pos;
        indexToPos(Eigen::Vector3i(x, y, z), pos);
        if (pos(2) > mp_.visualization_truncate_height_)
          continue;

        pt.x = pos(0);
        pt.y = pos(1);
        pt.z = pos(2);
        cloud.push_back(pt);
      }

  cloud.width = cloud.points.size();
  cloud.height = 1;
  cloud.is_dense = true;
  cloud.header.frame_id = mp_.frame_id_;
  sensor_msgs::PointCloud2 cloud_msg;

  pcl::toROSMsg(cloud, cloud_msg);
  map_inf_pub_.publish(cloud_msg);

  // ROS_INFO("pub map");
}

void GridMap::publishUnknown()//occupancy_buffer_中<clamp_min_log_的点存在点云中
{
  pcl::PointXYZ pt;
  pcl::PointCloud<pcl::PointXYZ> cloud;

  Eigen::Vector3i min_cut = md_.local_bound_min_;
  Eigen::Vector3i max_cut = md_.local_bound_max_;

  boundIndex(max_cut);
  boundIndex(min_cut);

  for (int x = min_cut(0); x <= max_cut(0); ++x)
    for (int y = min_cut(1); y <= max_cut(1); ++y)
      for (int z = min_cut(2); z <= max_cut(2); ++z)
      {

        if (md_.occupancy_buffer_[toAddress(x, y, z)] < mp_.clamp_min_log_ - 1e-3)
        {
          Eigen::Vector3d pos;
          indexToPos(Eigen::Vector3i(x, y, z), pos);
          if (pos(2) > mp_.visualization_truncate_height_)
            continue;

          pt.x = pos(0);
          pt.y = pos(1);
          pt.z = pos(2);
          cloud.push_back(pt);
        }
      }

  cloud.width = cloud.points.size();
  cloud.height = 1;
  cloud.is_dense = true;
  cloud.header.frame_id = mp_.frame_id_;

  sensor_msgs::PointCloud2 cloud_msg;
  pcl::toROSMsg(cloud, cloud_msg);
  unknown_pub_.publish(cloud_msg);
}

bool GridMap::odomValid() { return md_.has_odom_; }

bool GridMap::hasDepthObservation() { return md_.has_first_depth_; }

Eigen::Vector3d GridMap::getOrigin() { return mp_.map_origin_; }

// int GridMap::getVoxelNum() {
//   return mp_.map_voxel_num_[0] * mp_.map_voxel_num_[1] * mp_.map_voxel_num_[2];
// }

void GridMap::getRegion(Eigen::Vector3d &ori, Eigen::Vector3d &size)
{
  ori = mp_.map_origin_, size = mp_.map_size_;
}

void GridMap::depthOdomCallback(const sensor_msgs::ImageConstPtr &img,
                                const nav_msgs::OdometryConstPtr &odom)
{
  /* get pose */
  Eigen::Quaterniond body_q = Eigen::Quaterniond(odom->pose.pose.orientation.w,
                                                 odom->pose.pose.orientation.x,
                                                 odom->pose.pose.orientation.y,
                                                 odom->pose.pose.orientation.z);    
  Eigen::Matrix3d body_r_m = body_q.toRotationMatrix();   
  Eigen::Matrix4d body2world;
  body2world.block<3, 3>(0, 0) = body_r_m;
  body2world(0, 3) = odom->pose.pose.position.x;
  body2world(1, 3) = odom->pose.pose.position.y;
  body2world(2, 3) = odom->pose.pose.position.z;
  body2world(3, 3) = 1.0;
  
  Eigen::Matrix4d cam_T = body2world * md_.cam2body_;
  md_.camera_pos_(0) = cam_T(0, 3);
  md_.camera_pos_(1) = cam_T(1, 3);
  md_.camera_pos_(2) = cam_T(2, 3);
  md_.camera_q_ = Eigen::Quaterniond(cam_T.block<3, 3>(0, 0));

  /* get depth image */
  cv_bridge::CvImagePtr cv_ptr;
  cv_ptr = cv_bridge::toCvCopy(img, img->encoding);
  if (img->encoding == sensor_msgs::image_encodings::TYPE_32FC1)
  {
    (cv_ptr->image).convertTo(cv_ptr->image, CV_16UC1, mp_.k_depth_scaling_factor_);
  }
  cv_ptr->image.copyTo(md_.depth_image_);

  md_.occ_need_update_ = true;
}

// GridMap
