#include <Eigen/Eigen>
#include <cmath>
#include <iostream>
#include <plan_env/raycast.h>

int signum(int x) {
  return x == 0 ? 0 : x < 0 ? -1 : 1;
}//符号函数

double mod(double value, double modulus) {
  return fmod(fmod(value, modulus) + modulus, modulus);
}//取余数
double intbound(double s, double ds) {
  // Find the smallest positive t such that s+t*ds is an integer.
  //求最小的t，使s+t*ds为整数 
  if (ds < 0) {
    return intbound(-s, -ds);
  } else {
    s = mod(s, 1);
    // problem is now s+t*ds = 1
    return (1 - s) / ds;
  }
}

void Raycast(const Eigen::Vector3d& start, const Eigen::Vector3d& end, const Eigen::Vector3d& min,
             const Eigen::Vector3d& max, int& output_points_cnt, Eigen::Vector3d* output) {
  //    std::cout << start << ' ' << end << std::endl;
  // From "A Fast Voxel Traversal Algorithm for Ray Tracing"
  // by John Amanatides and Andrew Woo, 1987
  // <http://www.cse.yorku.ca/~amana/research/grid.pdf>
  // <http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.42.3443>
  // Extensions to the described algorithm:
  //   鈥?Imposed a distance limit.
  //   鈥?The face passed through to reach the current cube is provided to
  //     the callback.

  // The foundation of this algorithm is a parameterized representation of
  // the provided ray,
  //                    origin + t * direction,
  // except that t is not actually stored; rather, at any given point in the
  // traversal, we keep track of the *greater* t values which we would have
  // if we took a step sufficient to cross a cube boundary along that axis
  // (i.e. change the integer part of the coordinate) in the variables
  // tMaxX, tMaxY, and tMaxZ.

  // Cube containing origin point. 包含原点的立方体
  int x = (int)std::floor(start.x());
  int y = (int)std::floor(start.y());
  int z = (int)std::floor(start.z());//向下取整找起点体素角标 
  int endX = (int)std::floor(end.x());
  int endY = (int)std::floor(end.y());
  int endZ = (int)std::floor(end.z());//向下取整找终点体素角标  
  Eigen::Vector3d direction = (end - start);//方向向量 
  double maxDist = direction.squaredNorm();//最长距离：起点和终点之间的距离 

  // Break out direction vector.
  double dx = endX - x;
  double dy = endY - y;
  double dz = endZ - z;//起点终点体素角标相减 

  // Direction to increment x,y,z when stepping.
  int stepX = (int)signum((int)dx);
  int stepY = (int)signum((int)dy);
  int stepZ = (int)signum((int)dz);//步进时候增加的方向（1代表正方向 -1代表负方向 0代表在这个维度上不移动 

  // See description above. The initial values depend on the fractional
  // part of the origin.
  double tMaxX = intbound(start.x(), dx);
  double tMaxY = intbound(start.y(), dy);
  double tMaxZ = intbound(start.z(), dz);//以dx dy dz 为速度从起点开始，第一次和下一个x y z边界相交的时间 

  // The change in t when taking a step (always positive).
  double tDeltaX = ((double)stepX) / dx;
  double tDeltaY = ((double)stepY) / dy;
  double tDeltaZ = ((double)stepZ) / dz;//以dx dy dz 为速度，每前进一个体素所需要的的时间 
   

  // Avoids an infinite loop.避免无限循环 
  if (stepX == 0 && stepY == 0 && stepZ == 0) return;//如果步进方向全为0直接返回（起点终点在一个位置） 

  double dist = 0;
  while (true) {//xyz为体素角标 
    if (x >= min.x() && x < max.x() && y >= min.y() && y < max.y() && z >= min.z() && z < max.z()) {
      //如果体素在规定的最小最大体素范围内 
	  output[output_points_cnt](0) = x;
      output[output_points_cnt](1) = y;
      output[output_points_cnt](2) = z;//先存放起点所在的体素角标 

      output_points_cnt++;//指向下一个存储位置 
      dist = sqrt((x - start(0)) * (x - start(0)) + (y - start(1)) * (y - start(1)) +
                  (z - start(2)) * (z - start(2)));//到起点的距离 

      if (dist > maxDist) return;//如果到起点的距离超过了最大距离直接返回 

      /*            if (output_points_cnt > 1500) {
                      std::cerr << "Error, too many racyast voxels." <<
         std::endl;
                      throw std::out_of_range("Too many raycast voxels");
                  }*/
    }

    if (x == endX && y == endY && z == endZ) break;//如果到了终点所在的体素，停止 

    // tMaxX stores the t-value at which we cross a cube boundary along the
    // X axis, and similarly for Y and Z. Therefore, choosing the least tMax
    // chooses the closest cube boundary. Only the first case of the four
    // has been commented in detail.
    if (tMaxX < tMaxY) {//到x的时间比到y的时间短 
      if (tMaxX < tMaxZ) {//到x的时间比到z的时间短 
        // Update which cube we are now in.
        x += stepX;//x先到边界 
        // Adjust tMaxX to the next X-oriented boundary crossing.
        tMaxX += tDeltaX;//时间更新 
      } else {
        z += stepZ;//否则z先到边界 
        tMaxZ += tDeltaZ;//时间更新 
      }
    } else {
      if (tMaxY < tMaxZ) {
        y += stepY;
        tMaxY += tDeltaY;
      } else {
        z += stepZ;
        tMaxZ += tDeltaZ;//这几个同理，说白了就是看先到哪个边界，更新体素 
      }
    }
  }
}

void Raycast(const Eigen::Vector3d& start, const Eigen::Vector3d& end, const Eigen::Vector3d& min,
             const Eigen::Vector3d& max, std::vector<Eigen::Vector3d>* output) {
  //    std::cout << start << ' ' << end << std::endl;
  // From "A Fast Voxel Traversal Algorithm for Ray Tracing"
  // by John Amanatides and Andrew Woo, 1987
  // <http://www.cse.yorku.ca/~amana/research/grid.pdf>
  // <http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.42.3443>
  // Extensions to the described algorithm:
  //   鈥?Imposed a distance limit.
  //   鈥?The face passed through to reach the current cube is provided to
  //     the callback.

  // The foundation of this algorithm is a parameterized representation of
  // the provided ray,
  //                    origin + t * direction,
  // except that t is not actually stored; rather, at any given point in the
  // traversal, we keep track of the *greater* t values which we would have
  // if we took a step sufficient to cross a cube boundary along that axis
  // (i.e. change the integer part of the coordinate) in the variables
  // tMaxX, tMaxY, and tMaxZ.

  // Cube containing origin point.
  int x = (int)std::floor(start.x());
  int y = (int)std::floor(start.y());
  int z = (int)std::floor(start.z());
  int endX = (int)std::floor(end.x());
  int endY = (int)std::floor(end.y());
  int endZ = (int)std::floor(end.z());
  Eigen::Vector3d direction = (end - start);
  double maxDist = direction.squaredNorm();

  // Break out direction vector.
  double dx = endX - x;
  double dy = endY - y;
  double dz = endZ - z;

  // Direction to increment x,y,z when stepping.
  int stepX = (int)signum((int)dx);
  int stepY = (int)signum((int)dy);
  int stepZ = (int)signum((int)dz);

  // See description above. The initial values depend on the fractional
  // part of the origin.
  double tMaxX = intbound(start.x(), dx);
  double tMaxY = intbound(start.y(), dy);
  double tMaxZ = intbound(start.z(), dz);

  // The change in t when taking a step (always positive).
  double tDeltaX = ((double)stepX) / dx;
  double tDeltaY = ((double)stepY) / dy;
  double tDeltaZ = ((double)stepZ) / dz;

  output->clear();

  // Avoids an infinite loop.
  if (stepX == 0 && stepY == 0 && stepZ == 0) return;

  double dist = 0;
  while (true) {
    if (x >= min.x() && x < max.x() && y >= min.y() && y < max.y() && z >= min.z() && z < max.z()) {
      output->push_back(Eigen::Vector3d(x, y, z));

      dist = (Eigen::Vector3d(x, y, z) - start).squaredNorm();

      if (dist > maxDist) return;

      if (output->size() > 1500) {
        std::cerr << "Error, too many racyast voxels." << std::endl;
        throw std::out_of_range("Too many raycast voxels");//太多体素 
      }
    }

    if (x == endX && y == endY && z == endZ) break;

    // tMaxX stores the t-value at which we cross a cube boundary along the
    // X axis, and similarly for Y and Z. Therefore, choosing the least tMax
    // chooses the closest cube boundary. Only the first case of the four
    // has been commented in detail.
    if (tMaxX < tMaxY) {
      if (tMaxX < tMaxZ) {
        // Update which cube we are now in.
        x += stepX;
        // Adjust tMaxX to the next X-oriented boundary crossing.
        tMaxX += tDeltaX;
      } else {
        z += stepZ;
        tMaxZ += tDeltaZ;
      }
    } else {
      if (tMaxY < tMaxZ) {
        y += stepY;
        tMaxY += tDeltaY;
      } else {
        z += stepZ;
        tMaxZ += tDeltaZ;
      }
    }
  }
}
//设置输入，给定起点和终点的位置 返回是否可以遍历（如果起点和终点不是一个点就可以） 
bool RayCaster::setInput(const Eigen::Vector3d& start,
                         const Eigen::Vector3d& end /* , const Eigen::Vector3d& min,
                         const Eigen::Vector3d& max */) {
  start_ = start;
  end_ = end;
  // max_ = max;
  // min_ = min; 
  x_ = (int)std::floor(start_.x());
  y_ = (int)std::floor(start_.y());
  z_ = (int)std::floor(start_.z());
  endX_ = (int)std::floor(end_.x());
  endY_ = (int)std::floor(end_.y());
  endZ_ = (int)std::floor(end_.z());
  direction_ = (end_ - start_);
  maxDist_ = direction_.squaredNorm();

  // Break out direction vector.
  dx_ = endX_ - x_;
  dy_ = endY_ - y_;
  dz_ = endZ_ - z_;

  // Direction to increment x,y,z when stepping.
  stepX_ = (int)signum((int)dx_);
  stepY_ = (int)signum((int)dy_);
  stepZ_ = (int)signum((int)dz_);

  // See description above. The initial values depend on the fractional
  // part of the origin.
  tMaxX_ = intbound(start_.x(), dx_);
  tMaxY_ = intbound(start_.y(), dy_);
  tMaxZ_ = intbound(start_.z(), dz_);

  // The change in t when taking a step (always positive).
  tDeltaX_ = ((double)stepX_) / dx_;
  tDeltaY_ = ((double)stepY_) / dy_;
  tDeltaZ_ = ((double)stepZ_) / dz_;

  dist_ = 0;

  step_num_ = 0;

  // Avoids an infinite loop.
  if (stepX_ == 0 && stepY_ == 0 && stepZ_ == 0)
    return false;
  else
    return true;
}
//一步更新 
bool RayCaster::step(Eigen::Vector3d& ray_pt) {
  // if (x_ >= min_.x() && x_ < max_.x() && y_ >= min_.y() && y_ < max_.y() &&
  // z_ >= min_.z() && z_ <
  // max_.z())
  ray_pt = Eigen::Vector3d(x_, y_, z_);//当前的体素角标 

  // step_num_++;

  // dist_ = (Eigen::Vector3d(x_, y_, z_) - start_).squaredNorm();

  if (x_ == endX_ && y_ == endY_ && z_ == endZ_) {
    return false;//如果到了终点体素 
  }

  // if (dist_ > maxDist_)
  // {
  //   return false;
  // }

  // tMaxX stores the t-value at which we cross a cube boundary along the
  // X axis, and similarly for Y and Z. Therefore, choosing the least tMax
  // chooses the closest cube boundary. Only the first case of the four
  // has been commented in detail.
  if (tMaxX_ < tMaxY_) {
    if (tMaxX_ < tMaxZ_) {
      // Update which cube we are now in.
      x_ += stepX_;
      // Adjust tMaxX to the next X-oriented boundary crossing.
      tMaxX_ += tDeltaX_;
    } else {
      z_ += stepZ_;
      tMaxZ_ += tDeltaZ_;
    }
  } else {
    if (tMaxY_ < tMaxZ_) {
      y_ += stepY_;
      tMaxY_ += tDeltaY_;
    } else {
      z_ += stepZ_;
      tMaxZ_ += tDeltaZ_;
    }
  }

  return true;
}