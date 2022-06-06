#ifndef _UNIFORM_BSPLINE_H_
#define _UNIFORM_BSPLINE_H_

#include <Eigen/Eigen>
#include <algorithm>
#include <iostream>

using namespace std;

namespace ego_planner
{
  // An implementation of non-uniform B-spline with different dimensions
  // It also represents uniform B-spline which is a special case of non-uniform
  class UniformBspline
  {
  private:
    // control points for B-spline with different dimensions.
    // Each row represents one single control point
    // The dimension is determined by column number
    // e.g. B-spline with N points in 3D space -> Nx3 matrix
    Eigen::MatrixXd control_points_;

    int p_, n_, m_;     // p degree, n+1 control points, m = n+p+1
    Eigen::VectorXd u_; // knots vector
    double interval_;   // knot span \delta t

    Eigen::MatrixXd getDerivativeControlPoints();//获得求导后的控制点

    double limit_vel_, limit_acc_, limit_ratio_, feasibility_tolerance_; //物理约束 physical limits and time adjustment ratio

  public:
    UniformBspline() {}
    UniformBspline(const Eigen::MatrixXd &points, const int &order, const double &interval);
    ~UniformBspline();

    Eigen::MatrixXd get_control_points(void) { return control_points_; }//得到控制点

    // initialize as an uniform B-spline
    void setUniformBspline(const Eigen::MatrixXd &points, const int &order, const double &interval);

    // get / set basic bspline info

    void setKnot(const Eigen::VectorXd &knot);//设置knotspan
    Eigen::VectorXd getKnot();//取出knotspan
    Eigen::MatrixXd getControlPoint();//取出控制点
    double getInterval();//取出时间间隔
    bool getTimeSpan(double &um, double &um_p);//得到开始和结束的时间，返回是否时间间隔计算正确

    // compute position / derivative 计算位置和求导

    Eigen::VectorXd evaluateDeBoor(const double &u);   // use u \in [um, u_mp] u代表在整个knot span时间上的值
    //输入一个时间，计算有所影响的控制点(一共p_+1个)所对应的基函数（p+1阶次）在这个时间的值。

    inline Eigen::VectorXd evaluateDeBoorT(const double &t) { return evaluateDeBoor(t + u_(p_)); } // t代表在有效时间区间上的值use t \in [0, duration]
    
    UniformBspline getDerivative();//求导

    // 3D B-spline interpolation of points in point_set, with boundary vel&acc
    // constraints
    //输入k个插值点的位置和开头结尾的速度和加速度，一共k+4个约束。得到k+2个控制点位置
    static void parameterizeToBspline(const double &ts, const vector<Eigen::Vector3d> &point_set,
                                      const vector<Eigen::Vector3d> &start_end_derivative,
                                      Eigen::MatrixXd &ctrl_pts);

    /* check feasibility, adjust time 检查动力可行性，调整时间 */

    void setPhysicalLimits(const double &vel, const double &acc, const double &tolerance);//设置动力学约束
    bool checkFeasibility(double &ratio, bool show = false);//检查是否可行即是否超过了速度加速度限制
    void lengthenTime(const double &ratio);//时间扩展

    /* for performance evaluation 曲线效果评估 */

    double getTimeSum();//总时间
    double getLength(const double &res = 0.01);//总长度
    double getJerk();//得到曲线的jerk值
    void getMeanAndMaxVel(double &mean_v, double &max_v);//得到此时b样条曲线上对应的速度
    void getMeanAndMaxAcc(double &mean_a, double &max_a);//得到此时b样条曲线上对应的加速度

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  };
} // namespace ego_planner
#endif