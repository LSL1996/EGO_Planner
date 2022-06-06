#include "bspline_opt/uniform_bspline.h"
#include <ros/ros.h>

namespace ego_planner
{

  UniformBspline::UniformBspline(const Eigen::MatrixXd &points, const int &order,
                                 const double &interval)//得到均匀b样条曲线的knot span
  {
    setUniformBspline(points, order, interval);
  }

  UniformBspline::~UniformBspline() {}

  void UniformBspline::setUniformBspline(const Eigen::MatrixXd &points, const int &order,
                                         const double &interval)//输入控制点，样条曲线阶数，时间间隔span 计算knot span
  {
    control_points_ = points;//控制点
    p_ = order;//次数(阶数为p+1)
    interval_ = interval;//时间间隔

    n_ = points.cols() - 1;//控制点最大坐标（总个数-1，因为从0开始）
    m_ = n_ + p_ + 1;//knot span的最大坐标

    u_ = Eigen::VectorXd::Zero(m_ + 1);//得到knot_span
  //对于b样条曲线而言，有效的knot span为[up,un+1]有定义.所以分为三段定义
    for (int i = 0; i <= m_; ++i)
    {

      if (i <= p_)
      {
        u_(i) = double(-p_ + i) * interval_;
      }
      else if (i > p_ && i <= m_ - p_)
      {
        u_(i) = u_(i - 1) + interval_;
      }
      else if (i > m_ - p_)
      {
        u_(i) = u_(i - 1) + interval_;
      }
    }
  }

  void UniformBspline::setKnot(const Eigen::VectorXd &knot) { this->u_ = knot; }//得到knots vector

  Eigen::VectorXd UniformBspline::getKnot() { return this->u_; }//取出knots vector

  bool UniformBspline::getTimeSpan(double &um, double &um_p)//得到开始和结束的时间，返回是否时间间隔计算正确
  {
    if (p_ > u_.rows() || m_ - p_ > u_.rows())//发生错误
      return false;

    um = u_(p_);
    um_p = u_(m_ - p_);

    return true;
  }

  Eigen::MatrixXd UniformBspline::getControlPoint() { return control_points_; }//取出控制点

  Eigen::VectorXd UniformBspline::evaluateDeBoor(const double &u)//deboor递推公式计算p_+1阶基函数 返回这个时间b样条曲线所对应的位置
  //输入一个时间，得到它所在的knot span时间段。然后计算会对这段时间有所影响的控制点(一共p_+1个)所对应的基函数（p+1阶次）
  //
  {

    double ub = min(max(u_(p_), u), u_(m_ - p_));//把时间限制在有效定义的时间区间内

    // determine which [ui,ui+1] lay in
    //从有效时间段的起点开始遍历，看其在哪个时间段内, 在[uk,uk+1]内
    int k = p_;
    while (true)
    {
      if (u_(k + 1) >= ub)
        break;
      ++k;
    }

    /* deBoor's alg  deboor递推公式 */
    vector<Eigen::VectorXd> d;
    for (int i = 0; i <= p_; ++i)
    {
      d.push_back(control_points_.col(k - p_ + i));
      // cout << d[i].transpose() << endl;
    }
    //基函数deboor递推
    for (int r = 1; r <= p_; ++r)
    {
      for (int i = p_; i >= r; --i)
      {
        double alpha = (ub - u_[i + k - p_]) / (u_[i + 1 + k - r] - u_[i + k - p_]);
        // cout << "alpha: " << alpha << endl;
        d[i] = (1 - alpha) * d[i - 1] + alpha * d[i];
      }
    }

    return d[p_];//返回这个时间b样条曲线所对应的位置
  }

  // Eigen::VectorXd UniformBspline::evaluateDeBoorT(const double& t) {
  //   return evaluateDeBoor(t + u_(p_));
  // }

  Eigen::MatrixXd UniformBspline::getDerivativeControlPoints()//得到求导后的控制点
  {
    // The derivative of a b-spline is also a b-spline, its order become p_-1  求导后阶数也少了1
    // control point Qi = p_*(Pi+1-Pi)/(ui+p_+1-ui+1)
    Eigen::MatrixXd ctp(control_points_.rows(), control_points_.cols() - 1);//求导后的控制点个数会少一个点
    for (int i = 0; i < ctp.cols(); ++i)
    {
      ctp.col(i) =
          p_ * (control_points_.col(i + 1) - control_points_.col(i)) / (u_(i + p_ + 1) - u_(i + 1));
    }//计算具体值
    return ctp;//返回控制点导数矩阵
  }

  UniformBspline UniformBspline::getDerivative()
  {
    Eigen::MatrixXd ctp = getDerivativeControlPoints();//得到求导后的控制点
    UniformBspline derivative(ctp, p_ - 1, interval_);//求导后的控制点class

    /* cut the first and last knot */
    Eigen::VectorXd knot(u_.rows() - 2);//去除开头和结尾的knot
    knot = u_.segment(1, u_.rows() - 2);
    derivative.setKnot(knot);

    return derivative;//返回求导后的控制点class
  }

  double UniformBspline::getInterval() { return interval_; }//得到时间间隔delta

  void UniformBspline::setPhysicalLimits(const double &vel, const double &acc, const double &tolerance)//设置物理约束
  {
    limit_vel_ = vel;
    limit_acc_ = acc;
    limit_ratio_ = 1.1;
    feasibility_tolerance_ = tolerance;//放宽限制的大小
  }

  bool UniformBspline::checkFeasibility(double &ratio, bool show)//检查是否可行
  {
    bool fea = true;//先设置为true

    Eigen::MatrixXd P = control_points_;//控制点矩阵
    int dimension = control_points_.rows();//得到维度

    /* check vel feasibility and insert points 检查速度可行性*/
    double max_vel = -1.0;//最大速度 先设置为负数
    double enlarged_vel_lim = limit_vel_ * (1.0 + feasibility_tolerance_) + 1e-4;//放宽限制后的速度约束
    for (int i = 0; i < P.cols() - 1; ++i)//循环控制点
    {
      Eigen::VectorXd vel = p_ * (P.col(i + 1) - P.col(i)) / (u_(i + p_ + 1) - u_(i + 1));//计算速度控制点

      if (fabs(vel(0)) > enlarged_vel_lim || fabs(vel(1)) > enlarged_vel_lim ||
          fabs(vel(2)) > enlarged_vel_lim)//速度控制点的任何一个维度值超出了限制
      {

        if (show)//这个标志判断是不是要输出显示
          cout << "[Check]: Infeasible vel " << i << " :" << vel.transpose() << endl;
        fea = false;//超出了限制，设置为false

        for (int j = 0; j < dimension; ++j)//循环维度
        {
          max_vel = max(max_vel, fabs(vel(j)));//得到所有速度控制点所有维度的最大速度
        }
      }
    }

    /* acc feasibility检查加速度可行性 */
    //求解过程和上面速度类似
    double max_acc = -1.0;
    double enlarged_acc_lim = limit_acc_ * (1.0 + feasibility_tolerance_) + 1e-4;
    for (int i = 0; i < P.cols() - 2; ++i)
    {

      Eigen::VectorXd acc = p_ * (p_ - 1) *
                            ((P.col(i + 2) - P.col(i + 1)) / (u_(i + p_ + 2) - u_(i + 2)) -
                             (P.col(i + 1) - P.col(i)) / (u_(i + p_ + 1) - u_(i + 1))) /
                            (u_(i + p_ + 1) - u_(i + 2));

      if (fabs(acc(0)) > enlarged_acc_lim || fabs(acc(1)) > enlarged_acc_lim ||
          fabs(acc(2)) > enlarged_acc_lim)
      {

        if (show)
          cout << "[Check]: Infeasible acc " << i << " :" << acc.transpose() << endl;
        fea = false;//超出了限制，设置为false

        for (int j = 0; j < dimension; ++j)
        {
          max_acc = max(max_acc, fabs(acc(j)));//得到所有加速度控制点所有维度的最大加速度
        }
      }
    }

    ratio = max(max_vel / limit_vel_, sqrt(fabs(max_acc) / limit_acc_));//公式14

    return fea;//返回是否超出限制
  }

  void UniformBspline::lengthenTime(const double &ratio)//输入放大系数，延长时间
  {
    //问题：为什么前面和后面的点不放大时间？？？？？
    int num1 = 5;//起始的knotspan角标
    int num2 = getKnot().rows() - 1 - 5;//结束的knotspan角标

    double delta_t = (ratio - 1.0) * (u_(num2) - u_(num1));//这一段要考虑的时间 所放大的时间长度
    double t_inc = delta_t / double(num2 - num1);//每一段时间interval所放大的时间长度
    for (int i = num1 + 1; i <= num2; ++i)//循环要考虑的knotspan
      u_(i) += double(i - num1) * t_inc;
    for (int i = num2 + 1; i < u_.rows(); ++i)
      u_(i) += delta_t;//得到时间放大后的knotspan
  }

  // void UniformBspline::recomputeInit() {}

  void UniformBspline::parameterizeToBspline(const double &ts, const vector<Eigen::Vector3d> &point_set,
                                             const vector<Eigen::Vector3d> &start_end_derivative,
                                             Eigen::MatrixXd &ctrl_pts)
  //参考 https://blog.csdn.net/qq_37866732/article/details/118633062
  {
    if (ts <= 0)//时间间隔是负数
    {
      cout << "[B-spline]:time step error." << endl;
      return;
    }

    if (point_set.size() <= 3)//插值点个数太少
    {
      cout << "[B-spline]:point set have only " << point_set.size() << " points." << endl;
      return;
    }

    if (start_end_derivative.size() != 4)//导数错误 四个导数约束(两个速度两个加速度）
    {
      cout << "[B-spline]:derivatives error." << endl;
    }

    int K = point_set.size();//插值点的个数

    // write A
    Eigen::Vector3d prow(3), vrow(3), arow(3);
    prow << 1, 4, 1;//位置约束行
    vrow << -1, 0, 1;//速度约束行
    arow << 1, -2, 1;//加速度约束行

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(K + 4, K + 2);//行数为插值点的个数，列数为控制点的个数

    for (int i = 0; i < K; ++i)//书写计算矩阵A
    //.block 函数其实就是分块矩阵  参考:https://blog.csdn.net/zhangqian_shai/article/details/101756145
      A.block(i, i, 1, 3) = (1 / 6.0) * prow.transpose();//位置的约束 一共k个插值点位置

    A.block(K, 0, 1, 3) = (1 / 2.0 / ts) * vrow.transpose();//速度的约束 开头结尾两个速度约束
    A.block(K + 1, K - 1, 1, 3) = (1 / 2.0 / ts) * vrow.transpose();

    A.block(K + 2, 0, 1, 3) = (1 / ts / ts) * arow.transpose();//加速度的约束 开头加结尾两个速度约束
    A.block(K + 3, K - 1, 1, 3) = (1 / ts / ts) * arow.transpose();

    //cout << "A" << endl << A << endl << endl;

    // write b 书写向量b
    Eigen::VectorXd bx(K + 4), by(K + 4), bz(K + 4);
    for (int i = 0; i < K; ++i)//插值点xyz三个维度的位置
    {
      bx(i) = point_set[i](0);
      by(i) = point_set[i](1);
      bz(i) = point_set[i](2);
    }

    for (int i = 0; i < 4; ++i)//插值点起点和终点xyz三个维度的速度 加结尾
    {
      bx(K + i) = start_end_derivative[i](0);
      by(K + i) = start_end_derivative[i](1);
      bz(K + i) = start_end_derivative[i](2);
    }

    // solve Ax = b 求解方程
    Eigen::VectorXd px = A.colPivHouseholderQr().solve(bx);
    Eigen::VectorXd py = A.colPivHouseholderQr().solve(by);
    Eigen::VectorXd pz = A.colPivHouseholderQr().solve(bz);

    // convert to control pts 得到控制点的三维坐标矩阵
    ctrl_pts.resize(3, K + 2);
    ctrl_pts.row(0) = px.transpose();
    ctrl_pts.row(1) = py.transpose();
    ctrl_pts.row(2) = pz.transpose();

    // cout << "[B-spline]: parameterization ok." << endl;
  }

  double UniformBspline::getTimeSum()//得到knotspan上总的有效时间
  {
    double tm, tmp;
    if (getTimeSpan(tm, tmp))
      return tmp - tm;
    else
      return -1.0;
  }

  double UniformBspline::getLength(const double &res)//从当前时间res时b样条曲线的位置开始，计算其到b样条曲线结束时位置，b样条曲线的长度
  {
    double length = 0.0;
    double dur = getTimeSum();//knotspan上总的有效时间
    Eigen::VectorXd p_l = evaluateDeBoorT(0.0), p_n;//p_l :b样条曲线起点的位置(t=0)
    for (double t = res; t <= dur + 1e-4; t += res)
    {
      p_n = evaluateDeBoorT(t);//b样条曲线点的位置
      length += (p_n - p_l).norm();//b样条曲线上两点之间的长度
      p_l = p_n;
    }
    return length;
  }

  double UniformBspline::getJerk()//返回jerk控制点所求得的jerk值
  {
    UniformBspline jerk_traj = getDerivative().getDerivative().getDerivative();//求三次导数

    Eigen::VectorXd times = jerk_traj.getKnot();//得到knotspan
    Eigen::MatrixXd ctrl_pts = jerk_traj.getControlPoint();//得到jerk所对应的控制点
    int dimension = ctrl_pts.rows();//维度

    double jerk = 0.0;
    for (int i = 0; i < ctrl_pts.cols(); ++i)//循环控制点
    {
      for (int j = 0; j < dimension; ++j)//循环维度
      {
        jerk += (times(i + 1) - times(i)) * ctrl_pts(j, i) * ctrl_pts(j, i);//为什么要乘t？？？？？？？
      }
    }

    return jerk;
  }

  void UniformBspline::getMeanAndMaxVel(double &mean_v, double &max_v)//得到平均和最大速度
  {
    UniformBspline vel = getDerivative();//速度
    double tm, tmp;
    vel.getTimeSpan(tm, tmp);//起始和终止时间

    double max_vel = -1.0, mean_vel = 0.0;
    int num = 0;
    for (double t = tm; t <= tmp; t += 0.01)
    {
      Eigen::VectorXd vxd = vel.evaluateDeBoor(t);//得到此时b样条曲线上对应的速度
      double vn = vxd.norm();//取模

      mean_vel += vn;//加上现在的速度
      ++num;
      if (vn > max_vel)
      {
        max_vel = vn;//如果当前速度超过了最大速度，更新最大速度值
      }
    }

    mean_vel = mean_vel / double(num);//整段时间的平均速度
    mean_v = mean_vel;
    max_v = max_vel;//整段时间上最大的速度
  }

  void UniformBspline::getMeanAndMaxAcc(double &mean_a, double &max_a)//得到此时b样条曲线上对应的加速度
  {
    UniformBspline acc = getDerivative().getDerivative();
    double tm, tmp;
    acc.getTimeSpan(tm, tmp);

    double max_acc = -1.0, mean_acc = 0.0;
    int num = 0;
    for (double t = tm; t <= tmp; t += 0.01)
    {
      Eigen::VectorXd axd = acc.evaluateDeBoor(t);
      double an = axd.norm();

      mean_acc += an;
      ++num;
      if (an > max_acc)
      {
        max_acc = an;
      }
    }

    mean_acc = mean_acc / double(num);
    mean_a = mean_acc;
    max_a = max_acc;
  }
} // namespace ego_planner
