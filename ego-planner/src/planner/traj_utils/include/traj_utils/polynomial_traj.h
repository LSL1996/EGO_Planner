#ifndef _POLYNOMIAL_TRAJ_H
#define _POLYNOMIAL_TRAJ_H

#include <Eigen/Eigen>
#include <vector>

using std::vector;

class PolynomialTraj
{
private:
  vector<double> times;       //轨迹时间分配  time of each segment

  //每一段三个维度的多项式系数
  vector<vector<double>> cxs; // coefficient of x of each segment, from n-1 -> 0
  vector<vector<double>> cys; // coefficient of y of each segment
  vector<vector<double>> czs; // coefficient of z of each segment

  double time_sum;//总时间 
  int num_seg;//轨迹段数

  /* evaluation  评估*/
  vector<Eigen::Vector3d> traj_vec3d;//轨迹 存放着许多轨迹点
  double length;//长度

public:
  PolynomialTraj(/* args */)
  {
  }
  ~PolynomialTraj()
  {
  }

  void reset()//重置
  {
    times.clear(), cxs.clear(), cys.clear(), czs.clear();
    time_sum = 0.0, num_seg = 0;
  }//时间清零 系数清零

  void addSegment(vector<double> cx, vector<double> cy, vector<double> cz, double t)
  {
    cxs.push_back(cx), cys.push_back(cy), czs.push_back(cz), times.push_back(t);
  }//增加一段，存储三个维度多项式系数和这一段的时间

  void init()//初始化
  {
    num_seg = times.size();
    time_sum = 0.0;
    for (int i = 0; i < times.size(); ++i)
    {
      time_sum += times[i];
    }
  }//得到段数和总时间

  vector<double> getTimes()
  {
    return times;
  }//取得时间

  vector<vector<double>> getCoef(int axis)//取得不同维度多项式系数
  {
    switch (axis)
    {
    case 0:
      return cxs;
    case 1:
      return cys;
    case 2:
      return czs;
    default:
      std::cout << "\033[31mIllegal axis!\033[0m" << std::endl;
    }

    vector<vector<double>> empty;//否则返回空向量
    return empty;
  }

  Eigen::Vector3d evaluate(double t)//给定时间t返回位置
  {
    /* detetrmine segment num 看这个时间在哪个时间段中 times[idx]代表这一段的时间*/
    int idx = 0;
    while (times[idx] + 1e-4 < t)
    {
      t -= times[idx];
      ++idx;
    }

    /* evaluation */
    int order = cxs[idx].size();//这一段轨迹多项式系数的个数
    Eigen::VectorXd cx(order), cy(order), cz(order), tv(order);
    //循环多项式系数
    for (int i = 0; i < order; ++i)
    {
      cx(i) = cxs[idx][i], cy(i) = cys[idx][i], cz(i) = czs[idx][i];//存系数 
      tv(order - 1 - i) = std::pow(t, double(i));//存t的次方 c = [c5, c4, c3,...c0]' 所以要倒过来(order - 1 - i)存放
    }

    Eigen::Vector3d pt;
    pt(0) = tv.dot(cx), pt(1) = tv.dot(cy), pt(2) = tv.dot(cz);//得到时间t时候的位置，三个维度位置分开计算
    return pt;
  }

  Eigen::Vector3d evaluateVel(double t)//给定时间t返回速度
  {
    /* detetrmine segment num */
    int idx = 0;
    while (times[idx] + 1e-4 < t)
    {
      t -= times[idx];
      ++idx;
    }

    /* evaluation */
    int order = cxs[idx].size();
    Eigen::VectorXd vx(order - 1), vy(order - 1), vz(order - 1);

    /* coef of vel */
    for (int i = 0; i < order - 1; ++i)
    {
      vx(i) = double(i + 1) * cxs[idx][order - 2 - i];
      vy(i) = double(i + 1) * cys[idx][order - 2 - i];
      vz(i) = double(i + 1) * czs[idx][order - 2 - i];
    }
    double ts = t;
    Eigen::VectorXd tv(order - 1);
    for (int i = 0; i < order - 1; ++i)
      tv(i) = pow(ts, i);

    Eigen::Vector3d vel;
    vel(0) = tv.dot(vx), vel(1) = tv.dot(vy), vel(2) = tv.dot(vz);
    return vel;
  }

  Eigen::Vector3d evaluateAcc(double t)//给定时间t返回加速度
  {
    /* detetrmine segment num */
    int idx = 0;
    while (times[idx] + 1e-4 < t)
    {
      t -= times[idx];
      ++idx;
    }

    /* evaluation */
    int order = cxs[idx].size();
    Eigen::VectorXd ax(order - 2), ay(order - 2), az(order - 2);

    /* coef of vel */
    for (int i = 0; i < order - 2; ++i)
    {
      ax(i) = double((i + 2) * (i + 1)) * cxs[idx][order - 3 - i];
      ay(i) = double((i + 2) * (i + 1)) * cys[idx][order - 3 - i];
      az(i) = double((i + 2) * (i + 1)) * czs[idx][order - 3 - i];
    }
    double ts = t;
    Eigen::VectorXd tv(order - 2);
    for (int i = 0; i < order - 2; ++i)
      tv(i) = pow(ts, i);

    Eigen::Vector3d acc;
    acc(0) = tv.dot(ax), acc(1) = tv.dot(ay), acc(2) = tv.dot(az);
    return acc;
  }

  /* for evaluating traj, should be called in sequence!!! */
  double getTimeSum()
  {
    return this->time_sum;//取得总时间 
  }

  vector<Eigen::Vector3d> getTraj()//得到整个轨迹的位置
  {
    double eval_t = 0.0;
    traj_vec3d.clear();//轨迹清零
    while (eval_t < time_sum)//循环时间
    {
      Eigen::Vector3d pt = evaluate(eval_t);//得到这个时间的位置
      traj_vec3d.push_back(pt);//存放
      eval_t += 0.01;
    }
    return traj_vec3d;
  }

  double getLength()//得到轨迹总长度
  {
    length = 0.0;//初始长度为0

    Eigen::Vector3d p_l = traj_vec3d[0], p_n;//p_l:起点位置 
    for (int i = 1; i < traj_vec3d.size(); ++i)//循环轨迹点计算每一段的长度
    {
      p_n = traj_vec3d[i];
      length += (p_n - p_l).norm();
      p_l = p_n;
    }
    return length;//轨迹总长度
  }

  double getMeanVel()//平均速度
  {
    double mean_vel = length / time_sum;//轨迹长度/总时间
  }

  double getAccCost()//加速度cost 为什么要这么计算？？？？？
  {
    double cost = 0.0;
    int order = cxs[0].size();//起始段轨迹多项式系数的个数

    for (int s = 0; s < times.size(); ++s)//循环每一段
    {
      Eigen::Vector3d um;
      um(0) = 2 * cxs[s][order - 3], um(1) = 2 * cys[s][order - 3], um(2) = 2 * czs[s][order - 3];
      cost += um.squaredNorm() * times[s];//不太明白含义????????
    }

    return cost;
  }

  double getJerk()//计算三阶导数jerk    ppt327
  {
    double jerk = 0.0;

    /* evaluate jerk */
    for (int s = 0; s < times.size(); ++s)//循环每一段
    {
      Eigen::VectorXd cxv(cxs[s].size()), cyv(cys[s].size()), czv(czs[s].size());
      /* convert coefficient */
      int order = cxs[s].size();//多项式系数的个数
      for (int j = 0; j < order; ++j)//循环每一个系数
      {
        cxv(j) = cxs[s][order - 1 - j], cyv(j) = cys[s][order - 1 - j], czv(j) = czs[s][order - 1 - j];
      }// c = [c5, c4, c3,...c0]' 倒过来存放这一段多项式系数，即变为[c0, c1, c2,...c5]'
      double ts = times[s];//这一段的时间

      /* jerk matrix */
      Eigen::MatrixXd mat_jerk(order, order);
      mat_jerk.setZero();//初始化为0
      for (double i = 3; i < order; i += 1)
        for (double j = 3; j < order; j += 1)
        {
          mat_jerk(i, j) =
              i * (i - 1) * (i - 2) * j * (j - 1) * (j - 2) * pow(ts, i + j - 5) / (i + j - 5);
        }
        //ppt327所对应的M矩阵，只不过从四阶导变为了三阶导
      jerk += (cxv.transpose() * mat_jerk * cxv)(0, 0);//计算jerk
      jerk += (cyv.transpose() * mat_jerk * cyv)(0, 0);
      jerk += (czv.transpose() * mat_jerk * czv)(0, 0);
    }

    return jerk;
  }

  void getMeanAndMaxVel(double &mean_v, double &max_v)//平均速度最大速度
  {
    int num = 0;
    mean_v = 0.0, max_v = -1.0;
    for (int s = 0; s < times.size(); ++s)//循环每一段
    {
      int order = cxs[s].size();
      Eigen::VectorXd vx(order - 1), vy(order - 1), vz(order - 1);

      /* coef of vel */
      for (int i = 0; i < order - 1; ++i)//循环得到速度多项式的系数
      {
        vx(i) = double(i + 1) * cxs[s][order - 2 - i];
        vy(i) = double(i + 1) * cys[s][order - 2 - i];
        vz(i) = double(i + 1) * czs[s][order - 2 - i];
      }
      double ts = times[s];//这一段的时间

      double eval_t = 0.0;
      while (eval_t < ts)//循环这段轨迹，以时间t为参数
      {
        Eigen::VectorXd tv(order - 1);
        for (int i = 0; i < order - 1; ++i)//循环速度多项式的系数
          tv(i) = pow(ts, i);//存放t的次方
        Eigen::Vector3d vel;
        vel(0) = tv.dot(vx), vel(1) = tv.dot(vy), vel(2) = tv.dot(vz);//得到当前的速度
        double vn = vel.norm();//速度模值
        mean_v += vn;
        if (vn > max_v)
          max_v = vn;//寻找最大速度
        ++num;

        eval_t += 0.01;
      }
    }

    mean_v = mean_v / double(num);//这一段的平均速度（用每一个位置的速度和/点个数）
  }

  void getMeanAndMaxAcc(double &mean_a, double &max_a)//平均加速度最大加速度
  //计算同速度
  {
    int num = 0;
    mean_a = 0.0, max_a = -1.0;
    for (int s = 0; s < times.size(); ++s)
    {
      int order = cxs[s].size();
      Eigen::VectorXd ax(order - 2), ay(order - 2), az(order - 2);

      /* coef of acc */
      for (int i = 0; i < order - 2; ++i)
      {
        ax(i) = double((i + 2) * (i + 1)) * cxs[s][order - 3 - i];
        ay(i) = double((i + 2) * (i + 1)) * cys[s][order - 3 - i];
        az(i) = double((i + 2) * (i + 1)) * czs[s][order - 3 - i];
      }
      double ts = times[s];

      double eval_t = 0.0;
      while (eval_t < ts)
      {
        Eigen::VectorXd tv(order - 2);
        for (int i = 0; i < order - 2; ++i)
          tv(i) = pow(ts, i);
        Eigen::Vector3d acc;
        acc(0) = tv.dot(ax), acc(1) = tv.dot(ay), acc(2) = tv.dot(az);
        double an = acc.norm();
        mean_a += an;
        if (an > max_a)
          max_a = an;
        ++num;

        eval_t += 0.01;
      }
    }

    mean_a = mean_a / double(num);
  }

  static PolynomialTraj minSnapTraj(const Eigen::MatrixXd &Pos, const Eigen::Vector3d &start_vel,
                                    const Eigen::Vector3d &end_vel, const Eigen::Vector3d &start_acc,
                                    const Eigen::Vector3d &end_acc, const Eigen::VectorXd &Time);

  static PolynomialTraj one_segment_traj_gen(const Eigen::Vector3d &start_pt, const Eigen::Vector3d &start_vel, const Eigen::Vector3d &start_acc,
                                             const Eigen::Vector3d &end_pt, const Eigen::Vector3d &end_vel, const Eigen::Vector3d &end_acc,
                                             double t);
};

#endif