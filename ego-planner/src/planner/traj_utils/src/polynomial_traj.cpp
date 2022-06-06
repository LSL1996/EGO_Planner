#include <iostream>
#include <traj_utils/polynomial_traj.h>

PolynomialTraj PolynomialTraj::minSnapTraj(const Eigen::MatrixXd &Pos, const Eigen::Vector3d &start_vel,
                                           const Eigen::Vector3d &end_vel, const Eigen::Vector3d &start_acc,
                                           const Eigen::Vector3d &end_acc, const Eigen::VectorXd &Time)//minimumsnap轨迹
//输入位置、开始速度加速度、结束速度加速度和时间向量
{
  int seg_num = Time.size();//曲线段数
  Eigen::MatrixXd poly_coeff(seg_num, 3 * 6);//存放着所有参数，每段6个系数，三个维度
  Eigen::VectorXd Px(6 * seg_num), Py(6 * seg_num), Pz(6 * seg_num);//P代表每一段的多项式系数

  int num_f, num_p; // 固定变量和自由变量个数 number of fixed and free variables
  int num_d;        // 所有段的导数个数（每个点三个导数，每一段六个）number of all segments' derivatives

  const static auto Factorial = [](int x) {
    int fac = 1;
    for (int i = x; i > 0; i--)
      fac = fac * i;
    return fac;
  };//阶乘函数

  /* ---------- end point derivative ---------- */
  Eigen::VectorXd Dx = Eigen::VectorXd::Zero(seg_num * 6);
  Eigen::VectorXd Dy = Eigen::VectorXd::Zero(seg_num * 6);
  Eigen::VectorXd Dz = Eigen::VectorXd::Zero(seg_num * 6);
//每一段有6个未知的变量，起点和终点各3个 
//每一段存储顺序为：起点位置 终点位置 起点速度 终点速度 起点加速度 终点加速度
  for (int k = 0; k < seg_num; k++)//循环每一段
  {
    /* position to derivative */
    //每一段点的位置
    Dx(k * 6) = Pos(0, k);
    Dx(k * 6 + 1) = Pos(0, k + 1);
    Dy(k * 6) = Pos(1, k);
    Dy(k * 6 + 1) = Pos(1, k + 1);
    Dz(k * 6) = Pos(2, k);
    Dz(k * 6 + 1) = Pos(2, k + 1);
    //起点速度加速度 
    if (k == 0)
    {
      Dx(k * 6 + 2) = start_vel(0);
      Dy(k * 6 + 2) = start_vel(1);
      Dz(k * 6 + 2) = start_vel(2);

      Dx(k * 6 + 4) = start_acc(0);
      Dy(k * 6 + 4) = start_acc(1);
      Dz(k * 6 + 4) = start_acc(2);
    }
    //终点速度加速度
    else if (k == seg_num - 1)
    {
      Dx(k * 6 + 3) = end_vel(0);
      Dy(k * 6 + 3) = end_vel(1);
      Dz(k * 6 + 3) = end_vel(2);

      Dx(k * 6 + 5) = end_acc(0);
      Dy(k * 6 + 5) = end_acc(1);
      Dz(k * 6 + 5) = end_acc(2);
    }
  }

  /* ---------- Mapping Matrix A    A矩阵的计算 对应ppt里的M矩阵，只是换了行 ---------- */
  //因为每一段位置及其导数存储顺序为：起点位置 终点位置 起点速度 终点速度 起点加速度 终点加速度
  Eigen::MatrixXd Ab;
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(seg_num * 6, seg_num * 6);

  for (int k = 0; k < seg_num; k++)//循环每一段
  {
    Ab = Eigen::MatrixXd::Zero(6, 6);//每一段的小矩阵
    for (int i = 0; i < 3; i++)
    {
      Ab(2 * i, i) = Factorial(i);
      for (int j = i; j < 6; j++)
        Ab(2 * i + 1, j) = Factorial(j) / Factorial(j - i) * pow(Time(k), j - i);
    }
    A.block(k * 6, k * 6, 6, 6) = Ab;
  }

  /* ---------- Produce Selection Matrix C' C矩阵的构建---------- */
  Eigen::MatrixXd Ct, C;

  num_f = 2 * seg_num + 4; // 3 + 3 + (seg_num - 1) * 2 = 2m + 4固定变量个数 中间段的位置，(seg_num - 1) * 2代表每一段两个位置（起点终点）
  num_p = 2 * seg_num - 2; //(seg_num - 1) * 2 = 2m - 2自由变量个数(每一个点的速度和加速度是自由的，中间点的速度加速度是未知的，因为上一段终点和下一段起点速度加速度相同所以在这里按照一个自由变量来处理)
  //可以按照ppt第341页标出的图去理解

  num_d = 6 * seg_num;//所有段的变量数目
  Ct = Eigen::MatrixXd::Zero(num_d, num_f + num_p);//Ct矩阵的计算
  //行代表等式左侧所有变量，列代表右侧的固定变量和自由变量。
//左侧的存储顺序（对于每一段）为：起点位置 终点位置 起点速度 终点速度 起点加速度 终点加速度
//右侧的存储顺序：固定变量：起始点位置 速度 加速度 下一个点的位置 下一个点的位置…结束点的位置速度加速度
  
  //第一段存储
  Ct(0, 0) = 1;//起点位置 
  Ct(2, 1) = 1;//起点速度
  Ct(4, 2) = 1; //起点加速度
  Ct(1, 3) = 1;//第一段结束点的位置
  Ct(3, 2 * seg_num + 4) = 1;//第一段结束点的速度(列数靠后是因为已经在自由变量所对应的列了)
  Ct(5, 2 * seg_num + 5) = 1;//第一段结束点的加速度(列数靠后是因为已经在自由变量所对应的列了)

//最后一段存储
  Ct(6 * (seg_num - 1) + 0, 2 * seg_num + 0) = 1;//最后一段起点的位置
  Ct(6 * (seg_num - 1) + 1, 2 * seg_num + 1) = 1; // 最后一段终点的位置
  Ct(6 * (seg_num - 1) + 2, 4 * seg_num + 0) = 1;//最后一段起点的速度(列数靠后是因为已经在自由变量)
  Ct(6 * (seg_num - 1) + 3, 2 * seg_num + 2) = 1; // 最后一段终点的速度(列数靠前是因为在固定变量)
  Ct(6 * (seg_num - 1) + 4, 4 * seg_num + 1) = 1;//最后一段起点的加速度
  Ct(6 * (seg_num - 1) + 5, 2 * seg_num + 3) = 1; // 最后一段终点的加速度
//pad上画个图来理解


  for (int j = 2; j < seg_num; j++)//中间段
  {
    //每一段的位置
    Ct(6 * (j - 1) + 0, 2 + 2 * (j - 1) + 0) = 1;//2代表第一段起点的速度加速度，j-1代表在这段之前有几段，*2是因为每段有两个位置固定量  +0代表这段起点（+1就是终点）
    Ct(6 * (j - 1) + 1, 2 + 2 * (j - 1) + 1) = 1;
    //每一段的速度 
    Ct(6 * (j - 1) + 2, 2 * seg_num + 4 + 2 * (j - 2) + 0) = 1;// 2 * seg_num + 4代表固定变量个数  后面代表在自由变量中的位置
    Ct(6 * (j - 1) + 3, 2 * seg_num + 4 + 2 * (j - 1) + 0) = 1;
    //每一段的加速度
    Ct(6 * (j - 1) + 4, 2 * seg_num + 4 + 2 * (j - 2) + 1) = 1;
    Ct(6 * (j - 1) + 5, 2 * seg_num + 4 + 2 * (j - 1) + 1) = 1;
  }

  C = Ct.transpose();//转置

//从各个段的导数转换到固定变量和自由变量的形式
  Eigen::VectorXd Dx1 = C * Dx;
  Eigen::VectorXd Dy1 = C * Dy;
  Eigen::VectorXd Dz1 = C * Dz;

  /* ---------- minimum snap matrix  Q矩阵构建 ppt327---------- */
  Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(seg_num * 6, seg_num * 6);

  for (int k = 0; k < seg_num; k++)
  {
    for (int i = 3; i < 6; i++)
    {
      for (int j = 3; j < 6; j++)
      {
        Q(k * 6 + i, k * 6 + j) =
            i * (i - 1) * (i - 2) * j * (j - 1) * (j - 2) / (i + j - 5) * pow(Time(k), (i + j - 5));
      }
    }
  }

  /* ---------- R matrix ---------- */
  Eigen::MatrixXd R = C * A.transpose().inverse() * Q * A.inverse() * Ct;//R矩阵

  Eigen::VectorXd Dxf(2 * seg_num + 4), Dyf(2 * seg_num + 4), Dzf(2 * seg_num + 4);//固定变量

//取出固定变量
  Dxf = Dx1.segment(0, 2 * seg_num + 4);
  Dyf = Dy1.segment(0, 2 * seg_num + 4);
  Dzf = Dz1.segment(0, 2 * seg_num + 4);

//分块矩阵
  Eigen::MatrixXd Rff(2 * seg_num + 4, 2 * seg_num + 4);
  Eigen::MatrixXd Rfp(2 * seg_num + 4, 2 * seg_num - 2);
  Eigen::MatrixXd Rpf(2 * seg_num - 2, 2 * seg_num + 4);
  Eigen::MatrixXd Rpp(2 * seg_num - 2, 2 * seg_num - 2);

  Rff = R.block(0, 0, 2 * seg_num + 4, 2 * seg_num + 4);
  Rfp = R.block(0, 2 * seg_num + 4, 2 * seg_num + 4, 2 * seg_num - 2);
  Rpf = R.block(2 * seg_num + 4, 0, 2 * seg_num - 2, 2 * seg_num + 4);
  Rpp = R.block(2 * seg_num + 4, 2 * seg_num + 4, 2 * seg_num - 2, 2 * seg_num - 2);

  /* ---------- close form solution 闭式求解自由变量最优值---------- */

  Eigen::VectorXd Dxp(2 * seg_num - 2), Dyp(2 * seg_num - 2), Dzp(2 * seg_num - 2);
  Dxp = -(Rpp.inverse() * Rfp.transpose()) * Dxf;
  Dyp = -(Rpp.inverse() * Rfp.transpose()) * Dyf;
  Dzp = -(Rpp.inverse() * Rfp.transpose()) * Dzf;

  Dx1.segment(2 * seg_num + 4, 2 * seg_num - 2) = Dxp;
  Dy1.segment(2 * seg_num + 4, 2 * seg_num - 2) = Dyp;
  Dz1.segment(2 * seg_num + 4, 2 * seg_num - 2) = Dzp;

//转换到多项式系数
  Px = (A.inverse() * Ct) * Dx1;
  Py = (A.inverse() * Ct) * Dy1;
  Pz = (A.inverse() * Ct) * Dz1;

//存储多项式系数
  for (int i = 0; i < seg_num; i++)
  {
    poly_coeff.block(i, 0, 1, 6) = Px.segment(i * 6, 6).transpose();
    poly_coeff.block(i, 6, 1, 6) = Py.segment(i * 6, 6).transpose();
    poly_coeff.block(i, 12, 1, 6) = Pz.segment(i * 6, 6).transpose();
  }

  /* ---------- use polynomials ---------- */
  PolynomialTraj poly_traj;//
  for (int i = 0; i < poly_coeff.rows(); ++i)//循环每一行即每一段
  {
    vector<double> cx(6), cy(6), cz(6);
    for (int j = 0; j < 6; ++j)
    {
      cx[j] = poly_coeff(i, j), cy[j] = poly_coeff(i, j + 6), cz[j] = poly_coeff(i, j + 12);
    }
    //取个反向 从c = [c0, c1, c2,...c5]'   到c = [c5, c4, c3,...c0]'
    reverse(cx.begin(), cx.end());
    reverse(cy.begin(), cy.end());
    reverse(cz.begin(), cz.end());
    double ts = Time(i);//给出分配的时间
    poly_traj.addSegment(cx, cy, cz, ts);//这一段的轨迹
  }

  return poly_traj;//返回总的轨迹
}

PolynomialTraj PolynomialTraj::one_segment_traj_gen(const Eigen::Vector3d &start_pt, const Eigen::Vector3d &start_vel, const Eigen::Vector3d &start_acc,
                                                    const Eigen::Vector3d &end_pt, const Eigen::Vector3d &end_vel, const Eigen::Vector3d &end_acc,
                                                    double t)
//一段轨迹的获取
{
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(6, 6), Crow(1, 6);
  Eigen::VectorXd Bx(6), By(6), Bz(6);

//C矩阵(其实应该是M矩阵 ppt339页)
  C(0, 5) = 1;
  C(1, 4) = 1;
  C(2, 3) = 2;
  Crow << pow(t, 5), pow(t, 4), pow(t, 3), pow(t, 2), t, 1;
  C.row(3) = Crow;
  Crow << 5 * pow(t, 4), 4 * pow(t, 3), 3 * pow(t, 2), 2 * t, 1, 0;
  C.row(4) = Crow;
  Crow << 20 * pow(t, 3), 12 * pow(t, 2), 6 * t, 2, 0, 0;
  C.row(5) = Crow;
//起点终点的速度加速度约束
  Bx << start_pt(0), start_vel(0), start_acc(0), end_pt(0), end_vel(0), end_acc(0);
  By << start_pt(1), start_vel(1), start_acc(1), end_pt(1), end_vel(1), end_acc(1);
  Bz << start_pt(2), start_vel(2), start_acc(2), end_pt(2), end_vel(2), end_acc(2);
//求解出这段的多项式系数 C*p=b（C为M矩阵，b为导数约束，即d）参考ppt338
  Eigen::VectorXd Cofx = C.colPivHouseholderQr().solve(Bx);
  Eigen::VectorXd Cofy = C.colPivHouseholderQr().solve(By);
  Eigen::VectorXd Cofz = C.colPivHouseholderQr().solve(Bz);


//多项式系数
  vector<double> cx(6), cy(6), cz(6);
  for (int i = 0; i < 6; i++)
  {
    cx[i] = Cofx(i);
    cy[i] = Cofy(i);
    cz[i] = Cofz(i);
  }

  PolynomialTraj poly_traj;
  poly_traj.addSegment(cx, cy, cz, t);//把这一段加进去

  return poly_traj;
}
