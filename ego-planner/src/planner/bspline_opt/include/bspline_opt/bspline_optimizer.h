#ifndef _BSPLINE_OPTIMIZER_H_
#define _BSPLINE_OPTIMIZER_H_

#include <Eigen/Eigen>
#include <path_searching/dyn_a_star.h>
#include <bspline_opt/uniform_bspline.h>
#include <plan_env/grid_map.h>
#include <ros/ros.h>
#include "bspline_opt/lbfgs.hpp"

// Gradient and elasitc band optimization 梯度和弹性带优化 

// Input: 有符号距离场和点序列a signed distance field and a sequence of points
// Output: 优化的点序列the optimized sequence of points
// The format of points:控制点矩阵，N行3列 N x 3 matrix, each row is a point
namespace ego_planner
{

  class ControlPoints//控制点的类
  {
  public:
    double clearance;
    int size;//控制点个数
    Eigen::MatrixXd points;//一系列点（向量）组成的矩阵B-spline control points
    std::vector<std::vector<Eigen::Vector3d>> base_point; // 方向向量起点处的点（碰撞点）The point at the statrt of the direction vector (collision point)
    std::vector<std::vector<Eigen::Vector3d>> direction;  // 方向向量（归一化）Direction vector, must be normalized.
    std::vector<bool> flag_temp;                          //标志，使用到时会初始化它 A flag that used in many places. Initialize it everytime before using it.
    // std::vector<bool> occupancy;

    void resize(const int size_set)//设置矩阵和向量的大小
    {
      size = size_set;
      //初始全部清空
      base_point.clear();
      direction.clear();
      flag_temp.clear();
      // occupancy.clear();

      points.resize(3, size_set);//三行size列，即表示成一个一个点。三行为xyz，一列代表一个点
      base_point.resize(size);
      direction.resize(size);
      flag_temp.resize(size);
      // occupancy.resize(size);
    }
  };

  class BsplineOptimizer//B样条优化的类
  {

  public:
    BsplineOptimizer() {}
    ~BsplineOptimizer() {}

    /* main API */
    void setEnvironment(const GridMap::Ptr &env);//设置环境加载地图
    void setParam(ros::NodeHandle &nh);//设置参数
    //找不到定义和引用？？？？？
    Eigen::MatrixXd BsplineOptimizeTraj(const Eigen::MatrixXd &points, const double &ts,
                                        const int &cost_function, int max_num_id, int max_time_id);

    /* helper function */

    // required inputs
    void setControlPoints(const Eigen::MatrixXd &points);//初始化控制点，返回碰撞段A*算法搜索的路径
    void setBsplineInterval(const double &ts);//设置b样条曲线的knot span（时间间隔）
    void setCostFunction(const int &cost_function);//输入cost值
    void setTerminateCond(const int &max_num_id, const int &max_time_id);//限制，最大的个数id和时间id

    // optional inputs
    void setGuidePath(const vector<Eigen::Vector3d> &guide_pt);//引导路径
    void setWaypoints(const vector<Eigen::Vector3d> &waypts,
                      const vector<int> &waypt_idx); // 设置路径点N-2 constraints at most

    void optimize();//找不到定义和引用？？？？？

    Eigen::MatrixXd getControlPoints();//找不到定义和引用？？？？？

    AStar::Ptr a_star_;//A*路径搜索的指针
    std::vector<Eigen::Vector3d> ref_pts_;//refine（细化）点

    std::vector<std::vector<Eigen::Vector3d>> initControlPoints(Eigen::MatrixXd &init_points, bool flag_first_init = true);//初始化控制点，返回碰撞段A*算法搜索的路径
    bool BsplineOptimizeTrajRebound(Eigen::MatrixXd &optimal_points, double ts); //输入控制点和时间间隔，返回是否规划成功 must be called after initControlPoints()
    bool BsplineOptimizeTrajRefine(const Eigen::MatrixXd &init_points, const double ts, Eigen::MatrixXd &optimal_points);//返回轨迹是否安全

    inline int getOrder(void) { return order_; }//获得样条曲线次数

  private:
    GridMap::Ptr grid_map_;//栅格地图指针

    enum FORCE_STOP_OPTIMIZE_TYPE//强制停止的原因
    {
      DONT_STOP,
      STOP_FOR_REBOUND,
      STOP_FOR_ERROR
    } force_stop_type_;

    // main input主要的输入
    // Eigen::MatrixXd control_points_;     // B-spline control points, N x dim 曲线控制点
    double bspline_interval_; // B-spline knot span时间间隔
    Eigen::Vector3d end_pt_;  // end of the trajectory 终点
    // int             dim_;                // dimension of the B-spline b样条曲线的纬度
    //
    vector<Eigen::Vector3d> guide_pts_; //  几何引导路径点 geometric guiding path points, N-6
    vector<Eigen::Vector3d> waypoints_; //路径点约束 waypts constraints
    vector<int> waypt_idx_;             // 路径点约束的index waypts constraints index
                                        //
    int max_num_id_, max_time_id_;      // 停止标准 stopping criteria
    int cost_function_;                 // 损失函数值 used to determine objective function
    double start_time_;                 // 开始时间global time for moving obstacles

    /* optimization parameters 优化参数*/
    int order_;                    // b样条曲线阶数 bspline degree
    //优化量系数lamada
    double lambda1_;               // jerk smoothness weight
    double lambda2_, new_lambda2_; // distance weight
    double lambda3_;               // feasibility weight
    double lambda4_;               // curve fitting

    int a;
    //
    double dist0_;             // 安全距离safe distance
    double max_vel_, max_acc_; //动力学约束，最大速度加速度dynamic limits

    int variable_num_;              // 优化变量数目optimization variables
    int iter_num_;                  // 求解器的迭代次数iteration of the solver
    Eigen::VectorXd best_variable_; //最优变量
    double min_cost_;               //最小cost

    ControlPoints cps_;//控制点

    /* cost function损失函数方程 */
    /* calculate each part of cost function with control points q as input */
    //输入控制点q，计算损失函数的各个部分
    static double costFunction(const std::vector<double> &x, std::vector<double> &grad, void *func_data);
    void combineCost(const std::vector<double> &x, vector<double> &grad, double &cost);

    // q contains all control points
    void calcSmoothnessCost(const Eigen::MatrixXd &q, double &cost,
                            Eigen::MatrixXd &gradient, bool falg_use_jerk = true);//计算光滑项的cost
    void calcFeasibilityCost(const Eigen::MatrixXd &q, double &cost,
                             Eigen::MatrixXd &gradient);//计算可行性项的cost
    void calcDistanceCostRebound(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient, int iter_num, double smoothness_cost);
    void calcFitnessCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient);//计算碰撞项的cost
    bool check_collision_and_rebound(void);//检查是否碰撞和rebound（如果新的障碍物有效就为true，否则为false)

    static int earlyExit(void *func_data, const double *x, const double *g, const double fx, const double xnorm, const double gnorm, const double step, int n, int k, int ls);
    static double costFunctionRebound(void *func_data, const double *x, double *grad, const int n);
    static double costFunctionRefine(void *func_data, const double *x, double *grad, const int n);

    bool rebound_optimize();
    bool refine_optimize();
    void combineCostRebound(const double *x, double *grad, double &f_combine, const int n);
    void combineCostRefine(const double *x, double *grad, double &f_combine, const int n);

    /* for benckmark evaluation only */
  public:
    typedef unique_ptr<BsplineOptimizer> Ptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  };

} // namespace ego_planner
#endif