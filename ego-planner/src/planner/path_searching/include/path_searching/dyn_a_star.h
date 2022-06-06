#ifndef _DYN_A_STAR_H_
#define _DYN_A_STAR_H_

#include <iostream>
#include <ros/ros.h>
#include <ros/console.h>
#include <Eigen/Eigen>
#include <plan_env/grid_map.h>
#include <queue>

constexpr double inf = 1 >> 20;
struct GridNode;//节点结构体
typedef GridNode *GridNodePtr;//节点结构体指针

struct GridNode//定义GridNode结构体，代表每个节点
{
	enum enum_state{//枚举变量：1：openlist    2：closelist    3：未定义（还未被搜索到）
		OPENSET = 1,
		CLOSEDSET = 2,
		UNDEFINED = 3
	};

	int rounds{0}; // Distinguish every call
	enum enum_state state
	{
		UNDEFINED
	};//节点的初始状态全部未定义
	Eigen::Vector3i index;//三维向量，代表index

	double gScore{inf}, fScore{inf};//初始的g和f均为正无穷
	GridNodePtr cameFrom{NULL};//父节点的指针,初始为NULL
};

class NodeComparator//给定两个节点结构体指针，返回其所对应的结构体中谁的f大  node1大为true，反之为false
{
public:
	bool operator()(GridNodePtr node1, GridNodePtr node2)
	{
		return node1->fScore > node2->fScore;
	}
};

class AStar//定义Astar 的类,包括很多A*的函数
{
private://私密变量，不可以在类之外调用
	GridMap::Ptr grid_map_;//建立了一个地图GridMap

	inline void coord2gridIndexFast(const double x, const double y, const double z, int &id_x, int &id_y, int &id_z);
//启发式函数值获取函数
	double getDiagHeu(GridNodePtr node1, GridNodePtr node2);
	double getManhHeu(GridNodePtr node1, GridNodePtr node2);
	double getEuclHeu(GridNodePtr node1, GridNodePtr node2);
	inline double getHeu(GridNodePtr node1, GridNodePtr node2);
//改变起点终点到合理的位置
	bool ConvertToIndexAndAdjustStartEndPoints(const Eigen::Vector3d start_pt, const Eigen::Vector3d end_pt, Eigen::Vector3i &start_idx, Eigen::Vector3i &end_idx);
//节点index转换到实际位置
	inline Eigen::Vector3d Index2Coord(const Eigen::Vector3i &index) const;
	//给定位置计算index，返回是否在地图范围内
	inline bool Coord2Index(const Eigen::Vector3d &pt, Eigen::Vector3i &idx) const;

	//bool (*checkOccupancyPtr)( const Eigen::Vector3d &pos );
//输入位置，检查是否是障碍true or falsegetInflateOccupancy这个函数本来int返回为0（不是障碍物），1（是障碍物）-1（不在地图内）。bool处理后，是障碍物和不再地图内都算作障碍
	inline bool checkOccupancy(const Eigen::Vector3d &pos) { return (bool)grid_map_->getInflateOccupancy(pos); }
//回溯路径上的节点指针
	std::vector<GridNodePtr> retrievePath(GridNodePtr current);

	double step_size_, inv_step_size_;
	Eigen::Vector3d center_;//中心位置的实际位置
	Eigen::Vector3i CENTER_IDX_, POOL_SIZE_;//CENTER_IDX_:中心位置的index POOL_SIZE_:地图最大的index
	const double tie_breaker_ = 1.0 + 1.0 / 10000;//tie_breaker_防止路径的对称性

	std::vector<GridNodePtr> gridPath_;//路径回溯的节点指针向量

	GridNodePtr ***GridNodeMap_;//三维地图数组，存储节点结构体指针
	std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, NodeComparator> openSet_;// openSet_队列

	int rounds_{0};//节点是否被扩展的标志

public:
	typedef std::shared_ptr<AStar> Ptr;

	AStar(){};
	~AStar();
//初始化地图
	void initGridMap(GridMap::Ptr occ_map, const Eigen::Vector3i pool_size);
//Astar搜索是否成功
	bool AstarSearch(const double step_size, Eigen::Vector3d start_pt, Eigen::Vector3d end_pt);
//获取路径实际位置点
	std::vector<Eigen::Vector3d> getPath();
};
//加入 tie_breaker_后的启发式距离h
inline double AStar::getHeu(GridNodePtr node1, GridNodePtr node2)
{
	return tie_breaker_ * getDiagHeu(node1, node2);
}

inline Eigen::Vector3d AStar::Index2Coord(const Eigen::Vector3i &index) const
{
	return ((index - CENTER_IDX_).cast<double>() * step_size_) + center_;
};//从角标转化成实际位置

inline bool AStar::Coord2Index(const Eigen::Vector3d &pt, Eigen::Vector3i &idx) const//给定位置计算index，返回是否在地图范围内
{
	idx = ((pt - center_) * inv_step_size_ + Eigen::Vector3d(0.5, 0.5, 0.5)).cast<int>() + CENTER_IDX_;
	//疑问：为什么这么转化？center inv_step_size CENTER_IDX_是什么意思？.cast是什么意思？
	if (idx(0) < 0 || idx(0) >= POOL_SIZE_(0) || idx(1) < 0 || idx(1) >= POOL_SIZE_(1) || idx(2) < 0 || idx(2) >= POOL_SIZE_(2))
	{//超过了地图最大index或者小于最小index
		ROS_ERROR("Ran out of pool, index=%d %d %d", idx(0), idx(1), idx(2));
		return false;
	}

	return true;
};

#endif
