#include "path_searching/dyn_a_star.h"

using namespace std;
using namespace Eigen;

AStar::~AStar()
{
    for (int i = 0; i < POOL_SIZE_(0); i++)
        for (int j = 0; j < POOL_SIZE_(1); j++)
            for (int k = 0; k < POOL_SIZE_(2); k++)
                delete GridNodeMap_[i][j][k];//局部地图清零
}

void AStar::initGridMap(GridMap::Ptr occ_map, const Eigen::Vector3i pool_size)//输入地图和一个局部地图大小
{
    POOL_SIZE_ = pool_size;//局部地图大小
    CENTER_IDX_ = pool_size / 2;//中心位置

    GridNodeMap_ = new GridNodePtr **[POOL_SIZE_(0)];//三级指针建立三维数组 同深蓝学院第二节
    for (int i = 0; i < POOL_SIZE_(0); i++)
    {
        GridNodeMap_[i] = new GridNodePtr *[POOL_SIZE_(1)];
        for (int j = 0; j < POOL_SIZE_(1); j++)
        {
            GridNodeMap_[i][j] = new GridNodePtr[POOL_SIZE_(2)];
            for (int k = 0; k < POOL_SIZE_(2); k++)
            {
                GridNodeMap_[i][j][k] = new GridNode;
            }
        }
    }

    grid_map_ = occ_map;//地图
}

double AStar::getDiagHeu(GridNodePtr node1, GridNodePtr node2)//给定两个节点，返回Diag启发函数值
{
    double dx = abs(node1->index(0) - node2->index(0));
    double dy = abs(node1->index(1) - node2->index(1));
    double dz = abs(node1->index(2) - node2->index(2));

    double h = 0.0;
    int diag = min(min(dx, dy), dz);
    dx -= diag;
    dy -= diag;
    dz -= diag;

    if (dx == 0)
    {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dy, dz) + 1.0 * abs(dy - dz);
    }
    if (dy == 0)
    {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dx, dz) + 1.0 * abs(dx - dz);
    }
    if (dz == 0)
    {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dx, dy) + 1.0 * abs(dx - dy);
    }
    return h;
}

double AStar::getManhHeu(GridNodePtr node1, GridNodePtr node2)//给定两个节点，返回Manhadun启发函数值
{
    double dx = abs(node1->index(0) - node2->index(0));
    double dy = abs(node1->index(1) - node2->index(1));
    double dz = abs(node1->index(2) - node2->index(2));

    return dx + dy + dz;
}

double AStar::getEuclHeu(GridNodePtr node1, GridNodePtr node2)//给定两个节点，返回欧式距离启发函数值
{
    return (node2->index - node1->index).norm();//
}

vector<GridNodePtr> AStar::retrievePath(GridNodePtr current)//输入当前节点，回溯路径
{
    vector<GridNodePtr> path;//path为节点指针向量
    path.push_back(current);//加入当前节点指针

    while (current->cameFrom != NULL)//当父节点指针不为空的时候
    {
        current = current->cameFrom;
        path.push_back(current);//递推寻找父节点指针，加入节点指针向量
    }

    return path;
}

bool AStar::ConvertToIndexAndAdjustStartEndPoints(Vector3d start_pt, Vector3d end_pt, Vector3i &start_idx, Vector3i &end_idx)//输入开始点和结束点的实际位置和节点index
//改变起点终点到合理的位置
{
    if (!Coord2Index(start_pt, start_idx) || !Coord2Index(end_pt, end_idx))//如果起点或终点任何一个不在地图范围内
        return false;

    if (checkOccupancy(Index2Coord(start_idx)))//判断起点是否为障碍  超出地图的范围和是障碍物都算作是障碍，要改变起点
    {
        //ROS_WARN("Start point is insdide an obstacle.");
        do
        {
            start_pt = (start_pt - end_pt).normalized() * step_size_ + start_pt;//normalized参考以下文章。其实就是归一化https://zhidao.baidu.com/question/1513657298408540660.html
            //改变起点的位置，但是为什么这么改变？？
            if (!Coord2Index(start_pt, start_idx))//不在地图范围内 return false
                return false;
        } while (checkOccupancy(Index2Coord(start_idx)));//当起点为障碍物时循环
    }

    if (checkOccupancy(Index2Coord(end_idx)))//判断终点是否为障碍物
    {
        //ROS_WARN("End point is insdide an obstacle.");
        do
        {
            end_pt = (end_pt - start_pt).normalized() * step_size_ + end_pt;
            if (!Coord2Index(end_pt, end_idx))
                return false;
        } while (checkOccupancy(Index2Coord(end_idx)));
    }

    return true;
}

bool AStar::AstarSearch(const double step_size, Vector3d start_pt, Vector3d end_pt)
{
    ros::Time time_1 = ros::Time::now();//时间开始
    ++rounds_;//设置为1

    step_size_ = step_size;
    inv_step_size_ = 1 / step_size;
    center_ = (start_pt + end_pt) / 2;//中间实际位置点，（起点+终点）/2

    Vector3i start_idx, end_idx;
    if (!ConvertToIndexAndAdjustStartEndPoints(start_pt, end_pt, start_idx, end_idx))
    {
        ROS_ERROR("Unable to handle the initial or end point, force return!");
        return false;
    }

    // if ( start_pt(0) > -1 && start_pt(0) < 0 )
    //     cout << "start_pt=" << start_pt.transpose() << " end_pt=" << end_pt.transpose() << endl;

    GridNodePtr startPtr = GridNodeMap_[start_idx(0)][start_idx(1)][start_idx(2)];//起点节点
    GridNodePtr endPtr = GridNodeMap_[end_idx(0)][end_idx(1)][end_idx(2)];//终点节点

    std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, NodeComparator> empty;
    //优先队列的定义，参考https://blog.csdn.net/weixin_36888577/article/details/79937886
    //priority_queue<Type, Container, Functional>     
    //Type 就是数据类型，Container 就是容器类型（Container必须是用数组实现的容器，比如vector,deque等等，但不能用 list。STL里面默认用的是vector），Functional 就是比较的方式
    //此处的比较方式为比较节点的f值大小

    openSet_.swap(empty);//初始化openSet全空

    GridNodePtr neighborPtr = NULL;
    GridNodePtr current = NULL;

    //初始化起点节点
    startPtr->index = start_idx;
    startPtr->rounds = rounds_;
    startPtr->gScore = 0;
    startPtr->fScore = getHeu(startPtr, endPtr);
    startPtr->state = GridNode::OPENSET; //put start node in open set起点状态设置为open
    startPtr->cameFrom = NULL;
    openSet_.push(startPtr); //put start in open set 起点加到openSet_中

    endPtr->index = end_idx;

    double tentative_gScore;

    int num_iter = 0;//初始化迭代次数为0

    //主循环过程
    while (!openSet_.empty())//openSet_不为空
    {
        num_iter++;//迭代次数 +1
        current = openSet_.top();//弹出f最小的作为当前节点
        openSet_.pop();//把第一个删除，即将当前节点从openSet_中拿出   https://zouzhongliang.com/index.php/2019/07/21/priority_queueyouxianduiliepophanshuyunyongshili/
 
        // if ( num_iter < 10000 )
        //     cout << "current=" << current->index.transpose() << endl;

        if (current->index(0) == endPtr->index(0) && current->index(1) == endPtr->index(1) && current->index(2) == endPtr->index(2))//如果这个点就是终点
        {
            // ros::Time time_2 = ros::Time::now();
            // printf("\033[34mA star iter:%d, time:%.3f\033[0m\n",num_iter, (time_2 - time_1).toSec()*1000);
            // if((time_2 - time_1).toSec() > 0.1)
            //     ROS_WARN("Time consume in A star path finding is %f", (time_2 - time_1).toSec() );
            gridPath_ = retrievePath(current);//从当前点开始回溯所有节点
            return true;
        }
        current->state = GridNode::CLOSEDSET; //当前节点状态变为close

        for (int dx = -1; dx <= 1; dx++)//循环邻居节点
            for (int dy = -1; dy <= 1; dy++)
                for (int dz = -1; dz <= 1; dz++)
                {
                    if (dx == 0 && dy == 0 && dz == 0)//不考虑这个点本身
                        continue;

                    Vector3i neighborIdx;//邻居节点的index
                    neighborIdx(0) = (current->index)(0) + dx;
                    neighborIdx(1) = (current->index)(1) + dy;
                    neighborIdx(2) = (current->index)(2) + dz;

                    if (neighborIdx(0) < 1 || neighborIdx(0) >= POOL_SIZE_(0) - 1 || neighborIdx(1) < 1 || neighborIdx(1) >= POOL_SIZE_(1) - 1 || neighborIdx(2) < 1 || neighborIdx(2) >= POOL_SIZE_(2) - 1)
                    {
                        continue;//邻居节点超过了地图的范围
                    }//

                    neighborPtr = GridNodeMap_[neighborIdx(0)][neighborIdx(1)][neighborIdx(2)];//邻居节点结构体指针
                    neighborPtr->index = neighborIdx;//赋值index

                    bool flag_explored = neighborPtr->rounds == rounds_;//是否扩展的标志 为1(代表已经扩展)返回true，不相等(代表未扩展)返回false

                    if (flag_explored && neighborPtr->state == GridNode::CLOSEDSET)//如果已经扩展并且状态为close，代表在closeset中
                    {
                        continue; //in closed set.
                    }

                    neighborPtr->rounds = rounds_;//将邻居节点的round设置为1,代表已经扩展

                    if (checkOccupancy(Index2Coord(neighborPtr->index)))
                    {
                        continue;//如果邻居节点是障碍物
                    }

                    double static_cost = sqrt(dx * dx + dy * dy + dz * dz);//从当前点到邻居节点所走的距离delta g
                    tentative_gScore = current->gScore + static_cost;//从当前节点到邻居节点的g

                    if (!flag_explored)//如果节点未被扩展过
                    {
                        //discover a new node
                        neighborPtr->state = GridNode::OPENSET;//邻居节点状态设为open
                        neighborPtr->cameFrom = current;//当前节点为其父节点
                        neighborPtr->gScore = tentative_gScore;//邻居节点g值设置
                        neighborPtr->fScore = tentative_gScore + getHeu(neighborPtr, endPtr);//邻居节点f值设置
                        openSet_.push(neighborPtr); //邻居节点加入openSet_   put neighbor in open set and record it.
                    }
                    else if (tentative_gScore < neighborPtr->gScore)//如果已经在open中了而且经当前节点到达邻居节点的g比之前的小，就更新g值
                    { //in open set and need update
                        neighborPtr->cameFrom = current;//父节点改为当前节点
                        neighborPtr->gScore = tentative_gScore;//更新g值
                        neighborPtr->fScore = tentative_gScore + getHeu(neighborPtr, endPtr);//更新f值
                    }
                }
        ros::Time time_2 = ros::Time::now();//时间
        if ((time_2 - time_1).toSec() > 0.2)
        {
            ROS_WARN("Failed in A star path searching !!! 0.2 seconds time limit exceeded.");
            return false;//流程超过0.2s就报错
        }
    }

    ros::Time time_2 = ros::Time::now();//结束时间

    if ((time_2 - time_1).toSec() > 0.1)
        ROS_WARN("Time consume in A star path finding is %.3fs, iter=%d", (time_2 - time_1).toSec(), num_iter);

    return false;
}

vector<Vector3d> AStar::getPath()//获得路径
{
    vector<Vector3d> path;

    for (auto ptr : gridPath_)//循环节点向量geid_path
        path.push_back(Index2Coord(ptr->index));//将节点index转为实际位置后存储在path中

    reverse(path.begin(), path.end());//
    return path;//翻转
}
