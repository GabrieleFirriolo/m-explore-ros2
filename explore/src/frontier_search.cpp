#include <explore/costmap_tools.h>
#include <explore/frontier_search.h>

#include <geometry_msgs/msg/point.hpp>
#include <mutex>

#include "nav2_costmap_2d/cost_values.hpp"

namespace frontier_exploration
{
using nav2_costmap_2d::FREE_SPACE;
using nav2_costmap_2d::LETHAL_OBSTACLE;
using nav2_costmap_2d::NO_INFORMATION;

FrontierSearch::FrontierSearch(nav2_costmap_2d::Costmap2D* costmap,
                               double potential_scale, double gain_scale,
                               double min_frontier_size)
  : costmap_(costmap)
  , potential_scale_(potential_scale)
  , gain_scale_(gain_scale)
  , min_frontier_size_(min_frontier_size)
{
}

double euclideanDist(uint x1, uint y1, uint x2, uint y2) {
  return sqrt(pow((double(x1) - double(x2)), 2.0) + pow((double(y1) - double(y2)), 2.0));
}

int FrontierSearch::astarCost(uint startx, uint starty, uint goalx, uint goaly) {
  struct Node {
    uint x, y;
    int cost;
    double heuristic;
    bool operator<(const Node& other) const {
      return (cost + heuristic) > (other.cost + other.heuristic);
    }
  };

  std::priority_queue<Node> open_set;
  std::set<std::pair<uint, uint>> closed_set;

  open_set.push({startx, starty, 0, euclideanDist(startx, starty, goalx, goaly)});

  while (!open_set.empty()) {
    Node current = open_set.top();
    open_set.pop();

    if (current.x == goalx && current.y == goaly) {
      return current.cost;
    }

    if (closed_set.count({current.x, current.y})) {
      continue;
    }
    closed_set.insert({current.x, current.y});

    if (closed_set.size() > 20000) {
      return 10000;
    }

    uint idx = costmap_->getIndex(current.x, current.y);
    for (uint nbr : nhood4(idx, *costmap_)) {
      if (map_[nbr] != LETHAL_OBSTACLE) {
        uint nx, ny;
        costmap_->indexToCells(nbr, nx, ny);
        if (!closed_set.count({nx, ny})) {
          open_set.push({nx, ny, current.cost + 1, euclideanDist(nx, ny, goalx, goaly)});
        }
      }
    }
  }

  return -1; // Percorso non trovato
}

std::vector<Frontier>
FrontierSearch::searchFrom(geometry_msgs::msg::Point position)
{
  std::vector<Frontier> frontier_list;

  unsigned int mx, my;
  if (!costmap_->worldToMap(position.x, position.y, mx, my)) {
    RCLCPP_ERROR(rclcpp::get_logger("FrontierSearch"), "Robot out of costmap bounds");
    return frontier_list;
  }

  std::lock_guard<nav2_costmap_2d::Costmap2D::mutex_t> lock(*(costmap_->getMutex()));
  map_ = costmap_->getCharMap();
  size_x_ = costmap_->getSizeInCellsX();
  size_y_ = costmap_->getSizeInCellsY();

  std::vector<bool> frontier_flag(size_x_ * size_y_, false);
  std::vector<bool> visited_flag(size_x_ * size_y_, false);
  std::queue<unsigned int> bfs;

  unsigned int clear, pos = costmap_->getIndex(mx, my);
  if (nearestCell(clear, pos, FREE_SPACE, *costmap_)) {
    bfs.push(clear);
  } else {
    bfs.push(pos);
    RCLCPP_WARN(rclcpp::get_logger("FrontierSearch"), "Fallback to current position");
  }
  visited_flag[bfs.front()] = true;

  while (!bfs.empty()) {
    unsigned int idx = bfs.front();
    bfs.pop();

    for (unsigned nbr : nhood4(idx, *costmap_)) {
      if (map_[nbr] <= map_[idx] && !visited_flag[nbr]) {
        visited_flag[nbr] = true;
        bfs.push(nbr);
      } else if (isNewFrontierCell(nbr, frontier_flag)) {
        frontier_flag[nbr] = true;
        Frontier new_frontier = buildNewFrontier(nbr, pos, frontier_flag);
        if (new_frontier.size * costmap_->getResolution() >= min_frontier_size_) {
          frontier_list.push_back(new_frontier);
        }
      }
    }
  }

  // Calcolo distanza con A* e riordino frontiere
  std::vector<Frontier> backup = frontier_list;
  int count = 0;
  for (auto & frontier : frontier_list) {
    if (count++ < 15) {
      unsigned int gx, gy;
      if (costmap_->worldToMap(frontier.centroid.x, frontier.centroid.y, gx, gy)) {
        int dist = astarCost(mx, my, gx, gy);
        frontier.min_distance = (dist < 0) ? 10000 : static_cast<double>(dist);
      } else {
        frontier.min_distance = 10000;
      }
    } else {
      frontier.min_distance = 10000;
    }
    frontier.cost = frontierCost(frontier);
  }

  std::sort(frontier_list.begin(), frontier_list.end(),
            [](const Frontier& f1, const Frontier& f2) {
              return f1.cost < f2.cost;
            });

  if (!frontier_list.empty() && frontier_list[0].min_distance == 10000) {
    return backup;  // fallback se A* fallisce
  }

  return frontier_list;
}

Frontier FrontierSearch::buildNewFrontier(unsigned int initial_cell,
                                          unsigned int reference,
                                          std::vector<bool>& frontier_flag)
{
  // initialize frontier structure
  Frontier output;
  output.centroid.x = 0;
  output.centroid.y = 0;
  output.size = 1;
  output.min_distance = std::numeric_limits<double>::infinity();

  // record initial contact point for frontier
  unsigned int ix, iy;
  costmap_->indexToCells(initial_cell, ix, iy);
  costmap_->mapToWorld(ix, iy, output.initial.x, output.initial.y);

  // push initial gridcell onto queue
  std::queue<unsigned int> bfs;
  bfs.push(initial_cell);

  // cache reference position in world coords
  unsigned int rx, ry;
  double reference_x, reference_y;
  costmap_->indexToCells(reference, rx, ry);
  costmap_->mapToWorld(rx, ry, reference_x, reference_y);

  while (!bfs.empty()) {
    unsigned int idx = bfs.front();
    bfs.pop();

    // try adding cells in 8-connected neighborhood to frontier
    for (unsigned int nbr : nhood8(idx, *costmap_)) {
      // check if neighbour is a potential frontier cell
      if (isNewFrontierCell(nbr, frontier_flag)) {
        // mark cell as frontier
        frontier_flag[nbr] = true;
        unsigned int mx, my;
        double wx, wy;
        costmap_->indexToCells(nbr, mx, my);
        costmap_->mapToWorld(mx, my, wx, wy);

        geometry_msgs::msg::Point point;
        point.x = wx;
        point.y = wy;
        output.points.push_back(point);

        // update frontier size
        output.size++;

        // update centroid of frontier
        output.centroid.x += wx;
        output.centroid.y += wy;

        // determine frontier's distance from robot, going by closest gridcell
        // to robot
        double distance = sqrt(pow((double(reference_x) - double(wx)), 2.0) +
                               pow((double(reference_y) - double(wy)), 2.0));
        if (distance < output.min_distance) {
          output.min_distance = distance;
          output.middle.x = wx;
          output.middle.y = wy;
        }

        // add to queue for breadth first search
        bfs.push(nbr);
      }
    }
  }

  // average out frontier centroid
  output.centroid.x /= output.size;
  output.centroid.y /= output.size;
  return output;
}

bool FrontierSearch::isNewFrontierCell(unsigned int idx,
                                       const std::vector<bool>& frontier_flag)
{
  // check that cell is unknown and not already marked as frontier
  if (map_[idx] != NO_INFORMATION || frontier_flag[idx]) {
    return false;
  }

  // frontier cells should have at least one cell in 4-connected neighbourhood
  // that is free
  for (unsigned int nbr : nhood4(idx, *costmap_)) {
    if (map_[nbr] == FREE_SPACE) {
      return true;
    }
  }

  return false;
}

double FrontierSearch::frontierCost(const Frontier& frontier)
{
  return (potential_scale_ * frontier.min_distance *
          costmap_->getResolution()) -
         (gain_scale_ * frontier.size * costmap_->getResolution());
}
}  // namespace frontier_exploration
