#ifndef NODE_H
#define NODE_H

#include <Eigen/Dense>
#include <functional>
#include <iostream>
#include <numeric>
#include <vector>

enum class NodeState {ACTIVE, WALL, INLET, OUTLET};

struct d2q9Node {
  static const std::vector<double>          d2q9_weights;
  static const std::vector<Eigen::Vector2i> d2q9_directions;

  std::vector<double> f_density, temp_density;
  NodeState state;

  Eigen::Vector2d massCurrent();
  double density();
  Eigen::Vector2d velocity();
  static int reverse(int);

  d2q9Node() {
    f_density.resize(9, 0);
    temp_density.resize(9, 0);
  }
};

#endif
