#include "node.h"

const std::vector<double> d2q9Node::d2q9_weights = {1.0/36, 1.0/9, 1.0/36,
                                                    1.0/9 , 4.0/9, 1.0/9,
                                                    1.0/36, 1.0/9, 1.0/36};

const std::vector<Eigen::Vector2i> d2q9Node::d2q9_directions = {
  Eigen::Vector2i(-1, -1), Eigen::Vector2i(0,-1), Eigen::Vector2i(1,-1),
  Eigen::Vector2i(-1,0),   Eigen::Vector2i(0,0),  Eigen::Vector2i(1,0),
  Eigen::Vector2i(-1,1),   Eigen::Vector2i(0,1),  Eigen::Vector2i(1,1)
};

double d2q9Node::density() {
  double rho = std::accumulate(f_density.begin(), f_density.end(), 0);
  return rho;
};

Eigen::Vector2d d2q9Node::velocity() {
  Eigen::Vector2d vel(0,0);
  for(int i = 0; i < 9; ++i) {
    vel += d2q9_directions[i].cast<double>()*f_density[i];
  }

  return vel;
}


