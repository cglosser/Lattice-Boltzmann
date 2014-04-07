#include "node.h"

const std::vector<double> d2q9Node::d2q9_weights = {4.0/9,
                                                    1.0/9,  1.0/9,  1.0/9,  1.0/9,
                                                    1.0/36, 1.0/36, 1.0/36, 1.0/36};

const std::vector<Eigen::Vector2i> d2q9Node::d2q9_directions = {
  Eigen::Vector2i(0, 0), 
  Eigen::Vector2i(1,0), Eigen::Vector2i(0,1),  Eigen::Vector2i(-1,0),  Eigen::Vector2i(0,-1),
  Eigen::Vector2i(1,1), Eigen::Vector2i(-1,1), Eigen::Vector2i(-1,-1), Eigen::Vector2i(1,-1)
};

Eigen::Vector2d d2q9Node::massCurrent() {
  Eigen::Vector2d vel(0,0);
  for(int i = 0; i < 9; ++i) 
    vel += d2q9_directions[i].cast<double>()*f_density[i];

  return vel;
}

double d2q9Node::density() {
  double rho = std::accumulate(f_density.begin(), f_density.end(), 0.0);
  return rho;
}

Eigen::Vector2d d2q9Node::velocity() {
  return massCurrent()/density();
}

std::vector<double> d2q9Node::equilibrium() {
  std::vector<double> result(9, 0);

  for(int dir = 0; dir < 9; ++dir) {
    double e_dot_u = velocity().dot(d2q9_directions[dir].cast<double>());
    result[dir]    = density()*d2q9_weights[dir]*
      (1 + 3*e_dot_u + 4.5*e_dot_u*e_dot_u - 1.5*velocity().squaredNorm());
  }

  return result;
}

std::vector<double> d2q9Node::equilibrium(double rho, Eigen::Vector2d vel) {
  std::vector<double> result(9, 0);

  for(int dir = 0; dir < 9; ++dir) {
    double e_dot_u = vel.dot(d2q9_directions[dir].cast<double>());
    result[dir] = rho*d2q9_weights[dir]*(1 + 3*e_dot_u + 4.5*e_dot_u*e_dot_u - 1.5*vel.squaredNorm());
  }

  return result;
}

int d2q9Node::reverse(int dir) {
  if(dir == 0)
    return dir;
  else if (dir <= 4)
    return (dir + 1)%4 + 1;
  else 
    return (dir + 1)%4 + 5;
}
