#include "lattice.h"
using namespace std;

double d2q9_weights[] = {4.0/9, 1.0/9, 1.0/36, 1.0/9, 1.0/36, 1.0/9, 1.0/36, 1.0/9, 1.0/36};

Lattice::Lattice(const int x_size, const int y_size):
    XDIM(x_size), YDIM(y_size), NUM_WEIGHTS(9), 
    weight(d2q9_weights, d2q9_weights + NUM_WEIGHTS) {
  f_density.resize(boost::extents[XDIM][YDIM][NUM_WEIGHTS]);
  push_density.resize(boost::extents[XDIM][YDIM][NUM_WEIGHTS]);
  neighbors.resize(boost::extents[XDIM][YDIM][NUM_WEIGHTS]);

  for(int x = 0; x < XDIM; x++) {
    for(int y = 0; y < YDIM; y++) {
      for(int n = 0; n < NUM_WEIGHTS; n++) {
        f_density[x][y][n] = 1;
        if(x == XDIM/2 && y == YDIM/2 && n == 0) f_density[x][y][n] = 2;
      }
    }
  }

  buildNeighbors();

  return;
}

double Lattice::density(const int x, const int y) {
  double rho = 0;
  for(int n = 0; n < NUM_WEIGHTS; n++) rho += f_density[x][y][n];
  return rho;
}

void Lattice::update() {
  streamingUpdate();
  collisionUpdate();
  return;
}

void Lattice::print(ostream &os) {
  for(int y = YDIM - 1; y >= 0; y--) {
    for(int x = 0; x < XDIM; x++) {
      os << density(x, y);
      if (x < XDIM - 1) cout << ",";
    }
    os << endl;
  }
  return;
}

void Lattice::buildNeighbors() {
  for(int x = 0; x < XDIM; x++) {
    for(int y = 0; y < YDIM; y++) {
      for(int n = 0; n < NUM_WEIGHTS; n++) {
        Eigen::Vector2d dr(directionToSteps(n));
        int xprime = x + int(dr[0]), yprime = y + int(dr[1]), nprime = n;

        if(xprime < 0) {
          xprime += XDIM;
        } else if(xprime >= XDIM) {
          xprime -= XDIM;
        }

        if(yprime < 0) {
          xprime = x; yprime = y; nprime = n - 4;
        } else if(yprime >= YDIM) {
          xprime = x; yprime = y; nprime = n + 4;
        }

        neighbors[x][y][n] = &push_density[xprime][yprime][nprime];
      }
    }
  }
  return;
}

void Lattice::streamingUpdate() {
  for(int x = 0; x < XDIM; x++) {
    for(int y = 0; y < YDIM; y++) {
      for(int n = 0; n < NUM_WEIGHTS; n++) {
        *neighbors[x][y][n] = f_density[x][y][n];
      }
    }
  }
  f_density = push_density;
  return;
}

void Lattice::collisionUpdate() {
  for(int x = 0; x < XDIM; x++) {
    for(int y = 0; y < YDIM; y++) {
      Eigen::Vector2d macro_vel(0,0);
      double density = 0;

      for(int n = 0; n < NUM_WEIGHTS; n++) {
        macro_vel += f_density[x][y][n]*directionToSteps(n);
        density += f_density[x][y][n];
      }

      if(density != 0) {
        macro_vel /= density;
      } else {
        continue; // No fluid at node, nothing further to compute
      }

      for(int n = 0; n < NUM_WEIGHTS; n++) {
        double e_dot_u = directionToSteps(n).dot(macro_vel);
        double equilibrium = (1 + 3*e_dot_u + 9.0/2*pow(e_dot_u,2)
          - 3.0/2*macro_vel.squaredNorm())*weight[n]*density; // c = 1

        f_density[x][y][n] -= f_density[x][y][n] - equilibrium;
      }
    }
  }

  return;
}


Eigen::Vector2d directionToSteps(const int n) {
  if (n == 0) {
    return Eigen::Vector2d(0, 0);
  } else {
    int r = n - 1;
    double theta = r*atan(1.0); // equals 2*pi*r/8
    int dx = round(cos(theta)) + 0, dy = round(sin(theta)) + 0;
    return Eigen::Vector2d(dx, dy);
  }
}
