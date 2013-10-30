#include <boost/multi_array.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

typedef boost::multi_array<double, 3> Lattice3D;
typedef boost::multi_array<double, 2> Lattice2D;

void streamingUpdate(Lattice3D &, const int, const int, const int);
void collisionUpdate(const int, const int, const int, const vector<double> &,
    Lattice3D &);
Eigen::Vector2d dirToSteps(const int);
void printLattice(const Lattice3D &, const int, const int, const int);
double latticeDistance(int);
int mod(int, int);

int main() {
  // Define weightings for a D2Q9 lattice.
  double weights[] = {4.0/9, 1.0/9, 1.0/36,
                     1.0/9, 1.0/36, 1.0/9,
                      1.0/36, 1.0/9, 1.0/36};
  const int XDIM = 5, YDIM = 5, NUM_WEIGHTS = 9;
  const vector<double> weight(weights, weights + NUM_WEIGHTS);
  Lattice3D lboltz(boost::extents[XDIM][YDIM][NUM_WEIGHTS]);

  for(int y = YDIM - 1; y >= 0; y--) {
    for(int x = 0; x < XDIM; x++) {
      for(int w = 0; w < NUM_WEIGHTS; w++) {
        lboltz[x][y][w] = 0;
      }
    }
  }

  lboltz[XDIM/2][YDIM/2][1] = 1;
  printLattice(lboltz, XDIM, YDIM, 1);
  cout << endl;

  for(int t = 0; t < 10; t++) {
    streamingUpdate(lboltz, XDIM, YDIM, NUM_WEIGHTS);
    printLattice(lboltz, XDIM, YDIM, 1);
    cout << endl;
  }

  return 0;
}

void streamingUpdate(Lattice3D &lb, const int XDIM, const int YDIM, 
    const int NUM_WEIGHTS) {
  Lattice3D temp_lattice(boost::extents[XDIM][YDIM][NUM_WEIGHTS]);

  for(int x = 0; x < XDIM; x++) {
    for(int y = 0; y < YDIM; y++) {
      for(int w = 0; w < NUM_WEIGHTS; w++) {
        Eigen::Vector2d dr(dirToSteps(w));
        // Periodic boundary conditions, for the time being
        // off the edge
        int xprime = x + int(dr[0]), yprime = y + int(dr[0]);
        xprime = mod(xprime, XDIM); yprime = mod(yprime, YDIM);

        temp_lattice[xprime][yprime][w] = lb[x][y][w];
      }
    }
  }

  lb = temp_lattice;

  return;
}

void collisionUpdate(const int XDIM, const int YDIM, 
    const int NUM_WEIGHTS, const vector<double> &weight, Lattice3D &lb) {
  const double m = 1, c = 1, tau = 1;

  for(int x = 0; x < XDIM; x++) {
    for(int y = 0; y < YDIM; y++) {

      Eigen::Vector2d macro_vel(0,0);
      double density = 0;
      for(int w = 0; w < NUM_WEIGHTS; w++) {
        macro_vel += c*lb[x][y][w]*dirToSteps(w);
        density += lb[x][y][w];
      }
      if (density != 0) macro_vel /= density;

      for(int w = 0; w < NUM_WEIGHTS; w++) {
        double equilibrium = 0;

        double e_dot_u = dirToSteps(w).dot(macro_vel);
 
        equilibrium = 1 + 3/c*e_dot_u + 9/(2*pow(c,2))*pow(e_dot_u,2)
            - 3/(2*pow(c,2))*macro_vel.dot(macro_vel);
        equilibrium *= weight[w]*density;

        lb[x][y][w] = lb[x][y][w] - 1/tau*(lb[x][y][w] - equilibrium);
      }
    }
  }

}


/**
 * \brief Convert a direction index into displacements on the lattice
 *
 * \param[in]  w  Direction value to convert
 * \param[out] dx x shift on lattice
 * \param[out] dy y shift on lattice
 *
 * \return Nothing
 */
Eigen::Vector2d dirToSteps(const int w) {
  switch(w) {
    default:
    case 0:
      return Eigen::Vector2d(0, 0);
      break;
    case 1:
      return Eigen::Vector2d(1, 0);
      break;
    case 2:
      return Eigen::Vector2d(1, 1);
      break;
    case 3:
      return Eigen::Vector2d(0, 1);
      break;
    case 4:
      return Eigen::Vector2d(-1, 1);
      break;
    case 5:
      return Eigen::Vector2d(-1, 0);
      break;
    case 6:
      return Eigen::Vector2d(-1, -1);
      break;
    case 7:
      return Eigen::Vector2d(0, -1);
      break;
    case 8:
      return Eigen::Vector2d(1, -1);
      break;
  }
}

void printLattice(const Lattice3D &lb, const int XDIM, const int YDIM, const int w) {
  for(int y = YDIM - 1; y >= 0; y--) {
    for(int x = 0; x < XDIM; x++) {
      cout << lb[x][y][w] << " ";
    }
    cout << endl;
  }
}

double latticeDistance(int direction) {
    if (direction == 0) {
      return 0;
    } else if (1 <= direction && direction <= 4) {
      return 1;
    } else if (4 < direction && direction <= 8) {
      return sqrt(2.0);
    }
    return -1;
}

int mod(int a, int b) {
  int ret = a % b;
  if (ret < 0)
    ret += b;
  return ret;
}
