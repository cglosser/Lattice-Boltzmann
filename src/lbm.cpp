#include <boost/multi_array.hpp>
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

typedef boost::multi_array<double, 3> Lattice3D;
typedef boost::multi_array<double, 2> Lattice2D;

void streamingUpdate(Lattice3D &, const int, const int, const int);
void collisionUpdate(const int, const int, const int, const vector<double> &,
    Lattice3D &);
void dirToSteps(const int, int &, int &);
void printLattice(const Lattice3D &, const int, const int, const int);

int main() {
  // Define weightings for a D2Q9 lattice.
  double weights[] = {4.0/9,  1.0/9,  1.0/9, 
                      1.0/9,  1.0/9,  1.0/36, 
                      1.0/36, 1.0/36, 1.0/36}; 

  const int XDIM = 3, YDIM = 3, NUM_WEIGHTS = 9;
  const vector<double> weight(weights, weights + NUM_WEIGHTS);
  Lattice3D lboltz(boost::extents[XDIM][YDIM][NUM_WEIGHTS]);

  for(int y = YDIM - 1; y >= 0; y--) {
    for(int x = 0; x < XDIM; x++) {
      for(int w = 0; w < NUM_WEIGHTS; w++) {
        if(x == 1 && y == 1) {
          lboltz[x][y][w] = 0;
        }
      }
    }
  }

  return 0;
}

void streamingUpdate(Lattice3D &lb, const int XDIM, const int YDIM, 
    const int NUM_WEIGHTS) {
  Lattice3D temp_lattice(boost::extents[XDIM][YDIM][NUM_WEIGHTS]);

  for(int x = 0; x < XDIM; x++) {
    for(int y = 0; y < YDIM; y++) {
      for(int w = 0; w < NUM_WEIGHTS; w++) {
        int dx = 0, dy = 0;
        dirToSteps(w, dx, dy);
        // Free boundary conditions, for the time being - particles just fall
        // off the edge
        if ((x + dx) >= XDIM || (x + dx) < 0) continue;
        if ((y + dy) >= YDIM || (y + dy) < 0) continue;

        temp_lattice[x + dx][y + dy][w] = lb[x][y][w];
      }
    }
  }

  lb = temp_lattice;

  return;
}

void collisionUpdate(const int XDIM, const int YDIM, 
    const int NUM_WEIGHTS, const vector<double> &weight, Lattice3D &lb) {
  Lattice2D density(boost::extents[XDIM][YDIM]);
  Lattice3D equilibrium(boost::extents[XDIM][YDIM][NUM_WEIGHTS]);

  for(int x = 0; x < XDIM; x++) {
    for(int y = 0; y < YDIM; y++) {
      density[x][y] = 0;
      for(int w = 0; w < NUM_WEIGHTS; w++) {
        density[x][y] += lb[x][y][w];
      }
    }
  }

  for(int x = 0; x < XDIM; x++) {
    for(int y = 0; y < ydim; y++) {
      for(int w = 0; w < NUM_WEIGHTS; w++) {
        double fi = lb[x][y][w], rho = density[x][y], u_dot_u = 0;

        for(int w = 0; w < NUM_WEIGHTS; w++) {
          u_dot_u += pow(lb[x][y][w], 2)*pow(latticeDistance(w),2);
        }

        equilibrium[x][y][w] = rho*weight[w](1 + 
            3*fi/rho*pow(latticeDistance(w),2) +
            9*pow(fi,2)/(2*pow(rho,2))*pow(latticeDistance(w),4) -
            3/(2*pow(rho,2))*u_dot_u);
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
void dirToSteps(const int w, int &dx, int &dy) {
  switch(w) {
    case 0:
      dx = 0; dy = 0;
      break;
    case 1:
      dx = 1; dy = 0;
      break;
    case 2:
      dx = 1; dy = 1;
      break;
    case 3:
      dx = 0; dy = 1;
      break;
    case 4:
      dx = -1; dy = 1;
      break;
    case 5:
      dx = -1; dy = 0;
      break;
    case 6:
      dx = -1; dy = -1;
      break;
    case 7:
      dx = 0; dy = -1;
      break;
    case 8:
      dx = 1; dy = -1;
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


