#include "lattice.h"
using namespace std;

Lattice::Lattice(int x_size, int y_size) {
  XDIM = x_size; YDIM = y_size;
  f_density.resize(boost::extents[XDIM][YDIM][9]);
  push_density.resize(boost::extents[XDIM][YDIM][9]);
  neighbors.resize(boost::extents[XDIM][YDIM][9]);

  for(int x = 0; x < XDIM; x++) {
    for(int y = 0; y < YDIM; y++) {
      for(int n = 0; n < 9; n++) {
        f_density[x][y][n] = 1;
      }
    }
  }

  buildNeighbors();

  return;
}

void Lattice::buildNeighbors() {
  for(int x = 0; x < XDIM; x++) {
    for(int y = 0; y < YDIM; y++) {
      for(int n = 0; n < 9; n++) {
        int dx, dy;
        directionToSteps(n, dx, dy); // Set the values of dx & dy from n

        int xprime = x + dx, yprime = y + dy, nprime = n;
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


void directionToSteps(const int n, int &dx, int &dy) {
  if (n == 0) {
    dx = dy = 0;
  } else {
    int r = n - 1;
    double theta = r*atan(1.0); // equals 2*pi*r/8
    dx = round(cos(theta)); 
    dy = round(sin(theta));
  }
}
        



