#include <iostream>
#include <vector>
#include <boost/multi_array.hpp>

using namespace std;

typedef boost::multi_array<double, 3> Lattice;

void updateLattice(Lattice &, const int, const int, const int);
void dirToSteps(const int, int &, int &);

int main() {
  //double weights[] = {4.0/9,  1.0/36, 1.0/9, 
                      //1.0/36, 1.0/9,  1.0/36, 
                      //1.0/9,  1.0/36, 1.0/9}; 
  // Above are "Old" weights - check convention!

  /* Define weightings for a D2Q9 lattice. Follows the "mathematical"
   * convention for rotating around a circle, with weights[0] giving the
   * self-term. */
  double weights[] = {4.0/9,  1.0/9, 1.0/36, 
                      1.0/9, 1.0/36, 1.0/9, 
                      1.0/36, 1.0/9, 1.0/36};
  const int XDIM = 3, YDIM = 3, NUM_WEIGHTS = 9;
  const vector<double> weight(weights, weights + NUM_WEIGHTS);
  boost::multi_array<double, 3> lboltz(boost::extents[XDIM][YDIM][NUM_WEIGHTS]);

  int dir = 8;
  for(int y = YDIM - 1; y >= 0; y--) {
    for(int x = 0; x < XDIM; x++) {
      for(int w = 0; w < NUM_WEIGHTS; w++) {
        lboltz[x][y][w] = 0;
      }
    }
  }
  lboltz[0][2][dir] = 1;
  for(int y = YDIM - 1; y >= 0; y--) {
    for(int x = 0; x < XDIM; x++) {
      //for( int w = 0; w < NUM_WEIGHTS; w++ ) {
        cout << lboltz[x][y][dir] << " ";
      //}
    }
    cout << endl;
  }
  cout << endl;

  updateLattice(lboltz, XDIM, YDIM, NUM_WEIGHTS);

  for(int y = YDIM - 1; y >= 0; y--) {
    for( int x = 0; x < XDIM; x++ ) {
      //for( int w = 0; w < NUM_WEIGHTS; w++ ) {
        cout << lboltz[x][y][dir] << " ";
      //}
    }
    cout << endl;
  }
  cout << endl;

  updateLattice(lboltz, XDIM, YDIM, NUM_WEIGHTS);
  for(int y = YDIM - 1; y >= 0; y--) {
    for( int x = 0; x < XDIM; x++ ) {
      //for( int w = 0; w < NUM_WEIGHTS; w++ ) {
        cout << lboltz[x][y][dir] << " ";
      //}
    }
    cout << endl;
  }


  return 0;
}

void updateLattice(Lattice &lb, const int XDIM, const int YDIM, 
    const int NUM_WEIGHTS) {
  boost::multi_array<double, 3> 
    temp_lattice(boost::extents[XDIM][YDIM][NUM_WEIGHTS]);

  for(int x = 0; x < XDIM; x++) {
    for(int y = 0; y < YDIM; y++) {
      for(int w = 1; w < NUM_WEIGHTS; w++) {
        int dx, dy;
        dirToSteps(w, dx, dy);
        // Free boundary conditions, for the time being
        if ((x + dx) >= XDIM || (x + dx) < 0) continue;
        if ((y + dy) >= YDIM || (y + dy) < 0) continue;

        temp_lattice[x + dx][y + dy][w] = lb[x][y][w];
      }
    }
  }

  lb = temp_lattice;

  return;
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

