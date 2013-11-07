#include "lattice.h"
using namespace std;

int main() {
  Lattice myLat(20, 20);

  cout << myLat.neighbors[1][1][5] << endl;
  cout << &myLat.push_density[0][1][5] << endl;
  cout << endl;

  cout << myLat.neighbors[0][0][5] << endl;
  cout << &myLat.push_density[19][0][5] << endl;
  cout << endl;

  cout << myLat.neighbors[3][3][2] << endl;
  cout << &myLat.push_density[4][4][2] << endl;
  cout << endl;


  return 0;
}
