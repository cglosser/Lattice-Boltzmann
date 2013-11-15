#include "lattice.h"
using namespace std;

int main() {
  Lattice myLat(5,5);


  for(int s = 0; s < myLat.NUM_SITES; s++) {
    cout << s << endl;
    for(int n = 1; n < myLat.NUM_WEIGHTS; s++) {
      cout << myLat.neighbors[s][n] << "  ";
    }
    cout << endl;
  }
  return 0;
}
