#include "lattice.h"
using namespace std;

int main() {
  Lattice myLat(35, 135);

  for(int t = 0; t < 200; t++) {
    myLat.print(cout);
    myLat.update();
  }

  return 0;
}
