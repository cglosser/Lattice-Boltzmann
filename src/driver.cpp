#include "lattice.h"
using namespace std;

int main() {
  Lattice myLat(200,20);

  for(int time = 0; time < 1500; time++) {
    myLat.update();
    myLat.print(cout);
  }

  return 0;
}
