#include "lattice.h"
using namespace std;

int main() {
  Lattice myLat(5,5);
  myLat.print(cout);

  for(int time = 0; time < 500; time++) {
    myLat.update();
    myLat.print(cout);
  }

  return 0;
}
