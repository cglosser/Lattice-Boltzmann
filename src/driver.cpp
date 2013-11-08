#include "lattice.h"
using namespace std;

int main() {
  Lattice myLat(101,21);

  for(int time = 0; time < 1000; time++) {
    myLat.print(cout);
    myLat.update();
  }
  myLat.print(cout);

  return 0;
}
