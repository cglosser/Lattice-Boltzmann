#include "lattice.h"
using namespace std;

int main() {
  Lattice myLat(5,5);

  for(int time = 0; time < 3000; time++) {
    myLat.print(cout);
    myLat.update();
  }
  myLat.print(cout);

  return 0;
}
