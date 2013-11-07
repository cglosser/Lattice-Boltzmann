#include "lattice.h"
using namespace std;

int main() {
  Lattice myLat(5, 5);

  for(int i = 0; i < 10; i++) {
    myLat.print(cout);
    myLat.update();
  }


  return 0;
}
