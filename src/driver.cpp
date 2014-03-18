#include "lattice.h"
using namespace std;

int main() {
  int dim = 30;
  Lattice myLat(3*dim,dim);

  myLat.print(cout);
  for(int time = 0; time < 3000; time++) {
    myLat.update();
    myLat.print(cout);
  }

  return 0;
}
