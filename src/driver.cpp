#include "lattice.h"
using namespace std;

int main() {
  int dim = 30;
  Lattice myLat(8*dim,dim);

  myLat.print(cout);
  for(int time = 0; time < 5000; time++) {
    myLat.update();
    myLat.print(cout);
  }

  return 0;
}
