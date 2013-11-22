#include "lattice.h"
using namespace std;

int main() {
  int dim = 100;
  Lattice myLat(3*dim,dim);

  myLat.print(cout);
  for(int time = 0; time < 5000; time++) {
    myLat.update();
    if(time % 100 == 0) myLat.print(cout);
  }

  return 0;
}
