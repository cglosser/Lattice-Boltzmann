#include "lattice.h"
using namespace std;

int main() {
  Lattice myLat(50,50);

  for(int t = 0; t < 5000; t++) {
    myLat.print(cout);
    //cout << endl;

    myLat.update();
  }
  return 0;
}
