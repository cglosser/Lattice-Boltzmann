#include "node.h"
#include <iostream>

using namespace std;

int main() {

  d2q9Node myNode;

  for(int i = 0; i < 9; i++)
    myNode.f_density[i] = (10*i);
  

  cout << myNode.density() << endl;
  cout << myNode.velocity().transpose() << endl;

  for(int i = 0; i < 9; ++i) {
    cout << i << " " << d2q9Node::reverse(i) << endl;
  }

  return 0;
}
