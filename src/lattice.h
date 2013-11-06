#ifndef LATTICE_H
#define LATTICE_H

#include "boost/multi_array.hpp"
#include <iostream>
#include <cmath>

typedef boost::multi_array<double, 3> array3D;

class Lattice {
 public:

  Lattice(int, int);

 private:

  int XDIM, YDIM;
  array3D f_density, push_density;

};

void directionToSteps(const int, int &, int&);

#endif
