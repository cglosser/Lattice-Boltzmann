#ifndef LATTICE_H
#define LATTICE_H

#include "boost/multi_array.hpp"
#include <iostream>
#include <cmath>

typedef boost::multi_array<double, 3> array3D;

class Lattice {
 public:

  Lattice(int, int);
  double density(int, int);
  void update();
  void print(std::ostream &);

 private:

  const int XDIM, YDIM, NUM_WEIGHTS;
  const std::vector<double> weight;
  array3D f_density, push_density;
  boost::multi_array<double*, 3> neighbors;
  void buildNeighbors();
  void streamingUpdate();


};

void directionToSteps(const int, int &, int&);

#endif
