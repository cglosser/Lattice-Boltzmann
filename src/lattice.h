#ifndef LATTICE_H
#define LATTICE_H

#include "boost/multi_array.hpp"
#include <cmath>
#include <Eigen/Dense>
#include <iostream>

/**
 * @brief Define a 2d + velocities array for D2Q9 lattice geometries.
 */
typedef boost::multi_array<double, 3> array3D;

  /**
   * @brief   Standard two-dimensional Boltzmann lattice.
   * @details Maintains a record of the density at each site (for each velocity
   * direction) and each site's neighbors as separate arrays. Also evolves the
   * system in time through an update method. 
   * @todo Implement subclassing to allow D2Q7, D3Q* lattices
   * @todo Swap multi-dimensional arrays for vectors of density arrays
   */
class Lattice {
 public:

  /**
   * @brief Constructs a lattice object given x & y dimensions. Also builds an
   * appropriate neighbor list with reference entries. */
  Lattice(const int, const int);

  /** @brief Compute and return the particle density at site (x,y) as the sum
   * of the Boltzmann densities. */
  double density(const int, const int);

  /** @brief Perform a stream & collision update to evolve the system by one
   * timestep. */
  void update();

  /** 
   * @brief   Write the lattice to the specified output stream.  
   * @details Loops over nodes in row-major order, printing node particle
   * densities on the specified stream delimited by commas with newlines at the
   * end of a row. Generally requires reshaping/partitioning to produce a
   * visualization.
   * @todo Consider rewriting to overload stream output operators.
   */
  void print(std::ostream &os);

 private:

  const int XDIM, YDIM, NUM_WEIGHTS;

  /** @brief Directioinal weights for the equilibrium value of the D2Q9 lattice 
   */
  const std::vector<double> weight; 

  /** @brief Array containing the Boltzmann densities for every node in the
   * lattice. */
  array3D f_density;

  /** 
   * @brief   Array used to hold temporary Boltzmann density values during the
   * streaming update.
   * @details This array exists only to allow the simultaneous update during
   * streaming. It goes in scope at the instantiation of the class so that the
   * neighbor list can contain valid references to equivalent lattice
   * locations.
   */
  array3D  push_density;

  /** @brief 3d array of references to equivalent lattice locations in
   * push_density */
  boost::multi_array<double*, 3> neighbors;

  /**
   * @brief Generate an array of references to neighboring sites in a copy of
   * the canonical lattice.*/
  void buildNeighbors();

  /**
   * @brief Perform the streaming step of an update wherein lattice values move
   * along their velocity components to adjacent neighbors.
   * @details The values must stream simultaneously, in that they move to a
   * secondary lattice prior to getting copied back to the original.
   * @todo Implement detection of NULL neighbor references to handle
   * arbitrary-direction bouncebacks.
   */
  void streamingUpdate();

  /**
   * @brief   Perform the collision step of an update
   * @details Approximates the collisions of particles at each site on the
   * lattice, tending to push the system towards an equilibrium value. The form
   * of the equilibrium value is highly dependent on the lattice geometry.
   */
  void collisionUpdate();

};

/**
 * @brief   Convert a direction number to +-1 integer steps on a lattice.
 * @details Uses rounding of trigonometric functions to give unit steps for a specified direction. Because of the trigonometry, this function will be MUCH slower than necessary and should be avoided for anything except initialization.
 * @param[in] n direction
 * @param[out] dx step in the x direction
 * @param[out] dy step in the y direction
 * @todo Reimplement with a more efficient ordering. 1:NUM_STEPS with a
 * top-down looks more promising than using 0:NUM_STEPS-1 counter clockwise
 */
void directionToSteps(const int n, int &dx, int&dy);

#endif
