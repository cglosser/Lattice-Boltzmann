#ifndef LATTICE_H
#define LATTICE_H

#include "boost/multi_array.hpp"
#include <cmath>
#include <Eigen/Dense>
#include <iostream>

/**
 * @brief Define a 2d + velocities array for D2Q9 lattice geometries.
 */
typedef boost::multi_array<double, 2> array2D;

/** @brief Simplify a boost::multi_array range object. The directional
 * components of the lattice are most easily indexed through [1:NUM_WEIGHTS +
 * 1), so it is convenient for the multi_array to assume the same range.
 */
typedef boost::multi_array_types::extent_range range;

/**
 * @brief   Standard two-dimensional Boltzmann lattice.
 * @details Maintains a record of the density at each site (for each velocity
 * direction) and each site's neighbors as separate arrays. Also evolves the
 * system in time through an update method. 
 * @todo Implement subclassing to allow D2Q7, D3Q* lattices
 */
class Lattice {
 public:

  /** @brief Constructs a lattice object given x & y dimensions. Also builds an
   * appropriate neighbor list with integer entries corresponding to site
   * indices */
  Lattice(const int, const int);

  /** @brief Compute and return the particle density at site i as the sum
   * of the Boltzmann densities. */
  double density(const int);

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

  /**
   * @brief Convert an index to coordinates on the lattice.
   */
  Eigen::Vector2i idx2coord(const int);

  /**
   * @brief Convert lattice coordinates to a list index.
   */
  int coord2idx(const Eigen::Vector2i);

 private:

  enum state {INACTIVE, ACTIVE};

  const int XDIM,        ///< Length of lattice in the x dimension
            YDIM,        ///< Length of lattice in the y dimension
            NUM_SITES,   ///< Total number of sites
            NUM_WEIGHTS; ///< Number of discretized momentum directions

  /** @brief Directioinal weights for the equilibrium value of the D2Q9 lattice 
   */
  const std::vector<double> weight; 

  /** @brief Array containing the Boltzmann densities for every node in the
   * lattice. */
  array2D f_density;

  /** 
   * @brief   Array used to hold temporary Boltzmann density values during the
   * streaming update.
   * @details This array exists only to allow the simultaneous update during
   * streaming. It goes in scope at the instantiation of the class so that the
   * neighbor list can contain valid references to equivalent lattice
   * locations.
   */
  array2D  push_density;

  /** @brief Array of indices to neighboring sites. The second dimension runs
   * from [1, NUM_WEIGHTS].
   */
  boost::multi_array<int, 2> neighbors;

  boost::multi_array<state, 1> node_state;

  /** @brief Generate an array integers corresponding to neighbor nodes. By
   * default, assumes periodic boundary conditions in x & y -- rigid boundaries
   * must be added by deactivating edge nodes.
   */
  void buildNeighbors();

  /** @brief Set the state of all nodes in the lattice.
   * @todo Implement IO to read in lattice boundaries from a file.
   */
  void setStates();

  /** @brief Perform the streaming step of an update wherein lattice values
   * move along their velocity components to adjacent neighbors.
   * @details The values must stream simultaneously, in that they move to a
   * secondary lattice prior to getting copied back to the original.  @todo
   * Implement detection of inactive neighbors to handle arbitrary-direction
   * bouncebacks.
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
 * @details Converts an integer in the range [1, NUM_WEIGHTS] to unit steps on
 * the lattice. For the D2Q9 lattice, the directions are as follows: \n
 * <pre>
 * 1   2   3
 *   \ | /  
 * 4 - 5 - 6
 *   / | \
 * 7   8   9
 * </pre>
 * so that the reverse of any given direction becomes simply
 * (NUM_WEIGHTS + 1) - n
 * @throws  domain_error if the input does not correspond to one of the nearest
 * neighbor directions
 * @returns A two dimensional integer containing steps along x and y. This
 * vector must be cast to a double for general purpose computations!
 *
 * @todo Buld into Lattice class/subclasses (steps on the lattice are part of
 * the lattice geometry and thus this function should be a class member.)
 */
Eigen::Vector2i directionToSteps(const int n);

#endif
