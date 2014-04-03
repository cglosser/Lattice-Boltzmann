#ifndef LATTICE_H
#define LATTICE_H

#include "node.h"
#include <cmath>
#include <iostream>
#include <iterator>

typedef std::vector<d2q9Node> nodeArray;

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
  Lattice(int, int, double, double);

  /** @brief Perform a stream & collision update to evolve the system by one
   * timestep. */
  void update();

  double density(const unsigned, const unsigned);
  Eigen::Vector2d velocity(const unsigned, const unsigned);

  double poiseuilleX(int);
  double poiseuilleY(int);

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

  const int XDIM_,        ///< Length of lattice in the x dimension
            YDIM_,        ///< Length of lattice in the y dimension
            NUM_SITES_,   ///< Total number of sites
            NUM_WEIGHTS_; ///< Number of discretized momentum directions
  const double maxVelocity_, Reynolds_;
  double nu_, omega_;

  nodeArray sites_;
  std::vector<std::vector<int>> neighborhood_;

  void buildNeighbors_();

  /** @brief Set the state of all nodes in the lattice.
   *  @todo Implement IO to read in lattice boundaries from a file.
   */
  void setStates_();

  /** @brief Perform the streaming step of an update wherein lattice values
   * move along their velocity components to adjacent neighbors.
   * @details The values must stream simultaneously, in that they move to a
   * secondary lattice prior to getting copied back to the original.  @todo
   * Implement detection of inactive neighbors to handle arbitrary-direction
   * bouncebacks.
   */
  void streamingUpdate_();

  /**
   * @brief   Perform the collision step of an update
   * @details Approximates the collisions of particles at each site on the
   * lattice, tending to push the system towards an equilibrium value. The form
   * of the equilibrium value is highly dependent on the lattice geometry.
   */
  void collisionUpdate_();
  
  void inletOutlet_();
  void inletZou_(int);
  void outletZou_(int);
};

/**
 * @brief   Convert a direction number to +-1 integer steps on a lattice.
 * @details Converts an integer in the range [1, NUM_WEIGHTS] to unit steps on
 * the lattice. For the D2Q9 lattice, the directions are as follows: \n
 * <pre>
 *
 * yhat
 *  ^
 *  |
 *  |
 *  |
 *
 * 6   2   5
 *   \ | /  
 * 3 - 0 - 1
 *   / | \
 * 7   4   8  xhat--->
 * </pre>
 * @returns A two dimensional integer vector containing steps along x and y.
 * This vector must be cast to a double for general purpose computations!
 *
 * @todo Buld into Lattice class/subclasses (steps on the lattice are part of
 * the lattice geometry and thus this function should be a class member.)
 */
Eigen::Vector2i directionToSteps(const int n);

#endif
