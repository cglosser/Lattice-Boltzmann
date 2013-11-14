#include "lattice.h"
using namespace std;

#define tau 1

double d2q9_weights[] = {1.0/36, 1.0/9, 1.0/36, 1.0/9, 4.0/9, 1.0/9, 1.0/36, 1.0/9, 1.0/36};

Lattice::Lattice(const int x_size, const int y_size):
    XDIM(x_size), YDIM(y_size), NUM_SITES(XDIM*YDIM), NUM_WEIGHTS(9), 
    weight(d2q9_weights, d2q9_weights + NUM_WEIGHTS) {
  f_density.resize(
      boost::extents[range(-2, NUM_SITES)][range(1, NUM_WEIGHTS+1)]);
  push_density.resize(
      boost::extents[range(-2, NUM_SITES)][range(1, NUM_WEIGHTS+1)]);
  neighbors.resize(
      boost::extents[range(-2, NUM_SITES)][range(1, NUM_WEIGHTS+1)]);

  for(int site = 0; site < NUM_SITES; site++) {
    for(int n = 1; n <= NUM_WEIGHTS; n++) {
      f_density[site][n] = 1;
    }
  }
  f_density[NUM_SITES/2][6] = 3;

  buildNeighbors();

  return;
}

double Lattice::density(const int site) {
  double rho = 0;
  for(int n = 1; n <= NUM_WEIGHTS; n++) rho += f_density[site][n];
  return rho;
}

void Lattice::update() {
  streamingUpdate();
  collisionUpdate();
  return;
}

void Lattice::print(ostream &os) {
  for(int site = 0; site < NUM_SITES; site++) {
    os << density(site);
    if(site % XDIM == 0) {
      cout << endl;
    } else {
      cout << ",";
    }
  }
  return;
}

Eigen::Vector2i Lattice::idx2coord(const int idx) {
  int x = idx % XDIM, y = idx/XDIM;
  return Eigen::Vector2i(x, y);
}

int Lattice::coord2idx(Eigen::Vector2i r) {
  return r[0] + r[1]*XDIM;
}

void Lattice::buildNeighbors() {
  for(int site = 0; site < NUM_SITES; site++) {
    Eigen::Vector2i r(idx2coord(site));
    for(int n = 1; n <= NUM_WEIGHTS; n++) {
      Eigen::Vector2i dr = directionToSteps(n), rprime = r + dr;

      if(rprime[0] < 0) {
        rprime += Eigen::Vector2i(XDIM,0);
      } else if(rprime[0]>= XDIM) {
        rprime -= Eigen::Vector2i(XDIM,0);
      }

      neighbors[coord2idx(r)][n] = coord2idx(rprime); 
    }
  }

  for(int site = 0; site < NUM_SITES; site++) {
    cout << site << endl;
    for(int n = 1; n <= NUM_WEIGHTS; n++) {
      cout << neighbors[site][n] << "  ";
    }
    cout << endl << endl;
  }
  return;
}

void Lattice::streamingUpdate() {
  for(int site = 0; site < NUM_SITES; site++) {
    for(int n = 1; n <= NUM_WEIGHTS; n++) {
      push_density[neighbors[site][n]][n] = f_density[site][n];
    }
  }
  f_density = push_density;
  return;
}

void Lattice::collisionUpdate() {
  for(int site = 0; site < NUM_SITES; site++) {
    Eigen::Vector2d macro_vel(0,0);
    double density = 0;

    for(int n = 1; n <= NUM_WEIGHTS; n++) {
      macro_vel += f_density[site][n]*(directionToSteps(n).cast<double>());
      density += f_density[site][n];
    }

    if(density != 0) macro_vel /= density;

    for(int n = 1; n <= NUM_WEIGHTS; n++) {
      double e_dot_u = directionToSteps(n).cast<double>().dot(macro_vel);
      double equilibrium = (1 + 3*e_dot_u + (9.0/2)*pow(e_dot_u,2)
        - (3.0/2)*macro_vel.squaredNorm())*weight[n - 1]*density; // c = 1
               //weight is a std::vector here^, hence n - 1
        
      f_density[site][n] -= 1/tau*(f_density[site][n] - equilibrium);
    }
  }

  return;
}


Eigen::Vector2i directionToSteps(const int n) {
  if (n < 1 || n > 9 ) {
    throw domain_error("Direction out of bounds");
  }
  int count = 0;
  Eigen::Vector2i result;
  for(int dy = 1; dy >= -1; dy--) {
    for(int dx = -1; dx <= 1; dx++) {
      if (++count == n)
        result = Eigen::Vector2i(dx, dy);
    }
  }
  return result;
}
