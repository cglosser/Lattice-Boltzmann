#include "lattice.h"
using namespace std;

#define tau 1

Lattice::Lattice(const int x_size, const int y_size): XDIM_(x_size),
    YDIM_(y_size), NUM_SITES_(XDIM_*YDIM_), NUM_WEIGHTS_(9) { 
  sites_ = std::vector<d2q9Node>(XDIM_*YDIM_);

  buildNeighbors_();
  setStates_();

  return;
}

double Lattice::density(const unsigned x, const unsigned y) {
  
  return sites_[coord2idx(Eigen::Vector2i(x, y))].density();
}

Eigen::Vector2d Lattice::velocity(const unsigned x, const unsigned y) {
  return sites_[coord2idx(Eigen::Vector2i(x, y))].velocity();
}

void Lattice::update() {
  streamingUpdate_();
  collisionUpdate_();
  return;
}

void Lattice::buildNeighbors_() {
  neighborhood_.reserve(XDIM_*YDIM_);
  for(int i = 0; i < NUM_SITES_; ++i) {
    Eigen::Vector2i r = idx2coord(i);
    vector<int> nbd; nbd.reserve(9);

    for(auto dir : d2q9Node::d2q9_directions) {
      Eigen::Vector2i rprime = r + dir;
      if(rprime[0] < 0)           rprime[0] += XDIM_;
      else if(rprime[0] >= XDIM_) rprime[0] -= XDIM_;
      
      if(rprime[1] < 0)           rprime[1] += YDIM_;
      else if(rprime[1] >= YDIM_) rprime[1] -= YDIM_;
      
      nbd.push_back(coord2idx(rprime));
    }
    neighborhood_.push_back(nbd);
  }
}

void Lattice::print(ostream &os) {
  /** This really odd indexing is to have xhat and yhat correspond to
   * traditional cartesian directions, i.e., xhat is right, yhat is up.*/
  for(int y = YDIM_ - 1; y >= 0; --y) {
    for(int x = 0; x < XDIM_; ++x) {
      os << sites_[y*XDIM_ + x].density() << ", ";
    }
    os << std::endl;
  }
}

Eigen::Vector2i Lattice::idx2coord(const int idx) {
  int x = idx % XDIM_, y = idx/XDIM_;
  return Eigen::Vector2i(x, y);
}

int Lattice::coord2idx(Eigen::Vector2i r) {
  return r[0] + r[1]*XDIM_;
}

void Lattice::setStates_() {
  for(auto &it : sites_) {
    it.state = NodeState::ACTIVE;
    std::copy(d2q9Node::d2q9_weights.begin(), d2q9Node::d2q9_weights.end(),
        it.f_density.begin());
  }
  sites_[2].f_density[7] *= 1.1;
}

void Lattice::streamingUpdate_() {
  for(int idx = 0; idx < NUM_SITES_; ++idx) {
    d2q9Node &site = sites_[idx];
    for(int dir = 0; dir < NUM_WEIGHTS_; ++dir) {
      d2q9Node &neighbor = sites_[neighborhood_[idx][dir]];
      if(neighbor.state == NodeState::ACTIVE)
        neighbor.temp_density[dir] = site.f_density[dir];
      else if(neighbor.state == NodeState::INACTIVE)
        site.temp_density[(NUM_WEIGHTS_ - 1) - dir] = site.f_density[dir];
    }
  }

  //Finish off simultaneous update
  for(auto &site : sites_) site.f_density.swap(site.temp_density);

  return;
}

void Lattice::collisionUpdate_() {
  for(auto &site : sites_) {
    Eigen::Vector2d macro_vel(site.velocity());
    double density = site.density();
    for(int dir = 0; dir < NUM_WEIGHTS_; ++dir) {
      double e_dot_u =
        d2q9Node::d2q9_directions[dir].cast<double>().dot(macro_vel);
      double equilib = (1 + 3.0*e_dot_u + (9.0/2)*pow(e_dot_u,2) -
          (3.0/2)*macro_vel.squaredNorm())*d2q9Node::d2q9_weights[dir]*density;

      site.f_density[dir] -= 1.0/tau*(site.f_density[dir] - equilib);
    }
  }

  return;
}

Eigen::Vector2i directionToSteps(const int n) {
  return d2q9Node::d2q9_directions[n];
}
