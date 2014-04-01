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
      os << sites_[y*XDIM_ + x].density();
      if(x != XDIM_ - 1) os << ",";

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

  for(int idx = 0; idx < NUM_SITES_; ++idx) {
    Eigen::Vector2i r(idx2coord(idx));
    if(r[0] == 0 || r[0] == XDIM_ - 1) {// || r[1] == 0 || r[1] == YDIM_ - 1) {
      auto &site = sites_[idx];
      site.state = NodeState::INACTIVE;
      std::fill(site.f_density.begin(), site.f_density.end(), 0);
    }
  }
}

void Lattice::streamingUpdate_() {
  for(int idx = 0; idx < NUM_SITES_; ++idx) {
    d2q9Node &site = sites_[idx];
    if(site.state == NodeState::INACTIVE) continue;

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
  for(int idx = 0; idx < NUM_SITES_; ++idx) {
    d2q9Node &site     = sites_[idx];
    if(site.state == NodeState::INACTIVE) continue;
    vector<double> &fs = site.f_density;
    Eigen::Vector2i r  = idx2coord(idx);
    if(r[0] == 1) {
      double x_velocity = 0.006666667,
             rho0       = (fs[4] + fs[7] + fs[1] +2*(fs[3] + fs[0] + fs[6]))/(1 + x_velocity),
             ru         = rho0*x_velocity;
      fs[5] = fs[3] + (2.0/3)*ru;
      fs[8] = fs[0] + (1.0/6)*ru - 0.5*(fs[7] - fs[1]);
      fs[2] = fs[6] + (1.0/6)*ru - 0.5*(fs[1] - fs[7]);
    } else if (r[0] == XDIM_ - 2) {
      double rho0       = 2,
             x_velocity = -1 + (fs[4] + fs[7] + fs[1] + 2*(fs[5] + fs[8] + fs[2]))/rho0,
             ru         = rho0*x_velocity;
      fs[3] = fs[5] - (2.0/3)*ru;
      fs[0] = fs[8] + (1.0/6)*ru + 0.5*(fs[7] - fs[1]);
      fs[6] = fs[2] + (1.0/6)*ru + 0.5*(fs[1] - fs[7]);
    } else {
      double density = site.density();
      Eigen::Vector2d macro_vel(site.velocity());
      for(int dir = 0; dir < NUM_WEIGHTS_; ++dir) {
        double e_dot_u =
          d2q9Node::d2q9_directions[dir].cast<double>().dot(macro_vel);
        double equilib = (density + 3.0*e_dot_u + (9.0/2)*pow(e_dot_u,2) -
            (3.0/2)*macro_vel.squaredNorm())*d2q9Node::d2q9_weights[dir];

        site.f_density[dir] -= (site.f_density[dir] - equilib)/tau;
      }
    }
  }

  return;
}

Eigen::Vector2i directionToSteps(const int n) {
  return d2q9Node::d2q9_directions[n];
}
