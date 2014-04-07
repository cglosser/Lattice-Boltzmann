#include "lattice.h"
using namespace std;

#define tau 1.0

Lattice::Lattice(int x_size, int y_size, double mV, double Re): XDIM_(x_size),
    YDIM_(y_size), NUM_SITES_(XDIM_*YDIM_), NUM_WEIGHTS_(9), maxVelocity_(mV), Re_(Re) { 
  sites_ = std::vector<d2q9Node>(NUM_SITES_);
  nu_    = 2*maxVelocity_*YDIM_/Re_;
  omega_ = 1.0/(3*nu_ + 0.5);

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

double Lattice::poiseuilleX(int x) {
  double radius = 0.5*(XDIM_ - 1);
  return -maxVelocity_*(pow(std::abs(1 - x/radius), 2) - 1);
}

double Lattice::poiseuilleY(int y) {
  double radius = 0.5*(YDIM_ - 1);
  return -maxVelocity_*(pow(std::abs(1 - y/radius), 2) - 1);
}

void Lattice::update() {
  inletOutlet_();
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

void Lattice::setStates_() {
  for(int idx = 0; idx < NUM_SITES_; ++idx) {
    d2q9Node &site = sites_[idx];
    Eigen::Vector2i r(idx2coord(idx));
    Eigen::Vector2d u0(poiseuilleY(r[1]), 0);

    if(r[1] == 0 || r[1] == YDIM_ - 1) {
      site.state = NodeState::WALL;
      std::fill(site.f_density.begin(), site.f_density.end(), 0);
    } else {
      std::vector<double> equil_f = site.equilibrium(1, u0);
      std::copy(equil_f.begin(), equil_f.end(), site.f_density.begin());
    }
  }
}

void Lattice::print(ostream &os) {
  /** This really odd indexing is to have xhat and yhat correspond to
   * traditional cartesian directions, i.e., xhat is right, yhat is up.*/
  for(int y = YDIM_ - 1; y >= 0; --y) {
    for(int x = 0; x < XDIM_; ++x) {
      os << sites_[y*XDIM_ + x].velocity()[0];
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

void Lattice::streamingUpdate_() {
  for(int idx = 0; idx < NUM_SITES_; ++idx)
  {
    d2q9Node &site = sites_[idx];
    if(site.state == NodeState::WALL) continue;
    for(int dir = 0; dir < NUM_WEIGHTS_; ++dir) {
      d2q9Node &neighbor = sites_[neighborhood_[idx][dir]];
      if(neighbor.state != NodeState::WALL)
        neighbor.temp_density[dir] = site.f_density[dir];
      else 
        site.temp_density[site.reverse(dir)] = site.f_density[dir];
    }
  }

  //Finish off simultaneous update
  for(auto &site : sites_) {
    if(site.state == NodeState::ACTIVE)
      site.f_density.swap(site.temp_density);
  }

  return;
}

void Lattice::collisionUpdate_() {
  for(int idx = 0; idx < NUM_SITES_; ++idx) {
    d2q9Node &site = sites_[idx];
    if(site.state == NodeState::WALL) continue;
    //Eigen::Vector2i r = idx2coord(idx);
    double density = site.density();
    Eigen::Vector2d velocity = site.velocity();
      for(int dir = 0; dir < NUM_WEIGHTS_; ++dir)
      {
        double e_dot_u = site.d2q9_directions[dir].cast<double>().dot(velocity);
        double equilib = (density + 3.0*e_dot_u + (9.0/2)*pow(e_dot_u,2) -
            (3.0/2)*velocity.squaredNorm())*d2q9Node::d2q9_weights[dir];


        site.f_density[dir] += (equilib - site.f_density[dir])/tau;
      }

  }

  return;
}

void Lattice::inletOutlet_() {
  for(int idx = 0; idx < NUM_SITES_; ++idx) {
    if      (sites_[idx].state == NodeState::INLET)   inletZou_(idx);
    else if (sites_[idx].state == NodeState::OUTLET) outletZou_(idx);
  }
}

void Lattice::inletZou_(int idx)
{
  d2q9Node &site = sites_[idx];
  std::vector<double> &fi = site.f_density;
  Eigen::Vector2i r(idx2coord(idx));
  Eigen::Vector2d u0(poiseuilleY(r[1]), 0);
  //Assume a western wall
  double fint  = fi[0] + fi[2] + fi[4],
         fint2 = fi[3] + fi[6] + fi[7],
         rho   = (fint + 2*fint2)/(1 - u0[0]);
  
  double fdiff = 0.5*(fi[2] - fi[4]),
         rhoUx = rho*u0[0]/6.0,
         rhoUy = 0.5*rho*u0[1];

  fi[1] = fi[3] + 4*rhoUx;
  fi[5] = fi[7] - fdiff + rhoUx + rhoUy;
  fi[8] = fi[6] + fdiff + rhoUx - rhoUy;
}

void Lattice::outletZou_(int idx)
{
  d2q9Node &site = sites_[idx];
  std::vector<double> &fi = site.f_density;
  Eigen::Vector2i r(idx2coord(idx));
  Eigen::Vector2d u0(poiseuilleY(r[1]), 0);
  //Assume an eastern wall
  double fint  = fi[0] + fi[2] + fi[4],
         fint2 = fi[1] + fi[8] + fi[5],
         rho   = (fint + 2*fint2)/(1 - u0[0]);
  
  double fdiff = 0.5*(fi[2] - fi[4]),
         rhoUx = rho * u0[0]/6,
         rhoUy = 0.5*rho*u0[1];
  fi[3] = fi[1] - 4*rhoUx;
  fi[7] = fi[5] + fdiff - rhoUx - rhoUy;
  fi[6] = fi[8] - fdiff - rhoUx + rhoUy;
}

Eigen::Vector2i directionToSteps(const int n) {
  return d2q9Node::d2q9_directions[n];
}
