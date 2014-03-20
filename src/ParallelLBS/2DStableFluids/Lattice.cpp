#include "stdafx.h"
#include "Lattice.h"

#define tau 1

double d2q9_weights[] = {1.0/36, 1.0/9, 1.0/36, 1.0/9, 4.0/9, 1.0/9, 1.0/36, 1.0/9, 1.0/36};

Lattice::Lattice(const int x_size, const int y_size, const int threadCount = 1):
    XDIM(x_size), YDIM(y_size), NUM_SITES(XDIM*YDIM), NUM_WEIGHTS(9), 
    weight(d2q9_weights, d2q9_weights + NUM_WEIGHTS) {
  f_density.resize(
      boost::extents[range(-2, NUM_SITES)][range(1, NUM_WEIGHTS+1)]);
  push_density.resize(
      boost::extents[range(-2, NUM_SITES)][range(1, NUM_WEIGHTS+1)]);
  neighbors.resize(
      boost::extents[range(-2, NUM_SITES)][range(1, NUM_WEIGHTS+1)]);
  node_state.resize(
      boost::extents[range(-2, NUM_SITES)]);

  for(int site = 0; site < NUM_SITES; site++) {
    for(int n = 1; n <= NUM_WEIGHTS; n++) {
      f_density[site][n] = 2;
    }
  }
  //f_density[XDIM + 1][3] = 20;

  buildNeighbors();
  setStates();

  updateStates = vector<int>(threadCount - 1, 0);

  threads = vector<thread>();
  for (int i = 1; i < threadCount; i++)
  {
	  threads.push_back(thread(&Lattice::WorkerThread,
							  this,
							  i-1,
							  (i * NUM_SITES) / threadCount,
							  ((i + 1) * NUM_SITES) / threadCount));
  }

  return;
}

Lattice::~Lattice()
{
	done = true;

	for (unsigned i = 0; i < threads.size(); i++)
	{
		if (threads[i].joinable())
		{
			threads[i].join();
		}
	}
}

double Lattice::density(const int site) {
  double rho = 0;
  for(int n = 1; n <= NUM_WEIGHTS; n++) rho += f_density[site][n];
  return rho;
}

void Lattice::update() {
  preStreamingUpdate();
  streamingUpdate(0, NUM_SITES / (threads.size() + 1));
  postStreamingUpdate();
  preCollisionUpdate();
  collisionUpdate(0, NUM_SITES / (threads.size() + 1));
  postCollisionUpdate();
  return;
}

void Lattice::print(std::ostream &os) {
  for(int y = YDIM - 1; y >= 0; y--) {
    for(int x = 0; x < XDIM; x++) {
      int site = coord2idx(Eigen::Vector2i(x, y));
      os << density(site);
      if(x != XDIM - 1) os << ",";
    }
    os << std::endl;
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

void Lattice::setStates() {
  for(int site = 0; site < NUM_SITES; site++) {
    node_state[site] = ACTIVE;

    //Pipe conditions
    if(site < XDIM || site >= NUM_SITES - XDIM) {
      node_state[site] = INACTIVE;
    }

  }
}

void Lattice::SetDensity(int site, double newDensity)
{
	for (int i = 1; i <= NUM_WEIGHTS; i++)
	{
		f_density[site][i] = newDensity;
	}
}

void Lattice::Reset()
{
	for (int site = 0; site < NUM_SITES; site++) {
		for (int n = 1; n <= NUM_WEIGHTS; n++) {
			f_density[site][n] = 2;
		}
	}

	lastUpdateTime = std::chrono::system_clock::now();
}

void Lattice::SetLastUpdateTime()
{
	lastUpdateTime = std::chrono::system_clock::now();
}

double Lattice::LastTimeStep()
{
	return timeStep.count();
}

void Lattice::buildNeighbors() {
  for(int site = 0; site < NUM_SITES; site++) {
    Eigen::Vector2i r(idx2coord(site));

    for(int n = 1; n <= NUM_WEIGHTS; n++) {

      Eigen::Vector2i r = idx2coord(site), dr = directionToSteps(n),
                      rprime = r + dr;

      rprime[0] = (rprime[0] <  0    ? rprime[0] + XDIM : rprime[0]);
      rprime[0] = (rprime[0] >= XDIM ? rprime[0] - XDIM : rprime[0]);

      rprime[1] = (rprime[1] <  0    ? rprime[1] + XDIM : rprime[1]);
      rprime[1] = (rprime[1] >= XDIM ? rprime[1] - XDIM : rprime[1]);

      neighbors[site][n] = coord2idx(rprime);
    }
  }

  return;
}

void Lattice::preStreamingUpdate()
{
	for (unsigned i = 0; i < updateStates.size(); i++)
	{
		updateStates[i] = 1;
	}
}

void Lattice::streamingUpdate(int start, int end) {
  for(int site = start; site < end; site++) {
    if(node_state[site] == INACTIVE) continue;

    for(int n = 1; n <= NUM_WEIGHTS; n++) {
      int neighbor_idx = neighbors[site][n];

      if(node_state[neighbor_idx] == INACTIVE) {
        push_density[site][(NUM_WEIGHTS + 1) - n] = f_density[site][n];
      } else {
        push_density[neighbor_idx][n] = f_density[site][n];
      }
    }
  }
  return;
}

void Lattice::postStreamingUpdate()
{
	for (unsigned i = 0; i < updateStates.size(); i++)
	{
		while (updateStates[i] != 0){}
	}

	f_density = push_density;
}

void Lattice::preCollisionUpdate()
{
	for (unsigned i = 0; i < updateStates.size(); i++)
	{
		updateStates[i] = 2;
	}
}

void Lattice::collisionUpdate(int start, int end) {
  currentUpdateTime = std::chrono::system_clock::now();
  timeStep = currentUpdateTime - lastUpdateTime;
  for(int site = start; site < end; site++) {
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

void Lattice::postCollisionUpdate()
{
	for (unsigned i = 0; i < updateStates.size(); i++)
	{
		while (updateStates[i] != 0){}
	}

	lastUpdateTime = currentUpdateTime;
}

void Lattice::WorkerThread(int homeAddress, int start, int end)
{
	while (!done)
	{
		if (updateStates[homeAddress] == 0)
		{
			std::this_thread::yield();
			continue;
		}

		if (updateStates[homeAddress] == 1)
		{
			streamingUpdate(start, end);
		}

		else
		{
			collisionUpdate(start, end);
		}

		updateStates[homeAddress] = 0;
	}
}

Eigen::Vector2i directionToSteps(const int n) {
	if (n < 1 || n > 9) {
		throw "Direction out of bounds";
	}
	return Eigen::Vector2i((n-1) % 3 - 1, (9 - n) / 3 - 1);
}