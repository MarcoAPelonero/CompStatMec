#ifndef _PARAMS_
#define _PARAMS_
#include "./pvector.hpp"

class simpars
{
  using ntype=double;
public:
  int nx, ny, nz; /* nx*ny*nz particelle */
  double T, P; // temperature and pressure
  int Np; // numero di particelle
  long int maxadjstps, eqstps, adjstps, save_mgl_snapshot;
  long int savemeasure, outstps, totsteps; // Nsteps = simulations steps, outstps steps print something on current simulation status
  double rho, rc; // density
  int simtype; // simulation type (see below)
  int seed; // -1 means random
  pvector<double, 3> L; // box
  double sigma, epsilon, mass; // Lennard-Jones parameters
  double dt, deltra, vmax; // parameter of MC moves
  simpars()
    {
      simtype = 1; // 0 NTV, 1 NPT
      nx = 8;  // number of particles along each direction
      ny = 8;
      nz = 8;
      sigma=1.0;
      epsilon=1.0;
      rho = 0.5;
      rc = 4.0;
      seed=0;
      mass=1.0;
      adjstps = -200;
      maxadjstps = 2000;
      eqstps=500;
      totsteps = 20000;
      save_mgl_snapshot = 1000; 
      savemeasure=20;
      outstps=200;
      T = 2.0;
      P = 3.838; //se P*=\beta*P*v0, 1 < P* < 10 dove v0 è il volume di una particella 
      deltra = 0.2;
      vmax=10.0;
      dt = 0.002;
    }
};
#endif
