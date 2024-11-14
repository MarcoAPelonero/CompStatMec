#ifndef _PARTICLE_
#define _PARTICLE_
/*
 * MC 
 *
 * [ ] vij
 * [ ] store/restore methods
 * [ ] tra_move()
 *
 */
#include "./pvector.hpp"
#include "./params.hpp"
using ntype = double;
class particle
{
  pvector<ntype,3> rold; // to store particle's position 
protected:
  ntype vcut;
public:
  ntype sigma, epsilon, rc;
  pvector<ntype,3> r; // particle's position and velocity

  // vcut = vLJ(rc), i.e. the value of LJ potential at r=rc
  // this method should be called when initializing particles
  void set_vcut(void);
 
  void tra_move(pvector<ntype,3> delr)
    {
      r+=delr; 
    }
  void store()
    {
      rold = r;
    }
  void restore()
    {
      r = rold;
    }

  void set_sigma(ntype sig)
    {
      sigma = sig;
    }
  void set_epsilon(ntype eps)
    {
      epsilon = eps;
    }
  void set_rcut(ntype rcut)
    {
      rc = rcut;
    }
  particle()
    {
      sigma=1.0;
      epsilon=1.0;
      rc=2.5;
      vcut = 0.0;
    }
};

class particleLJ: public particle
{
public:
    ntype vij(particleLJ P, pvector<ntype,3> L)
    {
        ntype ene = 0.0;
        pvector<ntype, 3> dr = r - P.r;

        // Apply minimum image convention for periodic boundary conditions
        for (int i = 0; i < 3; ++i) {
            if (dr(i) > 0.5 * L(i)) dr(i) -= L(i);
            if (dr(i) < -0.5 * L(i)) dr(i) += L(i);
        }

        ntype r2 = dr.norm();
        if (r2*r2 < rc * rc) {
            ntype r6 = (sigma * sigma) / r2;
            r6 = r6 * r6 * r6;
            ene = 4 * epsilon * (r6 * r6 - r6) - vcut;
            return ene;
        }
        else {
          return 0.0;
        }
    }

  // vcut = vLJ(rc), i.e. the value of LJ potential at r=rc
  // this method should be called when initializing particles
  void set_vcut(void)
    {
      vcut = 4 * epsilon * (std::pow(sigma / rc, 12) - std::pow(sigma / rc, 6));
    }
};
// Now do the same for the soft sphere potential
class particleSS: public particle
{

  ntype vij(particleLJ P, pvector<ntype,3> L)
    {
      ntype ene=0.0;
      // Implement the periodic boundary conditions for the case 
      pvector<ntype, 3> dr = r - P.r;
      // WARNING: There is no softness parameter in the configuration file, it will be initialized to 1 by simply ignoring it

      for (int i = 0; i < 3; ++i) {
            if (dr(i) > 0.5 * L(i)) dr(i) -= L(i);
            if (dr(i) < -0.5 * L(i)) dr(i) += L(i);
        }

      ntype r2 = dr.norm();
      if (r2*r2 < rc * rc) {
            ntype r6 = epsilon*(sigma / r2);
            ene = epsilon * (sigma/r2) - vcut;
            return ene;
        }
      else {
          return 0.0;
      }
    }

  void set_vcut(void)
    {
      vcut = 4 * epsilon * (std::pow(sigma / rc, 12) - std::pow(sigma / rc, 6));
    }
};
#endif
