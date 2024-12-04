#ifndef _PARTICLE_
#define _PARTICLE_
/*
 * MC 
 *
 * [X] vij
 * [X] store/restore methods
 * [X] tra_move()
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
      ntype ene;
      pvector<ntype,3> Dr;
     
      Dr = r - P.r;
      // MINIMUM IMAGE CONVENTION 
      Dr = Dr - L.mulcw(rint(Dr.divcw(L))); // Dr - L*rint(Dr/L)
      ntype rsq, rn = Dr.norm();
      rsq = rn*rn; 
      if (rsq < rc*rc) // interaction potential cut-off
        ene = 4.0*epsilon*(pow(sigma/rn,12.0)-pow(sigma/rn,6));
      else
        ene = 0.0;
      return ene;
    }

  // vcut = vLJ(rc), i.e. the value of LJ potential at r=rc
  // this method should be called when initializing particles
  void set_vcut(void)
    {
      // WRITE THIS
    }
};

class particleSS: public particle
{

  ntype vij(particleLJ P, pvector<ntype,3> L)
    {
      ntype ene;
      pvector<ntype,3> Dr;
     
      Dr = r - P.r;
      // MINIMUM IMAGE CONVENTION 
      Dr = Dr - L.mulcw(rint(Dr.divcw(L))); // Dr - L*rint(Dr/L)
      ntype rsq, rn = Dr.norm();
      rsq = rn*rn; 
      if (rsq < rc*rc) // interaction potential cut-off
        ene = 4.0*epsilon*(pow(sigma/rn,12.0));
      else
        ene = 0.0;
      return ene;
    }

  // vcut = vLJ(rc), i.e. the value of LJ potential at r=rc
  // this method should be called when initializing particles
  void set_vcut(void)
    {
      // WRITE THIS
    }
};
#endif 
