#ifndef _SIMCLASS_
#define _SIMCLASS_
/*
 *  base class (common stuff for MC and MD) 
 *   [*]  calcenergyi()
 *   [*]  calctotenergy()
 *   [ ]  pbc() [ to apply pbc ]
 *   [ ]  prepare_initial_conf() 
 *   [ ]  init_rng()
 *   [ ]  run() method (main loop) 
 *   [ ]  save_mgl_snaphot() save configurations in mgl format
 *
 ** MC class
 *
 *  NTV ensemble
 *
 *   [ ]  trial_move()
 *   [ ]  acc_move()
 *   [ ]  move_NTV()
 *
 *   [ ]  calc_acceptance and adjust()
 *   [ ]  init_measures() truncate measure files
 *   [ ]  save_measures() append measures
 *
 *  optional exercises:
 *
 *   [ ]  save_snapshot() save configurations <--- mostrare Lab #7 o #8
 *   [ ]  save_restart()
 */
#include <vector>
#include <string>
#include <fstream>
#include "./params.hpp"
#include "./particle.hpp"
#include "./randnumgen.hpp"
#include <iomanip> // for setprecision()

using ntype=double;

template<typename particle_type>
class sim
{
  using simp = simpars;
protected:
  simp pars;  // parametri per composizione
  std::vector<particle_type> parts;

  // if opt=1 calculate energies only for i < j
  ntype calcenergyi(int i, int opt=0) 
    {
      int j;
      ntype enei=0.0;
      for (j=0; j < pars.Np; j++) // pars.Np è il numero totale di particelle
        {
          if (opt==1 && i >= j)
            continue;
          if (i==j)
            continue;
          // la classe particelle deve implementare un metodo vij per il calcolo dell'energia d'interazione
          enei += parts[i].vij(parts[j], pars.L); 
          // pars.L è un vettore con i lati del box 
        }
      return enei;
    }
  
  ntype totenergy(void)
    {
      ntype ene=0.0;
      for (auto i=0; i < parts.size(); i++)
        {
          ene+=calcenergyi(i, 1);
        }
      return ene;
    }

  void pbc(int i)
    {
      for (int j = 0; j < 3; ++j) {
            if (parts[i].r(j) > 0.5 * parts.L(i)) parts[i].r(j) -= parts.L(j);
            if (parts[i].r(j) < -0.5 * parts.L(i)) parts[i].r(j) += parts.L(j);
        }
    }

public:
  void set_sim_type(int type)
    {
      pars.simtype = type;
    }
 
  void init_rng(int n)
  {
    rng.seed(42); 
  }
  void prepare_initial_conf(void)
    {
      for (int i=0; i < pars.Np; i++)
        {
          for (int j=0; j < 3; j++)
            {
              parts[i].r(j) = pars.L(j)*rng.ranf();
            }
        }
    }

  void run(void) 
    {
      // intentionally void
    }; 
};

template<typename particle_type>
class mcsim: public sim<particle_type>
{
  // for calc_acceptance_and_adjust: total trial moves and accepted ones
  // for calculating acceptance rates.
  using bc=sim<particle_type>;
  using bc::calcenergyi, bc::parts, bc::pars, 
        bc::totenergy;
  
  // counters used to calculate acceptance rates
  long int tot_tra, tot_vol, tra_rej, vol_rej;

  void alpha(int i)
    pvector rand_orient; 
    rand_orient.random_orient();
    {
      for (int j = 0; j <3; j++){
        parts[i].r(j) += rand_orient.v(j) * deltra;
      }
    }
  ntype calcenergyi(int i, int opt=0) 
    {
      // calculate interaction energy for particle i
      // WRITE THIS
      return 0.0;
    }
 
  void acc(int i, ntype eno)
    {
      // acceptance of trial move according to Metropolis MC algorithm
      // assume kB=1 (Boltzmann constants) so that beta=1/T 
      // 1) calcolo la nuova energia 
      // 2) accetto la mossa con probabilità min{1taU)}
      // 3) se deltaU < 0 accetto la mossa
      //    altrimenti genero xi un numero a caso tra 0 e 1 
      //    e se xi < exp(-beta*deltaU) accetto altrimenti rifiuto
      // 4) accetta la mossa con criterio Metropolis MC

      // WRITE THIS
    }
  
  void move_NTV(int i)
    {
      // 1) calcolare l'energia della particella i-esima
      // 2) memorizzo la posizione della particella i
      // 3) trial move 
      // 4) acceptance per la particella i

    }
  void calc_acceptance_and_adjust(void)
    {
      // CALC ACCEPTANCE RATES AND 
      // ADJUST pars.deltra 
    }

  void init_measures(void)
    {
      // open files in writing mode to truncate them to 0
      // WRITE THIS
    }

 void save_measures(long int t)
    {
      // OPEN FILES AND ADD NEW MEASUREMENT

    }
 public:
  void run(void) 
    {
      // ciclo sui passi MC
      int i, t, ip;
      tot_tra = tot_vol = tra_rej = vol_rej = 0;
      init_measures();
      for (t = 0; t < pars.totsteps; t++)
        {
          // ATTEMPT TO MOVE ALL PARTICLES:
          // WRITE THIS

          if (t > 0 && pars.savemeasure > 0 && t % pars.savemeasure == 0)
            {
              save_measures(t);
            }
          if (t > 0 && t % pars.outstps == 0) 
            {
              // WRITE OUTPUT TO SCREEN 
            }
#if 0
          if (t > 0 && pars.save_mgl_snapshot > 0 && 
              t % pars.save_mgl_snapshot == 0)
            {
              save_mgl_snapshot(t);
            }
#endif
          if (t > 0 && pars.adjstps > 0 && t % pars.adjstps == 0 && 
              t < pars.maxadjstps)
            {
              calc_acceptance_and_adjust();
            }
        }
    } 
};
#endif
