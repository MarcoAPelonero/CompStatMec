#ifndef _SIMCLASS_
#define _SIMCLASS_
/*
 * X = done 
 * * = todo
 * E = exercise
 * 
 * OE = optional exercise
 *
 *  ===== base class (common stuff for MC and MD) =====  [ 12 and 19 Nov 2024 ] 
 *
 *   [X]  calcenergyi()
 *   [X]  calctotenergy()
 *   [X]  pbc() [ to apply pbc ]
 *   [X]  save_mgl_snaphot() save configurations in mgl format
 *   [X]  prepare_initial_conf() 
 *   [X]  init_rng()
 *   [X]  run() method (main loop) 
 *
 *  ===== MC class =====
 *
 *  NTV ensemble [ 12 and 19 Nov 2024 ]
 *
 *   [X]  alpha()
 *   [X]  acc()
 *   [X]  move_NTV()
 *
 *   [X]  calc_acceptance and adjust()
 *   [X]  init_measures() truncate measure files
 *   [X]  save_measures() append measures
 *
 *   [OE]  tail corrections
 *
 *  optional exercises:
 *
 *   // in "sim" base class
 *   [OE]  save_snapshot() save configurations 
 *   [OE]  save_restart()
 *   
 *   NPT ensemble [ 26/11/2024 ]
 *
 *   // in "mcsim" derived class
 *
 *   [*] init_measures() and save_measures(): add density saving
 *
 *   [*] modify run() method to attempt box moves too
 *   
 *   [ 26/11/2024 ] EXERCISES FOR NPT ensemble:
 *
 *   [E] restore_all_pars (restore positions of all particles)
 *
 *   [E] store_all_pars (store position of all particles)
 *
 *   [E] alpha_box() i.e. trial box move
 *   
 *   [E] acc_box() i.e. acceptance of trial box move
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
      auto Dr = parts[i].r;
      Dr = pars.L.mulcw(rint(Dr.divcw(pars.L))); // L*rint(Dr/L)
      parts[i].r = parts[i].r - Dr;
    }

  void save_mgl_snapshot(long int t) 
    {
       std::fstream f;
       std::string s;

       s = "cnf-" + std::to_string(t) + ".mgl";
       f.open(s, std::ios::out|std::ios::trunc);
       for (int i=0; i < pars.Np; i++)
         {
           f << parts[i].r(0) << " " << parts[i].r(1) << " " <<
             parts[i].r(2) << " @ " << pars.sigma*0.5 << "\n";
         }
       f.close();
    }


public:
  void prepare_initial_conf(void) 
    {
      // SC 
      int ix, iy, iz;
      int cc=0;
      pars.Np = pars.nx*pars.ny*pars.nz;
      parts.resize(pars.Np);
      ntype vcell = pow(pars.sigma,3.0);
      ntype rhomax = 1.0/vcell;
      ntype sf;
      sf = cbrt(rhomax/pars.rho);
      pars.L ={ntype(pars.nx), ntype(pars.ny), ntype(pars.nz)};
      ntype clen = sf*pars.sigma;
      pars.L *= clen;
      for (ix = 0; ix < pars.nx; ix++)
        for (iy = 0; iy < pars.ny; iy++)
          for (iz = 0; iz < pars.nz; iz++)
            {
              parts[cc].r = {ix*clen, iy*clen, iz*clen};
              parts[cc].r -= pars.L*0.5;
              parts[cc].set_sigma(pars.sigma);
              parts[cc].set_epsilon(pars.epsilon);
              parts[cc].set_rcut(pars.rc);
              cc++;
            }
      // ...or BCC or FCC lattice
    }

  void set_sim_type(int type)
    {
      pars.simtype = type;
    }
 
  void init_rng(int n) 
    {
      if (n < 0)
        rng.rseed();
      else
        rng.seed(n);
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
        bc::totenergy, bc::save_mgl_snapshot, bc::pbc;
  
  // counters used to calculate acceptance rates
  long int tot_tra, tot_vol, tra_rej, vol_rej;

  void alpha(int i)
    {
      // muovere la particella usando come max 
      // displacement pars.deltra
      // applico le pbc (tramite il metodo pbc)
      pvector<ntype,3> delr;
      delr ={pars.deltra*2.0*(rng.ranf()-0.5),  pars.deltra*2.0*(rng.ranf()-0.5),
           pars.deltra*2.0*(rng.ranf()-0.5)};
      parts[i].tra_move(delr);
      pbc(i);
    }

  void acc(int i, ntype eno)
    {
      // assume kB=1 (Boltzmann constants) so that beta=1/T 
      // calcolo la nuova energia 
      // accetto la mossa con probabilità min{1taU)}
      // se deltaU < 0 accetto la mossa
      // altrimenti genero xi un numero a caso tra 0 e 1 
      // e se xi < exp(-beta*deltaU) accetto altrimenti rifiuto
      // accetta la mossa con criterio Metropolis MC
      ntype enn = calcenergyi(i);
      ntype delu= enn-eno;
      ntype xi = rng.ranf();
      if (delu > 0.0 && xi >= exp(-delu/pars.T))
        {
          // reject move
          tra_rej++;
          parts[i].restore();
        }
    }
 
  void move_NTV(int i)
    {
      // 1) calcolare l'energia della particella i-esima
      ntype eno;
      eno = calcenergyi(i);
      // 2) memorizzo la posizione della particella i
      parts[i].store();
      // 3) trial move 
      alpha(i);
      // 4) acceptance per la particella i
      acc(i, eno);
    }
  void calc_acceptance_and_adjust(void)
    {
      // CALC ACCEPTANCE RATES AND 
      // ADJUST pars.deltra 
      ntype r_tra, r_vol;
      if (tot_tra > 0)
        {
          // rate di accettazione delle mosse
          r_tra = ((double)(tot_tra - tra_rej))/tot_tra;
          std::cout << "rate tra: " << r_tra << " deltra=" << pars.deltra << "\n";
          if (r_tra > 0.5)
            {
              pars.deltra *= 1.1;
            }
          else
            {
              pars.deltra /= 1.1;
            }
          tot_tra=tra_rej = 0;
        }
     
      // adjust maximum volume "displacement" in alpha_box
      // so that acceptance rate of box move is around 0.5
      if (tot_vol > 0) // <------------------------------------------------
        {
          r_vol = ((ntype)(tot_vol - vol_rej)) / tot_vol;
          std::cout << "rate vol: " << r_vol << " vmax=" << pars.vmax << "\n";
          if (r_vol > 0.5) 
            {
              pars.vmax *= 1.1;
           }
          else
            {
              pars.vmax /= 1.1;
            }
          tot_vol=vol_rej=0.0;
        }

    }

  void init_measures(void)
    {
      // open files in writing mode to truncate them to 0
      std::fstream f;
      f.open("energy.dat", std::ios::out|std::ios::trunc);
      f.close();
      if (pars.simtype==1) // simtype == 1 means an NPT simluations
        {
          // clear file used to store densities in NPT simulation 
          f.open("density.dat", std::ios::out|std::ios::trunc);
          f.close();
        }

    }

 void save_measures(long int t)
    {
      std::fstream f;
      f.open("energy.dat", std::ios::out|std::ios::app);
      // save potential energy per particle
      f << t << " " << totenergy()/pars.Np << "\n";
      f.close();
      if (pars.simtype==1) // 1 means NPT simulation, save density in this case
        {
          f.open("density.dat", std::ios::out|std::ios::app);
          f << t << " " << pars.Np/(pars.L(0)*pars.L(1)*pars.L(2)) << "\n";
          f.close();
        }

    }


 void restore_all_pars()
   {
    for (int i = 0; i < pars.nx * pars.ny * pars.nz; i++)
    {
        parts[i].restore();
    }
   }


 void store_all_pars()
   {
    std::vector<particle_type> stored_pars(pars.nx * pars.ny * pars.nz);
    for (int i = 0; i < pars.nx * pars.ny * pars.nz; i++)
    {
        stored_pars[i] = parts[i];
    }
   }
/*
 void alpha_box(ntype& DG, ntype& fact)
    {

      // trial box move
      // WRITE THIS
      //
      // 1) choose randomly a new volume
      //
      // 2) scale all particle positions and box size by factor "fact"
      // 
      // 3) calculate new total interaction energy
      //
      // 4) calculate \DeltaG (see pdf of lectures)
    }
              
  void acc_box(ntype DG, ntype fact)
    {
      // accept or reject box move
      // WRITE THIS
      //
      // 1) calculate e^{-\Delta G} (see pdf of lectures)
      // 2) generate a random numner \xi and check wheter to accept the box trial move      
      // 3) if move is rejected:
      //    i)   restore all particle positions thourgh method restore_all_pars()
      //    ii)  restore box size
      //    iii) update counter vol_rej of rejected box moves for calculating acceptance ratio 
    }

 
 void move_box(void)// <------------------------------------------------
    {
      ntype DG, fact; // DG is \Delta G as discussed during lectures (see pdf of lectures)
      store_all_pars(); // store all particle positions before attempting a box move
      alpha_box(DG, fact);
      acc_box(DG, fact);
    }
*/
void alpha_box(ntype& DG, ntype& fact)
{
    // 1) Choose randomly a new volume
    ntype new_volume = pars.L(0) * pars.L(1) * pars.L(2) * (1.0 + (rng.ranf() - 0.5) * pars.vmax);
    fact = cbrt(new_volume / (pars.L(0) * pars.L(1) * pars.L(2)));

    // 2) Scale all particle positions and box size by factor "fact"
    for (int i = 0; i < pars.Np; i++)
    {
        parts[i].r *= fact;
    }
    pars.L *= fact;

    // 3) Calculate new total interaction energy
    ntype new_energy = 0.0;
    for (int i = 0; i < pars.Np; i++)
    {
        new_energy += calcenergyi(i);
    }

    // 4) Calculate \DeltaG
    DG = new_energy - totenergy();
}

void acc_box(ntype DG, ntype fact)
{
    // 1) Calculate e^{-\Delta G}
    ntype acceptance_prob = exp(-DG / pars.T);

    // 2) Generate a random number \xi and check whether to accept the box trial move
    ntype xi = rng.ranf();
    if (xi >= acceptance_prob)
    {
        // 3) If move is rejected:
        // i) Restore all particle positions through method restore_all_pars()
        restore_all_pars();

        // ii) Restore box size
        pars.L /= fact;

        // iii) Update counter vol_rej of rejected box moves for calculating acceptance ratio
        vol_rej++;
    }
    else
    {
        // Update total energy if the move is accepted
        totenergy() += DG;
    }
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
          // ogni passo MC sono Np tentativi di mossa
          for (i=0; i < pars.Np; i++)
            {
              // choose between box move and particle move
              if (pars.simtype==1 && rng.ranf() < 1.0/(pars.Np+1))
                {
                  move_box();// <------------------------------------------------
                  tot_vol++; 
                }
              else 
                {
                  ip = rng.ranf()*pars.Np;
                  move_NTV(ip);
                  tot_tra++;
                }
            }
        
          if (t > 0 && pars.savemeasure > 0 && t % pars.savemeasure == 0)
            {
              save_measures(t);
            }

          if (t > 0 && t % pars.outstps == 0) 
            {
              std::cout << "Step #" << t << "\n";
              // per confrontarsi con il Johnson si deve calcolare l'energia interna ovvero sommare il contributo
              // di energia cinetica
              std::cout << "total energy per particle is " << totenergy()/pars.Np << "\n";
              std::cout << "box size: " << pars.L(0) << " " << pars.L(1) << " " << pars.L(2) << "\n";
            }

          if (t > 0 && pars.save_mgl_snapshot > 0 && 
              t % pars.save_mgl_snapshot == 0)
            {
              save_mgl_snapshot(t);
            }

          if (t > 0 && pars.adjstps > 0 && t % pars.adjstps == 0 && 
              t < pars.maxadjstps)
            {
              calc_acceptance_and_adjust();
            }
        }
    } 
};
#endif
