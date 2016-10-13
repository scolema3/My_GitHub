/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(eco/force/hcp,FixECOForceHCP)

#else

#ifndef LMP_FIX_ECO_FORCE_HCP_H
#define LMP_FIX_ECO_FORCE_HCP_H

#include "fix.h"
#include "my_page.h"

namespace LAMMPS_NS {

class FixECOForceHCP: public Fix {
public:
  FixECOForceHCP(class LAMMPS *, int, char **);
  ~FixECOForceHCP();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  double compute_scalar();
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  double memory_usage();

private:
  int me;
  int nlevels_respa;
  double rcut, scalenorm, cutsq;
  double u0;                      // potential value
  int sign;
  char *orifilename;  // file name for crystal orientations
  int kappa[22];
  double eta;
  double ori[8][3], invLatVec[22][3];
  double added_energy;
  double *forcePrefact; // force prefactor per atom
  int nmax;                        // expose 2 per-atom quantities
  double **order;                  // order param and normalized order param

  class NeighList *list;

  int pgsize;                      // size of neighbor page
  int oneatom;                     // max # of neighbors for one atom
  int *sht_num, **sht_first;        // short-range neighbor list
  MyPage<int> *ipage;              // neighbor list pages
  void Short_neigh();
};

}

#endif
#endif
