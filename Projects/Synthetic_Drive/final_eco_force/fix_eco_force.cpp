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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "mpi.h"
#include "fix_eco_force.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "output.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixECOForce::FixECOForce(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {
  MPI_Comm_rank(world, &me);

  if (narg != 7)
    error->all(FLERR, "Illegal fix eco/force command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;

  peratom_flag = 1;
  size_peratom_cols = 2;
  peratom_freq = 1;

  u0 = atof(arg[3]);
  if(u0>0) { sign = 1; } else { sign = -1; }
  eta = atof(arg[4]);
  rcut = atof(arg[5]);
  int n = strlen(arg[6]) + 1;
  orifilename = new char[n];
  strcpy(orifilename, arg[6]);
  // read reference orientations from file
  if (me == 0) {
    char line[IMGMAX];
    char *result;
    int count;

    FILE *infile = fopen(orifilename, "r");
    if (infile == NULL)
      error->one(FLERR, "Fix eco/force file open failed");
    for (int i = 0; i < 6; i++) {
      result = fgets(line, IMGMAX, infile);
      if (!result)
        error->one(FLERR, "Fix eco/force file read failed");
      count = sscanf(line, "%lg %lg %lg", &ori[i][0], &ori[i][1], &ori[i][2]);
      if (count != 3)
        error->one(FLERR, "Fix eco/force file read failed");
    }
    fclose(infile);

    // calculate inverse lattice vectors
    double vol = 0.5 / MY_PI * (ori[0][0] * (ori[1][1] * ori[2][2] - ori[1][2] * ori[2][1]) + ori[0][1] * (ori[1][2] * ori[2][0] - ori[1][0] * ori[2][2]) + ori[0][2] * (ori[1][0] * ori[2][1] - ori[1][1] * ori[2][0]));
    invLatVec[0][0] = (ori[1][1] * ori[2][2] - ori[1][2] * ori[2][1]) / vol;
    invLatVec[0][1] = (ori[1][2] * ori[2][0] - ori[1][0] * ori[2][2]) / vol;
    invLatVec[0][2] = (ori[1][0] * ori[2][1] - ori[1][1] * ori[2][0]) / vol;
    invLatVec[1][0] = (ori[2][1] * ori[0][2] - ori[2][2] * ori[0][1]) / vol;
    invLatVec[1][1] = (ori[2][2] * ori[0][0] - ori[2][0] * ori[0][2]) / vol;
    invLatVec[1][2] = (ori[2][0] * ori[0][1] - ori[2][1] * ori[0][0]) / vol;
    invLatVec[2][0] = (ori[0][1] * ori[1][2] - ori[0][2] * ori[1][1]) / vol;
    invLatVec[2][1] = (ori[0][2] * ori[1][0] - ori[0][0] * ori[1][2]) / vol;
    invLatVec[2][2] = (ori[0][0] * ori[1][1] - ori[0][1] * ori[1][0]) / vol;
    vol = 0.5 / MY_PI * (ori[3][0] * (ori[4][1] * ori[5][2] - ori[4][2] * ori[5][1]) + ori[3][1] * (ori[4][2] * ori[5][0] - ori[4][0] * ori[5][2]) + ori[3][2] * (ori[4][0] * ori[5][1] - ori[4][1] * ori[5][0]));
    invLatVec[3][0] = (ori[4][1] * ori[5][2] - ori[4][2] * ori[5][1]) / vol;
    invLatVec[3][1] = (ori[4][2] * ori[5][0] - ori[4][0] * ori[5][2]) / vol;
    invLatVec[3][2] = (ori[4][0] * ori[5][1] - ori[4][1] * ori[5][0]) / vol;
    invLatVec[4][0] = (ori[5][1] * ori[3][2] - ori[5][2] * ori[3][1]) / vol;
    invLatVec[4][1] = (ori[5][2] * ori[3][0] - ori[5][0] * ori[3][2]) / vol;
    invLatVec[4][2] = (ori[5][0] * ori[3][1] - ori[5][1] * ori[3][0]) / vol;
    invLatVec[5][0] = (ori[3][1] * ori[4][2] - ori[3][2] * ori[4][1]) / vol;
    invLatVec[5][1] = (ori[3][2] * ori[4][0] - ori[3][0] * ori[4][2]) / vol;
    invLatVec[5][2] = (ori[3][0] * ori[4][1] - ori[3][1] * ori[4][0]) / vol;
  }

  // initializations
  MPI_Bcast(&ori[0][0], 18, MPI_DOUBLE, 0, world);
  MPI_Bcast(&invLatVec[0][0], 18, MPI_DOUBLE, 0, world);
  kappa[0] = kappa[1] = kappa[2] = 1;
  kappa[3] = kappa[4] = kappa[5] = -1;

  comm_forward = 1;
  added_energy = 0.0;
  nmax = atom->nmax;
  forcePrefact = (double *) memory->smalloc(nmax * sizeof(double), "eco/force:forcePrefact");
  sht_first = (int **) memory->smalloc(nmax * sizeof(int *), "eco/force:sht_first");
  memory->create(order, nmax, 2, "eco/force:order");
  memory->create(sht_num, nmax, "eco/force:sht_num");
  array_atom = order;

  ipage = NULL;
  pgsize = oneatom = 0;

  // zero the array since a variable may access it before first run
  for (int i = 0; i < atom->nlocal; i++)
    order[i][0] = order[i][1] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixECOForce::~FixECOForce() {
  delete[] orifilename;
  memory->sfree(forcePrefact);
  memory->destroy(order);
  memory->destroy(sht_num);
  memory->sfree(sht_first);
  delete[] ipage;
}

/* ---------------------------------------------------------------------- */

int FixECOForce::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixECOForce::init() {
  cutsq = rcut * rcut;
  MPI_Comm_rank(world, &me);
  if (me == 0) { //calculates the normalization factor for force, depending on cutoff radius
    double dijx, dijy, dijz, dikx, diky, dikz, rsqij, rsqik, wij, wik, skalarprod;
    int layer = 4; //will produce wrong results for rcut>3*lattice constant
    int neigh = 0;
    scalenorm = 0;
    for (int ix = -layer; ix <= layer; ix++) {
      for (int iy = -layer; iy <= layer; iy++) {
        for (int iz = -layer; iz <= layer; iz++) {
          dijx = ix * ori[3][0] + iy * ori[4][0] + iz * ori[5][0];
          dijy = ix * ori[3][1] + iy * ori[4][1] + iz * ori[5][1];
          dijz = ix * ori[3][2] + iy * ori[4][2] + iz * ori[5][2];
          rsqij = dijx * dijx + dijy * dijy + dijz * dijz;
          if (rsqij != 0 && rsqij < cutsq) {
            neigh++;
            for (int ia = -layer; ia <= layer; ia++) {
              for (int ib = -layer; ib <= layer; ib++) {
                for (int ic = -layer; ic <= layer; ic++) {
                  dikx = ia * ori[3][0] + ib * ori[4][0] + ic * ori[5][0];
                  diky = ia * ori[3][1] + ib * ori[4][1] + ic * ori[5][1];
                  dikz = ia * ori[3][2] + ib * ori[4][2] + ic * ori[5][2];
                  rsqik = dikx * dikx + diky * diky + dikz * dikz;
                  if (rsqik != 0 && rsqik < cutsq) {
                    for (int alpha = 0; alpha < 3; alpha++) {
                      skalarprod = invLatVec[alpha][0] * (dikx - dijx) + invLatVec[alpha][1] * (diky - dijy) + invLatVec[alpha][2] * (dikz - dijz);
                      wij = rsqij / cutsq * (rsqij / cutsq - 2) + 1;
                      wik = rsqik / cutsq * (rsqik / cutsq - 2) + 1;
                      scalenorm += wij * wik * (1 - cos(skalarprod));
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    printf("fix/eco/force: cutoff=%f neighbors=%i\n", rcut, neigh);
  }
  if (rcut > force->pair->cutforce) {
    error->all(FLERR, "Pair cutoff too small: use 'hybrid/overlay' to set the pair cutoff larger than the fix eco/force cutoff");
  }
  if (2 * rcut > force->pair->cutforce and 2 * rcut > comm->cutghostuser and u0!=0) {
    error->all(FLERR, "Pair cutoff too small: set 'comm_modify cutoff' to larger than twice the fix eco/force cutoff");
  }
  //send normfactors to other processors
  MPI_Bcast(&scalenorm, 1, MPI_DOUBLE, 0, world);
  if (strstr(update->integrate_style, "respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // need a full neighbor list, built whenever re-neighboring occurs
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->newton = 2;
  neighbor->requests[irequest]->ghost = 1;

  int create = 0;
  if (ipage == NULL)
    create = 1;
  if (pgsize != neighbor->pgsize)
    create = 1;
  if (oneatom != neighbor->oneatom)
    create = 1;

  if (create) {
    delete[] ipage;
    pgsize = neighbor->pgsize;
    oneatom = neighbor->oneatom;

    int nmypage = comm->nthreads;
    ipage = new MyPage<int> [nmypage];
    for (int i = 0; i < nmypage; i++)
      ipage[i].init(oneatom, pgsize);
  }
}

/* ---------------------------------------------------------------------- */

void FixECOForce::init_list(int id, NeighList *ptr) {
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixECOForce::setup(int vflag) {
  if (strstr(update->integrate_style, "verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa - 1);
    post_force_respa(vflag, nlevels_respa - 1, 0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa - 1);
  }
}

/* ---------------------------------------------------------------------- */

void FixECOForce::post_force(int vflag) {
  int i, j, k, ii, jj, kk, inum, jnum, knum, alpha;
  int *ilist, *jlist, *klist, *numneigh, **firstneigh;
  double dkjx, dkjy, dkjz, dijx, dijy, dijz, djkx, djky, djkz, dikx, diky, dikz, rsqik, rsqij, rsqkj, wij, wik, wkj, wijdot, skalarprodjk, skalarprodik, omega, sumalpha, forcemag1, forcemag2, forcemag3;
  // set local ptrs
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  // insure nbr and order data structures are adequate size
  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->sfree(forcePrefact);
    memory->destroy(order);
    memory->destroy(sht_num);
    memory->sfree(sht_first);
    forcePrefact = (double *) memory->smalloc(nmax * sizeof(double), "eco/force:forcePrefact");
    sht_first = (int **) memory->smalloc(nmax * sizeof(int *), "eco/force:sht_first");
    memory->create(order, nmax, 2, "eco/force:order");
    memory->create(sht_num, nmax, "eco/force:sht_num");
    array_atom = order;
  }
  Short_neigh(); //build neighbor lists containing only atoms within cutoff distance
  added_energy = 0.0;
  for (ii = 0; ii < inum; ii++) {
    order[ii][0] = order[ii][1] = 0;
  }
  for (ii = 0; ii < inum; ii++) {  //potentials and precalculations
    i = ilist[ii];
    jlist = sht_first[i];
    jnum = sht_num[i];
    sumalpha = 0;
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      dijx = x[i][0] - x[j][0];
      dijy = x[i][1] - x[j][1];
      dijz = x[i][2] - x[j][2];
      rsqij = (dijx * dijx + dijy * dijy + dijz * dijz) / cutsq;
      wij = rsqij * (rsqij - 2) + 1;
      for (kk = 0; kk < jnum; kk++) {
        k = jlist[kk];
        k &= NEIGHMASK;
        dikx = x[i][0] - x[k][0];
        diky = x[i][1] - x[k][1];
        dikz = x[i][2] - x[k][2];
        rsqik = (dikx * dikx + diky * diky + dikz * dikz) / cutsq;
        wik = rsqik * (rsqik - 2) + 1;
        for (alpha = 0; alpha < 6; alpha++) {
          skalarprodjk = invLatVec[alpha][0] * (dikx - dijx) + invLatVec[alpha][1] * (diky - dijy) + invLatVec[alpha][2] * (dikz - dijz);
          sumalpha += kappa[alpha] * wij * wik * cos(skalarprodjk);
        }  //alpha
      }  //atom k
    }  //atom j
    sumalpha /= scalenorm;
    if (sumalpha > eta) {
      added_energy += u0 * 0.5;
      forcePrefact[i] = 0;
      order[i][1] += sign;
    } else if (sumalpha < -eta) {
      added_energy -= u0 * 0.5;
      forcePrefact[i] = 0;
      order[i][1] -= sign;
    } else {
      omega = MY_PI2 / eta * sumalpha;
      order[i][1] += sign*sin(omega);
      forcePrefact[i] = -u0 * MY_PI2 / eta * cos(omega) / scalenorm;
      //forcePrefact[i]=0;
      added_energy += u0 * 0.5 * sin(omega);
    }
    order[i][0] += sumalpha;
  }  //atom i
//communicate to acquire forcePrefact data for ghost atoms
  if(u0!=0){
  comm->forward_comm_fix(this);
  for (ii = 0; ii < inum; ii++) {  //forces
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    jlist = sht_first[i];
    jnum = sht_num[i];
    if (forcePrefact[i]) {
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        dijx = x[i][0] - x[j][0];
        dijy = x[i][1] - x[j][1];
        dijz = x[i][2] - x[j][2];
        rsqij = (dijx * dijx + dijy * dijy + dijz * dijz) / cutsq;
        wijdot = 4 / cutsq * (rsqij - 1);
        for (kk = 0; kk < jnum; kk++) {
          k = jlist[kk];
          k &= NEIGHMASK;
          dikx = x[i][0] - x[k][0];
          diky = x[i][1] - x[k][1];
          dikz = x[i][2] - x[k][2];
          rsqik = (dikx * dikx + diky * diky + dikz * dikz) / cutsq;
          wik = rsqik * (rsqik - 2) + 1;
          djkx = x[j][0] - x[k][0];
          djky = x[j][1] - x[k][1];
          djkz = x[j][2] - x[k][2];
          for (alpha = 0; alpha < 6; alpha++) {
            skalarprodjk = invLatVec[alpha][0] * djkx + invLatVec[alpha][1] * djky + invLatVec[alpha][2] * djkz;
            forcemag1 = kappa[alpha] * forcePrefact[i] * wijdot * wik * cos(skalarprodjk);
            f[i][0] += dijx * forcemag1;
            f[i][1] += dijy * forcemag1;
            f[i][2] += dijz * forcemag1;
          }  //alpha
        } // atom k
      }  //atom j
    } // speed optimization
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      if (forcePrefact[j]) {
        dijx = x[i][0] - x[j][0];
        dijy = x[i][1] - x[j][1];
        dijz = x[i][2] - x[j][2];
        rsqij = (dijx * dijx + dijy * dijy + dijz * dijz) / cutsq;
        //if (rsqij<1){//unwanted
        wij = rsqij * (rsqij - 2) + 1;
        wijdot = 4 / cutsq * (rsqij - 1);
        klist = sht_first[j];
        knum = sht_num[j];//wanted
        for (kk = 0; kk < knum; kk++) {
          k = klist[kk];
          k &= NEIGHMASK;
          dkjx = x[k][0] - x[j][0];
          dkjy = x[k][1] - x[j][1];
          dkjz = x[k][2] - x[j][2];
          rsqkj = (dkjx * dkjx + dkjy * dkjy + dkjz * dkjz) / cutsq;
          wkj = rsqkj * (rsqkj - 2) + 1;
          dikx = x[i][0] - x[k][0];
          diky = x[i][1] - x[k][1];
          dikz = x[i][2] - x[k][2];
          for (alpha = 0; alpha < 6; alpha++) {
            skalarprodik = invLatVec[alpha][0] * dikx + invLatVec[alpha][1] * diky + invLatVec[alpha][2] * dikz;
            forcemag2 = kappa[alpha] * forcePrefact[j] * wijdot * wkj * cos(skalarprodik);
            forcemag3 = -kappa[alpha] * forcePrefact[j] * wij * wkj * sin(skalarprodik);
            f[i][0] += dijx * forcemag2 + invLatVec[alpha][0] * forcemag3;
            f[i][1] += dijy * forcemag2 + invLatVec[alpha][1] * forcemag3;
            f[i][2] += dijz * forcemag2 + invLatVec[alpha][2] * forcemag3;
          }  //alpha
        } // atom k
      } // speed optimization
    }  //atom j
  }  //atom i
  }  //do not calculate for u0=0
}

/* ---------------------------------------------------------------------- */

void FixECOForce::post_force_respa(int vflag, int ilevel, int iloop) {
  if (ilevel == nlevels_respa - 1)
    post_force(vflag);
}

/* ---------------------------------------------------------------------- */

double FixECOForce::compute_scalar() {
  double added_energy_total;
  MPI_Allreduce(&added_energy, &added_energy_total, 1, MPI_DOUBLE, MPI_SUM, world);
  return added_energy_total;
}

/* ---------------------------------------------------------------------- */

int FixECOForce::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc) {
  int i, k;
  int m = 0;
  for (i = 0; i < n; i++) {
    k = list[i];
    buf[m++] = forcePrefact[k];
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

void FixECOForce::unpack_comm(int n, int first, double *buf) {

  int i;
  int last = first + n;
  int m = 0;
  for (i = first; i < last; i++) {
    forcePrefact[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void FixECOForce::Short_neigh() {
  int nj;
  int inum, jnum, i, j, ii, jj;
  int *neighptrj, *ilist, *jlist, *numneigh;
  int **firstneigh;
  double xtmp, ytmp, ztmp, rsq, delrj[3];

  double **x = atom->x;

  inum = list->inum + list->gnum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
// create neighbor list

  ipage->reset();

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    nj = 0;
    neighptrj = ipage->vget();

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delrj[0] = xtmp - x[j][0];
      delrj[1] = ytmp - x[j][1];
      delrj[2] = ztmp - x[j][2];
      rsq = delrj[0] * delrj[0] + delrj[1] * delrj[1] + delrj[2] * delrj[2];

      if (rsq < cutsq) {
        neighptrj[nj++] = j;
      }
    }

    sht_first[i] = neighptrj;
    sht_num[i] = nj;
    ipage->vgot(nj);
    if (ipage->status())
      error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
  }
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based arrays
 ------------------------------------------------------------------------- */

double FixECOForce::memory_usage() {
  double bytes = 3 * nmax * sizeof(double);
  for (int i = 0; i < comm->nthreads; i++)
    bytes += ipage[i].size();
  return bytes;
}
