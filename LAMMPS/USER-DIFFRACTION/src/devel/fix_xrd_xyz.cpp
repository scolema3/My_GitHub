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

/* ----------------------------------------------------------------------
   Contributing author: Pieter in 't Veld (SNL)
   Incorporating XRD: Shawn Coleman (Arkansas)
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "fix_xrd_xyz.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "compute_xrd.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "math.h"
#include "domain.h"

#include <iostream>
using namespace std;


using namespace LAMMPS_NS;
using namespace FixConst;

enum{COMPUTE};
enum{ONE,RUNNING,WINDOW};
enum{FIRST,MULTI};


#define INVOKED_ARRAY 4
/* ---------------------------------------------------------------------- */

FixXRDxyz::FixXRDxyz(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix xrd/xyz command");

  MPI_Comm_rank(world,&me);

  nevery = force->inumeric(FLERR,arg[3]);
  nrepeat = force->inumeric(FLERR,arg[4]);
  nfreq = force->inumeric(FLERR,arg[5]);

  global_freq = nfreq;

  nvalues = 0;
  int iarg = 6;
  while (iarg < narg) {
    if (strncmp(arg[iarg],"c_",2) == 0){
      nvalues++;
      iarg++;
    } else break;
  }
  if (nvalues != 1) error->all(FLERR,"Illegal fix xrd/xyz command");

  options(narg,arg);

  ids = NULL;
  nvalues = 0;
  iarg = 6;

  while (iarg < narg) {
    if (strncmp(arg[iarg],"c_",2) == 0 ) {

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) error->all(FLERR,"Illegal fix xrd/xyz command");

      n = strlen(suffix) + 1;
      ids = new char[n];
      strcpy(ids,suffix);
      delete [] suffix;

cout << "icompute is next " <<  endl;

      int icompute = modify->find_compute(ids);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix xrd/xyz does not exist");


      Compute *compute = modify->compute[icompute];

      // Check that specified compute is for XRD
      compute_xrd = (ComputeXRD*) modify->compute[icompute];
      if (strcmp(compute_xrd->style,"xrd") != 0)
        error->all(FLERR,"Fix xrd/xyz has invalid compute assigned");

      // Gather varialbes from specified compute_xrd
      double *xrd_var = compute_xrd->xrd_var;

      lambda     = xrd_var[0];
      Max2Theta  = xrd_var[1];
      Min2Theta  = xrd_var[2];
      c[0]       = xrd_var[3];
      c[1]       = xrd_var[4];
      c[2]       = xrd_var[5];
      double manual_double = xrd_var[6];
      manual = false;
      if (manual_double == 1) manual = true;

      // Standard error check for fix/xrd/xyz
      if (compute->array_flag == 0)
        error->all(FLERR,"Fix xrd/xyz compute does not calculate a array");
      if (2 != compute->size_array_cols)
          error->all(FLERR,"Fix xrd/xyz compute array is accessed out-of-range");
      if (compute->extarray != 0)
        error->all(FLERR,"Illegal fix xrd/xyz command");
      nrows = compute->size_array_rows;
      nvalues++;
      iarg++;
    } else break;
  }

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable

  if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix xrd/xyz command");
  if (nfreq % nevery || (nrepeat-1)*nevery >= nfreq)
    error->all(FLERR,"Illegal fix xrd/xyz command");
  if (ave != RUNNING && overwrite)
    error->all(FLERR,"Illegal fix xrd/xyz command");

  // allocate memory for averaging

  vector = vector_total = NULL;
  vector_list = NULL;

  if (ave == WINDOW)
    memory->create(vector_list,nwindow,nvalues,"xrd/xyz:vector_list");

  memory->create(vector,nrows,"xrd/xyz:vector");
  memory->create(vector_total,nrows,"xrd/xyz:vector_total");

  extlist = NULL;

  vector_flag = 1;
  size_vector = nrows;

  if (nOutput == 0) {

    Kmax = 2 * sin(Max2Theta) / lambda;
    // XRD specific paramaters needed
    int *periodicity = domain->periodicity;

    double *prd;
    double ave_inv = 0.0;
    prd = domain->prd;

    if (periodicity[0]){
      prd_inv[0] = 1 / prd[0];
      ave_inv += prd_inv[0];
    }
    if (periodicity[1]){
      prd_inv[1] = 1 / prd[1];
      ave_inv += prd_inv[1];
    }
    if (periodicity[2]){
      prd_inv[2] = 1 / prd[2];
      ave_inv += prd_inv[2];
    }

    // Using the average inverse dimensions for non-periodic direction
    ave_inv = ave_inv / (periodicity[0] + periodicity[1] + periodicity[2]);
    if (!periodicity[0]){
      prd_inv[0] = ave_inv;
    }
    if (!periodicity[1]){
      prd_inv[1] = ave_inv;
    }
    if (!periodicity[2]){
      prd_inv[2] = ave_inv;
    }

    // Use manual mapping of reciprocal lattice
    if (manual) {
      for (int i=0; i<3; i++) {
        prd_inv[i] = 1.0;
      }
    }

    // Find reciprocal spacing and integer dimensions
    for (int i=0; i<3; i++) {
      dK[i] = prd_inv[i]*c[i];
      Knmax[i] = ceil(Kmax / dK[i]);
      Knmin[i] = -Knmax[i];
    }


   // Finding dimensions for vtk files
    for (int i=0; i<3; i++) {
      if ( ( (Knmin[i] > 0) && (Knmax[i] > 0) ) || ( (Knmin[i] < 0) && (Knmax[i] < 0) ) ){
        Dim[i] = abs( (int) Knmin[i] ) + abs( (int) Knmax[i] );
      } else Dim[i] = abs( (int) Knmin[i] ) + abs( (int) Knmax[i] ) + 1;
    }
  }

  // initialization

  irepeat = 0;
  iwindow = window_limit = 0;
  norm = 0;

  for (int i = 0; i < nrows; i++)
     vector_total[i] = 0.0;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);

}

/* ---------------------------------------------------------------------- */

FixXRDxyz::~FixXRDxyz()
{
  delete [] extlist;
  memory->destroy(vector);
  memory->destroy(vector_list);
  memory->destroy(vector_total);
  if (fp && me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixXRDxyz::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixXRDxyz::init()
{
  // set current indices for all computes,fixes,variables


  int icompute = modify->find_compute(ids);
  if (icompute < 0)
    error->all(FLERR,"Compute ID for fix xrd/xyz does not exist");

  // need to reset nvalid if nvalid < ntimestep b/c minimize was performed

  if (nvalid < update->ntimestep) {
    irepeat = 0;
    nvalid = nextvalid();
    modify->addstep_compute_all(nvalid);
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixXRDxyz::setup(int vflag)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixXRDxyz::end_of_step()
{
  // skip if not step which requires doing something
  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;
  invoke_vector(ntimestep);
}

/* ---------------------------------------------------------------------- */

void FixXRDxyz::invoke_vector(bigint ntimestep)
{

cout << "invoke vector is next " <<  endl;

  // zero if first step
  int icompute = modify->find_compute(ids);
  if (icompute < 0)
    error->all(FLERR,"Compute ID for fix xrd/xyz does not exist");

  if (irepeat == 0)
    for (int i = 0; i < nrows; i++)
       vector[i] = 0.0;

  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  // invoke compute if not previously invoked

  Compute *compute = modify->compute[icompute];

  if (!(compute->invoked_flag & INVOKED_ARRAY)) {
    compute->compute_array();
    compute->invoked_flag |= INVOKED_ARRAY;
  }
  double **carray = compute->array;
  int icol = 1;
  for (int i = 0; i < nrows; i++)
    vector[i] = carray[i][icol];

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    modify->addstep_compute(nvalid);
    return;
  }

  irepeat = 0;
  nvalid = ntimestep+nfreq - (nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // average the final result for the Nfreq timestep

  double repeat = nrepeat;
  for ( int i = 0; i < nrows; i++)
    vector[i] /= repeat;

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, combine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    for (int i = 0; i < nrows; i++) vector_total[i] = vector[i];
    norm = 1;

  } else if (ave == RUNNING) {
    for (int i = 0; i < nrows; i++) vector_total[i] += vector[i];
    norm++;

  } else if (ave == WINDOW) {
    for (int i = 0; i < nrows; i++) {
      vector_total[i] += vector[i];
      if (window_limit) vector_total[i] -= vector_list[iwindow][i];
        vector_list[iwindow][i] = vector[i];
    }

    iwindow++;
    if (iwindow == nwindow) {
      iwindow = 0;
      window_limit = 1;
    }
    if (window_limit) norm = nwindow;
    else norm = iwindow;
  }

  // output result to file

  if (fp && me == 0) {
    if (nOutput > 0) {
      fclose(fp);

      char nName [128];
      sprintf(nName,"%s.%d.xyz",filename,nOutput);
      fp = fopen(nName,"w");

      if (fp == NULL) {
        char str[128];
        sprintf(str,"Cannot open fix xrd/xyz file %s",nName);
        error->one(FLERR,str);
      }
    }

    // Header information with relp spacing
    fprintf(fp,"#LAMBDA %g\n", lambda);
    fprintf(fp,"#ASPECT_RATIO %g %g %g\n", dK[0], dK[1], dK[2]);

    filepos = ftell(fp);

    if (overwrite) fseek(fp,filepos,SEEK_SET);

   // Finding the intersection of the reciprocal space and Ewald sphere
    int NROW1 = 0;
    int NROW2 = 0;
    double dinv2 = 0.0;
    double r = 0.0;
    double K[3];

if ( (Threshold > 0)) {
  for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
    for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
      for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
        K[0] = i * dK[0];
        K[1] = j * dK[1];
        K[2] = k * dK[2];
        dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
        if  (4 >= dinv2 * lambda * lambda ) {
       	  ang = asin(lambda * sqrt(dinv2) / 2);
          if ( (ang <= Max2Theta) & (ang >= Min2Theta) ) {
            NROW1++;
            if ( (vector_total[NROW1-1]/norm >= Threshold )){
            fprintf(fp,"%g %g %g %g\n",K[0],K[1],K[2],vector_total[NROW1-1]/norm);
            }
	        }
        }
      }
    }
  }
}
else {//begin normal
  for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
    for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
      for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
        K[0] = i * dK[0];
        K[1] = j * dK[1];
        K[2] = k * dK[2];
        dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
        if  (4 >= dinv2 * lambda * lambda ) {
       	  ang = asin(lambda * sqrt(dinv2) / 2);
          if ( (ang <= Max2Theta) & (ang >= Min2Theta) ) {
            NROW1++;
            fprintf(fp,"%g %g %g %g\n",K[0],K[1],K[2],vector_total[NROW1-1]/norm);
	        }
        }
      }
    }
  }
}//end normal

  }
  nOutput++;
}

/* ----------------------------------------------------------------------
   return Ith vector value
------------------------------------------------------------------------- */

double FixXRDxyz::compute_vector(int i)
{
  if (norm) {
    return vector_total[i]/norm;
  }
  return 0.0;
}

/* ----------------------------------------------------------------------
   parse optional args
------------------------------------------------------------------------- */

void FixXRDxyz::options(int narg, char **arg)
{
  // option defaults

  fp = NULL;
  ave = ONE;
  startstep = 0;
  overwrite = 0;
  Threshold = -1;

  // optional args
  int iarg = 6 + nvalues;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix xrd/xyz command");
      if (me == 0) {

         nOutput = 0;
         int n = strlen(arg[iarg+1]) + 1;
         filename = new char[n];
         strcpy(filename,arg[iarg+1]);

        char nName [128];
         sprintf(nName,"%s.%d.xyz",filename,nOutput);
         fp = fopen(nName,"w");

        if (fp == NULL) {
          char str[128];
          sprintf(str,"Cannot open fix xrd/xyz file %s",nName);
          error->one(FLERR,str);
        }
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix xrd/xyz command");
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else if (strcmp(arg[iarg+1],"window") == 0) ave = WINDOW;
      else error->all(FLERR,"Illegal fix xrd/xyz command");
      if (ave == WINDOW) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix xrd/xyz command");
        nwindow = force->inumeric(FLERR,arg[iarg+2]);
        if (nwindow <= 0) error->all(FLERR,"Illegal fix xrd/xyz command");
      }
      iarg += 2;
      if (ave == WINDOW) iarg++;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix xrd/xyz command");
      startstep = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"overwrite") == 0) {
      overwrite = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"threshold") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix xrd/xyz command");
      Threshold = atof(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix xrd/xyz command");
  }
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
   startstep is lower bound on nfreq multiple
------------------------------------------------------------------------- */

bigint FixXRDxyz::nextvalid()
{
  bigint nvalid = (update->ntimestep/nfreq)*nfreq + nfreq;
  while (nvalid < startstep) nvalid += nfreq;
  if (nvalid-nfreq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += nfreq;
  return nvalid;
}

/* ---------------------------------------------------------------------- */

void FixXRDxyz::reset_timestep(bigint ntimestep)
{
  if (ntimestep > nvalid) error->all(FLERR,"Fix xrd/xyz missed timestep");
}
