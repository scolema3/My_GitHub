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
   Contributing authors: Shawn Coleman & Douglas Spearot (Arkansas)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "math_const.h"
#include "compute_xrd.h"
#include "compute_xrd_consts.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "citeme.h"
#include "memory.h"
#include "error.h"
#include "stdio.h"
#include "string.h"

#include <iostream>

#ifdef _OPENMP
#include "omp.h"
#include "comm.h"
#endif

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;

static const char cite_compute_xrd_c[] =
  "compute_xrd command:\n\n"
  "@Article{Coleman13,\n"
  " author = {S. P. Coleman, D. E. Spearot, L. Capolungo},\n"
  " title = {Virtual diffraction analysis of Ni [010] symmetric tilt grain boundaries},\n"
  " journal = {Modelling and Simulation in Materials Science and Engineering},\n"
  " year =    2013,\n"
  " volume =  21,\n"
  " pages =   {055020}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

ComputeXRD::ComputeXRD(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_compute_xrd_c);

  MPI_Comm_rank(world,&me);
  int ntypes = atom->ntypes;
  int natoms = group->count(igroup);
  int dimension = domain->dimension;
  int *periodicity = domain->periodicity;
  int triclinic = domain->triclinic;

  // Checking errors
  if (dimension == 2) 
     error->all(FLERR,"Compute XRD does not work with 2d structures");
  if (narg < 4+ntypes) 
     error->all(FLERR,"Illegal Compute XRD Command");
  if (triclinic == 1) 
     error->all(FLERR,"Compute XRD does not work with triclinic structures");

  array_flag = 1;
  extarray = 0;

  // Store radiation wavelength
  lambda = atof(arg[3]);
  if (lambda <= 0)
    error->all(FLERR,"Compute SAED: Wavelength must be greater than zero");

  // Define atom types for atomic scattering factor coefficents
  int iarg = 4;  
  ztype = new int[ntypes];
  for (int i = 0; i < ntypes; i++){
    ztype[i] = XRDmaxType + 1;
  }
  for (int i=0; i<ntypes; i++) {   
    for(int j = 0; j < XRDmaxType; j++){
      if (strcasecmp(arg[iarg],XRDtypeList[j]) == 0) {
        ztype[i] = j;
       }
     }
    if ( ztype[i] == XRDmaxType + 1 )
        error->all(FLERR,"Compute XRD: Invalid ASF atom type");
    iarg++;
  }

  // Set defaults for optional args
  Min2Theta = 1;
  Max2Theta = 179;  
  radflag ==1;
  c[0] = 1; c[1] = 1; c[2] = 1;
  LP = 1;
  manual = false;
  echo = false;

  // Process optional args
  while (iarg < narg) {
    if (strcmp(arg[iarg],"2Theta") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal Compute XRD Command");
      Min2Theta = atof(arg[iarg+1]) / 2;
      Max2Theta = atof(arg[iarg+2]) / 2;
      if (Max2Theta > MY_PI ){
        Min2Theta = Min2Theta * MY_PI / 180;  // converting to radians if necessary
        Max2Theta = Max2Theta * MY_PI / 180;
        radflag = 0;
      }
      if (Min2Theta <= 0) 
        error->all(FLERR,"Minimum 2theta value must be greater than zero");
      if (Max2Theta > MY_PI ) 
        error->all(FLERR,"Maximum 2theta value must be less than 180 degrees");
      if (Max2Theta-Min2Theta <= 0) 
        error->all(FLERR,"Two-theta range must be greater than zero");
      iarg += 3;

    } else if (strcmp(arg[iarg],"c") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal Compute XRD Command");
      c[0] = atof(arg[iarg+1]);
      c[1] = atof(arg[iarg+2]);
      c[2] = atof(arg[iarg+3]);
      if (c[0] <= 0 || c[1] <= 0 || c[2] <= 0) 
        error->all(FLERR,"Compute XRD: c's must be greater than 0");  
      iarg += 4;
      
    } else if (strcmp(arg[iarg],"LP") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal Compute XRD Command");
      LP = atof(arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"echo") == 0) {
      echo = true;
      iarg += 1;

    } else if (strcmp(arg[iarg],"manual") == 0) {
      manual = true;
      iarg += 1;        
      
    } else error->all(FLERR,"Illegal Compute XRD Command");
  }
 
  Kmax = 2 * sin(Max2Theta) / lambda;
 
  // Calculating spacing between reciprical lattice points using the prd
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
  
  // Find reprical spacing and integer dimensions
  for (int i=0; i<3; i++) {
    dK[i] = prd_inv[i]*c[i];
    Knmax[i] = ceil(Kmax / dK[i]);
  } 
  
  // Finding the intersection of the reciprical space and Ewald sphere
  int nRows = 0;
  double dinv2= 0.0;
  double r = 0.0;
  double ang = 0.0;
  double K[3];
  
  // Procedure to determine how many rows are needed given the constraints on 2theta
  for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
    for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
      for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
        
        K[0] = i * dK[0];
        K[1] = j * dK[1];
        K[2] = k * dK[2];
        dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
        if  (4 >= dinv2 * lambda * lambda ) {
       	  ang = asin(lambda * sqrt(dinv2) / 2);
          if (ang <= Max2Theta & ang >= Min2Theta) {
          nRows++;
	        }
        }
      } 
    }
  } 

  size_array_rows = nRows;
  size_array_cols = 2;
  
  if (me == 0) {
    if (screen &&  echo)
      fprintf(screen,"-----\nCompute XRD id:%s, # of atoms:%d, # of relp:%d\n",id,natoms,nRows);
      fprintf(screen,"Reciprocal point spacing in k1,k2,k3 = %g %g %g\n-----\n",
              dK[0], dK[1], dK[2]);
  }  
   
  memory->create(array,size_array_rows,size_array_cols,"xrd:array");
  memory->create(store_tmp,3*size_array_rows,"xrd:store_tmp");
}

/* ---------------------------------------------------------------------- */

ComputeXRD::~ComputeXRD()
{
  memory->destroy(array);
  memory->destroy(store_tmp);
  delete ztype;
}

/* ---------------------------------------------------------------------- */

void ComputeXRD::init()
{

  int mmax = (2*Knmax[0]+1)*(2*Knmax[1]+1)*(2*Knmax[2]+1);
  double K[3];
  double dinv2 = 0.0;
  double ang =0.0;

  double convf = 360 / MY_PI;
  if (radflag ==1){
  convf = 1;
  }

  int n = 0;
    for (int m = 0; m < mmax; m++) {
      int k = m%(2*Knmax[2]+1);
      int j = (m%((2*Knmax[2]+1)*(2*Knmax[1]+1))-k)/(2*Knmax[2]+1);
      int i = (m-j*(2*Knmax[2]+1)-k)/((2*Knmax[2]+1)*(2*Knmax[1]+1))-Knmax[0];
      j = j-Knmax[1];
      k = k-Knmax[2];
      K[0] = i * dK[0];
      K[1] = j * dK[1];
      K[2] = k * dK[2];
      dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
      if  (4 >= dinv2 * lambda * lambda ) {
         ang = asin(lambda * sqrt(dinv2) / 2);
         if (ang <= Max2Theta & ang >= Min2Theta) {
            store_tmp[3*n] = k;
            store_tmp[3*n+1] = j;
            store_tmp[3*n+2] = i;
            array[n][0] = ang * convf;
            n++;
         }
      }
   }

  if (n != size_array_rows)  
     error->all(FLERR,"Compute XRD compute_array() rows mismatch");

}

/* ---------------------------------------------------------------------- */

void ComputeXRD::compute_array()
{
  invoked_array = update->ntimestep;

  if (me == 0 && echo) {
      if (screen)
        fprintf(screen,"-----\nComputing XRD intensities");
  }

  double t0 = MPI_Wtime();

  double *Fvec1 = new double[size_array_rows]; // Strct factor 
  double *Fvec2 = new double[size_array_rows]; // Strct factor  (imaginary)
  // -- Note: array rows correspond to different RELP
 
  ntypes = atom->ntypes;
  nlocal = atom->nlocal;
  int *type  = atom->type;
  int natoms = group->count(igroup);
  int *mask = atom->mask;

  double *x = new double [3*nlocal];

  for (int ii = 0; ii < nlocal; ii++) {
     x[3*ii+0] = atom->x[ii][0];
     x[3*ii+1] = atom->x[ii][1];
     x[3*ii+2] = atom->x[ii][2];
  }

// Setting up OMP
  int nthreads = 1;
#ifdef _OPENMP
  nthreads = comm->nthreads;

  if (me == 0 && echo) {
    if (screen)
      fprintf(screen," using %d OMP threads",nthreads);
  }
#endif

  if (me == 0 && echo) {
    if (screen)
      fprintf(screen,"\n");
  }
  int m = 0;
  double frac = 0.1;

#pragma omp parallel num_threads(nthreads)
  {
    double *f = new double[ntypes];    // atomic structure factor by type
    int typei = 0;
    double Fatom1 = 0.0;               // structure factor per atom
    double Fatom2 = 0.0;               // structure factor per atom (imaginary)
    double S = 0.0;                    // sin(theta)/lambda
    double K[3];
    double dinv2 = 0.0;
    double ang =0.0;
    double inners = 0.0;
    double lp = 0.0;

    if (LP == 1) {
#pragma omp for
      for (int n = 0; n < size_array_rows; n++) {
        int k = store_tmp[3*n+0];
        int j = store_tmp[3*n+1];
        int i = store_tmp[3*n+2];
        K[0] = i * dK[0];
        K[1] = j * dK[1];
        K[2] = k * dK[2];
        dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
        ang = asin(lambda * sqrt(dinv2) / 2);

        Fatom1 = 0.0;
        Fatom2 = 0.0;

        // Calculate the atomic structre factor by type
        S = sin(ang) / lambda;
        for (int ii = 0; ii < ntypes; ii++){
          f[ii] = 0;
          int C = 0;
          while (C < 8) {
            f[ii] += ASFXRD[ztype[ii]][C] * exp(-1 * ASFXRD[ztype[ii]][C+1] * S * S);
            C += 2;
          }
          f[ii] += ASFXRD[ztype[ii]][8];
        }

        // Evaluate the structure factor equation -- looping over all atoms
        for (int ii = 0; ii < nlocal; ii++){
          typei=type[ii]-1;
          if (mask[ii] & groupbit) {
            inners = 2 * MY_PI * (K[0] * x[3*ii+0] + K[1] * x[3*ii+1] +
                      K[2] * x[3*ii+2]);
            Fatom1 += f[typei] * cos(inners);
            Fatom2 += f[typei] * sin(inners);
          }
        }
        lp = (1 + cos( 2 * ang ) * cos( 2 * ang )) /
             ( cos( ang ) * sin ( ang ) * sin( ang ));
        Fvec1[n] = Fatom1 * lp;
        Fvec2[n] = Fatom2 * lp;

        // reporting progress of calculation
        if ( echo ) {
          #pragma omp critical
          {
            if ( m == round(frac * size_array_rows) ) {
              if (me == 0 && screen) fprintf(screen," %0.0f%% -",frac*100);
              frac += 0.1;
            }
            m++;
          }
        }      
      } // End of pragma omp for region

    } else {
#pragma omp for
      for (int n = 0; n < size_array_rows; n++) {
        int k = store_tmp[3*n+0];
        int j = store_tmp[3*n+1];
        int i = store_tmp[3*n+2];
        K[0] = i * dK[0];
        K[1] = j * dK[1];
        K[2] = k * dK[2];
        dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
        ang = asin(lambda * sqrt(dinv2) / 2);

        Fatom1 = 0.0;
        Fatom2 = 0.0;

        // Calculate the atomic structre factor by type
        S = sin(ang) / lambda;
        for (int ii = 0; ii < ntypes; ii++){
          f[ii] = 0;
          int C = 0;
          while (C < 8) {
            f[ii] += ASFXRD[ztype[ii]][C] * exp(-1 * ASFXRD[ztype[ii]][C+1] * S * S);
            C += 2;
          }
          f[ii] += ASFXRD[ztype[ii]][8];
        }

        // Evaluate the structure factor equation -- looping over all atoms
        for (int ii = 0; ii < nlocal; ii++){
          typei=type[ii]-1;
          if (mask[ii] & groupbit) {
            inners = 2 * MY_PI * (K[0] * x[3*ii+0] + K[1] * x[3*ii+1] +
                      K[2] * x[3*ii+2]);
            Fatom1 += f[typei] * cos(inners);
            Fatom2 += f[typei] * sin(inners);
          }
        }
        Fvec1[n] = Fatom1;
        Fvec2[n] = Fatom2;

        // reporting progress of calculation
        if ( echo ) {
          #pragma omp critical
          {
            if ( m == round(frac * size_array_rows) ) {
              if (me == 0 && screen) fprintf(screen," %0.0f%% -",frac*100 );
              frac += 0.1;
            }
            m++;
          }
        }
      } // End of pragma omp for region
    } // End of if LP=1 check 
    delete [] f;
  } // End of pragma omp parallel region

  double *scratch1 = new double[size_array_rows];
  double *scratch2 = new double[size_array_rows];

  // Sum intensity for each ang-hkl combination across processors
  MPI_Allreduce(Fvec1,scratch1,size_array_rows,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(Fvec2,scratch2,size_array_rows,MPI_DOUBLE,MPI_SUM,world);

  for (int i = 0; i < size_array_rows; i++) {
    array[i][1] = (scratch1[i] * scratch1[i] + scratch2[i] * scratch2[i]) / natoms;
  }

  double t2 = MPI_Wtime();

  // compute memory usage per processor
  double bytes = size_array_rows * size_array_cols * sizeof(double); //array
  bytes +=  4.0 * size_array_rows * sizeof(double); //Fvec1 & 2, scratch1 & 2
  bytes += ntypes * sizeof(double); // f
  bytes += 3.0 * nlocal * sizeof(double); // x
  bytes += 3.0 * size_array_rows * sizeof(int); // store_temp
  
  if (me == 0 && echo) {
    if (screen)
      fprintf(screen," 100%% \nTime ellapsed during compute_xrd = %0.2f sec using %0.2f Mbytes/processor\n-----\n", t2-t0, bytes/1024.0/1024.0);
  }

  delete [] scratch1;
  delete [] scratch2;
  delete [] Fvec1;
  delete [] Fvec2;
  delete [] x;
}

/* ----------------------------------------------------------------------
 memory usage of arrays
 ------------------------------------------------------------------------- */

double ComputeXRD::memory_usage()
{
  double bytes = size_array_rows * size_array_cols * sizeof(double); //array
  bytes +=  4.0 * size_array_rows * sizeof(double); //Fvec1 & 2, scratch1 & 2
  bytes += 3.0 * nlocal * sizeof(double); // x
  bytes += ntypes * sizeof(double); // f
  bytes += 3.0 * size_array_rows * sizeof(int); // store_temp
  
  return bytes;
}

