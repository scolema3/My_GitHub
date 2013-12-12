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
                       : Eric Homer (BYU)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "math_const.h"
#include "compute_saedmic.h"
#include "compute_saed_consts.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "group.h"
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

/* ---------------------------------------------------------------------- */

ComputeSAEDMIC::ComputeSAEDMIC(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  MPI_Comm_rank(world,&me);
  int ntypes = atom->ntypes;
  int natoms = group->count(igroup);
  int dimension = domain->dimension;
  int *periodicity = domain->periodicity;
  int triclinic = domain->triclinic;

  // Checking errors specific to the compute
  if (dimension == 2) 
    error->all(FLERR,"Compute SAED does not work with 2d structures");
  if (narg < 4) 
    error->all(FLERR,"Illegal Compute SAED Command");
  if (!periodicity[0] && !periodicity[1] && !periodicity[2])
    error->all(FLERR,"Compute SAED must have at least one periodic boundary");
  if (triclinic == 1) 
    error->all(FLERR,"Compute SAED does not work with triclinic structures"); 
  
  vector_flag = 1;
  extvector = 0;
 
  // store the wavelength
  lambda = atof(arg[3]);
  if (lambda <= 0)
    error->all(FLERR,"Compute SAED: Wavelength must be greater than zero");
  
  //store the z number of the different atom types
  ztype = new int[ntypes];
  int iarg = 4;
  for (int i=0; i<ntypes; i++) {
    ztype[i] = atoi(arg[iarg]) - 1; //subtract one to index zero-based array
    if ( (ztype[i] < 0) || (ztype[i] > (ASFmaxZ - 1) ) )
      error->all(FLERR,"Compute SAED: Invalid atomic number for atom type");
    iarg++;
  }

  // Set defaults for optional args
  Kmax = 1.70;
  Zone[0] = 1; Zone[1] = 0; Zone[2] = 0;
  c[0] = 1; c[1] = 1; c[2] = 1;
  dR_Ewald = 0.01 / 2; 

  // Process optional args
   
  while (iarg < narg) {
    
    if (strcmp(arg[iarg],"Kmax") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal Compute SAED Command");
      Kmax = atof(arg[iarg+1]);
      if (Kmax / 2 <= 0 || Kmax / 2 > 6) 
        error->all(FLERR,"Compute SAED: |K|max/2 must be between 0 and 6 ");
      iarg += 2;
        
    } else if (strcmp(arg[iarg],"Zone") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal Compute SAED Command");
      Zone[0] = atof(arg[iarg+1]);
      Zone[1] = atof(arg[iarg+2]);
      Zone[2] = atof(arg[iarg+3]);
      iarg += 4;
      
    } else if (strcmp(arg[iarg],"c") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal Compute SAED Command");
      c[0] = atof(arg[iarg+1]);
      c[1] = atof(arg[iarg+2]);
      c[2] = atof(arg[iarg+3]);
      if (c[0] <= 0 || c[1] <= 0 || c[2] <= 0) 
        error->all(FLERR,"Compute SAED: dKs must be greater than 0");  
      iarg += 4;
      
    } else if (strcmp(arg[iarg],"dR_Ewald") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal Compute SAED Command");
      dR_Ewald = atof(arg[iarg+1]);
      if (dR_Ewald < 0) 
        error->all(FLERR,"Compute SAED: dR_Ewald slice must be greater than 0"); 
      iarg += 2;
    } else if (strcmp(arg[iarg],"echo") == 0) {
      echo = true;
      iarg += 1;
    } else if (strcmp(arg[iarg],"manual") == 0) {
      manual = true;
      iarg += 1;        
                
    } else error->all(FLERR,"Illegal Compute SAED Command");
  }
  
  // Zone flag to capture entire recrocal space volume
  if (  (Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0) ){
  } else {
      R_Ewald = (1 / lambda); 
      double Rnorm = R_Ewald/ sqrt(Zone[0] * Zone[0] + 
                     Zone[1] * Zone[1] +  Zone[2]* Zone[2]);
      Zone[0] = Zone[0] * Rnorm; 
      Zone[1] = Zone[1] * Rnorm; 
      Zone[2] = Zone[2] * Rnorm; 
  }

  // Procedure to determine how many rows are needed given the constraints on 2theta
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
  int n = 0;
  double dinv2 = 0.0;
  double r = 0.0;
  double K[3];
  
  // Zone flag to capture entire recrocal space volume
  if ( (Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0) ){
    for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
      for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
        for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
          K[0] = i * dK[0];
          K[1] = j * dK[1];
          K[2] = k * dK[2];
          dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
          if (dinv2 < Kmax * Kmax) n++;
        } 
      }
    } 
  } else {
    for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
      for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
        for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
          K[0] = i * dK[0];
          K[1] = j * dK[1];
          K[2] = k * dK[2];
          dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
          if (dinv2 < Kmax * Kmax) {
            r=0.0;
            for (int m=0; m<3; m++)
              r += pow(K[m] - Zone[m],2.0);
              r = sqrt(r);          
            if  ( (r >  (R_Ewald - dR_Ewald) ) && (r < (R_Ewald + dR_Ewald) ) ){
              n++;
            }
          }
        } 
      }
    } 
  }


  if (me == 0) {
    if (screen &&  echo)
      fprintf(screen,"-----\nCompute SAED id:%s, # of atoms:%d, # of relp:%d\n",id,natoms,n);
      fprintf(screen,"Reciprocal point spacing in k1,k2,k3 = %g %g %g\n-----\n",
              dK[0], dK[1], dK[2]);
  }  

  nRows = n;
  size_vector = n;
  memory->create(vector,size_vector,"saedmic:vector");
//  store_tmp = new int[3 * nRows];
  memory->create(store_tmp,3*size_vector,"saed:store_tmp");
}

/* ---------------------------------------------------------------------- */

ComputeSAEDMIC::~ComputeSAEDMIC()
{
  memory->destroy(vector);
  delete ztype;
 // delete [] store_tmp;
  memory->destroy(store_tmp);
}

/* ---------------------------------------------------------------------- */

void ComputeSAEDMIC::init()
{


  double dinv2, r;
  double K[3];
  int n = 0;
  
  // Zone 0 0 0 flag to capture entire recrocal space volume
  if ( (Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0) ){
    for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
      for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
        for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
          K[0] = i * dK[0];
          K[1] = j * dK[1];
          K[2] = k * dK[2];
          dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
          if (dinv2 < Kmax * Kmax) {
            store_tmp[3*n] = i;
            store_tmp[3*n+1] = j;
            store_tmp[3*n+2] = k;
            n++;
          }
        } 
      } 
    }
  } else {
    for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
      for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
        for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
          K[0] = i * dK[0];
          K[1] = j * dK[1];
          K[2] = k * dK[2];
          dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
          if (dinv2 < Kmax * Kmax) {
            r=0.0;
            for (int m=0; m<3; m++)
              r += pow(K[m] - Zone[m],2.0);
            r = sqrt(r);          
            if  ( (r >  (R_Ewald - dR_Ewald) ) && (r < (R_Ewald + dR_Ewald) ) ) {
              store_tmp[3*n] = i;
              store_tmp[3*n+1] = j;
              store_tmp[3*n+2] = k;
              n++;
            }
          }
        }
      } 
    } 
  }

  if (n != nRows)  error->all(FLERR,"Compute SAED Nrows inconsistent");

}

/* ---------------------------------------------------------------------- */

void ComputeSAEDMIC::compute_vector()
{
  invoked_vector = update->ntimestep;

  if (me == 0 && echo) {
      if (screen)
        fprintf(screen,"-----\nComputing SAED intensities");
  }

  double t0 = MPI_Wtime();
  double *Fvec1 = new double[nRows]; // Strct factor 
  double *Fvec2 = new double[nRows]; // Strct factor  (imaginary)
  // -- Note, vector entries correspond to different RELP


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

 // determining paramater set to use based on maximum S = sin(theta)/lambda
  double Smax = Kmax / 2;
  
  int offset = 0;                 // offset the ASF matrix for appropriate value
  if (Smax <= 2) offset = 0;
  if (Smax > 2)  offset = 10;     

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
    double r = 0.0; 
    double dinv2 = 0.0;
    double ang =0.0;
    double inners = 0.0;

#pragma omp for
    for (int n = 0; n < nRows; n++) {
      int i = store_tmp[3*n+0];
      int j = store_tmp[3*n+1];
      int k = store_tmp[3*n+2];
      K[0] = i * dK[0];
      K[1] = j * dK[1];
      K[2] = k * dK[2];
      dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
      ang = asin(lambda * sqrt(dinv2) / 2);

      Fatom1 = 0.0;
      Fatom2 = 0.0;

      // Calculate the atomic structre factor by type
      S = sin(ang) / lambda;
      // determining paramater set to use based on S = sin(theta)/lambda <> 2
      for (int ii = 0; ii < ntypes; ii++){
        f[ii] = 0;
        for (int C = 0; C < 5; C++){
          int D = C + offset;
          f[ii] += ASF[ztype[ii]][D] * exp(-1*ASF[ztype[ii]][5+D]*S*S);
        }
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
          if ( m == round(frac * nRows) ) {
            if (me == 0 && screen) fprintf(screen," %0.0f%% -",frac*100);
            frac += 0.1;
          }
          m++;
        }
      }      
    } // End of pragma omp for region
    delete [] f;
  }

  double *scratch1 = new double[nRows];
  double *scratch2 = new double[nRows];

  // Sum intensity for each ang-hkl combination across processors
  MPI_Allreduce(Fvec1,scratch1,nRows,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(Fvec2,scratch2,nRows,MPI_DOUBLE,MPI_SUM,world);


// #pragma omp parallel for
  for (int i = 0; i < nRows; i++) {
    vector[i] = (scratch1[i] * scratch1[i] + scratch2[i] * scratch2[i]) / natoms;
  }


  double t2 = MPI_Wtime();

  // compute memory usage per processor
  double bytes = nRows * sizeof(double); //vector
  bytes +=  4.0 * nRows * sizeof(double); //Fvec1 & 2, scratch1 & 2
  bytes += ntypes * sizeof(double); // f
  bytes += 3.0 * nlocal * sizeof(double); // x
  bytes += 3.0 * nRows * sizeof(int); // store_temp
  
  if (me == 0 && echo) {
    if (screen)
      fprintf(screen," 100%% \nTime ellapsed during compute_saedmic = %0.2f sec using %0.2f Mbytes/processor\n-----\n", t2-t0,  bytes/1024.0/1024.0);
  }


  delete [] x;  
  delete [] scratch1;
  delete [] scratch2;
  delete [] Fvec1;
  delete [] Fvec2;
}

/* ----------------------------------------------------------------------
 memory usage of arrays
 ------------------------------------------------------------------------- */

double ComputeSAEDMIC::memory_usage()
{
  double bytes = nRows * sizeof(double); //vector
  bytes +=  4.0 * nRows * sizeof(double); //Fvec1 & 2, scratch1 & 2
  bytes += ntypes * sizeof(double); // f
  bytes += 3.0 * nlocal * sizeof(double); // x
  bytes += 3.0 * nRows * sizeof(int); // store_temp
  return bytes;
}
