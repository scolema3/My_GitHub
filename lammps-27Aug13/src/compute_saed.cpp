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
#include "compute_saed.h"
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

using namespace LAMMPS_NS;
using namespace std;

/* ---------------------------------------------------------------------- */

ComputeSAED::ComputeSAED(LAMMPS *lmp, int narg, char **arg) :
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
  
  // Find reprical spacing and integer dimensions
  for (int i=0; i<3; i++) {
    dK[i] = prd_inv[i]*c[i];
    Knmax[i] = ceil(Kmax / dK[i]);
  }
  
  // Finding the intersection of the reciprical space and Ewald sphere
  nRows = 0;
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
          if (dinv2 < Kmax * Kmax) nRows++;
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
              nRows++;
            }
          }
        } 
      }
    } 
  }
 

  if (me == 0) {
    if (logfile)
      fprintf(logfile,"Compute SAED id:%s, # of atoms:%d, # of relp:%d\n",id,natoms,nRows);
      fprintf(logfile,"Reciprocal point spacing in k1,k2,k3 = %g %g %g\n",
              dK[0], dK[1], dK[2]);
  }

  size_vector = nRows;
  memory->create(vector,size_vector,"saed:vector");
}

/* ---------------------------------------------------------------------- */

ComputeSAED::~ComputeSAED()
{
  memory->destroy(vector);
  delete ztype;
}

/* ---------------------------------------------------------------------- */

void ComputeSAED::init()
{
  // Finding the intersection of the reciprical space and Ewald sphere
  double K_dist, dinv2, r;
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
            K_dist = sqrt(dinv2);
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
              K_dist = sqrt(dinv2);
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

void ComputeSAED::compute_vector()
{
  invoked_vector = update->ntimestep;
 
  double **x = atom->x;
  int ntypes = atom->ntypes;
  int *type  = atom->type;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int natoms = group->count(igroup);
   	
  int typei = 0;
  double PI = 4.0*atan(1.0);
  double angle=0.0;                     // Diffraction angle (radians)
  double S = 0.0;                       // sin(theta)/lambda -- used to 
                                        // determine atomic structure factor
  f = new double[ntypes];               // atomic structure factor by type

  double Fatom1 = 0.0;                  // structure factor per atom
  double Fatom2 = 0.0;                  // structure factor per atom (imaginary)

  double* Fvec1  = NULL;
  double* Fvec2  = NULL;
  Fvec1 = new double[nRows];            // Structure factor vector
  Fvec2 = new double[nRows];            // Structure factor vector (imaginary)

  double* scratch1  = NULL;
  double* scratch2  = NULL;
  scratch1 = new double[nRows];         // Scratch vectors used in mpi_allreduce
  scratch2 = new double[nRows];

  double K_dist, dinv2, r;
  double K[3];

  // determining paramater set to use based on maximum S = sin(theta)/lambda
  double Smax = Kmax / 2;
  
  int offset = 0;                 // offset the ASF matrix for appropriate value
  if (Smax <= 2) offset = 0;
  if (Smax > 2)  offset = 10;     

  double frac = 0.1;
  if (me == 0 && echo) {
    if (logfile) fprintf(logfile,"Beginning SAED compute\n");
  }

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
            K_dist = sqrt(dinv2);
            angle =  asin(lambda * K_dist / 2); // Angle - formally array[n][3]
            Fatom1 = 0.0;
            Fatom2 = 0.0;

            // Calculate the atomic structre factor by type
            S = sin(angle) / lambda;

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
              if (mask[ii] & groupbit) {
              typei=type[ii]-1;
              Fatom1 += f[typei] * cos(2 * PI * (K[0] * x[ii][0] +
                             K[1] * x[ii][1] + K[2] * x[ii][2]));
              Fatom2 += f[typei] * sin(2 * PI * (K[0] * x[ii][0] +
                             K[1] * x[ii][1] + K[2] * x[ii][2]));
              }
            }

            Fvec1[n] = Fatom1;
            Fvec2[n] = Fatom2;

            // Reporting progress of computation progress
            if ( n >= (frac * nRows) ) {
              frac += 0.1;
              if (me == 0 && echo) {
                if (logfile) fprintf(logfile," .");
              }
            }

            n++;

          } // Kmax 
        } // l-loop
      } //k-loop
    } //h-loop
   
  } else {
  // Special loop for specific zone axis

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

            // Determine if RELP is within specified zone axis selection
            if  ( (r >  (R_Ewald - dR_Ewald) ) && (r < (R_Ewald + dR_Ewald) ) ) {
              K_dist = sqrt(dinv2);
              angle =  asin(lambda * K_dist / 2); // Angle - formally array[n][3]
              Fatom1 = 0.0;
              Fatom2 = 0.0;

              // Calculate the atomic structre factor by type
              S = sin(angle) / lambda;

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
                if (mask[ii] & groupbit) {
                typei=type[ii]-1;
                Fatom1 += f[typei] * cos(2 * PI * (K[0] * x[ii][0] +
                             K[1] * x[ii][1] + K[2] * x[ii][2]));
                Fatom2 += f[typei] * sin(2 * PI * (K[0] * x[ii][0] +
                             K[1] * x[ii][1] + K[2] * x[ii][2]));
                }
              }

              Fvec1[n] = Fatom1;
              Fvec2[n] = Fatom2;

              // Reporting progress of computation progress
              if ( n >= (frac * nRows) ) {
                frac += 0.1;
                if (me == 0 && echo) {
                  if (logfile) fprintf(logfile," .");
                }
              }

            n++;

            } // Zone
          } // Kmax
        } // l-loop
      } //k-loop
    } //h-loop
  }

  // Sum intensity for each ang-hkl combination across processors
  MPI_Allreduce(Fvec1,scratch1,nRows,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(Fvec2,scratch2,nRows,MPI_DOUBLE,MPI_SUM,world);

  for (int i = 0; i < nRows; i++) {
    vector[i] = (scratch1[i] * scratch1[i] + scratch2[i] * scratch2[i]) / natoms;
  }
  
  delete [] scratch1;
  delete [] scratch2;
  delete [] Fvec1;
  delete [] Fvec2;

  if (me == 0 && echo) {
    if (logfile)
      fprintf(logfile,"\nFinished SAED compute array\n");
  }
}

/* ----------------------------------------------------------------------
 memory usage of arrays
 ------------------------------------------------------------------------- */

double ComputeSAED::memory_usage()
{
  double bytes = nRows * sizeof(double); // vector
  bytes +=  4.0 * nRows * sizeof(double); //Fvec1 & 2, scratch1 & 2
  return bytes;
}
