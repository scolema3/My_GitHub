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
   Contributing authors: Shawn Coleman (Arkansas)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "compute_xrdmic.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "stdio.h"
#include "string.h"
#include "omp.h"

#ifdef ENABLE_MIC
#include "offload.h"
#endif

#include <iostream>

#define __NUM_THREADS 240

using namespace LAMMPS_NS;
using namespace std;

/* ---------------------------------------------------------------------- */

ComputeXRDMIC::ComputeXRDMIC(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  MPI_Comm_rank(world,&me);
  ntypes = atom->ntypes;
  int natoms = group->count(igroup);
  int dimension = domain->dimension;
  int *periodicity = domain->periodicity;
  int triclinic = domain->triclinic;

  // Checking errors
  if (dimension == 2) 
     error->all(FLERR,"Compute XRD does not work with 2d structures");
  if (narg < 4) 
     error->all(FLERR,"Illegal Compute XRD Command");
  if (!periodicity[0] && !periodicity[1] && !periodicity[2])
     error->all(FLERR,"Compute XRD must have at least one periodic boundary");
  if (triclinic == 1) 
     error->all(FLERR,"Compute XRD does not work with triclinic structures");

  array_flag = 1;
  extarray = 0;

  // Grabbing information from a file to determine atomic scattering factors
  int n = strlen(arg[3]) + 1;
  ASFfilename = new char[n];
  strcpy(ASFfilename,arg[3]);  

  char line[512];
  char *result;
  int count;

  // Reading in the coeffiecients for analytical approximation to the scattering factor
  // Found in International Tables For Crystallography V. C (Table 6.11.1.4)
  FILE *infile = fopen(ASFfilename,"r");
  if (infile == NULL) error->all(FLERR,"Compute XRD - ASF file failed to open");
  
  ASF = new double[ntypes * 9];

  for (int i = 0; i < ntypes; i++) {
    result = fgets(line,512,infile);
    if (!result) error->all(FLERR,"Compute XRD - ASF file does not contain enough rows");
    count = sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf", &ASF[i*9+0], &ASF[i*9+1], 
    	&ASF[i*9+2], &ASF[i*9+3], &ASF[i*9+4], &ASF[i*9+5], &ASF[i*9+6], &ASF[i*9+8], &ASF[i*9+8]);
    if (count != 9) error->all(FLERR,"Compute XRD - ASF file does not contain enough columns");
  }
  fclose(infile);
  
  // Set defaults for optional args
  lambda = 1.541838;
  Min2Theta = 0.00872664626;
  Max2Theta = 1.56206968;  
  c[0] = 1; c[1] = 1; c[2] = 1;
  LP = 1;
  
  
  double PI = 4.0*atan(1.0);
   // Process optional args
   
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"lambda") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal Compute XRD Command");
      lambda = atof(arg[iarg+1]);
      if (lambda < 0) 
        error->all(FLERR,"Compute XRD: Wavelength must be greater than zero");
      iarg += 2;
      
    } else if (strcmp(arg[iarg],"2Theta") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal Compute XRD Command");
      Min2Theta = atof(arg[iarg+1]) / 2;
      Max2Theta = atof(arg[iarg+2]) / 2;
      if (Max2Theta > PI){
        Min2Theta = Min2Theta * PI / 180;  // converting to radians if necessary
        Max2Theta = Max2Theta * PI / 180;
      }
      if (Min2Theta < 0) 
        error->all(FLERR,"Minimum 2theta value must be greater than zero");
      if (Max2Theta > PI) 
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
    if (screen)
      fprintf(screen,"Reciprocal point spacing in k1,k2,k3 = %g %g %g\n",
              dK[0], dK[1], dK[2]);
    if (logfile)
      fprintf(logfile,"Reciprocal point spacing in k1,k2,k3 = %g %g %g\n",
              dK[0], dK[1], dK[2]);
  }

  if (me == 0 && echo) {
      if (screen)
        fprintf(screen,"Compute XRD id:%s, # of atoms:%d, # of relp:%d\n",id,natoms,nRows);
      if (logfile)
        fprintf(logfile,"Compute XRD id:%s, # of atoms:%d, # of relp:%d\n",id,natoms,nRows);
  }  

  my_mpi_rank = 0;
  num_mpi_procs = 1;

  size_array_rows_loc = size_array_rows/num_mpi_procs;
  size_array_rows_mod = size_array_rows%num_mpi_procs;
  if (my_mpi_rank < size_array_rows_mod) {
     size_array_rows_loc++;
  }

  if (my_mpi_rank == 0) {
     cout<<"size_array_rows    ="<<size_array_rows<<"\n";
     cout<<"size_array_cols    ="<<size_array_cols<<"\n";
     cout<<"size_array_rows_loc="<<size_array_rows_loc<<"\n";
  }
   
  memory->create(array,size_array_rows_loc,size_array_cols,"xrd:array");
  array_loc = array[0];
}

/* ---------------------------------------------------------------------- */

ComputeXRDMIC::~ComputeXRDMIC()
{
  fprintf(screen," Compute XRD Complete 3.\n");
  memory->destroy(array);
  delete [] ASFfilename;

}

/* ---------------------------------------------------------------------- */

void ComputeXRDMIC::init()
{

}

/* ---------------------------------------------------------------------- */

void ComputeXRDMIC::compute_array()
{
  invoked_array = update->ntimestep;

  Fvec1 = new double[size_array_rows_loc]; // Structure factor vector 
  Fvec2 = new double[size_array_rows_loc]; // Structure factor vector (imaginary)
                                           // -- Vector entries correspond to 
                                           //    diffraction angles
 
  ntypes = atom->ntypes;
  nlocal = atom->nlocal;
  type  = atom->type;
  int natoms = group->count(igroup);
  mask = atom->mask;

  f = new double[ntypes];        	// atomic structure factor by type

  if (me == 0 && echo) {
      if (screen)
        fprintf(screen,"Computing XRD intensities");
      if (logfile)
        fprintf(screen,"Computing XRD intensities");
  }
 
  double t0 = MPI_Wtime();

  double *x = new double [3*nlocal];
  for (int ii = 0; ii < nlocal; ii++) {
      x[ii+0] = atom->x[ii][0];
      x[ii+1] = atom->x[ii][1];
      x[ii+2] = atom->x[ii][2];
   }

char signal_var;

#ifdef ENABLE_MIC && __INTEL_OFFLOAD

#pragma offload target(mic:0) \
   in(lambda) \
   in(size_array_rows_loc) \
   in(size_array_cols) \
   in(nlocal) \
   in(ntypes) \
   in(groupbit) \
   in(x : length(3*nlocal)) \
   in(ASF : length(9*ntypes)) \
   in(array_loc : length(size_array_rows_loc*size_array_cols)) \
   in(mask : length(nlocal)) \ 
   in(type : length(nlocal)) \
   out(Fvec1 : length(size_array_rows_loc)) \
   out(Fvec2 : length(size_array_rows_loc)) \
   out(f : length(ntypes)) \
   signal(&signal_var)

#endif
{  // Start of the MIC region

#ifdef __MIC__
   printf("After pragma using offload from coprocessor, ntypes=%d\n",ntypes);
#else
   printf("After pragma using offload from processor, ntypes=%d\n",ntypes);
#endif
  double dinv2 = 0.0;
  double ang =0.0;
  int n = 0;
  int m = 0;

  double K[3];
  double S = 0.0;                       // sin(theta)/lambda -- used to 
                                        // determine atomic structure factor
  double PI = 4.0*atan(1.0);
  int typei = 0;

  double Fatom1 = 0.0;                  // structure factor per atom
  double Fatom2 = 0.0;                  // structure factor per atom (imaginary)

  int starting_row;
  if (my_mpi_rank < size_array_rows_mod) {
     starting_row = my_mpi_rank*size_array_rows_loc;
  }
  else {
     starting_row = my_mpi_rank*size_array_rows_loc+size_array_rows_mod;
  }

#pragma omp parallel num_threads(__NUM_THREADS)
   {
#pragma omp for
  // looping over all incident angles
    for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
      for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
        for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
        
          K[0] = i * dK[0];
          K[1] = j * dK[1];
          K[2] = k * dK[2];
          dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
          if  (4 >= dinv2 * lambda * lambda && n < starting_row+size_array_rows_loc) {
       	    ang = asin(lambda * sqrt(dinv2) / 2);
            if (ang <= Max2Theta & ang >= Min2Theta) {
              if (n >= starting_row) {
              
                array_loc[i*size_array_cols+0] = ang;
                
                Fatom1 = 0.0;
                Fatom2 = 0.0;
  
                // Calculate the atomic structre factor by type	
                S = sin(ang) / lambda;    	
                for (int ii = 0; ii < ntypes; ii++){
                  f[ii] = 0;
                  int C = ii * 9;
                  int D = C+8;
                  while (C < D) {
                    f[ii] += ASF[C] * exp( -1 * ASF[C+1] * S * S);
                    C += 2;
                  }
                  f[ii] += ASF[D];
                }
  
                // Evaluate the structure factor equation -- looping over all atoms
                for (int ii = 0; ii < nlocal; ii++){
                  if (mask[ii] & groupbit) {
                    typei=type[ii]-1;       
                    Fatom1 += f[typei] * cos(2 * PI * (K[0] * x[3*ii] + 
                                                       K[1] * x[3*ii+1] +
                                                       K[2] * x[3*ii+2]));
                    Fatom2 += f[typei] * sin(2 * PI * (K[0] * x[3*ii] + 
                                                       K[1] * x[3*ii+1] +
                                                       K[2] * x[3*ii+2]));
                  }
                }
               
                Fvec1[i] = Fatom1;
                Fvec2[i] = Fatom2;              
                m++;
              }
              n++;
            }
          }
        } 
      } 
    }
  }
}  // End of MIC region

#ifdef ENABLE_MIC && __INTEL_OFFLOAD

  // Start of CPU concurrent activity region   
  double frac=0.1;
  for (int i = 0; i < size_array_rows_loc; i++) {     
     if ( i >= (frac * size_array_rows_loc) ) {
        frac += 0.1;
        if (me == 0 && echo) {
           if (screen)
           fprintf(screen," .");
           if (logfile)
           fprintf(screen," .");
        }
     }    
  } // End of CPU concurrent activity region

#pragma offload_wait target(mic:0) wait(&signal_var)

#endif

  double t1 = MPI_Wtime();
#ifdef ENABLE_MIC && __INTEL_OFFLOAD
  printf("Time ellapsed in offload region on MIC =%f\n",t1-t0);
#else
  printf("Time ellapsed in offload region on CPU =%f\n",t1-t0);
#endif

  if (me == 0 && echo) {
     if (screen)
        fprintf(screen," Compute XRD Complete.\n");
     if (logfile)
        fprintf(screen," Compute XRD Complete.\n");
  }

  double *scratch1 = new double[size_array_rows_loc];
  double *scratch2 = new double[size_array_rows_loc];
  
  // Sum intensity for each ang-hkl combination across processors
  MPI_Allreduce(Fvec1,scratch1,size_array_rows_loc,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(Fvec2,scratch2,size_array_rows_loc,MPI_DOUBLE,MPI_SUM,world);


  if (LP == 1) {
    for (int i = 0; i < size_array_rows_loc; i++) {
       array_loc[i*size_array_cols+1] = (scratch1[i] * scratch1[i] + scratch2[i] * scratch2[i]) 
                                         * (1 + cos(2*array_loc[i*size_array_cols+0]) 
                                         * cos(2*array_loc[i*size_array_cols+0])) 
                                         / (cos(array_loc[i*size_array_cols+0]) * sin(array_loc[i*size_array_cols+0]) * sin(array_loc[i*size_array_cols+0])) / natoms;
    }
  } else {
    for (int i = 0; i < size_array_rows_loc; i++) {
       array_loc[i*size_array_cols+1] = (scratch1[i] * scratch1[i] + scratch2[i] * scratch2[i])  / natoms;
    }
  }
  
  
  for (int i = 0; i < size_array_rows_loc; i++) {
     array[i][0] = array_loc[i*size_array_cols+0];
     array[i][1] = array_loc[i*size_array_cols+1];
  }

  delete [] scratch1;
  delete [] scratch2;
  delete [] Fvec1;
  delete [] Fvec2;
}

/* ----------------------------------------------------------------------
 memory usage of arrays
 ------------------------------------------------------------------------- */

double ComputeXRDMIC::memory_usage()
{
  double bytes = size_array_rows_loc * size_array_cols * sizeof(double); //array
  bytes +=  4.0 * size_array_rows_loc * sizeof(double); //Fvec1 & 2, scratch1 & 2
 
  return bytes;
}
