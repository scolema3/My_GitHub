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

#include "offload.h"

#include <iostream>

#ifndef __NUM_CPU_THREADS
#ifndef OMP_NUM_THREADS
#define __NUM_CPU_THREADS 4
#else
#define __NUM_CPU_THREADS OMP_NUM_THREADS
#endif
#endif

#ifndef __NUM_MIC_THREADS
#ifndef MIC_OMP_NUM_THREADS
#define __NUM_MIC_THREADS 60
#else
#define __NUM_MIC_THREADS MIC_OMP_NUM_THREADS
#endif
#endif

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

  double *Fvec1 = new double[size_array_rows_loc]; // Structure factor vector 
  double *Fvec2 = new double[size_array_rows_loc]; // Structure factor vector (imaginary)
                                           // -- Vector entries correspond to 
                                           //    diffraction angles
 
  ntypes = atom->ntypes;
  nlocal = atom->nlocal;
  int *type  = atom->type;
  int natoms = group->count(igroup);
  int *mask = atom->mask;

  double *x = new double [3*nlocal];

  char sig0;
  double *array_mic = array_loc;
  double *ASF_mic = ASF;
  double *store_tmp = new double[4*size_array_rows_loc];

#pragma offload target(mic:0) \
   nocopy(ASF_mic : length(9*ntypes) alloc_if(1) free_if(0)) \
   nocopy(array_mic : length(size_array_rows_loc*size_array_cols) alloc_if(1) free_if(0)) \
   nocopy(mask : length(nlocal) alloc_if(1) free_if(0)) \ 
   nocopy(x : length(3*nlocal) alloc_if(1) free_if(0)) \
   nocopy(type : length(nlocal) alloc_if(1) free_if(0)) \
   nocopy(Fvec1 : length(size_array_rows_loc) alloc_if(1) free_if(0)) \
   nocopy(Fvec2 : length(size_array_rows_loc) alloc_if(1) free_if(0)) \
   nocopy(store_tmp : length(4*size_array_rows_loc) alloc_if(1) free_if(0)) signal(&sig0)
{}

  if (me == 0 && echo) {
      if (screen)
        fprintf(screen,"Computing XRD intensities\n");
      if (logfile)
        fprintf(screen,"Computing XRD intensities\n");
  }
 
  double t0 = MPI_Wtime();

  for (int ii = 0; ii < nlocal; ii++) {
     x[ii+0] = atom->x[ii][0];
     x[ii+1] = atom->x[ii][1];
     x[ii+2] = atom->x[ii][2];
  }

   int starting_row;
   if (my_mpi_rank < size_array_rows_mod) {
      starting_row = my_mpi_rank*size_array_rows_loc;
   }
   else {
      starting_row = my_mpi_rank*size_array_rows_loc+size_array_rows_mod;
   }

// ========================================
// The following code is moved from the openMP/MIC
// region to allow for openMP to run properly
// ========================================
   double K[3];
   double dinv2 = 0.0;
   double ang =0.0;
   int mmax = (2*Knmax[0]+1)*(2*Knmax[1]+1)*(2*Knmax[2]+1);
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
      if  (4 >= dinv2 * lambda * lambda && n < starting_row+size_array_rows_loc) {
         ang = asin(lambda * sqrt(dinv2) / 2);
         if (ang <= Max2Theta & ang >= Min2Theta) {
            store_tmp[4*n] = ang;
            store_tmp[4*n+1] = K[0];
            store_tmp[4*n+2] = K[1];
            store_tmp[4*n+3] = K[2];
            n++;
         }
      }
   }
   int nmax = n;
// ========================================


#pragma offload_wait target(mic:0) wait(&sig0)

char signal_var;
#pragma offload target(mic:0) \
   in(x : length(3*nlocal) alloc_if(0) free_if(1)) \
   in(ASF_mic : length(9*ntypes) alloc_if(0) free_if(1)) \
   in(mask : length(nlocal) alloc_if(0) free_if(1)) \ 
   in(store_tmp : length(4*nmax) alloc_if(0) free_if(1)) \
   in(type : length(nlocal) alloc_if(0) free_if(1)) \
   out(array_mic : length(nmax*size_array_cols) alloc_if(0) free_if(1)) \
   out(Fvec1 : length(nmax) alloc_if(0) free_if(1)) \
   out(Fvec2 : length(nmax) alloc_if(0) free_if(1)) \
   signal(&signal_var)
{  // Start of the MIC region


#ifdef __MIC__
   printf("After pragma using offload from coprocessor, ntypes=%d\n",ntypes);
   omp_set_num_threads(__NUM_MIC_THREADS);
#else
   printf("After pragma using offload from processor, ntypes=%d\n",ntypes);
   omp_set_num_threads(__NUM_CPU_THREADS);
#endif
   int num_threads = omp_get_max_threads();
   printf("No. of OMP Threads: %d\n",num_threads);

   double PI = 4.0*atan(1.0);

#pragma omp parallel num_threads(num_threads) // private(Fatom1, Fatom2, typei)
   {
      double *f = new double[ntypes];    // atomic structure factor by type
      int typei = 0;
      double Fatom1 = 0.0;               // structure factor per atom
      double Fatom2 = 0.0;               // structure factor per atom (imaginary)
      double K[3];
      double ang =0.0;
      double S = 0.0;                    // sin(theta)/lambda -- used to 
                                         // determine atomic structure factor

#pragma omp for
      for (int n=starting_row; n<nmax; n++) {
 
         ang = store_tmp[4*n];
         K[0] = store_tmp[4*n+1];
         K[1] = store_tmp[4*n+2];
         K[2] = store_tmp[4*n+3];

         array_mic[n*size_array_cols+0] = ang;
                
         Fatom1 = 0.0;
         Fatom2 = 0.0;

         // Calculate the atomic structre factor by type	
         S = sin(ang) / lambda;    	
         for (int ii = 0; ii < ntypes; ii++){
            f[ii] = 0;
            int C = ii * 9;
            int D = C+8;
            while (C < D) {
               f[ii] += ASF_mic[C] * exp( -1 * ASF_mic[C+1] * S * S);
               C += 2;
            }
            f[ii] += ASF_mic[D];
         }
  
         // Evaluate the structure factor equation -- looping over all atoms
         for (int ii = 0; ii < nlocal; ii++){
            typei=type[ii]-1;  // for some strange reason, this line needs to be here instead of below the following if-statement, otherwise MIC gives junk
            if (mask[ii] & groupbit) {
               Fatom1 += f[typei] * cos(2 * PI * (K[0] * x[3*ii] + 
                                                  K[1] * x[3*ii+1] +
                                                  K[2] * x[3*ii+2]));
               Fatom2 += f[typei] * sin(2 * PI * (K[0] * x[3*ii] + 
                                                  K[1] * x[3*ii+1] +
                                                  K[2] * x[3*ii+2]));
            }
         }
            
         Fvec1[n] = Fatom1;
         Fvec2[n] = Fatom2;              
      } // End of pragma omp for region
      delete [] f;
   } // End of pragma omp parallel region
}  // End of MIC region

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
  fprintf(screen," .\n");
#pragma offload_wait target(mic:0) wait(&signal_var)

  double t1 = MPI_Wtime();
  printf("Time ellapsed during the offload region =%f\n",t1-t0);

  if (me == 0 && echo) {
     if (screen)
        fprintf(screen," Compute XRD Complete.\n");
     if (logfile)
        fprintf(screen," Compute XRD Complete.\n");
  }

// I think size_array_rows_loc here should be replaced by nmax.... -Yang W
  double *scratch1 = new double[size_array_rows_loc];
  double *scratch2 = new double[size_array_rows_loc];
  
  // Sum intensity for each ang-hkl combination across processors
// I think size_array_rows_loc here should be replaced by nmax.... -Yang W
  MPI_Allreduce(Fvec1,scratch1,size_array_rows_loc,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(Fvec2,scratch2,size_array_rows_loc,MPI_DOUBLE,MPI_SUM,world);


  omp_set_num_threads(__NUM_CPU_THREADS);
  int num_threads = omp_get_max_threads();
#pragma omp parallel num_threads(num_threads)
  { // Begin of pragma omp parallel block

// I think size_array_rows_loc here should be replaced by nmax.... -Yang W
     if (LP == 1) {
#pragma omp for
        for (int i = 0; i < size_array_rows_loc; i++) {
           array_loc[i*size_array_cols+1] = (scratch1[i] * scratch1[i] + scratch2[i] * scratch2[i]) 
                               * (1 + cos(2*array_loc[i*size_array_cols+0]) 
                               * cos(2*array_loc[i*size_array_cols+0])) 
                               / (cos(array_loc[i*size_array_cols+0]) * sin(array_loc[i*size_array_cols+0]) * sin(array_loc[i*size_array_cols+0])) / natoms;
        }
     } else {
#pragma omp for
        for (int i = 0; i < size_array_rows_loc; i++) {
           array_loc[i*size_array_cols+1] = (scratch1[i] * scratch1[i] + scratch2[i] * scratch2[i])  / natoms;
        }
     }
  
#pragma omp for
     for (int i = 0; i < size_array_rows_loc; i++) {
        array[i][0] = array_loc[i*size_array_cols+0];
        array[i][1] = array_loc[i*size_array_cols+1];
     }

  } // End of pragma omp parallel block

  delete [] scratch1;
  delete [] scratch2;
  delete [] Fvec1;
  delete [] Fvec2;
  delete [] x;
  delete [] store_tmp;
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
