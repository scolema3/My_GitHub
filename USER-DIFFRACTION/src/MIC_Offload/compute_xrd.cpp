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

// ==========================================================
// The following lines are added ............................
//     Purpose: To implement many-thread calculations on MIC
// ----------------------------------------------------------
#ifdef ENABLE_MIC
#include "offload.h"
#endif

#ifndef __NUM_CPU_THREADS
#ifndef OMP_NUM_THREADS
#define __NUM_CPU_THREADS 4
#else
#define __NUM_CPU_THREADS OMP_NUM_THREADS
#endif
#endif

#ifndef __NUM_MIC_THREADS
#ifndef MIC_OMP_NUM_THREADS
#define __NUM_MIC_THREADS 240
#else
#define __NUM_MIC_THREADS MIC_OMP_NUM_THREADS
#endif
#endif
// ==========================================================

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

  // Checking errors specific to the compute
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

    } else if (strcmp(arg[iarg],"ratio") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal Compute XRD Command");
      ratio = atof(arg[iarg+1]);
      if (ratio < 0) 
        error->all(FLERR,"Ratio value must be greater than or equal to zero");
      if (ratio > 1) 
        error->all(FLERR,"Ratio value must be less than or equal to one");
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
        fprintf(screen,"-----\nComputing XRD intensities for %s\n",id);
  }

  double t0 = MPI_Wtime();
  double *Fvec = new double[2*size_array_rows]; // Strct factor  (real & imaginary)
  // -- Note: array rows correspond to different RELP
 
  ntypes = atom->ntypes;
  int nlocal = atom->nlocal;
  int *type  = atom->type;
  int natoms = group->count(igroup);
  int *mask = atom->mask;

  nlocalgroup = 0;
  for (int ii = 0; ii < nlocal; ii++) {
    if (mask[ii] & groupbit) {
     nlocalgroup++;
    }
  }

  double *xlocal = new double [3*nlocalgroup];
  int *typelocal = new int [nlocalgroup];

  int *store_omp = store_tmp;
  int *ztype_omp = ztype;
  
  int size_array_rows_MIC = 0;
  int size_array_rows_CPU = size_array_rows;

  int nthreadMIC = 1;
  int nthreadCPU = 1; 

#ifdef ENABLE_MIC
   char sig0;
   float MIC2CPU_ratio = ratio;
   size_array_rows_MIC = MIC2CPU_ratio*size_array_rows;
   size_array_rows_CPU = size_array_rows-size_array_rows_MIC;

   nthreadMIC = 240;

#pragma offload target(mic:0) \
   nocopy(ztype_omp : length(ntypes) alloc_if(1) free_if(0)) \
   nocopy(xlocal : length(3*nlocalgroup) alloc_if(1) free_if(0)) \
   nocopy(typelocal : length(nlocalgroup) alloc_if(1) free_if(0)) \
   nocopy(store_omp : length(3*size_array_rows_MIC) alloc_if(1) free_if(0)) \
   nocopy(Fvec : length(2*size_array_rows_MIC) alloc_if(1) free_if(0)) signal(&sig0)
{ }
#endif
// ==========================================================

  nlocalgroup = 0;
  for (int ii = 0; ii < nlocal; ii++) {
    if (mask[ii] & groupbit) {
     xlocal[3*nlocalgroup] = atom->x[ii][0];
     xlocal[3*nlocalgroup+1] = atom->x[ii][1];
     xlocal[3*nlocalgroup+2] = atom->x[ii][2];
     typelocal[nlocalgroup]=type[ii];
     nlocalgroup++;
    }
  } 
  

// Setting up OMP
#ifdef _OPENMP
  nthreadCPU = comm->nthreads;
  if (me == 0 && echo) {
    if (screen)
      fprintf(screen," - using %d OMP threads on CPU",nthreadCPU);
  }
#endif
  if (me == 0 && echo) {
    if (screen)
#ifdef ENABLE_MIC
      fprintf(screen," and using %d OMP threads on MIC",nthreadMIC);
#endif
      fprintf(screen,"\n");
  }

// ==========================================================
// The following lines are added ............................
//     Purpose: To offload calculations to MIC
// ----------------------------------------------------------
#ifdef ENABLE_MIC
  if (me == 0 && echo) {
    printf(" - size_array_rows_MIC = %d  (Ratio=%0.2f) \n",size_array_rows_MIC,ratio);
    printf(" - size_array_rows_CPU = %d  (Ratio=%0.2f) \n",size_array_rows_CPU,1-ratio);
  }
#endif
  
#ifdef ENABLE_MIC
#pragma offload_wait target(mic:0) wait(&sig0)
char signal_var;
#pragma offload target(mic:0) \
   in(ztype_omp : length(ntypes) alloc_if(0) free_if(1)) \
   in(xlocal : length(3*nlocalgroup) alloc_if(0) free_if(1)) \
   in(typelocal : length(nlocalgroup) alloc_if(0) free_if(1)) \
   in(store_omp : length(3*size_array_rows_MIC) alloc_if(0) free_if(1)) \
   out(Fvec : length(2*size_array_rows_MIC) alloc_if(0) free_if(1)) \
   signal(&signal_var)
{  // Start of the MIC region
#pragma omp parallel num_threads(nthreadMIC)
  {
    double *f = new double[ntypes];    // atomic structure factor by type
    int typei = 0;
    double Fatom1 = 0.0;               // structure factor per atom
    double Fatom2 = 0.0;               // structure factor per atom (imaginary)

    double K[3];
    double dinv2 = 0.0;
    double dinv  = 0.0;
    double SinTheta_lambda  = 0.0;     // sin(theta)/lambda
    double SinTheta = 0.0;
    double ang = 0.0;
    double Cos2Theta = 0.0;
    double CosTheta = 0.0;

    double inners = 0.0;
    double lp = 0.0;

    if (me == 0 && echo && omp_get_thread_num() == 0) {
       printf("Inside the parallel region on MIC, there are %d openMP thread(s) available for XRD\n",omp_get_max_threads());
     printf(" nthreadsCPU = %d\n",nthreadCPU);
     printf(" nthreadsMIC = %d\n",nthreadMIC);
     printf(" LP = %d\n",LP);
    }

    if (LP == 1) {
#pragma omp for
      for (int n = 0; n < size_array_rows_MIC; n++) {
        int k = store_omp[3*n];
        int j = store_omp[3*n+1];
        int i = store_omp[3*n+2];
        K[0] = i * dK[0];
        K[1] = j * dK[1];
        K[2] = k * dK[2];

        dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
        dinv = sqrt(dinv2);
        SinTheta_lambda = 0.5*dinv;
        SinTheta = SinTheta_lambda * lambda; 
        ang = asin( SinTheta );
        Cos2Theta = cos( 2 * ang);
        CosTheta = cos( ang ); 

        Fatom1 = 0.0;
        Fatom2 = 0.0;

        // Calculate the atomic structre factor by type
        for (int ii = 0; ii < ntypes; ii++){
          f[ii] = 0;
          for (int C = 0; C < 8 ; C+=2){
            f[ii] += ASFXRD[ztype_omp[ii]][C] * exp(-1 * ASFXRD[ztype_omp[ii]][C+1] * SinTheta_lambda * SinTheta_lambda );
          }
          f[ii] += ASFXRD[ztype_omp[ii]][8];
        }

        // Evaluate the structure factor equation -- looping over all atoms
        for (int ii = 0; ii < nlocalgroup; ii++){
            inners = 2 * MY_PI * (K[0] * xlocal[3*ii] + K[1] * xlocal[3*ii+1] +
                      K[2] * xlocal[3*ii+2]);
            typei=typelocal[ii]-1;
            Fatom1 += f[typei] * cos(inners);
            Fatom2 += f[typei] * sin(inners);
        }
        lp = (1 + Cos2Theta * Cos2Theta) /
             ( CosTheta * SinTheta * SinTheta);

        Fvec[2*n] = Fatom1 * lp;
        Fvec[2*n+1] = Fatom2 * lp;
      } // End of pragma omp for region

    } else {
#pragma omp for
      for (int n = 0; n < size_array_rows_MIC; n++) {
        int k = store_omp[3*n];
        int j = store_omp[3*n+1];
        int i = store_omp[3*n+2];
        K[0] = i * dK[0];
        K[1] = j * dK[1];
        K[2] = k * dK[2];

        dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
        dinv = sqrt(dinv2);
        SinTheta_lambda = 0.5*dinv;

        Fatom1 = 0.0;
        Fatom2 = 0.0;

        // Calculate the atomic structre factor by type
        for (int ii = 0; ii < ntypes; ii++){
          f[ii] = 0;
          for (int C = 0; C < 8 ; C+=2){
            f[ii] += ASFXRD[ztype_omp[ii]][C] * exp(-1 * ASFXRD[ztype_omp[ii]][C+1] * SinTheta_lambda * SinTheta_lambda );
          }
          f[ii] += ASFXRD[ztype_omp[ii]][8];
        }

        // Evaluate the structure factor equation -- looping over all atoms
        for (int ii = 0; ii < nlocalgroup; ii++){
            inners = 2 * MY_PI * (K[0] * xlocal[3*ii] + K[1] * xlocal[3*ii+1] +
                      K[2] * xlocal[3*ii+2]);
            typei=typelocal[ii]-1;
            Fatom1 += f[typei] * cos(inners);
            Fatom2 += f[typei] * sin(inners);
        }
        Fvec[2*n] = Fatom1;
        Fvec[2*n+1] = Fatom2;
      } // End of pragma omp for region
    } // End of if LP=1 check 
    delete [] f;
  } // End of pragma omp parallel region 
// ==========================================================
// The following lines are added ............................
//    Purpose: A bracket is needed to close the MIC region,
//             and a wait signal is needed to make sure that
//             MIC is done the calculation.
// ----------------------------------------------------------
}  // End of MIC region
#endif



  // Start of CPU concurrent activity region
  double tCPU0 = MPI_Wtime();
  int m = 0;
  double frac = 0.1;
#pragma omp parallel num_threads(nthreadCPU)
  {
    double *f = new double[ntypes];    // atomic structure factor by type
    int typei = 0;
    double Fatom1 = 0.0;               // structure factor per atom
    double Fatom2 = 0.0;               // structure factor per atom (imaginary)
    double S = 0.0;                    // sin(theta)/lambda
    double K[3];
    double dinv2 = 0.0;
    double dinv  = 0.0;
    double SinTheta_lambda  = 0.0;     // sin(theta)/lambda
    double SinTheta = 0.0;
    double ang = 0.0;
    double Cos2Theta = 0.0;
    double CosTheta = 0.0;

    double *inners = new double[nlocalgroup];
    double lp = 0.0;

    if (me == 0 && echo && omp_get_thread_num() == 0) {
       printf("Inside the CPU parallel region, there are %d openMP thread(s) available\n",omp_get_max_threads());
    }
    if (LP == 1) {
#pragma omp for
      for (int n = size_array_rows_MIC; n < size_array_rows; n++) {
        int k = store_omp[3*n];
        int j = store_omp[3*n+1];
        int i = store_omp[3*n+2];
        K[0] = i * dK[0];
        K[1] = j * dK[1];
        K[2] = k * dK[2];

        dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
        dinv = sqrt(dinv2);
        SinTheta_lambda = 0.5*dinv;
        SinTheta = SinTheta_lambda * lambda; 
        ang = asin( SinTheta );
        Cos2Theta = cos( 2 * ang);
        CosTheta = cos( ang ); 

        Fatom1 = 0.0;
        Fatom2 = 0.0;

        // Calculate the atomic structre factor by type
        for (int ii = 0; ii < ntypes; ii++){
          f[ii] = 0;
          for (int C = 0; C < 8 ; C+=2){
            f[ii] += ASFXRD[ztype_omp[ii]][C] * exp(-1 * ASFXRD[ztype_omp[ii]][C+1] * SinTheta_lambda * SinTheta_lambda );
          }
          f[ii] += ASFXRD[ztype_omp[ii]][8];
        }


        // Evaluate the structure factor equation -- looping over all atoms
        for (int ii = 0; ii < nlocalgroup; ii++){
            inners[ii] = 2 * MY_PI * (K[0] * xlocal[3*ii] + K[1] * xlocal[3*ii+1] +
                      K[2] * xlocal[3*ii+2]);
        }
        for (int ii = 0; ii < nlocalgroup; ii++){
            typei=typelocal[ii]-1;
            Fatom1 += f[typei] * cos(inners[ii]);
            Fatom2 += f[typei] * sin(inners[ii]);
        }
        lp = (1 + Cos2Theta * Cos2Theta) /
             ( CosTheta * SinTheta * SinTheta);
        Fvec[2*n] = Fatom1 * lp;
        Fvec[2*n+1] = Fatom2 * lp;

        // reporting progress of calculation
        if ( echo ) {
          #pragma omp critical
          {
            if ( m == round(frac * size_array_rows_CPU) ) {
              if (me == 0 && screen) fprintf(screen," %0.0f%% -",frac*100);
              frac += 0.1;
            }
            m++;
          }
        }      
      } // End of pragma omp for region

    } else {
#pragma omp for
      for (int n = size_array_rows_MIC; n < size_array_rows; n++) {
        int k = store_omp[3*n];
        int j = store_omp[3*n+1];
        int i = store_omp[3*n+2];
        K[0] = i * dK[0];
        K[1] = j * dK[1];
        K[2] = k * dK[2];

        dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
        dinv = sqrt(dinv2);
        SinTheta_lambda = 0.5*dinv;

        Fatom1 = 0.0;
        Fatom2 = 0.0;


        // Calculate the atomic structre factor by type
        for (int ii = 0; ii < ntypes; ii++){
          f[ii] = 0;
          for (int C = 0; C < 8 ; C+=2){
            f[ii] += ASFXRD[ztype_omp[ii]][C] * exp(-1 * ASFXRD[ztype_omp[ii]][C+1] * SinTheta_lambda * SinTheta_lambda );
          }
          f[ii] += ASFXRD[ztype_omp[ii]][8];
        }

        // Evaluate the structure factor equation -- looping over all atoms
        for (int ii = 0; ii < nlocalgroup; ii++){
            inners[ii] = 2 * MY_PI * (K[0] * xlocal[3*ii] + K[1] * xlocal[3*ii+1] +
                      K[2] * xlocal[3*ii+2]);
        }
        for (int ii = 0; ii < nlocalgroup; ii++){
            typei=typelocal[ii]-1;
            Fatom1 += f[typei] * cos(inners[ii]);
            Fatom2 += f[typei] * sin(inners[ii]);
        }
        Fvec[2*n] = Fatom1;
        Fvec[2*n+1] = Fatom2;


        // reporting progress of calculation
        if ( echo ) {
          #pragma omp critical
          {
            if ( m == round(frac * size_array_rows_CPU) ) {
              if (me == 0 && screen) fprintf(screen," %0.0f%% -",frac*100 );
              frac += 0.1;
            }
            m++;
          }
        }
      } // End of pragma omp for region
    } // End of if LP=1 check 
    delete [] f;
    delete [] inners;
  } // End of pragma omp parallel region on CPU
  double tCPU1 = MPI_Wtime();

// ==========================================================
#ifdef ENABLE_MIC
#pragma offload_wait target(mic:0) wait(&signal_var)
#endif

// ==========================================================
// The following lines are modified..........................
//    Purpose: MPI_Allreduce function is called only once.
// ----------------------------------------------------------
  double *scratch = new double[2*size_array_rows];

  // Sum intensity for each ang-hkl combination across processors
  MPI_Allreduce(Fvec,scratch,2*size_array_rows,MPI_DOUBLE,MPI_SUM,world);

#pragma omp parallel for
  for (int i = 0; i < size_array_rows; i++) {
    array[i][1] = (scratch[2*i] * scratch[2*i] + scratch[2*i+1] * scratch[2*i+1]) / natoms;
  }

  delete [] scratch;
  delete [] Fvec;
  delete [] xlocal;
  delete [] typelocal;
// ==========================================================

  double t2 = MPI_Wtime();

  // compute memory usage per processor
  double bytes = size_array_rows * size_array_cols * sizeof(double); //array
  bytes +=  4.0 * size_array_rows * sizeof(double); //Fvec1 & 2, scratch1 & 2
  bytes += ntypes * sizeof(double); // f
  bytes += 3.0 * nlocalgroup * sizeof(double); // xlocal
  bytes += nlocalgroup * sizeof(int); // x
  bytes += 3.0 * size_array_rows * sizeof(int); // store_temp
  

  if (me == 0 && echo) {
    if (screen)
      fprintf(screen," 100%%\nTime ellapsed during compute_xrd = %0.2f sec using %0.2f Mbytes/processor", t2-t0, bytes/1024.0/1024.0);
      fprintf(screen," \n -time ellapsed within CPU loop = %0.2f sec", tCPU1-tCPU0);
      fprintf(screen," \n -time waiting for MIC to finish = %0.2f sec\n-----\n", (t2-t0)-(tCPU1-tCPU0));
  }
}

/* ----------------------------------------------------------------------
 memory usage of arrays
 ------------------------------------------------------------------------- */

double ComputeXRD::memory_usage()
{
  double bytes = size_array_rows * size_array_cols * sizeof(double); //array
  bytes +=  4.0 * size_array_rows * sizeof(double); //Fvec1 & 2, scratch1 & 2
  bytes += 3.0 * nlocalgroup * sizeof(double); // x
  bytes += ntypes * sizeof(double); // f
  bytes += 3.0 * size_array_rows * sizeof(int); // store_temp
  
  return bytes;
}

