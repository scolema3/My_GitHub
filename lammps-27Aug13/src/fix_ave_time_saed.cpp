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
   Incorporating SAED: Shawn Coleman (Arkansas)
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "fix_ave_time_saed.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
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
enum{SCALAR,VECTOR};

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4

/* ---------------------------------------------------------------------- */

FixAveTimeSAED::FixAveTimeSAED(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix ave/time/saed command");

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

  options(narg,arg);

  // parse values until one isn't recognized
  // if mode = VECTOR and value is a global array:
  //   expand it as if columns listed one by one
  //   adjust nvalues accordingly via maxvalues

  which = argindex = value2index = offcol = NULL;
  ids = NULL;
  int maxvalues = nvalues;
  allocate_values(maxvalues);
  nvalues = 0;

  iarg = 6;
  while (iarg < narg) {
    if (strncmp(arg[iarg],"c_",2) == 0 ) {
      if (arg[iarg][0] == 'c') which[nvalues] = COMPUTE;
      else break;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Illegal fix ave/time command");
        argindex[nvalues] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[nvalues] = 0;

      n = strlen(suffix) + 1;
      ids[nvalues] = new char[n];
      strcpy(ids[nvalues],suffix);
      delete [] suffix;

      if (mode == VECTOR && which[nvalues] == COMPUTE &&
          argindex[nvalues] == 0) {
        int icompute = modify->find_compute(ids[nvalues]);
        if (icompute < 0)
          error->all(FLERR,"Compute ID for fix ave/time does not exist");
        if (modify->compute[icompute]->array_flag) {
          int ncols = modify->compute[icompute]->size_array_cols;
          maxvalues += ncols-1;
          allocate_values(maxvalues);
          argindex[nvalues] = 1;
          for (int icol = 1; icol < ncols; icol++) {
            which[nvalues+icol] = which[nvalues];
            argindex[nvalues+icol] = icol+1;
            n = strlen(ids[nvalues]) + 1;
            ids[nvalues+icol] = new char[n];
            strcpy(ids[nvalues+icol],ids[nvalues]);
          }
          nvalues += ncols-1;
        }
      }
      nvalues++;
      iarg++;
    } else break;

  }

  // set off columns now that nvalues is finalized

  for (int i = 0; i < nvalues; i++) offcol[i] = 0;
  for (int i = 0; i < noff; i++) {
    if (offlist[i] < 1 || offlist[i] > nvalues)
      error->all(FLERR,"Invalid fix ave/time off column");
    offcol[offlist[i]-1] = 1;
  }

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable

  if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix ave/time command");
  if (nfreq % nevery || (nrepeat-1)*nevery >= nfreq)
    error->all(FLERR,"Illegal fix ave/time command");
  if (ave != RUNNING && overwrite)
    error->all(FLERR,"Illegal fix ave/time command");

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE && mode == VECTOR) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/time does not exist");
      if (argindex[i] == 0 && modify->compute[icompute]->vector_flag == 0)
        error->all(FLERR,"Fix ave/time compute does not calculate a vector");
      if (argindex[i] && modify->compute[icompute]->array_flag == 0)
        error->all(FLERR,"Fix ave/time compute does not calculate an array");
      if (argindex[i] &&
          argindex[i] > modify->compute[icompute]->size_array_cols)
        error->all(FLERR,"Fix ave/time compute array is accessed out-of-range");
    } 
  }

  // if VECTOR mode, check that all columns are same length
  // nrows = # of rows in output array

  if (mode == VECTOR) {
    int length;
    for (int i = 0; i < nvalues; i++) {
      if (which[i] == COMPUTE) {
        int icompute = modify->find_compute(ids[i]);
        if (argindex[i] == 0) length = modify->compute[icompute]->size_vector;
        else length = modify->compute[icompute]->size_array_rows;
      } 
      if (i == 0) nrows = length;
      else if (length != nrows)
        error->all(FLERR,"Fix ave/time columns are inconsistent lengths");
    }

    column = new double[nrows];
  } else column = NULL;


 // SAED specific paramaters needed

  int *periodicity = domain->periodicity;
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

int ncount=0;
  // Find integer dimensions of the reciprocal lattice box bounds
  if ( (Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0) ){
    for (int i=0; i<3; i++) {
      dK[i] = prd_inv[i]*c[i];
      Knmax[i] = ceil(Kmax / dK[i]);
      Knmin[i] = -Knmax[i];
    } 
  } else {
  
  ncount++;
  
    for (int i=0; i<3; i++) {
      Knmax[i] = -10000;
      Knmin[i] =  10000;
    }       
    double dinv2 = 0.0;
    double r = 0.0;
    double K[3];
    int Ksearch[3];
    for (int i=0; i<3; i++) {
      dK[i] = prd_inv[i]*c[i];
      Ksearch[i] = ceil(Kmax / dK[i]);
    } 
cout << Ksearch[0] << " , " << Ksearch[1] << " , "   << Ksearch[2] << " , "   << endl;
    
    for (int k = -Ksearch[2]; k <= Ksearch[2]; k++) {
      for (int j = -Ksearch[1]; j <= Ksearch[1]; j++) {
        for (int i = -Ksearch[0]; i <= Ksearch[0]; i++) {
          K[0] = i * dK[0];
          K[1] = j * dK[1];
          K[2] = k * dK[2];
          dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
          if (dinv2 < Kmax * Kmax) {

            r=0.0;
            for (int m=0; m<3; m++) r += pow(K[m] - Zone[m],2.0);
            r = sqrt(r);     
cout << r << endl;                   
            if  ( (r >  (R_Ewald - dR_Ewald) ) && (r < (R_Ewald + dR_Ewald) ) ){
            ncount++;
              for (int ii=0; ii<3; ii++) {
                if ( Ksearch[ii] > Knmax[ii] ) Knmax[ii] = Ksearch[ii];
                if ( Ksearch[ii] < Knmin[ii] ) Knmin[ii] = Ksearch[ii];
                ncount++;
              }
            }
          }
        } 
      }
    } 
  }

cout << " Ncount " << ncount << endl;

  // Finding dimensions for vtk files
  int Dim[3];
  for (int i=0; i<3; i++) {
    if ( ( (Knmin[i] > 0) && (Knmax[i] > 0) ) || ( (Knmin[i] < 0) && (Knmax[i] < 0) ) ){
      Dim[i] = abs( (int) Knmin[i] ) + abs( (int) Knmax[i] );
    } else Dim[i] = abs( (int) Knmin[i] ) + abs( (int) Knmax[i] ) + 1;
  }

  if (fp && me == 0) {
    if (mode == VECTOR) {
      for (int i = 0; i < nvalues; i++) {
        if (which[i] == COMPUTE) {
          fprintf(fp,"# vtk DataFile Version 3.0 c_%s \n",ids[i]);
          fprintf(fp,"Image data set \n");
          fprintf(fp,"ASCII \n");
          fprintf(fp,"DATASET STRUCTURED_POINTS \n");
          fprintf(fp,"DIMENSIONS %d %d %d \n", Dim[0],  Dim[1], Dim[2]);
          fprintf(fp,"ASPECT_RATIO %g %g %g \n", dK[0], dK[1], dK[2]);
          fprintf(fp,"ORIGIN %g %g %g \n", Knmin[0] * dK[0],  Knmin[1] * dK[1], Knmin[2] * dK[2]);
          fprintf(fp,"POINT_DATA %d \n ",  Dim[0] *  Dim[1] * Dim[2] );
          fprintf(fp,"SCALARS intensity float \n");
          fprintf(fp,"LOOKUP_TABLE default \n");
        } else break;
        if (argindex[i]) fprintf(fp,"[%d]",argindex[i]);
      }
    }
    filepos = ftell(fp);
  }

  // allocate memory for averaging

  vector = vector_total = NULL;
  vector_list = NULL;
  array = array_total = NULL;
  array_list = NULL;

  memory->create(array,nrows,nvalues,"ave/time:array");
  memory->create(array_total,nrows,nvalues,"ave/time:array_total");

  // this fix produces either a global scalar or vector or array
  // SCALAR mode produces either a scalar or vector
  // VECTOR mode produces either a vector or array
  // intensive/extensive flags set by compute,fix,variable that produces value

  extlist = NULL;

  if (mode == SCALAR) {
  } else {
    if (nvalues == 1) {
      vector_flag = 1;
      size_vector = nrows;
      if (which[0] == COMPUTE) {
        Compute *compute = modify->compute[modify->find_compute(ids[0])];
        if (argindex[0] == 0) {
          extvector = compute->extvector;
          if (extvector == -1) {
            extlist = new int[nrows];
            for (int i = 0; i < nrows; i++) extlist[i] = compute->extlist[i];
          }
        } else extvector = compute->extarray;
      }

    } else {
      array_flag = 1;
      size_array_rows = nrows;
      size_array_cols = nvalues;
      int value;
      for (int i = 0; i < nvalues; i++) {
        if (which[i] == COMPUTE) {
          Compute *compute = modify->compute[modify->find_compute(ids[i])];
          if (argindex[i] == 0) value = compute->extvector;
          else value = compute->extarray;
        } 
        if (value == -1)
          error->all(FLERR,"Fix ave/time cannot set output array "
                     "intensive/extensive from these inputs");
        if (i == 0) extarray = value;
        else if (value != extarray)
          error->all(FLERR,"Fix ave/time cannot set output array "
                     "intensive/extensive from these inputs");
      }
    }
  }

  // initializations
  // set vector_total/array_total to zero since it accumulates

  irepeat = 0;
  iwindow = window_limit = 0;
  norm = 0;

  for (int i = 0; i < nrows; i++)
    for (int j = 0; j < nvalues; j++) array_total[i][j] = 0.0;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);

}

/* ---------------------------------------------------------------------- */

FixAveTimeSAED::~FixAveTimeSAED()
{
  memory->destroy(which);
  memory->destroy(argindex);
  memory->destroy(value2index);
  memory->destroy(offcol);
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  memory->sfree(ids);

  delete [] extlist;

  if (fp && me == 0) fclose(fp);

  delete [] vector;
  delete [] vector_total;
  delete [] column;
  memory->destroy(array);
  memory->destroy(array_total);
  memory->destroy(array_list);
}

/* ---------------------------------------------------------------------- */

int FixAveTimeSAED::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveTimeSAED::init()
{
  // set current indices for all computes,fixes,variables

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/time does not exist");
      value2index[i] = icompute;

    } 
  }

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

void FixAveTimeSAED::setup(int vflag)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveTimeSAED::end_of_step()
{
  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  invoke_vector(ntimestep);
}

/* ---------------------------------------------------------------------- */


void FixAveTimeSAED::invoke_vector(bigint ntimestep)
{
  int i,j,m;


  // zero if first step

  if (irepeat == 0)
    for (i = 0; i < nrows; i++)
      for (j = 0; j < nvalues; j++) array[i][j] = 0.0;

  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  for (j = 0; j < nvalues; j++) {
    m = value2index[j];

    // invoke compute if not previously invoked

    if (which[j] == COMPUTE) {
      Compute *compute = modify->compute[m];

      if (argindex[j] == 0) {
        if (!(compute->invoked_flag & INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= INVOKED_VECTOR;
        }
        double *cvector = compute->vector;
        for (i = 0; i < nrows; i++)
          column[i] = cvector[i];

      } else {
        if (!(compute->invoked_flag & INVOKED_ARRAY)) {
          compute->compute_array();
          compute->invoked_flag |= INVOKED_ARRAY;
        }
        double **carray = compute->array;
        int icol = argindex[j]-1;
        for (i = 0; i < nrows; i++)
          column[i] = carray[i][icol];
      }
    }

    // add columns of values to array or just set directly if offcol is set

    if (offcol[j]) {
      for (i = 0; i < nrows; i++)
        array[i][j] = column[i];
    } else {
      for (i = 0; i < nrows; i++)
        array[i][j] += column[i];
    }
  }

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
  for (i = 0; i < nrows; i++)
    for (j = 0; j < nvalues; j++)
      if (offcol[j] == 0) array[i][j] /= repeat;

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, combine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    for (i = 0; i < nrows; i++)
      for (j = 0; j < nvalues; j++) array_total[i][j] = array[i][j];
    norm = 1;

  } else if (ave == RUNNING) {
    for (i = 0; i < nrows; i++)
      for (j = 0; j < nvalues; j++) array_total[i][j] += array[i][j];
    norm++;

  } else if (ave == WINDOW) {
    for (i = 0; i < nrows; i++)
      for (j = 0; j < nvalues; j++) {
        array_total[i][j] += array[i][j];
        if (window_limit) array_total[i][j] -= array_list[iwindow][i][j];
        array_list[iwindow][i][j] = array[i][j];
      }

    iwindow++;
    if (iwindow == nwindow) {
      iwindow = 0;
      window_limit = 1;
    }
    if (window_limit) norm = nwindow;
    else norm = iwindow;
  }

  // insure any columns with offcol set are effectively set to last value

  for (i = 0; i < nrows; i++)
    for (j = 0; j < nvalues; j++)
      if (offcol[j]) array_total[i][j] = norm*array[i][j];

  // output result to file


  if (fp && me == 0) {
    if (overwrite) fseek(fp,filepos,SEEK_SET);

      // Finding the intersection of the reciprical space and Ewald sphere
      int NROW1 = 0;
      int NROW2 = 0;
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
              if (dinv2 < Kmax * Kmax) {
                 fprintf(fp,"%g\n",array_total[NROW1][0]/norm);
                 fflush(fp);
                 NROW1++;
                 NROW2++;
              } else {
              fprintf(fp,"%d\n",-1);
              fflush(fp);
              NROW2++;
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
                if  ( (r >  (R_Ewald - dR_Ewald) ) && (r < (R_Ewald + dR_Ewald) ) ){
                 fprintf(fp,"%g\n",array_total[NROW1][0]/norm);
                 fflush(fp);
                 NROW2++;
                 NROW1++;
                } else {
                  fprintf(fp,"%d\n",-1);
                  fflush(fp);
                  NROW2++;
                }
              } else {
              fprintf(fp,"%d\n",-1);
              fflush(fp);
              NROW2++;
             }
            }
          }
        }
      }
cout << "Nrows = " << NROW1 << " " << NROW2 << endl;
    }
}

/* ----------------------------------------------------------------------
   return scalar value
------------------------------------------------------------------------- */

double FixAveTimeSAED::compute_scalar()
{
  if (norm) return vector_total[0]/norm;
  return 0.0;
}

/* ----------------------------------------------------------------------
   return Ith vector value
------------------------------------------------------------------------- */

double FixAveTimeSAED::compute_vector(int i)
{
  if (norm) {
    if (mode == SCALAR) return vector_total[i]/norm;
    if (mode == VECTOR) return array_total[i][0];
  }
  return 0.0;
}

/* ----------------------------------------------------------------------
   return I,J array value
------------------------------------------------------------------------- */

double FixAveTimeSAED::compute_array(int i, int j)
{
  if (norm) return array_total[i][j]/norm;
  return 0.0;
}

/* ----------------------------------------------------------------------
   parse optional args
------------------------------------------------------------------------- */

void FixAveTimeSAED::options(int narg, char **arg)
{
  // option defaults

  fp = NULL;
  ave = ONE;
  startstep = 0;
  mode = SCALAR;
  noff = 0;
  offlist = NULL;
  overwrite = 0;
  title1 = NULL;
  title2 = NULL;
  title3 = NULL;

  // optional args

  int iarg = 6 + nvalues;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      if (me == 0) {
        fp = fopen(arg[iarg+1],"w");
        if (fp == NULL) {
          char str[128];
          sprintf(str,"Cannot open fix ave/time file %s",arg[iarg+1]);
          error->one(FLERR,str);
        }
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else if (strcmp(arg[iarg+1],"window") == 0) ave = WINDOW;
      else error->all(FLERR,"Illegal fix ave/time command");
      if (ave == WINDOW) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix ave/time command");
        nwindow = force->inumeric(FLERR,arg[iarg+2]);
        if (nwindow <= 0) error->all(FLERR,"Illegal fix ave/time command");
      }
      iarg += 2;
      if (ave == WINDOW) iarg++;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      startstep = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"mode") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      if (strcmp(arg[iarg+1],"scalar") == 0) mode = SCALAR;
      else if (strcmp(arg[iarg+1],"vector") == 0) mode = VECTOR;
      else error->all(FLERR,"Illegal fix ave/time command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"off") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      memory->grow(offlist,noff+1,"ave/time:offlist");
      offlist[noff++] = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"overwrite") == 0) {
      overwrite = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"title1") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      delete [] title1;
      int n = strlen(arg[iarg+1]) + 1;
      title1 = new char[n];
      strcpy(title1,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title2") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      delete [] title2;
      int n = strlen(arg[iarg+1]) + 1;
      title2 = new char[n];
      strcpy(title2,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title3") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      delete [] title3;
      int n = strlen(arg[iarg+1]) + 1;
      title3 = new char[n];
      strcpy(title3,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"Kmax") == 0) {
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
    } else if (strcmp(arg[iarg],"lambda") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal Compute SAED Command");
      lambda = atof(arg[iarg+1]);
      if (lambda <= 0 )
        error->all(FLERR,"Compute SAED: lambda must be greater than 0 ");
      iarg += 2;
    } else if (strcmp(arg[iarg],"dR_Ewald") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal Compute SAED Command");
      dR_Ewald = atof(arg[iarg+1]);
      if (dR_Ewald < 0)
        error->all(FLERR,"Compute SAED: dR_Ewald slice must be greater than 0");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix ave/time command");
  }
}

/* ----------------------------------------------------------------------
   reallocate vectors for each input value, of length N
------------------------------------------------------------------------- */

void FixAveTimeSAED::allocate_values(int n)
{
  memory->grow(which,n,"ave/time:which");
  memory->grow(argindex,n,"ave/time:argindex");
  memory->grow(value2index,n,"ave/time:value2index");
  memory->grow(offcol,n,"ave/time:offcol");
  ids = (char **) memory->srealloc(ids,n*sizeof(char *),"ave/time:ids");
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
   startstep is lower bound on nfreq multiple
------------------------------------------------------------------------- */

bigint FixAveTimeSAED::nextvalid()
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

void FixAveTimeSAED::reset_timestep(bigint ntimestep)
{
  if (ntimestep > nvalid) error->all(FLERR,"Fix ave/time missed timestep");
}
