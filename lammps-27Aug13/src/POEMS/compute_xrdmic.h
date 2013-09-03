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

#ifdef COMPUTE_CLASS

ComputeStyle(xrdmic,ComputeXRDMIC)

#else

#ifndef LMP_COMPUTE_XRDMIC_H
#define LMP_COMPUTE_XRDMIC_H

#include "stdio.h"
#include "compute.h"

namespace LAMMPS_NS {

class ComputeXRDMIC : public Compute {
 public:
  ComputeXRDMIC(class LAMMPS *, int, char **);
  ~ComputeXRDMIC();
  void init();
  void compute_array();
  double memory_usage();

 private:
  int     me;
  char 	  *ASFfilename;      // Filename containing atomic scattering factor info
  double  Min2Theta;         // Minimum 2theta value (input in 2theta rad)
  double  Max2Theta;         // Maximum 2theta value (input in 2theta rad)
  double  Kmax;              // Maximum reciprocal distance to explore  
  double  c[3];              // Resolution parameters for reciprocal space explored
  int     Knmax[3];          // maximum integer value for K points in each dimension
  double  dK[3];             // Parameters controlling resolution of reciprocal space explored
  double  prd_inv[3];        // Inverse spacing of unit cell
  int     LP;                // Switch to turn on Lorentz-Polarization factor 1=on
  bool    echo;              // echo compute_array progress
	  
  int my_mpi_rank;
  int num_mpi_procs;
  int size_array_rows_mod;

#ifdef ENABLE_MIC && __INTEL_OFFLOAD

#pragma offload_attribute(push, target(mic))
  int size_array_rows_loc;
  int ntypes;
  int nlocal;
  int *type;
  int *mask;
  double *x;
  double *array_loc;

  double *ASF;
  double lambda;
  double *Fvec1;
  double *Fvec2;
  double *f;

  int size_array_cols;
  int groupbit;
#pragma offload_attribute(pop)

#else

  int size_array_rows_loc;
  int ntypes;
  int nlocal;
  int *type;
  int *mask;
  double *x;
  double *array_loc;

  double *ASF;	             // Analytic atomic scattering factor parameters
  double lambda;            // Radiation wavelenght (distance units)
  double *Fvec1;
  double *Fvec2;
  double *f;

#endif
  
};

}

#endif
#endif
