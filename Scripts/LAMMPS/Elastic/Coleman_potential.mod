# NOTE: This script can be modified for different pair styles.  

# It is invoked multiple times, before deformation occurs in the 
#  various directions (after resetting the simulation).

######################## Choose Potential ########################

pair_style      reax/c NULL
pair_coeff      * * ffield.goddard.txt B C B C B
fix             QEQ all qeq/reax 1 0.0 10.0 1.0e-6 reax/c

#pair_style      sw
#pair_coeff      * * Si.sw Si

######################## Neighbor Options ########################

neighbor         1.0 nsq
neigh_modify     once no every 1 delay 0 check yes

####################### Minimization Style #######################

min_style	     ${MinStyle}
min_modify	     dmax ${dmax} line quadratic

######################### Output Options #########################

thermo	         ${NThermo}
thermo_style     custom ${ThermoStyle}
thermo_modify    norm no
