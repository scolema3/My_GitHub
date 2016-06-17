#!/bin/sh 
# Script to run GB translation/overlap optimization using data  
# files created by MATLAB script

Version=0.037

help_ani()
{
  echo "
        Usage: $0 Parameter_File 
        Version: $Version
        
        The Parmater_File must include:
        
        Potential=\"
        
        ( Multiple line input setting up LAMMPS potential -
          bounded by \" \"'s - example commands include: 
          pair_style, pair_coeff, neighbor, neigh_modify 
              
          Note: \\\${Species} writes species list by type
        )

        \"

        ########### Optional Values Listed Below ############

        ( Note: \$PWD variable lists current directory path ) 
        
        
        # Data Location & Names
        Data_Dir=\$PWD/Data
        Work_Dir=\$PWD
        
        # GB File Names
        Find_Prefix=
        Find_Suffix=.data        

        # Execute $MPIrun nMPI LAMMPS < INPUT.in > OUTPUT.out
        nMPI=12
        LAMMPS=lmp_kim
        MPIrun=\"mpirun -np\"

        # Minimization
        etol=1e-20
        ftol=1e-20
        maxiter=10000
        maxeval=10000

        # Optimization
        Nx=1            (Number of x-translations)
        Nz=1            (Number of z-translations)
        Noverlap=1      (Number of overlap distances)
        DeleteRegion=1  (1-upper 2-lower region for deletion)

        # LAMMPS Style
        Units=metal
        AtomStyle=atomic

        # Bulk Files
        Bulk_Prefix=Bulk
        Bulk_Suffix=.data
        Bulk_Name=bulk
        Bulk_Input=\${Bulk_Name}.in
        Bulk_Output=\${Bulk_Name}.out
        "
  exit 1
}

# KIM Example
# Potential="
#  pair_style  kim KIMvirial EAM_Dynamo_Mishin_Farkas_Al__MO_651801486679_001
# pair_coeff  * * \${Species}
# neighbor    5.0 bin
# neigh_modify    delay 0 every 1 check yes"

# ReaxFF/c Example
# Potential="
# pair_style reax/c NULL
# pair_coeff * * ffield.reax \${Species}
# fix        0 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c"

# EAM Example
# Potential="
# pair_style eam/fs 
# pair_coeff * * /home/scoleman/Research/GB_Generation_Code/Matlab/Fe3CFe/CaseStudy/Fe-C_Hepburn_Ackland.eam.fs \${Species} 
# "

## Default Parameter Values

PWD=`pwd`
PWD0=$PWD

# Data
Data_Dir=$PWD/Data
Work_Dir=$PWD
Bulk_Prefix=Bulk
Bulk_Suffix=.Data

FullPeriod=0 # Delete_soon

# Execute $MPIrun nMPI LAMMPS < INPUT.in > OUTPUT.out
nMPI=12
LAMMPS=lmp_kim
MPIrun="mpirun -np"

# GB file names
Find_Prefix=
Find_Suffix=.data

# Minimization
etol=1e-10
ftol=1e-10
maxiter=10000
maxeval=10000

# Optimization
Nx=1
Nz=1
Noverlap=1
DeleteRegion=1

# Style
Units=metal
AtomStyle=atomic

# Bulk Files
Bulk_Dir=bulk
Bulk_Name=bulk
Bulk_Input=${Bulk_Name}.in
Bulk_Output=${Bulk_Name}.out


#################################################################
#################################################################
################### Load Parameters & Setup #####################
#################################################################
#################################################################

echo " Runing $0  "
echo " Version: $Version" 
if [ $# -ne 1 ]; then
  help_ani
else 
  Parameter_File="$1"
fi

# Checking for files
if [ -f $PWD0/$Parameter_File ]; then
  echo " "
  echo "Using parameters from: $Parameter_File"
  source $PWD0/$Parameter_File
else
   echo "File $Parameter_File does not exist."
   exit 1
fi

NGBs=`ls -l ${Data_Dir}/${Find_Prefix}*${Find_Suffix} | wc -l`
if [ $NGBs -eq "0" ]; then 
    echo " No GB files found in ${Data_Dir} ";  exit 1
else
    echo " Found $NGBs grain boundary structure files "
fi

NBulk=`ls -l ${Data_Dir}/${Bulk_Prefix}*${Bulk_Suffix} | wc -l`
if [ $NBulk -eq "0" ]; then 
    echo " No Bulk files found in ${Data_Dir} ";  exit 1
else
    echo " Found $NBulk bulk structures files "
fi

# Adjusting for atom_style and units
if [ "${AtomStyle}" == "charge" ];then
  DumpSpecial="q"
elif [ "${AtomStyle}" == "atomic" ];then
  DumpSpecial=""
else
  echo "AtomStyle \"$AtomStyle\" not recognized. Choose charge or atomic."
  exit
fi

if [ "${Units}" == "real" ];then
  eConv=6.9477E-21*1000 # mJ
elif [ "${Units}" == "metal" ];then
  eConv=1.60217733E-19*1000 # mJ
else
  echo "Units \"$Units\" not recognized. Choose metal or charge."
  exit
fi


mkdir -p $Work_Dir
mkdir -p $Work_Dir/$Bulk_Dir/
mkdir -p $Work_Dir/min_gb

#################################################################
#################################################################
#################################################################
####################### Bulk Simulation #########################
#################################################################
#################################################################
#################################################################
lammps_bulk()
{

count=0
for i in ${Data_Dir}/${Bulk_Prefix}*${Bulk_Suffix}
do

FullPath=`echo ${i} | sed 's./.\\/.g'` # full path 
Index=${FullPath%%$Bulk_Suffix}
Index=${Index#*$Bulk_Prefix}     

# Test if previously computed bulk properties and ask to reuse
NBulkDat=`ls -l ${Work_Dir}/$Bulk_Dir/${Bulk_Name}${Index}.dat 2>/dev/null | wc -l`
if [ $NBulkDat -eq "1" ]; then 
  echo " "
  echo " Bulk minimization already completed. "
  echo -n " Would you like to reuse this data? Type (yes/no) and press [ENTER]: "
  read yno
    
  case $yno in
    [yY] | [yY][Ee][Ss] )
      echo " Reusing computed values..."; CompBulk=0 ;;
    [nN] | [n|N][O|o] )
      echo " Computing new values..."; CompBulk=1 ;;
    * )    
      echo " Invalid Answer"; exit ;;
  esac   
  
  else
    CompBulk=1
fi

if [ $CompBulk -eq "1" ]; then

# Obtaining values from datafile
nTypes=`grep "atom types" $FullPath`
nTypes=${nTypes%atom*}

Types=`grep " Types:" $FullPath`
Types=${Types#*:}
IFS=', ' read -a Types <<< $Types

Species=`grep Species: $FullPath`
Species=${Species#*:}
IFS=', ' read -a Species <<< $Species

Masses=`grep Masses: $FullPath`
Masses=${Masses#*:}
IFS=', ' read -a Masses <<< $Masses


# Writing LAMMPS Bulk Input File
cat > $Work_Dir/$Bulk_Dir/$Bulk_Input${Index} << _EOF_
log ${Bulk_Name}${Index}.log

# Comptute bulk properties (a0 and cohesive energy) for grain 
# boundary energy calculations

#################################################################
#    LAMMPS Input File for Grain Boundary Optimization
#    Shawn P. Coleman, July 2015
#################################################################

#################################################################
# --------------- Setup Simulation Variables ------------------ #
#################################################################

variable     Species      string   "${Species[*]}"

#################################################################
# ------------------- Initialize simulation ------------------- #
#################################################################
clear

units         ${Units} 
dimension     3 
boundary      p p p 
atom_style    ${AtomStyle} 

read_data     ${FullPath}
atom_modify   sort 0 0.0

${Potential}

#################################################################
# ------------------------ Mimization ------------------------- #
#################################################################

## Box Relax
fix        1 all box/relax aniso 0.0 vmax 0.001
minimize   ${etol} ${ftol} ${maxiter} ${maxeval} 
unfix      1

variable   natoms equal "count(all)" 
variable   teng   equal pe
variable   length equal "vol^(1/3)"
variable   ecoh   equal "v_teng/v_natoms"

print      "Total energy (eV) = \${teng}" file ${Bulk_Name}${Index}.dat
print      "Number of atoms = \${natoms}" append ${Bulk_Name}${Index}.dat
print      "Lattice constant (Angstoms) = \${length}" append ${Bulk_Name}${Index}.dat
print      "Cohesive energy (eV) = \${ecoh}" append ${Bulk_Name}${Index}.dat

dump       1 all custom 1 ${Bulk_Name}${Index}.dump id type x y z ${DumpSpecial}
run        0
_EOF_

# Run bulk simulation
cd ${Work_Dir}/$Bulk_Dir
$MPIrun $nMPI $LAMMPS < $Bulk_Input${Index}  >& $Bulk_Output
cd ${PWD0}
fi

# Read data from bulk minimization 
cd ${Work_Dir}/$Bulk_Dir
a0=`grep "Lattice constant" ${Bulk_Name}${Index}.dat`
a0=${a0#* = }
IFS=', ' read -a a0 <<< $a0

E=`grep "Cohesive energy" ${Bulk_Name}${Index}.dat`
E=${E#*=}
Ecoh[$count]=${E#*=}

rm -f log.lammps
cd ${PWD0}

unset Name
unset Types
unset Species
unset Masses

count=$(( $count + 1))
done

if [ ${#Ecoh[@]} -eq "1" ];then
Ecoh[1]=${Ecoh[0]} 
fi
}


#################################################################
#################################################################
#################################################################
####################### GB Simulations ##########################
#################################################################
#################################################################
#################################################################
lammps_gb()
{

# LAMMPS GBE computes (0=slab 1=fully periodic)
GBE[0]="(c_GBpe-(${Ecoh[0]}*\${natomsGB1}+${Ecoh[1]}*\${natomsGB2}))/\${gbarea}"
GBE[1]="(pe-(${Ecoh[0]}*\${natoms1}+${Ecoh[1]}*\${natoms2}))/(2*\${gbarea})"
GBE[0]=`(echo -e "${GBE[0]}" | tr -d '[[:space:]]')`
GBE[1]=`(echo -e "${GBE[1]}" | tr -d '[[:space:]]')`

# Start looping
T0="$(date +%s)"
echo "Running in $Work_Dir" > Progress.txt
echo " " >> Progress.txt
echo " -- Looping over $NGBs Structures " >> Progress.txt
echo " " >> Progress.txt


# Creating LAMMPS master script that will be used for each GB ##
cat > $Work_Dir/master_gb.in << _EOF_
log SED_NAME_SED.log

#################################################################
#    LAMMPS Input File for Grain Boundary Optimization
#    Shawn P. Coleman, July 2015
#################################################################

#################################################################
# --------------- Setup Simulation Variables ------------------ #
#################################################################

## General
variable     Name         string   SED_NAME_SED
variable     Species      string  "SED_SPECIES_SED"

## GB Properties
variable     xlen         equal    SED_LX_SED
variable     zlen         equal    SED_LZ_SED

variable     Dx            equal    \${xlen}/${Nx}
variable     Dz            equal    \${zlen}/${Nz} 

variable     NeighDist     equal    sqrt(3)/2*${a0}
variable     OverlapStart  equal    \${NeighDist}*0.1
variable     0verlapStop   equal    \${NeighDist}
variable     Doverlap      equal    (\${0verlapStop}-\${OverlapStart})/${Noverlap}

#################################################################
# --------------- Define loops for simulation ----------------- #
#################################################################

variable     counter      equal    0 
variable     atomprev     equal    0
variable     optdx        equal    0
variable     optdz        equal    0
variable     optoverlap   equal    0
variable     optnatoms    equal    0
variable     gbeprev      equal    100000000
variable     gbeRprev     equal    100000000

print        "Count dx dz overlapdist natoms gbe" file \${Name}.dat

## x-Translations 
label        loopa 
variable     a            loop     ${Nx} 
variable     tx           equal    (\${a}-1)*\${Dx}
 
## z-Translations  
label        loopb 
variable     b            loop     ${Nz} 
variable     tz           equal    (\${b}-1)*\${Dz}

## Delete Region (upper/lower, 1/2) 
label        loopd
variable     d            loop     ${DeleteRegion} 

## Overlap Distances
label        loopc 
variable     c           loop     ${Noverlap} 
variable     OverlapDist equal    0+\${Doverlap}*(\${c}-1) 


variable      counter     equal    \${counter}+1 
#################################################################
# ------------------- Initialize simulation ------------------- #
#################################################################
clear
print          "Counter: \${counter}" 
print          "Overlap Distance: \${OverlapDist}"

units           ${Units}
dimension       3 
boundary        p p p 
atom_style      ${AtomStyle}

read_data       SED_DATA_FILE_SED
atom_modify     sort 0 0.0

group           upper type SED_LAT1TYPES_SED
group           lower type SED_LAT2TYPES_SED
group           gb    type SED_GBTYPES_SED
group           gb1   intersect upper gb
group           gb2   intersect lower gb

${Potential}

balance         0.9 shift y 50 1.0

## Translations
displace_atoms  upper move \${tx} 0 \${tz} units box 

## Deletions 
if              "\$d == 1" then "delete_atoms overlap \${OverlapDist} lower upper" 
if              "\$d == 2" then "delete_atoms overlap \${OverlapDist} upper lower" 

## Check Uniqueness
variable        natoms   equal    count(all)
print           "Previous: \${atomprev}, Present: \${natoms}" 
if              "\${atomprev} == \${natoms}" then "jump SED_NAME_SED.in loopend"

## Computes 
compute         PE     all pe/atom
compute         GBpe   gb reduce sum c_PE

#################################################################
# ------------------------ Mimization ------------------------- #
#################################################################

## Constant Volume
reset_timestep  0 
thermo          10 
thermo_style    custom step pe lx ly lz press pxx pyy pzz c_GBpe
minimize        ${etol} ${ftol} ${maxiter} ${maxeval} 

## Box Relax
reset_timestep   0 
thermo           10 
thermo_style     custom step pe lx ly lz press pxx pyy pzz c_GBpe

fix              1 all box/relax aniso 0 vmax 0.001
minimize         ${etol} ${ftol} ${maxiter} ${maxeval}
unfix            1 

## Grain Boundary Energy Computation
variable        natoms1      equal  count(upper)
variable        natoms2      equal  count(lower)
variable        natomsGB     equal  count(gb)
variable        natomsGB1    equal  count(gb1)
variable        natomsGB2    equal  count(gb2)
variable        esum         equal  ${Ecoh}*\${natomsGB}
variable        gbarea       equal  lx*lz

variable        gbe          equal  SED_GBE_SED
variable        gbemJm2      equal  \${gbe}*$eConv*1e20
variable        gbernd       equal  round(\${gbemJm2})

## Determine Global Minimum GBE
print           "Previous: \${gbeprev}, Present: \${gbemJm2}" 
if              "\${gbemJm2} < \${gbeprev}" then &
                "variable     optdx        equal   \${tx}" &
                "variable     optdz        equal   \${tz}" &
                "variable     optoverlap   equal   \${OverlapDist}" &
                "variable     optnatoms    equal   \${natoms}" &
                "variable     gbeprev      equal   \${gbemJm2}" &
                "variable     gbeRprev     equal   \${gbernd}"

## Report Findings
reset_timestep  0 
dump            1 all custom 1 \${Name}_GBE_\${gbernd}.dump id type x y z c_PE ${DumpSpecial}
run             0 

print          "GB energy is \${gbemJm2} mJ/m^2" 
print          "\${counter} \${tx} \${tz} \${OverlapDist} \${natoms} \${gbernd}" append \${Name}.dat
 
variable        atomprev     equal  \${natoms}


#################################################################
# ----------------- Close loops for simulation ---------------- #
#################################################################
label           loopend 
next            c 
jump            SED_NAME_SED.in loopc 
variable        c delete 
next            d 
jump            SED_NAME_SED.in loopd 
variable        d delete 
next            b 
jump            SED_NAME_SED.in loopb 
variable        b delete 
next            a 
jump            SED_NAME_SED.in loopa 

print          'Optimal: \${optdx} \${optdz} \${optoverlap} \${optnatoms} \${gbeRprev} \${gbeprev}'
print          "All done"
_EOF_




# Looping over all data files here...
for i in ${Data_Dir}/${Find_Prefix}*${Find_Suffix}
do
    Tloop1="$(date +%s)"

    # Naming convention for files
    FullPath=`echo ${i} | sed 's./.\\/.g'` # full path 
    JustName=`echo ${i}  | sed 's/.*\///'` # datafile name only
    Name=${JustName%%$Find_Suffix}         # remove suffix 
    Name=${Name#*$Find_Prefix}             # remove prefix 
    mkdir -p $Work_Dir/$Name/

    # Obtaining values from datafile
    nTypes=`grep "atom types" $FullPath`
    nTypes=${nTypes%atom*}
  
    Types=`grep " Types:" $FullPath`
    Types=${Types#*:}
    IFS=', ' read -a Types <<< $Types

    Species=`grep Species: $FullPath`
    Species=${Species#*:}
    IFS=', ' read -a Species <<< $Species

    Masses=`grep Masses: $FullPath`
    Masses=${Masses#*:}
    IFS=', ' read -a Masses <<< $Masses

    Lat1Types=`grep Lat1Types: $FullPath`
    Lat1Types=${Lat1Types#*:}
    IFS=', ' read -a Lat1Types <<< $Lat1Types

    Lat2Types=`grep Lat2Types: $FullPath`
    Lat2Types=${Lat2Types#*:}
    IFS=', ' read -a Lat2Types <<< $Lat2Types

    GBTypes=`grep GBTypes: $FullPath`
    GBTypes=${GBTypes#*:}
    IFS=', ' read -a GBTypes <<< $GBTypes
    
    FullPeriod=`grep FullyPeriodic: $FullPath`
    FullPeriod=${FullPeriod#*:}
    FullPeriod=${FullPeriod%%(*}
    IFS=', ' read -a FullPeriod <<< $FullPeriod
    
    gx=`grep gx: $FullPath`
    gx=${gx#*:}
    IFS=', ' read -a gx <<< $gx

    gy=`grep gy: $FullPath`
    gy=${gy#*:}
    IFS=', ' read -a gy <<< $gy

    gz=`grep gz: $FullPath`
    gz=${gz#*:}
    IFS=', ' read -a gz <<< $gz

    Xdim=`grep xlo $FullPath`
    Xdim=${Xdim%xlo*}
    IFS=', ' read -a Xdim <<< $Xdim
    Xdim=`echo "scale=6;${Xdim[1]}-${Xdim[0]}" | bc`

    Zdim=`grep zlo $FullPath`
    Zdim=${Zdim%zlo*}
    IFS=', ' read -a Zdim <<< $Zdim
    Zdim=`echo "scale=6;${Zdim[1]}-${Zdim[0]}" | bc`

    # Assigning values in input that are specific for each GB
    cp  $Work_Dir/master_gb.in $Work_Dir/$Name/$Name.in
    sed -i "s+SED_DATA_FILE_SED+${FullPath}+g" $Work_Dir/$Name/$Name.in
    sed -i "s+SED_NAME_SED+${Name}+g" $Work_Dir/$Name/$Name.in
    sed -i "s+SED_NTYPES_SED+${Ntypes}+g" $Work_Dir/$Name/$Name.in
    sed -i "s+SED_LAT1TYPES_SED+${Lat1Types[*]}+g" $Work_Dir/$Name/$Name.in
    sed -i "s+SED_LAT2TYPES_SED+${Lat2Types[*]}+g" $Work_Dir/$Name/$Name.in
    sed -i "s+SED_GBTYPES_SED+${GBTypes[*]}+g" $Work_Dir/$Name/$Name.in
    sed -i "s+SED_SPECIES_SED+${Species[*]}+g" $Work_Dir/$Name/$Name.in  
    sed -i "s+SED_LX_SED+${gx}+g" $Work_Dir/$Name/$Name.in
    sed -i "s+SED_LZ_SED+${gz}+g" $Work_Dir/$Name/$Name.in   
    sed -i "s%SED_GBE_SED%${GBE[$FullPeriod]}%g" $Work_Dir/$Name/$Name.in
    
    cd $Work_Dir/$Name/
    $MPIrun $nMPI $LAMMPS < $Work_Dir/$Name/$Name.in >& $Name.out
    rm log.lammps

    OptCond=`grep Optimal: $Name.out`
    OptCond=${OptCond#*:}
    IFS=', ' read -a OptCond <<< $OptCond

    cp ${Name}_GBE_${OptCond[4]}.dump $Work_Dir/min_gb/

    Tloop2="$(($(date +%s)-Tloop1))"
    T1="$(($(date +%s)-T0))"
    printf "GB $Name: Loop Time: %02d:%02d:%02d:%02d " "$(($Tloop2/86400))" "$(($Tloop2/3600%24))" "$(($Tloop2/60%60))" "$(($Tloop2%60))" >> $PWD0/Progress.txt
    printf "Total Time: %02d:%02d:%02d:%02d  GBE: %8.8f\n" "$(($T1/86400))" "$(($T1/3600%24))" "$(($T1/60%60))" "$(($T1%60))" "${OptCond[5]}" >> $PWD0/Progress.txt
   
    unset FullPath
    unset JustName
    unset Name
    unset Types
    unset Species
    unset Masses
    unset Lat1Types
    unset Lat2Types
    unset GBTypes
    
done
rm $Work_Dir/${Find_Prefix}master_gb.in
}


lammps_bulk
lammps_gb


