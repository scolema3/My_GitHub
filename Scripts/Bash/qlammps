#!/bin/sh 

# This script works with a LAMMPS_Base to compute diffraction patterns for all files within the current dir containing the same prefix OR suffix.  The script is designed to work with both LAMMPS dump and data files.

Version=0.18

printf '\n  Running Version %s of %s\n' ${Version} $0

# help guide
help_ani()
{
  printf '
  Usage: Qlammps (-Optional arguments) LAMMPS_INPUT

  The last argument is the LAMMPS input deck

  Optional arguments:

    -A account id for time charging
       default: %s 
    -q name of queue i.e. d=debug s=standard p=phi
       default: %s
    -n number of nodes requested (%s cores)
       default: %s
    -P control the number of processors (MPI, OpenMP, #)
       default: %s       
    -t wall-time for computation in dd:hh:mm
       default: %s
    -M email option y or n
       default: %s
    -E email address for to send correspondence
       default: %s
  ' ${ACCOUNT} ${QUEUE} ${NCPUS} ${NODES} ${PPN_SET} ${Time} ${MAIL} ${EMAIL}
  exit 1
}


#############  General Default Values ############# 
nInput=$# 
isdef=0
PWD=`pwd`

# LAMMPS
Application_Name=LAMMPS
LAMMPS=~/../scoleman/Public/bin/lmp_diff   # mic versions included in phi queue options

# Naming 
Name="${@:${#@}}"    # use last input
InputPrefix=""
InputSuffix=".in"


# Job
QUEUE="s"
NODES=4
Time_Debug="00:10"                 # 10 min
Time_Standard="10:00"              # 10 hr
ACCOUNT="ARLAP38753387"
EMAIL=""
MAIL="n"

# Defualt by user
if [ $USER == "scoleman" ]; then 
  EMAIL="scolema3@gmail.com"
  MAIL="y"
elif [ $USER == "ehernan"]; then 
 EMAIL="efrain.hernandez18.ctr@mail.mil"
 MAIL="y"
fi

# PPN_SET=OMP for using max OpenMP Threads
# PPN_SET=MPI for using max MPI processes
# PPN_SET=# for directly specifing # of processors per node
PPN_SET=M

# System
Host=`hostname`
case ${Host} in
  ls*)        CLUSTER='Spirit';;    # workstaiton, debugging
  shepard*)   CLUSTER='Shepard';;
  armstrong*) CLUSTER='Armstrong';;
  garnet*)    CLUSTER='Garnet';;
  k*)         CLUSTER='Kilrain';;
  h*)         CLUSTER='Haise';;
  excalibur*) CLUSTER='Excalibur';;
  spirit*)    CLUSTER='Spirit';;
  gordon*)    CLUSTER='Gordon';;
  conrad*)    CLUSTER='Conrad';;
  thunder*)   CLUSTER='Thunder';;    
  topaz*)     CLUSTER='Topaz';;    
  copper*)    CLUSTER='Copper';;
  *)  "Host not recognized! "; exit 0 ;;
esac

# Optional inputs flags
while (($#)); do
  case "$1" in
    -A) ACCOUNT="$2";;
    -n) NODES="$2";;
    -t) TIME="$2";;
    -M) MAIL="$2";;
    -q) QUEUE="$2";;
    -P) PPN_SET="$2"; shift;;
    -E) EMAIL="$2";;
    -*)  echo "${@:1} is not a valid input"; exit 1;;
    *)  ;;
  esac
  shift
done

#################################################################
#################################################################
##         Defaults for Specific DOD-HPCMP Machines            ##
#################################################################
#################################################################

############# Excalibur Default Values #############
if [ $CLUSTER == "Excalibur" ]; then

  ## Architecture and Environment
  NCPUS=32
  if [ ${PPN_SET:0:1} == "M" ]; then
    PPN=${NCPUS}
  elif [ ${PPN_SET:0:1} == "O" ]; then
    PPN=2
  else
    PPN=${PPN_SET}
  fi
  nOMP=$(( $NCPUS / $PPN ))
  nMPI=$(( $PPN * $NODES ))
  Setup_Env="export OMP_NUM_THREADS=${nOMP}"

  ## Specific Queue Options ##
  if [ ${QUEUE:0:1} == "d" ]; then
    QUEUE="debug"
    Time=${Time_Debug}
    if [ $NODES -gt 16 ]; then
      NODES=16
      echo "Exceded node max: -n set back to $NODES"
    fi
  elif [ ${QUEUE:0:1} == "s" ]; then
    QUEUE="standard"
    Time=${Time_Standard}
  else
    echo "Unknown queue requested for $CLUSTER"
    exit 1   
  fi
  
  ## PBS and Submission ##
  PBS_extra="#PBS -l place=scatter:excl"
  MPI_launch="aprun -n ${nMPI} -N ${PPN} -d ${nOMP}"
  Submit_Command=qsub


############# Spirit Default Values #############
elif [ $CLUSTER == "Spirit" ]; then
  
  ## Architecture and Environment
  NCPUS=16
  if [ ${PPN_SET:0:1} == "M" ]; then
    PPN=${NCPUS}
  elif [ ${PPN_SET:0:1} == "O" ]; then
    PPN=2
  else
    PPN=${PPN_SET}
  fi
  nOMP=$(( $NCPUS / $PPN ))
  nMPI=$(( $PPN * $NODES ))
  Setup_Env="export OMP_NUM_THREADS=${nOMP}"

  ## Specific Queue Options ##    
  if [ ${QUEUE:0:1} == "d" ]; then
    QUEUE="debug"
    Time=${Time_Debug}
    if [ $NODES -gt 46 ]; then
      NODES=46
    fi
  elif [ ${QUEUE:0:1} == "s" ]; then
    QUEUE="standard"
    Time=${Time_Standard}
    if [ $NODES -gt 2295 ]; then
      NODES=2295
      echo "Exceded node max: -n set back to $NODES"
    fi
  else
    echo "Unknown queue requested for $CLUSTER"
    exit 1      
  fi

  ## PBS and Submission ##  
  PBS_extra=""  
  MPI_launch="mpiexec_mpt -n ${nMPI} omplace"
  Submit_Command=qsub


############# Garnet or Copper Default Values #############
elif [ $CLUSTER == "Garnet" ] || [ $CLUSTER == "Copper" ] ; then
  
  ## Architecture and Environment
  NCPUS=32
  if [ ${PPN_SET:0:1} == "M" ]; then
    PPN=${NCPUS}
  elif [ ${PPN_SET:0:1} == "O" ]; then
    PPN=2
  else
    PPN=${PPN_SET}
  fi
  nOMP=$(( $NCPUS / $PPN ))
  nMPI=$(( $PPN * $NODES ))
  Setup_Env="export OMP_NUM_THREADS=${nOMP}"

  ## Specific Queue Options ##    
  if [ ${QUEUE:0:1} == "d" ]; then
    QUEUE="debug"
    Time=${Time_Debug}
    if [ $NODES -gt 256 ]; then
      NODES=256
    fi
  elif [ ${QUEUE:0:1} == "s" ]; then
    QUEUE="standard_sm"
    Time=${Time_Standard}
    if [ $NODES -gt 160 ]; then
      QUEUE="standard_lg"
    elif [ $NODES -gt 3200 ]; then
      NODES=3200
      echo "Exceded node max: -n set back to $NODES"
    fi
  else
    echo "Unknown queue requested for $CLUSTER"
    exit 1      
  fi

  ## PBS and Submission ##  
  PBS_extra=""  
  MPI_launch="aprun -n ${nMPI} -N ${PPN} -d ${nOMP}"
  Submit_Command=qsub


############# Shepard or Armstrong Default Values #############
elif [ $CLUSTER == "Shepard" ] || [ $CLUSTER == "Armstrong" ] ; then
  
  ## Architecture and Environment
  NCPUS=24
  if [ ${PPN_SET:0:1} == "M" ]; then
    PPN=${NCPUS}
  elif [ ${PPN_SET:0:1} == "O" ]; then
    PPN=1
  else
    PPN=${PPN_SET}
  fi
  nOMP=$(( $NCPUS / $PPN ))
  nMPI=$(( $PPN * $NODES ))
  Setup_Env="export OMP_NUM_THREADS=${nOMP}"
  
  ## Specific Queue Options ## 
  if [ ${QUEUE:0:1} == "d" ]; then
    QUEUE="debug"
    Time=${Time_Debug}
    if [ $NODES -gt 150 ]; then
      NODES=150
      echo "Exceded node max: -n set back to $NODES"
    fi
  elif [ ${QUEUE:0:1} == "s" ]; then
    QUEUE="standard"
    Time=${Time_Standard}
    if [ $NODES -gt 480 ]; then
      NODES=480
      echo "Exceded node max: -n set back to $NODES"
    fi
  elif [ ${QUEUE:0:1} == "p" ]; then
    QUEUE="phi"
    LAMMPS=~/../scoleman/Public/bin/lmp_diff_mic
    NCPUS=10
    MIC=236
    if [ ${PPN_SET:0:1} == "M" ]; then
     PPN=${NCPUS}
    elif [ ${PPN_SET:0:1} == "O" ]; then
     PPN=1
    else
     PPN=${PPN_SET}
    fi
    nOMP=$(( $NCPUS / $PPN ))
    nMPI=$(( $PPN * $NODES ))   
    Time=${Time_Standard}
    Setup_Env="
module swap PrgEnv-cray PrgEnv-intel 
module unload libsci atp
source /opt/intel/composer_xe_2013_sp1.2.144/bin/compilervars.sh intel64
        
export OMP_NUM_THREADS=${nOMP}  
export MIC_ENV_PREFIX=MIC
export MIC_OMP_NUM_THREADS=${MIC}
export MIC_KMP_AFFINITY=\"granularity=fine,compact\"
export MIC_KMP_PLACE_THREADS=\"59c,4t\" "   
    if [ $NODES -gt 124 ]; then
      NODES=124
      echo "Exceded node max: -n set back to $NODES"
    fi    
  fi

  ## PBS and Submission ##  
  PBS_extra=""  
  MPI_launch="aprun -n ${nMPI} -N ${PPN} -d ${nOMP} -cc none"
  Submit_Command=qsub



############# Conrad or Gordon Default Values #############
elif [ $CLUSTER == "Conrad" ] || [ $CLUSTER == "Gordon" ] ; then
  
  ## Architecture and Environment
  NCPUS=32
  if [ ${PPN_SET:0:1} == "M" ]; then
    PPN=${NCPUS}
  elif [ ${PPN_SET:0:1} == "O" ]; then
    PPN=1
  else
    PPN=${PPN_SET}
  fi
  nOMP=$(( $NCPUS / $PPN ))
  nMPI=$(( $PPN * $NODES ))
  Setup_Env="export OMP_NUM_THREADS=${nOMP}"
  
  ## Specific Queue Options ## 
  if [ ${QUEUE:0:1} == "d" ]; then
    QUEUE="debug"
    Time=${Time_Debug}
    if [ $NODES -gt 48 ]; then
      NODES=48
      echo "Exceded node max: -n set back to $NODES"
    fi
  elif [ ${QUEUE:0:1} == "s" ]; then
    QUEUE="standard"
    Time=${Time_Standard}
    if [ $NODES -gt 70 ]; then
      NODES=70
      echo "Exceded node max: -n set back to $NODES"
    fi
  elif [ ${QUEUE:0:1} == "p" ]; then
    QUEUE="phi"
    LAMMPS=~/../scoleman/Public/bin/lmp_diff_mic
    NCPUS=10
    MIC=236
    if [ ${PPN_SET:0:1} == "M" ]; then
     PPN=${NCPUS}
    elif [ ${PPN_SET:0:1} == "O" ]; then
     PPN=1
    else
     PPN=${PPN_SET}
    fi
    nOMP=$(( $NCPUS / $PPN ))
    nMPI=$(( $PPN * $NODES ))   
    Time=${Time_Standard}
    Setup_Env="
export OMP_NUM_THREADS=${nOMP}  
export MIC_ENV_PREFIX=MIC
export MIC_OMP_NUM_THREADS=${MIC}
export MIC_KMP_AFFINITY=\"granularity=fine,compact\"
export MIC_KMP_PLACE_THREADS=\"59c,4t\" "   
    if [ $NODES -gt 168 ]; then
      NODES=168
      echo "Exceded node max: -n set back to $NODES"
    fi    
  fi

  ## PBS and Submission ##  
  PBS_extra=""  
  MPI_launch="aprun -n ${nMPI} -N ${PPN} -d ${nOMP} -cc none"
  Submit_Command=qsub



############# Haise or Kilrain Default Values #############
elif [ $CLUSTER == "Haise"  ] || [ $CLUSTER == "Kilrain"  ]; then

  ## Architecture and Environment
  NCPUS=16
  MIC=236
  if [ ${PPN_SET:0:1} == "M" ]; then
    PPN=${NCPUS}
  elif [ ${PPN_SET:0:1} == "O" ]; then
    PPN=1
  else
    PPN=${PPN_SET}
  fi
  nOMP=$(( $NCPUS / $PPN ))
  nMPI=$(( $PPN * $NODES ))
  Setup_Env="export OMP_NUM_THREADS=${nOMP}"
  
  ## Specific Queue Options ##     
  if [ ${QUEUE:0:1} == "d" ]; then
    QUEUE="debug"
    Time=${Time_Debug}
    if [ $NODES -gt 64 ]; then
      NODES=64
      echo "Exceded node max: -n set back to $NODES"
    fi
  elif [ ${QUEUE:0:1} == "s" ]; then
    QUEUE="standard"
    Time=${Time_Standard}
    if [ $NODES -gt 256 ]; then
      NODES=256
      echo "Exceded node max: -n set back to $NODES"
    fi   
  elif [ ${QUEUE:0:1} == "p" ]; then
    QUEUE="phi"
    LAMMPS=~/../scoleman/Public/bin/lmp_diff_mic
    Time=${Time_Standard}
    Setup_Env="
module swap PrgEnv-cray PrgEnv-intel 
    
export OMP_NUM_THREADS=${nOMP}  
export MIC_ENV_PREFIX=MIC
export MIC_OMP_NUM_THREADS=${MIC}
export MIC_KMP_AFFINITY=\"granularity=fine,compact\"
export MIC_KMP_PLACE_THREADS=\"59c,4t\" "        
    if [ $NODES -gt 12 ]; then
      NODES=12
      echo "Exceded node max: -n set back to $NODES"
   fi
  else
    echo "Unknown queue requested for $CLUSTER"
    exit 1   
  fi

  ## PBS and Submission ##   
  PBS_extra=""  
  MPI_launch="mpirun -np ${nMPI} -ppn ${PPN}"
  Submit_Command=qsub



############# Topaz Default Values #############
elif [ $CLUSTER == "Topaz" ]; then
  
  ## Architecture and Environment
  NCPUS=36
  if [ ${PPN_SET:0:1} == "M" ]; then
    PPN=${NCPUS}
  elif [ ${PPN_SET:0:1} == "O" ]; then
    PPN=2
  else
    PPN=${PPN_SET}
  fi
  nOMP=$(( $NCPUS / $PPN ))
  nMPI=$(( $PPN * $NODES ))
  Setup_Env="export OMP_NUM_THREADS=${nOMP}"

  ## Specific Queue Options ##    
  if [ ${QUEUE:0:1} == "d" ]; then
    QUEUE="debug"
    Time=${Time_Debug}
    if [ $NODES -gt 180 ]; then
      NODES=180
    fi
  elif [ ${QUEUE:0:1} == "s" ]; then
    QUEUE="standard"
    Time=${Time_Standard}
    if [ $NODES -gt 1728 ]; then
      NODES=1728
      echo "Exceded node max: -n set back to $NODES"
    fi
  else
    echo "Unknown queue requested for $CLUSTER"
    exit 1      
  fi

  ## PBS and Submission ##  
  PBS_extra=""  
  MPI_launch="mpiexec_mpt -n ${nMPI} omplace"
  Submit_Command=qsub



############# Thunder Default Values #############
elif [ $CLUSTER == "Thunder" ]; then
  
  ## Architecture and Environment
  NCPUS=36
  if [ ${PPN_SET:0:1} == "M" ]; then
    PPN=${NCPUS}
  elif [ ${PPN_SET:0:1} == "O" ]; then
    PPN=2
  else
    PPN=${PPN_SET}
  fi
  nOMP=$(( $NCPUS / $PPN ))
  nMPI=$(( $PPN * $NODES ))
  Setup_Env="export OMP_NUM_THREADS=${nOMP}"

  ## Specific Queue Options ##    
  if [ ${QUEUE:0:1} == "d" ]; then
    QUEUE="debug"
    Time=${Time_Debug}
    if [ $NODES -gt 65 ]; then
      NODES=65
    fi
  elif [ ${QUEUE:0:1} == "s" ]; then
    QUEUE="standard"
    Time=${Time_Standard}
    if [ $NODES -gt 1608 ]; then
      NODES=1608
      echo "Exceded node max: -n set back to $NODES"
    fi
  else
    echo "Unknown queue requested for $CLUSTER"
    exit 1      
  fi

  ## PBS and Submission ##  
  PBS_extra=""  
  MPI_launch="mpiexec_mpt -n ${nMPI} omplace"
  Submit_Command=qsub


fi


#################################################################
#################################################################

### Error Checks

# If no inputs--display help
if [ $nInput -lt "1" ]; then
  help_ani
fi

### Naming formulation - if like file already exists increase suffix by 1

Title=`head -1 $Name | awk '{print $NF}'`

FLAG=0
test=$Title
count=0
while [  $FLAG == 0 ]; do
  if [ -f "run.$test" ]; then
    let count=count+1
    test=${Title}_$count
    else
    FLAG=1
    Title=$test
  fi
done

LammpsInput=${InputPrefix}${Title}${InputSuffix}

### Create specific file for new input script
cp ./$Name ./${LammpsInput}

### PBS mail commands
if [ $MAIL == "y" ]; then
PBS_mail="#PBS -m abe
#PBS -M $EMAIL"
else
PBS_mail="#"
fi

### Use TIME inputed on command line over default values
if [ -n "$TIME" ]; then
  Time=$TIME
fi

### Format time to be hh:mm:ss (at minimum)
Ncolon=`echo "${Time}" | grep -o ":" | wc -l`
if [ $Ncolon == 0 ]; then
  Time="00:$Time"
fi
PBS_time_suffix=":00"

Submit_File="run.${Title}"

### Create a specific runscript for each input script
cat > ${Submit_File} << _EOF_
#!/bin/bash
#PBS -N ${Title:0:9}
#PBS -l application=${Application_Name}
#PBS -l select=${NODES}:ncpus=${NCPUS}:mpiprocs=${PPN}
#PBS -l walltime=${Time}${PBS_time_suffix}
#PBS -q ${QUEUE}
#PBS -A ${ACCOUNT}
#PBS -j oe
${PBS_mail}
${PBS_extra}

echo "This job was created by Version ${Version} of $0"

cd ${PWD}

${Setup_Env}

T0="\$(date +%s)"
${MPI_launch} ${LAMMPS} < ${LammpsInput}

T1="\$((\$(date +%s)-T0))"
printf "Compute Time : %02d:%02d:%02d:%02d\n" "\$((\$T1/86400))" "\$((\$T1/3600%24))" "\$((\$T1/60%60))" "\$((\$T1%60))"""

_EOF_

echo "  Submitting $Title on $CLUSTER"

${Submit_Command} ${Submit_File}

echo " "

