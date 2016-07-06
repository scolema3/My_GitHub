variable    tmp             equal ${NSteps}+1
variable    LoadLoopStep    loop ${tmp}
variable    tmp             equal ${NSteps}*2+1
variable    UnLoadLoopStep  loop ${tmp}

#################################################################
# -----  Setup Variables for Computing Elastic Constants ------ #
#################################################################
 
## Loading+Unloading
variable Avepxx              equal 0
variable CovPxx              equal 0
variable VarPxx              equal 0
variable Avepyy              equal 0
variable CovPyy              equal 0
variable VarPyy              equal 0
variable Avepzz              equal 0
variable CovPzz              equal 0
variable VarPzz              equal 0
variable Avepxy              equal 0
variable CovPxy              equal 0
variable VarPxy              equal 0
variable Avepxz              equal 0
variable CovPxz              equal 0
variable VarPxz              equal 0
variable Avepyz              equal 0
variable CovPyz              equal 0
variable VarPyz              equal 0
variable AveStrain           equal 0

## Loading
variable AveLpxx             equal 0
variable CovLPxx             equal 0
variable VarLPxx             equal 0
variable AveLpyy             equal 0
variable CovLPyy             equal 0
variable VarLPyy             equal 0
variable AveLpzz             equal 0
variable CovLPzz             equal 0
variable VarLPzz             equal 0
variable AveLpxy             equal 0
variable CovLPxy             equal 0
variable VarLPxy             equal 0
variable AveLpxz             equal 0
variable CovLPxz             equal 0
variable VarLPxz             equal 0
variable AveLpyz             equal 0
variable CovLPyz             equal 0
variable VarLPyz             equal 0
variable AveLStrain          equal 0

## Unloading
variable AveUpxx             equal 0
variable CovUPxx             equal 0
variable VarUPxx             equal 0
variable AveUpyy             equal 0
variable CovUPyy             equal 0
variable VarUPyy             equal 0
variable AveUpzz             equal 0
variable CovUPzz             equal 0
variable VarUPzz             equal 0
variable AveUpxy             equal 0
variable CovUPxy             equal 0
variable VarUPxy             equal 0
variable AveUpxz             equal 0
variable CovUPxz             equal 0
variable VarUPxz             equal 0
variable AveUpyz             equal 0
variable CovUPyz             equal 0
variable VarUPyz             equal 0
variable AveUStrain          equal 0


#################################################################
# -----------  Clear and Set Reference State Data ------------- #
#################################################################

if "${dir} == 1" then &
   "variable len0 equal ${lx0}"  
if "${dir} == 2" then &
   "variable len0 equal ${ly0}"  
if "${dir} == 3" then &
   "variable len0 equal ${lz0}"  
if "${dir} == 4" then &
   "variable len0 equal ${lz0}"  
if "${dir} == 5" then &
   "variable len0 equal ${lz0}"  
if "${dir} == 6" then &
   "variable len0 equal ${ly0}"  

clear
box tilt large
read_restart ${Name}_restart.equil
include Coleman_potential.mod

label           LoadLoop
#################################################################
# --------------------------  Load ---------------------------- #
#################################################################

if "${LoadLoopStep} == 1" then &
   "# Intial Conditions & Structure : " &
   "variable     up            equal    0" &
   "variable     strain        equal    0" &
   "print '% Load_Step Strain pxx pyy pzz pxy pxz pyz' append ${Name}.m " &
   "print 'Load(${dir}).Data = [' append ${Name}.m " &
elif "${LoadLoopStep} > 1"  &
   "variable     up            equal    ${LExpand}/${NSteps}"

variable delta equal ${up}*${len0}

#################################################################################################### SHOULD THIS BE BASED ON THE ORIGINAL XY XZ AND YZ VALUES --- TESTING ON EXCALIBUR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#################################################################################################### SHOULD THIS BE BASED ON THE ORIGINAL XY XZ AND YZ VALUES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
variable deltaxy equal ${up}*xy
variable deltaxz equal ${up}*xz
variable deltayz equal ${up}*yz
#################################################################################################### SHOULD THIS BE BASED ON THE ORIGINAL XY XZ AND YZ VALUES --- TESTING ON EXCALIBUR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#################################################################################################### SHOULD THIS BE BASED ON THE ORIGINAL XY XZ AND YZ VALUES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if "${dir} == 1" then &
   "change_box all x delta 0 ${delta} xy delta ${deltaxy} xz delta ${deltaxz} remap units box"
if "${dir} == 2" then &
   "change_box all y delta 0 ${delta} yz delta ${deltayz} remap units box"
if "${dir} == 3" then &
   "change_box all z delta 0 ${delta} remap units box"
if "${dir} == 4" then &
   "change_box all yz delta ${delta} remap units box"
if "${dir} == 5" then &
   "change_box all xz delta ${delta} remap units box"
if "${dir} == 6" then &
   "change_box all xy delta ${delta} remap units box"

minimize        ${etol} ${ftol} ${maxiter} ${maxeval}

if "${DumpFilesBool} == 1" then &
   "dump            1 all custom 1 ${Name}_Load_${dir}_${LoadLoopStep}.dump id type x y z ${OptDumpData}" &
   "run             0" &
   "undump          1"

###################### Obtain stress tensor #####################

variable strain              equal ${strain}+${up}

variable Avepxx              equal ${Avepxx}+${tmp1}/(3*${NSteps}+2)
variable CovPxx              equal ${CovPxx}+${tmp1}*${strain}
variable VarPxx              equal ${VarPxx}+${strain}*${strain}
variable AveLpxx             equal ${AveLpxx}+${tmp1}/(${NSteps}+1)
variable CovLPxx             equal ${CovLPxx}+${tmp1}*${strain}
variable VarLPxx             equal ${VarLPxx}+${strain}*${strain}

variable Avepyy              equal ${Avepyy}+${tmp2}/(3*${NSteps}+2)
variable CovPyy              equal ${CovPyy}+${tmp2}*${strain}
variable VarPyy              equal ${VarPyy}+${strain}*${strain}
variable AveLpyy             equal ${AveLpyy}+${tmp2}/(${NSteps}+1)
variable CovLPyy             equal ${CovLPyy}+${tmp2}*${strain}
variable VarLPyy             equal ${VarLPyy}+${strain}*${strain}

variable Avepzz              equal ${Avepzz}+${tmp3}/(3*${NSteps}+2)
variable CovPzz              equal ${CovPzz}+${tmp3}*${strain}
variable VarPzz              equal ${VarPzz}+${strain}*${strain}
variable AveLpzz             equal ${AveLpzz}+${tmp3}/(${NSteps}+1)
variable CovLPzz             equal ${CovLPzz}+${tmp3}*${strain}
variable VarLPzz             equal ${VarLPzz}+${strain}*${strain}

variable Avepxy              equal ${Avepxy}+${tmp4}/(3*${NSteps}+2)
variable CovPxy              equal ${CovPxy}+${tmp4}*${strain}
variable VarPxy              equal ${VarPxy}+${strain}*${strain}
variable AveLpxy             equal ${AveLpxy}+${tmp4}/(${NSteps}+1)
variable CovLPxy             equal ${CovLPxy}+${tmp4}*${strain}
variable VarLPxy             equal ${VarLPxy}+${strain}*${strain}

variable Avepxz              equal ${Avepxz}+${tmp5}/(3*${NSteps}+2)
variable CovPxz              equal ${CovPxz}+${tmp5}*${strain}
variable VarPxz              equal ${VarPxz}+${strain}*${strain}
variable AveLpxz             equal ${AveLpxz}+${tmp5}/(${NSteps}+1)
variable CovLPxz             equal ${CovLPxz}+${tmp5}*${strain}
variable VarLPxz             equal ${VarLPxz}+${strain}*${strain}

variable Avepyz              equal ${Avepyz}+${tmp6}/(3*${NSteps}+2)
variable CovPyz              equal ${CovPyz}+${tmp6}*${strain}
variable VarPyz              equal ${VarPyz}+${strain}*${strain}
variable AveLpyz             equal ${AveLpyz}+${tmp6}/(${NSteps}+1)
variable CovLPyz             equal ${CovLPyz}+${tmp6}*${strain}
variable VarLPyz             equal ${VarLPyz}+${strain}*${strain}


variable AveStrain              equal ${AveStrain}+${strain}/(3*${NSteps}+2)
variable AveLStrain             equal ${AveLStrain}+${strain}/(${NSteps}+1)

print "${LoadLoopStep} ${strain} ${tmp1} ${tmp2} ${tmp3} ${tmp4} ${tmp5} ${tmp6}" append ${Name}.m 

next            LoadLoopStep
jump            SELF LoadLoop


print "];" append ${Name}.m 




label           UnLoadLoop
#################################################################
# -------------------------  UnLoad --------------------------- #
#################################################################

if "${UnLoadLoopStep} == 1" then &
   "# Intial Structure: " &
   "variable     up            equal    0" &
   "print '% UnLoad_Step Strain pxx pyy pzz pxy pxz pyz' append ${Name}.m " &
   "print 'UnLoad(${dir}).Data = [' append ${Name}.m " &
elif "${UnLoadLoopStep} > 1"  &
   "variable     up            equal    -${LExpand}/${NSteps}"

variable delta equal ${up}*${len0}
variable deltaxy equal ${up}*xy
variable deltaxz equal ${up}*xz
variable deltayz equal ${up}*yz

if "${dir} == 1" then &
   "change_box all x delta 0 ${delta} xy delta ${deltaxy} xz delta ${deltaxz} remap units box"
if "${dir} == 2" then &
   "change_box all y delta 0 ${delta} yz delta ${deltayz} remap units box"
if "${dir} == 3" then &
   "change_box all z delta 0 ${delta} remap units box"
if "${dir} == 4" then &
   "change_box all yz delta ${delta} remap units box"
if "${dir} == 5" then &
   "change_box all xz delta ${delta} remap units box"
if "${dir} == 6" then &
   "change_box all xy delta ${delta} remap units box"

minimize        ${etol} ${ftol} ${maxiter} ${maxeval}
if "${DumpFilesBool} == 1" then &
   "dump            1 all custom 1 ${Name}_UnLoad_${dir}_${LoadLoopStep}.dump id type x y z ${OptDumpData}" &
   "run             0" &
   "undump          1"

###################### Obtain stress tensor #####################

variable strain              equal ${strain}+${up}
 
variable Avepxx              equal ${Avepxx}+${tmp1}/(3*${NSteps}+2)
variable CovPxx              equal ${CovPxx}+${tmp1}*${strain}
variable VarPxx              equal ${VarPxx}+${strain}*${strain}
variable AveUpxx             equal ${AveUpxx}+${tmp1}/(${NSteps}+1)
variable CovUPxx             equal ${CovUPxx}+${tmp1}*${strain}
variable VarUPxx             equal ${VarUPxx}+${strain}*${strain}

variable Avepyy              equal ${Avepyy}+${tmp2}/(3*${NSteps}+2)
variable CovPyy              equal ${CovPyy}+${tmp2}*${strain}
variable VarPyy              equal ${VarPyy}+${strain}*${strain}
variable AveUpyy             equal ${AveUpyy}+${tmp2}/(${NSteps}+1)
variable CovUPyy             equal ${CovUPyy}+${tmp2}*${strain}
variable VarUPyy             equal ${VarUPyy}+${strain}*${strain}

variable Avepzz              equal ${Avepzz}+${tmp3}/(3*${NSteps}+2)
variable CovPzz              equal ${CovPzz}+${tmp3}*${strain}
variable VarPzz              equal ${VarPzz}+${strain}*${strain}
variable AveUpzz             equal ${AveUpzz}+${tmp3}/(${NSteps}+1)
variable CovUPzz             equal ${CovUPzz}+${tmp3}*${strain}
variable VarUPzz             equal ${VarUPzz}+${strain}*${strain}

variable Avepxy              equal ${Avepxy}+${tmp4}/(3*${NSteps}+2)
variable CovPxy              equal ${CovPxy}+${tmp4}*${strain}
variable VarPxy              equal ${VarPxy}+${strain}*${strain}
variable AveUpxy             equal ${AveUpxy}+${tmp4}/(${NSteps}+1)
variable CovUPxy             equal ${CovUPxy}+${tmp4}*${strain}
variable VarUPxy             equal ${VarUPxy}+${strain}*${strain}

variable Avepxz              equal ${Avepxz}+${tmp5}/(3*${NSteps}+2)
variable CovPxz              equal ${CovPxz}+${tmp5}*${strain}
variable VarPxz              equal ${VarPxz}+${strain}*${strain}
variable AveUpxz             equal ${AveUpxz}+${tmp5}/(${NSteps}+1)
variable CovUPxz             equal ${CovUPxz}+${tmp5}*${strain}
variable VarUPxz             equal ${VarUPxz}+${strain}*${strain}

variable Avepyz              equal ${Avepyz}+${tmp6}/(3*${NSteps}+2)
variable CovPyz              equal ${CovPyz}+${tmp6}*${strain}
variable VarPyz              equal ${VarPyz}+${strain}*${strain}
variable AveUpyz             equal ${AveUpyz}+${tmp6}/(${NSteps}+1)
variable CovUPyz             equal ${CovUPyz}+${tmp6}*${strain}
variable VarUPyz             equal ${VarUPyz}+${strain}*${strain}

print "${UnLoadLoopStep} ${strain} ${tmp1} ${tmp2} ${tmp3} ${tmp4} ${tmp5} ${tmp6}" append ${Name}.m 

next            UnLoadLoopStep
jump            SELF UnLoadLoop

print "];" append ${Name}.m



#################################################################
# --------------  Compute Eleasic Constants ------------------- #
#################################################################

# NOTE: Elastic constants are computed from pressure data using 
#       a linear regression WITHOUT intercept, i.e. y = mx

# Loading
variable LC1${dir} equal -${CovLPxx}/${VarLPxx}
variable LC2${dir} equal -${CovLPyy}/${VarLPyy}
variable LC3${dir} equal -${CovLPzz}/${VarLPzz}
variable LC4${dir} equal -${CovLPxy}/${VarLPxy}
variable LC5${dir} equal -${CovLPxz}/${VarLPxz}
variable LC6${dir} equal -${CovLPyz}/${VarLPyz}

# Unloading
variable UC1${dir} equal -${CovUPxx}/${VarUPxx}
variable UC2${dir} equal -${CovUPyy}/${VarUPyy}
variable UC3${dir} equal -${CovUPzz}/${VarUPzz}
variable UC4${dir} equal -${CovUPxy}/${VarUPxy}
variable UC5${dir} equal -${CovUPxz}/${VarUPxz}
variable UC6${dir} equal -${CovUPyz}/${VarUPyz}

# Loading+Unloading
variable C1${dir} equal -${CovPxx}/${VarPxx}
variable C2${dir} equal -${CovPyy}/${VarPyy}
variable C3${dir} equal -${CovPzz}/${VarPzz}
variable C4${dir} equal -${CovPxy}/${VarPxy}
variable C5${dir} equal -${CovPxz}/${VarPxz}
variable C6${dir} equal -${CovPyz}/${VarPyz}

