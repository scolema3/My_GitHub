import numpy
import math
import os
#from submit_maker import *
import time

#NOTE: ThermalData.py for the given material and crystaltype must have previously been run


#Inputs:

miller1x=numpy.array([2, -1, -1, 0])
miller1y=numpy.array([0, 1, -1, 1])
miller1z=numpy.array([1, 0, -1, 2])
miller2x=numpy.array([-2, 1, 1, 0])
miller2y=numpy.array([0, 1, -1, 1])
miller2z=numpy.array([1, 0, 1, -2])

datafile='mgtwin_thick.apf'



#Crystalography
material='Mg'
crystaltype='hcp'       #bcc, fcc, or hcp
a=3.21                  #Lattice constant approximation


#Temperatures to be analyzed
Temps=[0]

#Driving energies to be considered (+/- u0/2 added to each crystal for ECO, u0 added to one crystal for orient)
u0=[.005, .010, .015, .020, .025];
u0=[.005]
#Computational machine
machine='topaz'

#Parameters
rcut=3.9                  #ECO cutoff radius
eta=.25                    #ECO cutoff parmater (.25 in most cases)
cutlo=.4                   #orient cutoff parameters (between 0 and 1)
cuthi=.6
timesteps=3000

#Potential
pair_style='eam/fs'
pair_coeff='* * Mg.eam.fs Mg'
mass=24.3

#End of Inputs


os.system('mkdir MobilityData')
os.system('mkdir MobilityData/'+material+'_'+crystaltype)
Orientation='Largetest'
os.system('mkdir MobilityData/'+material+'_'+crystaltype+'/'+Orientation)
os.system('mkdir MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/InputScripts')
os.system('mkdir MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Dumps')
os.system('mkdir MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Outputfiles')
os.system('mkdir MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Submitfiles')

cutoff=2*rcut

#Create array of temperature dependent lattice constants
f = open('ThermalData/'+material+'_'+crystaltype+'/LatticeConst','r')
f.readline()
for i in range(0,21):
    line=f.readline()
    data=line.split(' ')
    if i == 0:
        TempData=numpy.array([float(data[0]),float(data[1]),float(data[3])])
    else:
        TempData=numpy.vstack([TempData,[float(data[0]),float(data[1]),float(data[3])]])
f.close()






for i in range(0,1):



#Find Particles on free surface
    f = open(datafile,'r')
    f.readline()
    line2=f.readline()
    data=line2.split(' ')
    N=int(data[0])
    f.readline()
    line4 = f.readline()
    xdata=line4.split(' ')
    Lx=float(xdata[1])-float(xdata[0])
    line5 = f.readline()
    ydata=line5.split(' ')
    Ly=float(ydata[1])-float(ydata[0])
    line6 = f.readline()
    zdata=line6.split(' ')
    Lz=float(zdata[1])-float(zdata[0])
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()

#Use far side of box as starting points
    ymin=float(zdata[1])
    ymax=float(zdata[0])

#Find min and max values of particles in y-direction
    for j in range(0,N):
    	line = f.readline()
    	data=line.split(' ')
    	if float(data[4])<ymin:
        	ymin=float(data[4])
    	if float(data[4])>ymax:
        	ymax=float(data[4])
    f.close()

#Calculate regions for vy=0,fy=0 (Necessary to keep box stationary)
    yhi=ymax-rcut
    ylo=ymin+rcut



#For temperatures to be analyzed, perform linear interpolation to find a and c
    for j in range(0,len(Temps)):
    	temp=Temps[j]
    	for i in range(0,21):
        	if temp<TempData[i][0]:
            		T1=TempData[i-1][0]
            		T2=TempData[i][0]
            		a1=TempData[i-1][1]
            		a2=TempData[i][1]
            		c1=TempData[i-1][2]
            		c2=TempData[i][2]
            		break
    	r=(temp-T1)/(T2-T1)
    	a=(1-r)*a1+r*a2
    	c=(1-r)*c1+r*c2

	ori1x=miller1x[0]*numpy.matrix([a, 0, 0])+miller1x[1]*numpy.matrix([-a/2, a*math.sqrt(3)/2, 0])+miller1x[2]*numpy.matrix([-a/2, -a*math.sqrt(3)/2, 0])+miller1x[3]*numpy.matrix([0, 0, c])
	ori1y=miller1y[0]*numpy.matrix([a, 0, 0])+miller1y[1]*numpy.matrix([-a/2, a*math.sqrt(3)/2, 0])+miller1y[2]*numpy.matrix([-a/2, -a*math.sqrt(3)/2, 0])+miller1y[3]*numpy.matrix([0, 0, c])
	ori2x=miller2x[0]*numpy.matrix([a, 0, 0])+miller2x[1]*numpy.matrix([-a/2, a*math.sqrt(3)/2, 0])+miller2x[2]*numpy.matrix([-a/2, -a*math.sqrt(3)/2, 0])+miller2x[3]*numpy.matrix([0, 0, c])
	ori2y=miller2y[0]*numpy.matrix([a, 0, 0])+miller2y[1]*numpy.matrix([-a/2, a*math.sqrt(3)/2, 0])+miller2y[2]*numpy.matrix([-a/2, -a*math.sqrt(3)/2, 0])+miller2y[3]*numpy.matrix([0, 0, c])
	#ori2y=ori1y	
	ori1z=numpy.cross(ori1x,ori1y)
	ori2z=numpy.cross(ori2x,ori2y)
	ori1x=ori1x.T
	ori1y=ori1y.T
	ori1z=ori1z.T
	ori2x=ori2x.T
	ori2y=ori2y.T
	ori2z=ori2z.T


	#Create Rotation Matrices
    	ori1x=ori1x/numpy.linalg.norm(ori1x)
    	ori1y=ori1y/numpy.linalg.norm(ori1y)
    	ori1z=ori1z/numpy.linalg.norm(ori1z)
    	ori2x=ori2x/numpy.linalg.norm(ori2x)
    	ori2y=ori2y/numpy.linalg.norm(ori2y)
    	ori2z=ori2z/numpy.linalg.norm(ori2z)

    	x=numpy.matrix([1,0,0])
    	y=numpy.matrix([0,1,0])
    	z=numpy.matrix([0,0,1])



    	R1=numpy.matrix([[float(numpy.dot(x,ori1x)),float(numpy.dot(y,ori1x)),float(numpy.dot(z,ori1x))],[float(numpy.dot(x,ori1y)),float(numpy.dot(y,ori1y)),float(numpy.dot(z,ori1y))],[float(numpy.dot(x,ori1z)),float(numpy.dot(y,ori1z)),float(numpy.dot(z,ori1z))]])
    	R2=numpy.matrix([[float(numpy.dot(x,ori2x)),float(numpy.dot(y,ori2x)),float(numpy.dot(z,ori2x))],[float(numpy.dot(x,ori2y)),float(numpy.dot(y,ori2y)),float(numpy.dot(z,ori2y))],[float(numpy.dot(x,ori2z)),float(numpy.dot(y,ori2z)),float(numpy.dot(z,ori2z))]])
	print(R1)
	print(R2)

    #Define basis vectors
    	a1=numpy.matrix([[a], [0],[0]])
    	a2=numpy.matrix([[-a/2.], [a*math.sqrt(3)/2.],[0]])
    	a3=numpy.matrix([[0], [a*math.sqrt(3)/3.],[c/2.0]])
        a4=numpy.matrix([[0], [a*math.sqrt(3)/3.],[-c/2.]])

    #Calculate Rotated Basis
    	b1_1=R1*a1
    	b2_1=R1*a2
    	b3_1=R1*a3
        b4_1=R1*a4
        b1_2=R2*a1
    	b2_2=R2*a2
    	b3_2=R2*a3
        b4_2=R2*a4

    #Write Orientation file (ECO)
    	f=open('MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/InputScripts/eco_'+repr(temp)+'K.ori','w')
    	f.write(str(b1_1[0]).strip('[]')+' '+str(b1_1[1]).strip('[]')+' '+str(b1_1[2]).strip('[]')+'\n')
    	f.write(str(b2_1[0]).strip('[]')+' '+str(b2_1[1]).strip('[]')+' '+str(b2_1[2]).strip('[]')+'\n')
    	f.write(str(b3_1[0]).strip('[]')+' '+str(b3_1[1]).strip('[]')+' '+str(b3_1[2]).strip('[]')+'\n')
    	f.write(str(b1_2[0]).strip('[]')+' '+str(b1_2[1]).strip('[]')+' '+str(b1_2[2]).strip('[]')+'\n')
    	f.write(str(b2_2[0]).strip('[]')+' '+str(b2_2[1]).strip('[]')+' '+str(b2_2[2]).strip('[]')+'\n')
    	f.write(str(b3_2[0]).strip('[]')+' '+str(b3_2[1]).strip('[]')+' '+str(b3_2[2]).strip('[]')+'\n')
    	f.write(str(b4_1[0]).strip('[]')+' '+str(b4_1[1]).strip('[]')+' '+str(b4_1[2]).strip('[]')+'\n')
        f.write(str(b4_2[0]).strip('[]')+' '+str(b4_2[1]).strip('[]')+' '+str(b4_2[2]).strip('[]')+'\n')
        f.close()


#Define additional vectors for orient
        a5=numpy.matrix([[-a], [0],[0]])
        a6=numpy.matrix([[a/2.], [a*math.sqrt(3)/2.],[0]])
        a7=numpy.matrix([[a/2.], [-a*math.sqrt(3)/2.],[0]])
        a8=numpy.matrix([[-a/2.], [-a*math.sqrt(3)/2.],[0]])
        a9=numpy.matrix([[-a/2.], [a*math.sqrt(3)/6.],[-c/2.]])
        a10=numpy.matrix([[a/2.], [a*math.sqrt(3)/6.],[c/2.]])
        a11=numpy.matrix([[-a/2.], [a*math.sqrt(3)/6.],[c/2.]])
        a12=numpy.matrix([[0], [-a*math.sqrt(3)/3.],[-c/2.0]])
        
        b5_1=R1*a5
        b6_1=R1*a6
        b7_1=R1*a7
        b8_1=R1*a8
        b9_1=R1*a9
        b10_1=R1*a10
        b11_1=R1*a11
        b12_1=R1*a12
        b5_2=R2*a5
        b6_2=R2*a6
        b7_2=R2*a7
        b8_2=R2*a8
        b9_2=R2*a9
        b10_2=R2*a10
        b11_2=R2*a11
        b12_2=R2*a12

#Write Orientation file (orient)
        f=open('orient1.file','w')
        f.write(str(b1_1[0]).strip('[]')+' '+str(b1_1[1]).strip('[]')+' '+str(b1_1[2]).strip('[]')+'\n')
        f.write(str(b2_1[0]).strip('[]')+' '+str(b2_1[1]).strip('[]')+' '+str(b2_1[2]).strip('[]')+'\n')
        f.write(str(b3_1[0]).strip('[]')+' '+str(b3_1[1]).strip('[]')+' '+str(b3_1[2]).strip('[]')+'\n')
        f.write(str(b4_1[0]).strip('[]')+' '+str(b4_1[1]).strip('[]')+' '+str(b4_1[2]).strip('[]')+'\n')
        f.write(str(b5_1[0]).strip('[]')+' '+str(b5_1[1]).strip('[]')+' '+str(b5_1[2]).strip('[]')+'\n')
        f.write(str(b6_1[0]).strip('[]')+' '+str(b6_1[1]).strip('[]')+' '+str(b6_1[2]).strip('[]')+'\n')
        f.write(str(b7_1[0]).strip('[]')+' '+str(b7_1[1]).strip('[]')+' '+str(b7_1[2]).strip('[]')+'\n')
        f.write(str(b8_1[0]).strip('[]')+' '+str(b8_1[1]).strip('[]')+' '+str(b8_1[2]).strip('[]')+'\n')
        f.write(str(b9_1[0]).strip('[]')+' '+str(b9_1[1]).strip('[]')+' '+str(b9_1[2]).strip('[]')+'\n')
        f.write(str(b10_1[0]).strip('[]')+' '+str(b10_1[1]).strip('[]')+' '+str(b10_1[2]).strip('[]')+'\n')
        f.write(str(b11_1[0]).strip('[]')+' '+str(b11_1[1]).strip('[]')+' '+str(b11_1[2]).strip('[]')+'\n')
        f.write(str(b12_1[0]).strip('[]')+' '+str(b12_1[1]).strip('[]')+' '+str(b12_1[2]).strip('[]')+'\n')
        f.close()
        f=open('orient2.file','w')
        f.write(str(b1_2[0]).strip('[]')+' '+str(b1_2[1]).strip('[]')+' '+str(b1_2[2]).strip('[]')+'\n')
        f.write(str(b2_2[0]).strip('[]')+' '+str(b2_2[1]).strip('[]')+' '+str(b2_2[2]).strip('[]')+'\n')
        f.write(str(b3_2[0]).strip('[]')+' '+str(b3_2[1]).strip('[]')+' '+str(b3_2[2]).strip('[]')+'\n')
        f.write(str(b4_2[0]).strip('[]')+' '+str(b4_2[1]).strip('[]')+' '+str(b4_2[2]).strip('[]')+'\n')
        f.write(str(b5_2[0]).strip('[]')+' '+str(b5_2[1]).strip('[]')+' '+str(b5_2[2]).strip('[]')+'\n')
        f.write(str(b6_2[0]).strip('[]')+' '+str(b6_2[1]).strip('[]')+' '+str(b6_2[2]).strip('[]')+'\n')
        f.write(str(b7_2[0]).strip('[]')+' '+str(b7_2[1]).strip('[]')+' '+str(b7_2[2]).strip('[]')+'\n')
        f.write(str(b8_2[0]).strip('[]')+' '+str(b8_2[1]).strip('[]')+' '+str(b8_2[2]).strip('[]')+'\n')
        f.write(str(b9_2[0]).strip('[]')+' '+str(b9_2[1]).strip('[]')+' '+str(b9_2[2]).strip('[]')+'\n')
        f.write(str(b10_2[0]).strip('[]')+' '+str(b10_2[1]).strip('[]')+' '+str(b10_2[2]).strip('[]')+'\n')
        f.write(str(b11_2[0]).strip('[]')+' '+str(b11_2[1]).strip('[]')+' '+str(b11_2[2]).strip('[]')+'\n')
        f.write(str(b12_2[0]).strip('[]')+' '+str(b12_2[1]).strip('[]')+' '+str(b12_2[2]).strip('[]')+'\n')
        f.close()



    #Calculate necessary size of simulation

        repx=int(math.ceil(50/Lx))
        repz=int(math.ceil(50/Lz))


        	#Write ECO 0K Input file
        f=open('MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/InputScripts/eco_0K.in','w')
        f.write('#ECO 0K input \n')
        f.write('units metal \n')
    	f.write('read_data '+datafile+' \n')
    	f.write('replicate '+repr(repx)+' 1 '+repr(repz)+' \n')
    	f.write('reset_timestep 0 \n')
    	f.write('pair_style '+pair_style+'\n')
    	f.write('pair_coeff '+pair_coeff+' \n')
    	f.write('comm_modify cutoff '+repr(cutoff)+'\n')
    	f.write('thermo 1000 \n')
        f.write('fix eco all eco/force/hcp 1 '+repr(eta)+' '+repr(rcut)+' MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/InputScripts/eco_'+repr(temp)+'K.ori \n')
        f.write('thermo_style custom step temp etotal press f_eco \n')
        f.write('dump save all custom 1000 MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Dumps/eco_0K.* id x y z f_eco[1] f_eco[2] f_eco[3] f_eco[4] f_eco[5] fx fy fz vx vy vz \n')
        f.write('dump_modify save sort id \n')
        f.write('run 0')
        f.close()

        os.system('mpiexec -np 12 ~/Desktop/lammps-30Jul16/src/lmp_mpi <MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/InputScripts/eco_0K.in>MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Outputfiles/eco_0K.out')

        f=open('MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/InputScripts/orient_0K.in','w')
        f.write('#Orient input '+repr(temp)+'K \n')
        f.write('units metal \n')
    	f.write('read_data '+datafile+' \n')
    	f.write('replicate '+repr(repx)+' 1 '+repr(repz)+' \n')
        f.write('reset_timestep 0 \n')
        f.write('pair_style '+pair_style+'\n')
        f.write('pair_coeff '+pair_coeff+' \n')
        f.write('comm_modify cutoff '+repr(cutoff)+'\n')
        f.write('thermo 100 \n')
        f.write('fix orient all orient/'+crystaltype+' 0 1 '+repr(a)+' 1 '+repr(cutlo)+' '+repr(cuthi)+' orient1.file orient2.file \n')
        f.write('thermo_style custom step temp etotal press f_orient \n')
        f.write('dump save all custom 1000 MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Dumps/orient_0K.* id x y z f_orient[1] f_orient[2] fx fy fz vx vy vz \n')
        f.write('dump_modify save sort id \n')
        f.write('run 0')
        f.close()

        #os.system('mpiexec -np 12 ~/Desktop/lammps-30Jul16/src/lmp_mpi <MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/InputScripts/orient_0K.in>MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Outputfiles/orient_0K.out')

    #Write restart file at temperature to be analyzed
    	f=open('MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/InputScripts/temp_equil_'+repr(temp)+'.in','w')
    	f.write('#Thermal Equilibration Input File \n')
    	f.write('dimension 3 \n')
    	f.write('boundary p p p \n')
    	f.write('units metal \n')
    	f.write('atom_style atomic \n')
    	f.write('lattice fcc '+repr(a)+' origin .01 .01 .01 \n')
    	f.write('region box block 0 1 0 1 0 1 units lattice \n')
    	f.write('create_box 1 box \n')
    	f.write('read_data '+datafile+' \n')
    	f.write('replicate '+repr(repx)+' 1 '+repr(repz)+' \n')
    	f.write('reset_timestep 0 \n')
    	f.write('pair_style '+pair_style+'\n')
    	f.write('pair_coeff '+pair_coeff+' \n')
    	f.write('comm_modify cutoff '+repr(cutoff)+'\n')
    	f.write('comm_style tiled \n')
    	f.write('balance 0.9 rcb \n')
    	f.write('timestep 0.001 \n')
    	f.write('thermo 1000 \n')
    	f.write('velocity all create 5 1 \n')
    	f.write('fix integrator all npt iso 0 0 1 temp 1 '+repr(temp)+' .1 \n')
    	f.write('thermo_style custom step temp etotal press \n')
    	f.write('run 20000 \n')
    	f.write('unfix integrator \n')
    	f.write('fix integrator all npt iso 0 0 1 temp '+repr(temp)+' '+repr(temp)+' .1 \n')
    	f.write('run 5000 \n')
    	f.write('write_restart MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Dumps/'+repr(temp)+'K.restart \n')
    	f.write('run 0')
    	f.close()


    	#write_submit_topaz('MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Submitfiles/submit_topaz_restartmaker'+repr(temp)+'K.bash',material+'_'+repr(temp)+'K_restart',72,'00:59:00','debug','/p/work1/mcg84/DrivingForce/test','MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/InputScripts/temp_equil_'+repr(temp)+'.in','MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Outputfiles/temp_equil_'+repr(temp)+'.out')
    	#os.system('qsub MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Submitfiles/submit_topaz_restartmaker'+repr(temp)+'K.bash')

    	for k in range(0,len(u0)):
        	u=format(u0[k], '.5f')
        	#Write ECO Input file
        	f=open('MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/InputScripts/eco_'+repr(temp)+'K_'+u+'eV.in','w')
        	f.write('#ECO input '+repr(temp)+'K_'+u+'eV \n')
        	f.write('read_restart MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Dumps/'+repr(temp)+'K.restart \n')
        	f.write('reset_timestep 0 \n')
        	f.write('pair_style '+pair_style+'\n')
    		f.write('pair_coeff '+pair_coeff+' \n')
            f.write('comm_modify cutoff '+repr(cutoff)+'\n')
    		f.write('thermo 100 \n')
            f.write('fix integrator all nvt temp '+repr(temp)+' '+repr(temp)+' .1 \n')
        	f.write('region ends block INF INF '+repr(ylo)+' '+repr(yhi)+' INF INF units box side out \n')
        	f.write('group end region ends \n')
        	f.write('fix eco all eco/force/hcp '+u+' '+repr(eta)+' '+repr(rcut)+' MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/InputScripts/eco_'+repr(temp)+'K.ori \n')
        	f.write('fix edges end setforce NULL 0.0 NULL \n')
        	f.write('velocity end set NULL 0.0 NULL \n')
        	f.write('comm_style tiled \n')
        	f.write('balance 0.9 rcb \n')
        	f.write('thermo_style custom step temp etotal press f_eco \n')
        	f.write('dump save all custom 1000 MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Dumps/eco_'+repr(temp)+'K_'+u+'eV.* id x y z f_eco[1] f_eco[2] f_eco[3] f_eco[4] f_eco[5] fx fy fz vx vy vz \n')
        	f.write('dump_modify save sort id \n')
        	f.write('run '+repr(timesteps))
        	f.close()
        	#write_submit_topaz('MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Submitfiles/submit_topaz_'+repr(temp)+'K_'+u+'eV.bash',material+'_'+repr(temp)+'K_'+u+'eV',36,'00:50:00','debug','/p/work1/mcg84/DrivingForce/test','MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/InputScripts/eco_'+repr(temp)+'K_'+u+'eV.in','MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Outputfiles/'+repr(temp)+'K_'+u+'eV.out')

        	#while not os.path.exists('MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Dumps/'+repr(temp)+'K.restart'):
            		#time.sleep(10)
            		#print('Waiting for restart file')
        	#os.system('qsub MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Submitfiles/submit_topaz_'+repr(temp)+'K_'+u+'eV.bash')



        #Write input file (orient)
  	        f=open('MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/InputScripts/orient_'+repr(temp)+'K_'+u+'eV.in','w')
 	        f.write('#Orient input '+repr(temp)+'K_'+u+'eV \n')
  	        f.write('read_restart MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Dumps/'+repr(temp)+'K.restart \n')
 	        f.write('reset_timestep 0 \n')
 	        f.write('pair_style '+pair_style+'\n')
 	        f.write('pair_coeff '+pair_coeff+' \n')
 	        f.write('comm_modify cutoff '+repr(cutoff)+'\n')
 	        f.write('thermo 100 \n')
 	        f.write('fix integrator all nvt temp '+repr(temp)+' '+repr(temp)+' .1 \n')
        	f.write('region ends block INF INF '+repr(ylo)+' '+repr(yhi)+' INF INF units box side out \n')
        	f.write('group end region ends \n')
            f.write('fix orient all orient/'+crystaltype+' 0 1 '+repr(a)+' '+u+' '+repr(cutlo)+' '+repr(cuthi)+' orient1.file orient2.file \n')
        	f.write('fix edges end setforce NULL 0.0 NULL \n')
        	f.write('velocity end set NULL 0.0 NULL \n')
        	f.write('balance 1.1 shift y 5 1.02 \n')
        	f.write('thermo_style custom step temp etotal press f_orient \n')
        	f.write('dump save all custom 1000 MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Dumps/orient_'+repr(temp)+'K_'+u+'eV.* id x y z f_orient[1] f_orient[2] fx fy fz vx vy vz \n')
        	f.write('dump_modify save sort id \n')
        	f.write('run '+repr(timesteps))
         	f.close()

        	#write_submit_excalibur('MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Submitfiles/submit_excal_orient_'+repr(temp)+'K_'+u+'eV.bash',material+'_o_'+repr(temp)+'K_'+u,64,'07:20:00','standard','/work/mcg84/test','MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/InputScripts/orient_'+repr(temp)+'K_'+u+'eV.in','MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Outputfiles/'+repr(temp)+'K_'+u+'eV_orient.out')
        	#os.system('qsub MobilityData/'+material+'_'+crystaltype+'/'+Orientation+'/Submitfiles/submit_excal_orient_'+repr(temp)+'K_'+u+'eV.bash')



