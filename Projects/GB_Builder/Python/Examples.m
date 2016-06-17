%% Four examples for the grain boundary builder script that generates 
% initial structures in 1) LAMMPS datafile 2) LAMMPS dumpfile 3) CarFile
% formats.

addpath('L:\Coleman\Documents\MATLAB\')


%% Example 1: Simplest Metal Test

clear;clc;close all

Lattice=[1 0 0  
         0 1 0
         0 0 1]*4.05;

Basis=[1 0.0 0.0 0.0 
       1 0.5 0.5 0.0
       1 0.0 0.5 0.5
       1 0.5 0.0 0.5];        

Species={'Al'};
Masses=[26.9815];  
Axis=[0 0 1];
Direction='tilt';
Rots=GB_GrainRot(Lattice,Basis,Species,Masses,Axis,Direction);
GBorientations=GB_Orientations(Rots);
GB_Construct(GBorientations)


%% Example 2: More Detailed Metal Test

clear;clc; close all

%%%%%%%%%%%%
Lattice1=[1 0 0  
          0 1 0
          0 0 1]*4.05;

Basis1=[1 0.0 0.0 0.0 
        1 0.5 0.5 0.0
        1 0.0 0.5 0.5
        1 0.5 0.0 0.5];        

Species1={'Al'};
Masses1=[26.9815];  
Axis1=[0 0 1];
Direction1='x';
Rot1_Options.MaxMiller=3;
Rot1_Options.Dis_Tol=0.01;
Rot1_Options.Deg_Tol=0.01;
Rot1_Options.Symmetry=true;
Rot1_Options.Sym_Tol=0.01;
Rot1_Options.Verbose=false;

GrainRot1=GB_GrainRot(Lattice1,Basis1,Species1,Masses1,Axis1,Direction1,Rot1_Options);

%%%%%%%%%%%%
Lattice2=[1 0 0  
          0 1 0
          0 0 1]*4.05;

Basis2=[1 0.0 0.0 0.0 
        2 0.5 0.5 0.0
        1 0.0 0.5 0.5
        1 0.5 0.0 0.5];        

Species2={'Al1', 'Al2'};
Masses2=[26 23];  
Axis2=[0 0 1];
Direction2='x';

Rot2_Options.MaxMiller=3;
Rot2_Options.Dis_Tol=0.01;
Rot2_Options.Deg_Tol=0.01;
Rot2_Options.Symmetry=true;
Rot2_Options.Sym_Tol=0.01;
Rot2_Options.Verbose=false;

GrainRot2=GB_GrainRot(Lattice2,Basis2,Species2,Masses2,Axis2,Direction2,Rot2_Options);

%%%%%%%%%%%%
Orient_Options.Strain_Tol=0.1;
Orient_Options.MaxArea=100000;
Orient_Options.Verbose=false;

GB=GB_Orientations(GrainRot1,GrainRot2,Orient_Options);

%%%%%%%%%%%%
Const_Options.GBSort='Area';
Const_Options.nGBs=4;
Const_Options.WriteStyle={'dump', 'data'};
Const_Options.Overlap_tol=1;
Const_Options.NormSlabDim=20;
Const_Options.Vacuum=20;
Const_Options.GBregion=10;
Const_Options.FullyPeriodic=false;
Const_Options.Archive=false;
Const_Options.Suffix='_1_cut';
Const_Options.Verbose=false;
Const_Options.Plane_Shift=[0 0];
Const_Options.Translate=[0 0 0 0 0 0];

GB_Construct(GB,Const_Options)


%% Example 3: Detailed Ceramic Test
clear; clc

a=5.4205; b=5.4205; c=12.3066;
alpha=90.0*pi/180; beta=90.0*pi/180 ;  gamma=120*pi/180;


gamma_star=acos( (cos(alpha)*cos(beta)-cos(gamma)) / (sin(alpha)*sin(beta)) );
Lattice=[a*sin(beta) 0 a*cos(beta)
           -b*sin(alpha)*cos(gamma_star)  b*sin(alpha)*sin(gamma_star)  b*cos(alpha)
            0 0 c];
       
Basis =[2   0.024015*2   0.024067*2	  0.645517
        2   0.190680*2   0.357410*2   0.312178
        2   0.357435*2   0.190657*2   0.739377
        2   0.190776*2   0.357337*2   0.072714        
        2   0.357352*2   0.190744*2   0.978843
        2   0.024105*2   0.024004*2   0.406045
        1   0.349283*2   0.436655*2   0.050429
        1   0.111446*2   0.436597*2   0.050437
        1   0.111386*2   0.198763*2   0.050447
        1   0.468964*2   0.079030*2   0.136153
        1   0.134057*2   0.078995*2   0.136164
        1   0.469013*2   0.413948*2   0.136183
        1   0.412501*2   0.301030*2   0.244264
        1   0.246977*2   0.135501*2   0.244273
        1   0.412536*2   0.135452*2   0.244285
        1   0.031952*2   0.277952*2   0.330585
        1   0.270027*2   0.278011*2   0.330598
        1   0.270084*2   0.016088*2   0.330600
        1   0.182609*2   0.103321*2   0.383759
        1   0.444773*2   0.103263*2   0.383773
        1   0.444715*2   0.365434*2   0.383791
        1   0.302294*2   0.245689*2   0.469488
        1   0.467387*2   0.245656*2   0.469499
        1   0.302337*2   0.080604*2   0.469514
        1   0.245817*2   0.467670*2   0.577589
        1   0.080295*2   0.302139*2   0.577611
        1   0.245857*2   0.302099*2   0.577618
        1   0.365290*2   0.444606*2   0.663908
        1   0.103355*2   0.444660*2   0.663930
        1   0.103412*2   0.182739*2   0.663942
        1   0.278108*2   0.269918*2   0.717099
        1   0.015940*2   0.269977*2   0.717107
        1   0.278050*2   0.032087*2   0.717116
        1   0.135642*2   0.412359*2   0.802822
        1   0.300731*2   0.412318*2   0.802823
        1   0.135686*2   0.247272*2   0.802849
        1   0.079167*2   0.134348*2   0.910926
        1   0.413647*2   0.468823*2   0.910937
        1   0.079203*2   0.468774*2   0.910953
        1   0.198625*2   0.111285*2   0.997253
        1   0.436697*2   0.111342*2   0.997262
        1   0.436754*2   0.349417*2   0.997271 ];
      
Species={'B','O'};
Masses=[10.811 15.994];  
Axis=[0 0 1];
Direction='x';
Rot_Options.MaxMiller=5;
Rot_Options.Dis_Tol=0.1;
Rot_Options.Deg_Tol=0.1;
Rot_Options.Symmetry=true;
Rot_Options.Sym_Tol=0.5;
Rot_Options.Verbose=true;

GrainRot=GB_GrainRot(Lattice,Basis,Species,Masses,Axis,Direction,Rot_Options);

Orient_Options.Strain_Tol=0.1;
Orient_Options.MaxArea=100000;
Orient_Options.Verbose=true;

GB=GB_Orientations(GrainRot,GrainRot,Orient_Options);


Const_Options.GBSort='Area';
Const_Options.nGBs=1;
Const_Options.WriteStyle={'dump', 'data'};
Const_Options.Overlap_tol=1;
Const_Options.NormSlabDim=[20 20];
Const_Options.Vacuum=20;
Const_Options.GBregion=10;
Const_Options.FullyPeriodic=false;
Const_Options.Suffix='_1_cut';
Const_Options.Verbose=true;
Const_Options.Plane_Shift=[0 0];
Const_Options.Translate=[1 0 0 0 0 0];


Const_Options.Stoich=true;
Const_Options.ForceStoich=true;

GB_Construct(GB,Const_Options)


%% Example 4: Alumina test

clear; clc

Lattice=[4.758 0 0
         0 8.241097742412718 0
         0 0 12.991];
       
Basis =[
        1	0	0.33333	0.0187
        1	0	0.33333	0.3145
        1	0	0.33333	0.5187
        1	0	0.33333	0.8145
        1	0	0.66666	0.6854
        1	0	0.66666	0.9812
        1	0	0.66666	0.1854
        1	0	0.66666	0.4812
        1	0.5	0.83333	0.0187
        1	0.5	0.83333	0.3145
        1	0.5	0.83333	0.5187
        1	0.5	0.83333	0.8145
        1	0	0	0.3521
        1	0	0	0.6479
        1	0	0	0.8521
        1	0	0	0.1479
        1	0.5	0.16666	0.6854
        1	0.5	0.16666	0.9812
        1	0.5	0.16666	0.1854
        1	0.5	0.16666	0.4812
        1	0.5	0.5	0.3521
        1	0.5	0.5	0.6479
        1	0.5	0.5	0.8521
        1	0.5	0.5	0.1479
        2	0.19355	0.5	0.25
        2	0.15325	0.4866	0.9166
        2	0.15325	0.1801	0.9166
        2	0.3065	0.33335	0.4166
        2	0.15325	0.8199	0.5833
        2	0.15325	0.5134	0.5833
        2	0.3065	0.66665	0.0833
        2	0.1533	0.84675	0.25
        2	0.3468	0.65325	0.75
        2	0.65325	0.9866	0.9166
        2	0.65325	0.6801	0.9166
        2	0.1935	0.83335	0.9166
        2	0.34675	0.6801	0.4166
        2	0.34675	0.9866	0.4166
        2	0.8065	0.83335	0.4166
        2	0.1533	0.15325	0.25
        2	0.69355	0	0.25
        2	0.3468	0.34675	0.75
        2	0.30655	0	0.75
        2	0.65325	0.3199	0.5833
        2	0.65325	0.0134	0.5833
        2	0.1935	0.16665	0.5833
        2	0.34675	0.0134	0.0833
        2	0.34675	0.3199	0.0833
        2	0.8065	0.16665	0.0833
        2	0.6533	0.65325	0.25
        2	0.6533	0.34675	0.25
        2	0.8468	0.84675	0.75
        2	0.8468	0.15325	0.75
        2	0.80655	0.5	0.75
        2	0.6935	0.33335	0.9166
        2	0.84675	0.1801	0.4166
        2	0.84675	0.4866	0.4166
        2	0.6935	0.66665	0.5833
        2	0.84675	0.5134	0.0833
        2	0.84675	0.8199	0.0833];
    
    
Species={'Al','O'};
Masses=[26.9815 15.994];  
Axis=[0 0 1];
Direction='x';
Rot_Options.MaxMiller=5;
Rot_Options.Dis_Tol=0.1;
Rot_Options.Deg_Tol=0.1;
Rot_Options.Symmetry=true;
Rot_Options.Sym_Tol=1.5;
Rot_Options.Verbose=true;

GrainRot=GB_GrainRot(Lattice,Basis,Species,Masses,Axis,Direction,Rot_Options);

Orient_Options.Strain_Tol=0.1;
Orient_Options.MaxArea=100000;
Orient_Options.Verbose=true;

GB=GB_Orientations(GrainRot,GrainRot,Orient_Options);


Const_Options.GBSort='Area';
Const_Options.nGBs=5;
Const_Options.WriteStyle={'dump', 'data'};
Const_Options.Overlap_tol=1;
Const_Options.NormSlabDim=[20 20];
Const_Options.Vacuum=20;
Const_Options.GBregion=10;
Const_Options.FullyPeriodic=false;
Const_Options.Suffix='_1_cut';
Const_Options.Verbose=true;
Const_Options.Plane_Shift=[0 0];
Const_Options.Translate=[1 0 0 0 0 0];


Const_Options.Stoich=true;
Const_Options.ForceStoich=true;

GB_Construct(GB,Const_Options)
    