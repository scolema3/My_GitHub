%% Three examples for the grain boundary builder script
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

Rot1_Options.MaxMiller=5;

Rots=GB_GrainRot(Lattice,Basis,Species,Masses,Axis,Direction,Rot1_Options);


GBorientations=GB_Orientations(Rots);

Const_Options.Angle_Range=[36.8 36.9];
Const_Options.Verbose=false;
GB_Construct(GBorientations,Const_Options)


%% Example 2: Detailed Metal Test


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


%% Example 2: Detailed Ceramic Test
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



%% Example 2: Detailed Ceramic Test B12-CCC
clear; clc

Lattice=[
5.17399979       0.00000000       0.00000000 
2.10415506       4.72639990       0.00000000 
2.10415506       1.36688995       4.52477980];
       
Basis =[
1   0.927845765340347   0.409863168018095   0.432652214368531
1   0.614782337539415   0.273745664185089   0.620750207557062
1   0.290590631028928   0.597961514632117   0.620728725848714
1   0.426672211645590   0.911082182284684   0.432618488970447
1   0.813548866818212   0.796772449167830   0.318323048560286
1   0.949639603494063   0.608648287635083   0.631422152300096
1   0.813558577842273   0.295528615696991   0.819535969463089
1   0.426684675892762   0.409841707201193   0.933825796340410
1   0.312384942712063   0.796747767053407   0.819501404245130
1   0.625447534409605   0.932862947513067   0.631401422009531
1   0.614773641211052   0.597976490420801   0.296542143332588
1   0.625454563242427   0.608632342113871   0.955609088424590
2   0.120130543395555   0.103265833620477   0.126066311558410
2   0.230320331695361   0.213526980336710   0.236314416891624
2   0.009937088314977  -0.006998401457806   0.01581681389224  
];

      
Species={'B','C'};
Masses=[10.811 12.0107];  
Axis=[0 0 1];
Direction='x';
Rot_Options.MaxMiller=15;
Rot_Options.Dis_Tol=0.5;
Rot_Options.Deg_Tol=0.5;
Rot_Options.Symmetry=true;
Rot_Options.Sym_Tol=0.5;
Rot_Options.Verbose=true;
GrainRot=GB_GrainRot(Lattice,Basis,Species,Masses,Axis,Direction,Rot_Options);



Orient_Options.Strain_Tol=0.1;
Orient_Options.MaxArea=100000;
Orient_Options.Verbose=true;

GB=GB_Orientations(GrainRot,GrainRot,Orient_Options);

%%
clc
Const_Options.Angle_Range=[46.9969 46.9971];
Const_Options.AtomStyle='charge';

Const_Options.GBSort='Area';
Const_Options.nGBs=1;
Const_Options.WriteStyle={'dump', 'data'};
Const_Options.Overlap_tol=0;
Const_Options.NormSlabDim=[20 20];
Const_Options.Vacuum=20;
Const_Options.GBregion=10;
Const_Options.FullyPeriodic=true;
Const_Options.Suffix='';
Const_Options.Verbose=true;
Const_Options.Plane_Shift=[0 0];
Const_Options.Translate=[0 0 0 0 .5 0];

Const_Options.Stoich=true;
Const_Options.ForceStoich=false;

GB_Construct(GB,Const_Options)
