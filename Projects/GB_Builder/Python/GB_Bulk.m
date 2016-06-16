function GB_Bulk(GBorientations,Write_Style,AtomStyle,Archive,Dir_Base,Verbose)

Lattice1=GBorientations(1).Info.Lattice1;
Basis1=GBorientations(1).Info.Basis1;
Species1=GBorientations(1).Info.Species1;
Masses1=GBorientations(1).Info.Masses1;
Axis1=GBorientations(1).Info.Axis1;
Direction1=GBorientations(1).Info.Direction1;
MaxMiller1=GBorientations(1).Info.MaxMiller1;
Dis_Tol1=GBorientations(1).Info.Dis_Tol1;
Deg_Tol1=GBorientations(1).Info.Deg_Tol1;
Symmetry1=GBorientations(1).Info.Symmetry1;
Sym_Tol1=GBorientations(1).Info.Sym_Tol1;

Lattice2=GBorientations(1).Info.Lattice2;
Basis2=GBorientations(1).Info.Basis2;
Species2=GBorientations(1).Info.Species2;
Masses2=GBorientations(1).Info.Masses2;
Axis2=GBorientations(1).Info.Axis2;
Direction2=GBorientations(1).Info.Direction2;
MaxMiller2=GBorientations(1).Info.MaxMiller2;
Dis_Tol2=GBorientations(1).Info.Dis_Tol2;
Deg_Tol2=GBorientations(1).Info.Deg_Tol2;
Symmetry2=GBorientations(1).Info.Symmetry2;
Sym_Tol2=GBorientations(1).Info.Sym_Tol2;


%% Find minimum bulk volume
MinVol1=1e10;
MinVol2=1e10;
for i = 1:length(GBorientations)
  if GBorientations(i).Lat(2).volume < MinVol1
    MinVol1=GBorientations(i).Lat(2).volume;
    Orient1=GBorientations(i).Lat(2);
    Lx1=GBorientations(i).Lat(2).norms(1);
    Ly1=GBorientations(i).Lat(2).norms(2);
    Lz1=GBorientations(i).Lat(2).norms(3);
  end
  if GBorientations(i).Lat(2).volume < MinVol2
    MinVol2=GBorientations(i).Lat(2).volume;
    Orient2=GBorientations(i).Lat(2);
    Lx2=GBorientations(i).Lat(2).norms(1);
    Ly2=GBorientations(i).Lat(2).norms(2);
    Lz2=GBorientations(i).Lat(2).norms(3);
  end  
 
end
555555555
%% Determine if same bulk structure
ComLat=norm(Lattice1-Lattice2)<10e-9;
ComBasL=length(Basis1)==length(Basis2);
ComSpecL=length(Species1)==length(Species2);
ComMassL=length(Masses1)==length(Masses2);
ComBasis=0; ComSpec=0; ComMass=0;
if ComBasL
  ComBasis=norm(Basis1-Basis2)<10e-9;
end
if ComSpecL
   ComSpec=strcmp(Species1,Species2);
end
if ComMassL
  ComMass=norm(Masses1-Masses2)<10e-9;
end


if ComLat & ComBasis & ComSpec & ComMass
  % Lattice 1 and Lattice 2 are the same material
  if MinVol1<=MinVol2
  AtomData=GB_FillRegion(Lattice1,Basis1,Orient1,Lx1,Ly1,Lz1);
  Corners=[0 Lx1 0 Ly1 0 Lz1];
  else
  AtomData=GB_FillRegion(Lattice2,Basis2,Orient2,Lx2,Ly2,Lz2);
  Corners=[0 Lx2 0 Ly1 0 Lz2];  
  end

  
  AtomData=GB_WrapPBC(AtomData,Corners(2:2:6),0.5);
  nAtoms=length(AtomData);
  AtomIDs=[1:nAtoms]';
  Header.Head1={'# Bulk: 1'};
  Header.Head2={''};
  Header.Head3={'# Bulk: 1'};
  GB_WriteFiles('Bulk',Write_Style,AtomData,Species1,Masses1,Corners,Header,AtomStyle,Archive,Dir_Base,Verbose)
  
else

  AtomData1=GB_FillRegion(Lattice1,Basis1,Orient1,Lx1,Ly1,Lz1);
  Corners1=[0 Lx1 0 Ly1 0 Lz1];
  
  AtomData1=GB_WrapPBC(AtomData1,Corners1(2:2:6),0.5);
  nAtoms1=length(AtomData1);
  AtomIDs1=[1:nAtoms1]';
  Header1.Head1={'# Bulk: 2\n# Lat 1'};
  Header1.Head2={''};
  Header1.Head3={'# Bulk: 2 | Lat 1'};
  GB_WriteFiles('Bulk1',Write_Style,AtomData1,Species1,Masses1,Corners1,Header1,AtomStyle,Archive,Dir_Base,Verbose)
  
  
  AtomData2=GB_FillRegion(Lattice2,Basis2,Orient2,Lx2,Ly2,Lz2);
  Corners2=[0 Lx2 0 Ly2 0 Lz2]; 

  AtomData2=GB_WrapPBC(AtomData2,Corners2(2:2:6),0.5);
  nAtoms2=length(AtomData2);
  AtomIDs2=[1:nAtoms2]';
  Header2.Head1={'# Bulk 2\n# Lat 2'};
  Header2.Head2={''};
  Header2.Head3={'# Bulk: 2 | Lat 1'};
  GB_WriteFiles('Bulk2',Write_Style,AtomData2,Species2,Masses2,Corners2,Header2,AtomStyle,Archive,Dir_Base,Verbose)
  
  
end

end

