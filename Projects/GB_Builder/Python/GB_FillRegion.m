function [AtomData]=GB_FillRegion(Lattice,Basis,Orientation,Lx,Ly,Lz,TypeOffset,Shift) %accepts shift but never receives it?

% Filling a volume region with atoms


%% Setup variables
w = whos;
for a = 1:length(w)
vars.(w(a).name) = eval(w(a).name);
end
vars_pass=GB_Vars(vars,'GB_FillRegion');

varNames=fieldnames(vars_pass);
for  i=1:length(varNames)
    eval([varNames{i} '=' 'vars_pass.(varNames{i});']);
end

%% Simulation dimensions in periodic directions - defining the 
%  corners of each lattice.

LatCorners=[0 0   0  Lx Lx Lx 0  Lx
            0 0   Ly 0  Ly 0  Ly Ly
            0 Lz  0  0  0  Lz Lz Lz]';

% Rotation matrix for each lattice.
Rot=normr(Orientation.vectors)';
Mins=min(LatCorners);
Maxs=max(LatCorners);

% Finding the basis coordinates from the simulation corners. 
% Which determine the bounds needed for translating the Basis.
BoundsLat=[];
for i = 1:8 
  BoundsLat=[BoundsLat;round(LatCorners(i,:)/(Lattice*Rot))];
end

BoundsLat=[min(BoundsLat)' max(BoundsLat)'];

clear A1 LatCorners i

%% Create the extended basis in order be more computationally
%  efficent.  NOTE: Block, can be optimized.

Block=25;
R = fullfact([Block Block Block])-1;
nR = length(R);
[nB,~]=size(Basis);   
Expand=zeros(nR*nB,4);

%% Replicating each basis atom
Tmp=ones(nR,1);
for m = 1:nB, 
  Expand((m-1)*nR+1:m*nR,:)=[Tmp*Basis(m,1)+TypeOffset ...
                             Tmp*Basis(m,2:4)+R];
end
clear R m Tmp
  
buff=0.1;
AtomData=[];
count=0;

%% Construct Lattice
for x=floor(BoundsLat(1,1)/Block):ceil(BoundsLat(1,2)/Block)
  for y=floor(BoundsLat(2,1)/Block):ceil(BoundsLat(2,2)/Block)
    for z=floor(BoundsLat(3,1)/Block):ceil(BoundsLat(3,2)/Block)
      % Translating the ExpandedBasis in 3D
      Shift=Expand(:,2:4)+...
            ones(nB*nR,1)*[x*(Block) y*(Block) z*(Block)];

      % Transform basis into the expanded lattices for each grain
      ExpandLattice=[Expand(:,1) Shift*Lattice*Rot];

      % Check boundary conditions
      Select=ExpandLattice(:,2)>=-buff+Mins(1) & ...
          ExpandLattice(:,2)<=buff+Maxs(1) & ...
          ExpandLattice(:,3)>=-buff+Mins(2) & ...
          ExpandLattice(:,3)<=buff+Maxs(2) & ...
          ExpandLattice(:,4)>=-buff+Mins(3) & ...
          ExpandLattice(:,4)<=buff+Maxs(3);
      nSelect=sum(Select);
      if nSelect==0
      count=count+1;
      end
      AtomData=[AtomData; ExpandLattice(Select,:)];
      clear Select nSelect Shift ExpandLattice
    end
  end
end

clear id x y z nR nB count Orientation Expand buff Block Rot
clear BoundsLat
end
