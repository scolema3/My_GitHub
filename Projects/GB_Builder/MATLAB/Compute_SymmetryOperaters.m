function [ActiveSymmetry,NSymmetry] = Compute_SymmetryOperaters(Lattice,Basis,Symmetry,Orientation,Sym_Tol)

if exist('Symmetry')==0
  Symmetry=true;
end
if exist('Orientation')==0
  Orientation=eye(3);
end
if exist('Tol')==0
  Sym_Tol=0.001;
end

if Symmetry==true
  
  Rot=normr(Orientation);

  Rx=@(theta)[ 1           0           0
               0          cos(theta) -sin(theta)
               0          sin(theta)  cos(theta)];
  Ry=@(theta)[ cos(theta) 0           sin(theta)
               0          1           0
              -sin(theta) 0           cos(theta)];
  Rz=@(theta)[ cos(theta) -sin(theta)  0
               sin(theta)  cos(theta)  0
               0           0           1];

  % Symmetry operations to test
  e(:,:,1)=[1 0 0; 0 1 0; 0 0 1]; 
  e(:,:,2)=[1 0 0;0 -1 0;0 0 -1]; 
  e(:,:,3)=[1 0 0;0 0 -1;0 1 0]; 
  e(:,:,4)=[1 0 0;0 0 1; 0 -1 0];
  e(:,:,5)=[-1 0 0;0 1 0;0 0 -1]; 
  e(:,:,6)=[-1 0 0;0 -1 0;0 0 1]; 
  e(:,:,7)=[-1 0 0;0 0 -1;0 -1 0]; 
  e(:,:,8)=[-1 0 0;0 0 1;0 1 0];
  e(:,:,9)=[0 1 0;-1 0 0;0 0 1]; 
  e(:,:,10)=[0 1 0;0 0 -1;-1 0 0]; 
  e(:,:,11)=[0 1 0;1 0 0;0 0 -1]; 
  e(:,:,12)=[0 1 0;0 0 1;1 0 0];
  e(:,:,13)=[0 -1 0;1 0 0;0 0 1]; 
  e(:,:,14)=[0 -1 0;0 0 -1;1 0 0]; 
  e(:,:,15)=[0 -1 0;-1 0 0;0 0 -1]; 
  e(:,:,16)=[0 -1 0;0 0 1;-1 0 0];
  e(:,:,17)=[0 0 1;0 1 0;-1 0 0]; 
  e(:,:,18)=[0 0 1;1 0 0;0 1 0]; 
  e(:,:,19)=[0 0 1;0 -1 0;1 0 0]; 
  e(:,:,20)=[0 0 1;-1 0 0;0 -1 0];
  e(:,:,21)=[0 0 -1;0 1 0;1 0 0]; 
  e(:,:,22)=[0 0 -1;-1 0 0;0 1 0]; 
  e(:,:,23)=[0 0 -1;0 -1 0;-1 0 0]; 
  e(:,:,24)=[0 0 -1;1 0 0;0 -1 0];
  e(:,:,25)=Rx(2*pi/3);
  e(:,:,26)=Ry(2*pi/3);
  e(:,:,27)=Rz(2*pi/3);
  e(:,:,28)=Rx(-2*pi/3);
  e(:,:,29)=Ry(-2*pi/3);
  e(:,:,30)=Rz(-2*pi/3);
  e(:,:,31)=Rx(2*pi/6);
  e(:,:,32)=Ry(2*pi/6);
  e(:,:,33)=Rz(2*pi/6);
  e(:,:,34)=Rx(-2*pi/6);
  e(:,:,35)=Ry(-2*pi/6);
  e(:,:,36)=Rz(-2*pi/6);

  [~,~,nSym]=size(e);

  % Create expanded crytal for reference
  MaxN=2;
  n = 2*MaxN + 1;
  R = fullfact([n n n])-MaxN-1;
  nR = length(R);
  [nB,~]=size(Basis);   
  Expand=zeros(nR*nB,4);

  % Replicating each basis atom
  Tmp1=ones(nR,1);
  for m = 1:nB, 
    Expand((m-1)*nR+1:m*nR,:)=[Tmp1*Basis(m,1) Tmp1*Basis(m,2:4)+R];
  end
  Ref=[Expand(:,1) Expand(:,2:end)*Lattice*Rot];

  [nRef,~]=size(Ref);

  Ref2=zeros(nRef,nSym*4);

  % Create test by applying symmetry opperation to crystal
  for n=1:nSym
    Test(:,(n-1)*4+1:4*n)=[Basis(:,1) Basis(:,2:4)*Lattice*Rot*e(:,:,n)];
  end

  Ref2(:,1:4:nSym*4)=Ref(:,1)*ones(1,nSym);
  Ref2(:,2:4:nSym*4)=Ref(:,2)*ones(1,nSym);
  Ref2(:,3:4:nSym*4)=Ref(:,3)*ones(1,nSym);
  Ref2(:,4:4:nSym*4)=Ref(:,4)*ones(1,nSym);

  for i=1:nB
    SymBasis(i,:)=sum(reshape(sum(reshape((Ref2-ones(nRef,1)*Test(i,:)).^2',4,[])),nSym,[])'<=Sym_Tol);
  end
  ActSym=find(sum(SymBasis)==nB);


  ActiveSymmetry=e(:,:,ActSym);
  [~,~,NSymmetry]=size(ActiveSymmetry);

else
  ActiveSymmetry=eye(3);
  NSymmetry=1;
end