function GrainRot=GB_GrainRot2(Lattice,Basis,Species,Masses,Axis1,Direction1,Axis2,Direction2,MaxInt,Rot_Options)

% Mapping lattice onto a sets of orthogoal rotations (orientations)

%% Setup variables


%% Check if same rotation has already been calculated



if [strcmpi(Direction1,'x')]
    DirIndex1=1;
elseif [strcmpi(Direction1,'y')]
    DirIndex1=2;
elseif [strcmpi(Direction1,'z')] 
    DirIndex1=3;
else
  error(['Direction1: (x, y, or z) must be set.'])
end

if [strcmpi(Direction2,'x')]
    DirIndex2=1;
elseif [strcmpi(Direction2,'y')]
    DirIndex2=2;
elseif [strcmpi(Direction2,'z')] 
    DirIndex2=3;
else
  error(['Direction2: (x, y, or z) must be set.'])
end

DirIndex3=setdiff([1,2,3],[DirIndex1,DirIndex2]);


  OffPlane(1,1)=0;
  OffPlane(2,1)=0;
  OffPlane(3,1)=0;
 
  Vectors(DirIndex1,:)=reshape(Axis1,1,[]);        
  Vectors(DirIndex2,:)=reshape(Axis2,1,[]);
  
  % Determining 3rd axis in closest integer form
  A=cross([Lattice*Axis1']',[Lattice*Axis2']')/Lattice;
  [~,id1]=min(abs(A));
  [~,id2]=max(abs(A));
  id3=setdiff([1,2,3],[id1,id2]);
  
  if A(id1)==0
    tmpid1=id1;
    id1=id3;
    id3=tmpid1;
  end
  
  MinError=1;
  for TolExp=1:7
    Tol=1*10^-TolExp;
    [R1,R2]=rat(A(id1)/A(id2),Tol);
    err1=A(id1)/A(id2)-R1/R2;
    R3=round(A(id3)/A(id2)*R2);
    err2=round(A(id3)/A(id2)*R2)/R2-A(id3)/A(id2);
    err=err1^2+err2^2;
        if  abs(R2)<=MaxInt && err<MinError
            MinError=err;
            Axis3(id1)=R1;
            Axis3(id2)=R2;
            Axis3(id3)=R3;
        end

  end
  

  Vectors(DirIndex3,:)=reshape(Axis3,1,[]);
  Vectors=Vectors*Lattice;
  Norms(1,1)=norm(Vectors(1,:));
  Norms(2,1)=norm(Vectors(2,:));
  Norms(3,1)=norm(Vectors(3,:));


  GrainRot(1).vectors=Vectors;
  GrainRot(1).offplane=0;
  GrainRot(1).angle=0;
  GrainRot(1).norms=Norms;
  GrainRot(1).volume=Norms(1)*Norms(2)*Norms(3);
  GrainRot(1).Sym(:,:,1)=eye(3);

GrainRot(1).Info.Lattice=Lattice;
GrainRot(1).Info.Basis=Basis;
GrainRot(1).Info.Species=Species;
GrainRot(1).Info.Masses=Masses;
GrainRot(1).Info.Axis=Axis1;
GrainRot(1).Info.Direction=Direction1;
GrainRot(1).Info.MaxMiller=MaxInt;
GrainRot(1).Info.Dis_Tol=0;
GrainRot(1).Info.Deg_Tol=0;
GrainRot(1).Info.Symmetry=0;
GrainRot(1).Info.Sym_Tol=0;


end