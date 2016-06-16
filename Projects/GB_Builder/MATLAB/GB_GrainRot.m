function GrainRot=GB_GrainRot(Lattice,Basis,Species,Masses,Axis,Direction,Rot_Options)

% Mapping lattice onto a sets of orthogoal rotations (orientations)

%% Setup variables

% Read optional paramaters stored in Options structure
if exist('Rot_Options')==1
  Opts=fieldnames(Rot_Options);
  for  i=1:length(Opts)
    eval([Opts{i} '=' 'Rot_Options.(Opts{i});']);
  end
end

w = whos;
for a = 1:length(w)
vars.(w(a).name) = eval(w(a).name);
end
vars_pass=GB_Vars(vars,'GB_GrainRot');

varNames=fieldnames(vars_pass);
for  i=1:length(varNames)
    eval([varNames{i} '=' 'vars_pass.(varNames{i});']);
end

%% Check if same rotation has already been calculated

Compute=true;

if exist(['Rot_' Species{:} '_' num2str(Axis(1)) '_' ...
          num2str(Axis(2)) '_' num2str(Axis(3)) '_' ...
          Direction '_' num2str(MaxMiller) '.mat'],'file')==2
  % Save original variables
  varNames=fieldnames(vars_pass);
  for  i=1:length(varNames)
    eval([varNames{i} '1=' 'vars_pass.(varNames{i});'])
  end
  
  % Load saved varialbes
  load(['Rot_' Species{:} '_' num2str(Axis(1)) '_' ...
         num2str(Axis(2)) '_' num2str(Axis(3)) '_' ...
         Direction '_' num2str(MaxMiller) '.mat'])
  Info=fieldnames(GrainRot(1).Info);
  for  i=1:length(Info)
    eval([Info{i} '=' 'GrainRot(1).Info.(Info{i});']);
  end
  
  Compute=false;
  % Compare original and saved variables  
  for  i=1:length(varNames)
    if isequal(eval([varNames{i} '1']),eval(varNames{i}))~=1
          Compute=true;
          eval( [varNames{i} '=' [varNames{i} '1'] ';'] );
    end
  end
end

if Compute==true
clear GrainRot

%% Searching over all lattice directions upto +/- maxMiller
%  Checking that each are the shortest possible and unique

n = 2*MaxMiller + 1;
IndexSearch = fullfact([n n n])-MaxMiller-1;
LatticeIndSearch=IndexSearch*Lattice;

% Ordering according to vector magnitude
[~, ids]=sort(sum([LatticeIndSearch.^2]')');
LatIndSearchSort=LatticeIndSearch(ids,:);

count=1;  % Discard origin (starting at 2)
LatIndSearchKeep(count,:)=LatIndSearchSort(2,:);
LatIndSearchKeepNormal(count,:)=normr(LatIndSearchSort(2,:));

% Testing uniqueness of vectors - i.e, 111 vs {111,222,333,etc}
for i=3:length(LatIndSearchSort)
  [Test,~]=min(sum([ones(count,1)*...
                     normr(LatIndSearchSort(i,:))-...
                     LatIndSearchKeepNormal]'.^2));
  if Test>=1e-10  
    count=count+1;
    LatIndSearchKeep(count,:)=LatIndSearchSort(i,:);
    LatIndSearchKeepNormal(count,:)=normr(LatIndSearchSort(i,:));
  end
end

LatIndex=LatIndSearchKeep;


%% Find where the expanded lattice intersects with the plane
%  * note, the plane is centered at the origin with the vector 
%    normal provided by the user

AxisLat=Axis*Lattice;

dPlane=abs(LatIndex(:,1)*AxisLat(1)...
           +LatIndex(:,2)*AxisLat(2)...
           +LatIndex(:,3)*AxisLat(3))/norm(AxisLat);
Select=dPlane<=Dis_Tol;

if sum(Select)==0 
  error('No points found on plane. Try increasing dis_tol.')
end

PlaneData=[LatIndex(Select,:) dPlane(Select)];
[~, index]=sortrows(abs(PlaneData)); 
PlaneData=PlaneData(index,:);    


if Verbose==true
fprintf('\nFound %d points within %0.2f from the {%d %d %d}'...
    ,sum(Select), max(dPlane(Select)),...
    Axis(1), Axis(2), Axis(3))
fprintf(' recentered on [0 0 0].\n')
end

%% Find the orthogonal vectors that lie on the chosen plane
%  * note, keeping the shortest unique orthogonal vector pairs.

dot_tol=cos(pi/180*(90-Deg_Tol)); % Convert deg_tol for dot prod
                          
% Computing dot prod of all vector pairs (v1,v2) lying in plane 
[L,~]=size(PlaneData);
A=[];
for i=1:L; A=[A; i*ones(L-i,1) [1+i:L]'];end
norm1=sqrt(sum(PlaneData(A(:,1),1:3).^2')');
norm2=sqrt(sum(PlaneData(A(:,2),1:3).^2')');
dotpro=dot(PlaneData(A(:,1),1:3),PlaneData(A(:,2),1:3),2)./...
       (norm1.*norm2);
OrthVecs=find(abs(dotpro)<=dot_tol);


%% Combining orthogonal vectors to create orthogoal rotations
%  * note, checks for uniquess in place taking
%    advantage of GB_SymOp.m 

if Verbose==true
    fprintf('Combining orthogonal vecotrs:')
end
  
V0=AxisLat;
norm0=norm(V0);
Unique=true;
Sym_Vec=[];
Progress=0.1;
count=0;
for K=1:length(OrthVecs)
  
  if Verbose==true
    if K/length(OrthVecs) >= Progress
      fprintf(' %3.0f%%',Progress*100)
      Progress=Progress+0.1;
    end
  end
    
  V1=PlaneData(A(OrthVecs(K),1),1:3);
  V2=PlaneData(A(OrthVecs(K),2),1:3);
  norm1=norm(V1);
  norm2=norm(V2);
  dotpro=dot(V1,V2)./(norm1*norm2);

  OffPlane(DirIdex(1),1)=0;
  OffPlane(DirIdex(2),1)=PlaneData(A(OrthVecs(K),1),4);
  OffPlane(DirIdex(3),1)=PlaneData(A(OrthVecs(K),2),4);
  
  Vectors(DirIdex(1),:)=V0;        
  Vectors(DirIdex(2),:)=V1;
  Vectors(DirIdex(3),:)=V2;
  
  Norms(DirIdex(1),1)=norm0;
  Norms(DirIdex(2),1)=norm1;
  Norms(DirIdex(3),1)=norm2;

  % Force right handedness
  Vectors_RH=Vectors;
  if det(Vectors)<=0
    Vectors_RH=RH*Vectors_RH;
    OffPlane=RH*OffPlane;
    Norms=RH*Norms;
  end

  if count==0               % Save initial vector pair data
    count=count+1;
    GrainRot(count).vectors=Vectors_RH;
    GrainRot(count).offplane=OffPlane;
    GrainRot(count).angle=acos(dotpro)*180/pi;
    GrainRot(count).norms=Norms;
    GrainRot(count).volume=Norms(1)*Norms(2)*Norms(3);
    
    % Identify symmetric orientations 
    [GrainRot(count).Sym,NSym] = GB_SymOp(Lattice,Basis,...
                                    Symmetry,Vectors_RH,Sym_Tol);
    for S=1:NSym
      Sym_Vec=[Sym_Vec; GrainRot(count).Sym(:,:,S)*Vectors_RH];
    end
    
  else
    
    % Check for uniqueness - using symmetry
    nSymTot=length(Sym_Vec)/3;
    Compare=zeros(nSymTot*3,3);
    Compare(1:3:end,:)=ones(nSymTot,1)*Vectors_RH(1,:);
    Compare(2:3:end,:)=ones(nSymTot,1)*Vectors_RH(2,:);
    Compare(3:3:end,:)=ones(nSymTot,1)*Vectors_RH(3,:);
    if sum(sum(reshape(sum((Sym_Vec-Compare).^2'),3,[]))<=1e-2)>0
      Unique=false;
    end
    
    % Save unique orthogonal rotations
    if Unique==true;
      count=count+1;
      GrainRot(count).vectors=Vectors_RH;
      GrainRot(count).offplane=OffPlane;
      GrainRot(count).angle=acos(dotpro)*180/pi;
      GrainRot(count).norms=Norms;
      GrainRot(count).volume=Norms(1)*Norms(2)*Norms(3);

      % Identify symmetric orientations 
      [GrainRot(count).Sym,NSym] = GB_SymOp(Lattice,Basis,...
                                    Symmetry,Vectors_RH,Sym_Tol);
      for S=1:NSym
        Sym_Vec=[Sym_Vec; GrainRot(count).Sym(:,:,S)*Vectors_RH];
      end
    else
      Unique=true;
    end
  end

end

nOrientations=count;

if nOrientations==0 
  error(['No orthogonal rotation vectors found.' ...
         ' Try increasing dis_tol, deg_tol, or maxMiller.'])
end

if Verbose==true
fprintf('\nFound %d ~orthogonal vectors with the [%d %d %d]' ...
        ,nOrientations,Axis(1), Axis(2), Axis(3))
fprintf(' oriented along the %c-direction.\n'...
        ,Direction)
      
fprintf('- Note vector angles are between %0.2f and %0.2f.\n' ...
        ,min([GrainRot.angle]),max([GrainRot.angle]))      
end

GrainRot(1).Info.Lattice=Lattice;
GrainRot(1).Info.Basis=Basis;
GrainRot(1).Info.Species=Species;
GrainRot(1).Info.Masses=Masses;
GrainRot(1).Info.Axis=Axis;
GrainRot(1).Info.Direction=Direction;
GrainRot(1).Info.MaxMiller=MaxMiller;
GrainRot(1).Info.Dis_Tol=Dis_Tol;
GrainRot(1).Info.Deg_Tol=Deg_Tol;
GrainRot(1).Info.Symmetry=Symmetry;
GrainRot(1).Info.Sym_Tol=Sym_Tol;

save(['Rot_' Species{:} '_' num2str(Axis(1)) '_' ...
      num2str(Axis(2)) '_' num2str(Axis(3)) '_' Direction ...
      '_' num2str(MaxMiller) '.mat'],'GrainRot')

else
  
fprintf('Reusing previously computed grain rotations\n')

end


end