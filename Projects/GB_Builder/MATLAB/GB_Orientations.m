function [GBorientations,id_Area,id_Angle]=GB_Orientations(GrainRot1,GrainRot2,Orient_Options)

% Combing two grain rotations into 

if nargin==1
  GrainRot2=GrainRot1;
end

%% Setup variables

% Read information stored in GrainRot1/2 structures
Info=fieldnames(GrainRot1(1).Info);
for  i=1:length(Info)
	eval([Info{i} '1=' 'GrainRot1(1).Info.(Info{i});']);
end

Info=fieldnames(GrainRot2(1).Info);
for  i=1:length(Info)
	eval([Info{i} '2=' 'GrainRot2(1).Info.(Info{i});']);
end

% Read optional paramaters stored in Options structure
if exist('Orient_Options')==1
  Opts=fieldnames(Orient_Options);
  for  i=1:length(Opts)
    eval([Opts{i} '=' 'Orient_Options.(Opts{i});']);
  end
end

clear Info

w = whos;
for a = 1:length(w)
  vars.(w(a).name) = eval(w(a).name);
end
vars_pass=GB_Vars(vars,'GB_Orientations');

varNames=fieldnames(vars_pass);
for  i=1:length(varNames)
    eval([varNames{i} '=' 'vars_pass.(varNames{i});'])
end


%% Check if same rotation has already been calculated

Compute=true;

if exist(['GBorient_' Species1{:} '_' num2str(Axis1(1)) '_' ...
      num2str(Axis1(2)) '_' num2str(Axis1(3)) '_' Direction1 ...
      '_' num2str(MaxMiller1) '_' Species2{:} '_' ...
      num2str(Axis2(1)) '_' num2str(Axis2(2)) '_' ...
      num2str(Axis2(3)) '_' Direction2 '_' ...
      num2str(MaxMiller2) '.mat'],'file')==2
  % Save original variables
  varNames=fieldnames(vars_pass);
  for  i=1:length(varNames)
    eval([varNames{i} '1=' 'vars_pass.(varNames{i});'])
  end
  
  % Load saved varialbes
  load(['GBorient_' Species1{:} '_' num2str(Axis1(1)) '_' ...
      num2str(Axis1(2)) '_' num2str(Axis1(3)) '_' Direction1 ...
      '_' num2str(MaxMiller1) '_' Species2{:} '_' ...
      num2str(Axis2(1)) '_' num2str(Axis2(2)) '_' ...
      num2str(Axis2(3)) '_' Direction2 '_' ...
      num2str(MaxMiller2) '.mat'])
  Info=fieldnames(GBorientations(1).Info);
  for  i=1:length(Info)
    eval([Info{i} '=' 'GBorientations(1).Info.(Info{i});']);
  end
  
  Compute=false;
  % Compare original and saved variables  
  for  i=1:length(varNames)
    if isequal(eval([varNames{i} '1']),eval(varNames{i}))~=1
          Compute=true;
          eval( [varNames{i} '=' [varNames{i} '1'] ';'])
    end
  end
end


if Compute==true;
clear GBorientations
%%

RotAxis1=find((sum((GrainRot1(1).vectors-GrainRot1(length(GrainRot1)).vectors)==0,2))==3);
RotAxis2=find((sum((GrainRot2(1).vectors-GrainRot2(length(GrainRot2)).vectors)==0,2))==3);
RotAxis1=RotAxis1(1);
RotAxis2=RotAxis2(1);

if RotAxis1==RotAxis2 && RotAxis1==2
  GBType='twist';
  SymmT=[0 0 1;0 1 0; 1 0 0];
  AngleAxis=1;
elseif RotAxis1==RotAxis2 && RotAxis1~=2
  GBType='tilt';
  if RotAxis1==1;
    SymmT=[1 0 0;0 0 1; 0 1 0];
  else
    SymmT=[0 1 0; 1 0 0; 0 0 1];
  end
  AngleAxis=2;
else
  GBType='mixed';
  SymmT=eye(3);
end
tmp1=GrainRot1(1).vectors/Lattice1;
tmp2=GrainRot1(1).vectors/Lattice2;
Directions=[tmp1(RotAxis1,:) tmp2(RotAxis2,:)];

if strcmp(GBType,'twist') || strcmp(GBType,'tilt');
  if Verbose==true
    fprintf('\nDetermining [%d %d %d] / [%d %d %d] %s grain boundary geometries:\n',Directions,GBType)   
  end
end

Istart=1;
Istop=length(GrainRot1);
Jstart=@(i)1;
Jstop=length(GrainRot2);
TotK=((length(GrainRot1))*length(GrainRot2));

% Preallocating variables
n1=zeros(3,1);
n2=zeros(3,1);
dim1=zeros(3,1);
dim2=zeros(3,1);
Strain=zeros(3,1);

K=0;
Progress=0.1;
count=1;
for i=Istart:Istop
  for j=Jstart(i):Jstop
    K=K+1;
    if Verbose==true
      if K/TotK>= Progress
        fprintf(' %3.0f%%',Progress*100)
        Progress=Progress+0.1;
      end
    end

    Lat1=GrainRot1(i);
    Lat2=GrainRot2(j);

    if i<=j % Symmetric Boundaries
      Lat2.vectors=SymmT*Lat2.vectors;
      Lat2.norms=SymmT*Lat2.norms;
      Lat2.offplane=SymmT*Lat2.norms;
    end

    propflag=1;    % Flag to calculate properties
    savedata=1;    % Flag to save data
    % Determine periodic repeating distance for the GB  
    % simulation (GB normal is always in the y-direction)
    for k=[1 3];
      n1(k,1)=1;
      n2(k,1)=1;  
      dim1(k,1)=n1(k)*norm(Lat1.vectors(k,:));
      dim2(k,1)=n2(k)*norm(Lat2.vectors(k,:));
      Strain(k,1)=(dim1(k)-dim2(k))/dim1(k);

      while abs(Strain(k)) >= Strain_Tol 
        % Determine number of repeating units
        dim1(k,1)=n1(k)*norm(Lat1.vectors(k,:));
        dim2(k,1)=n2(k)*norm(Lat2.vectors(k,:));
        Strain(k,1)=(dim1(k)-dim2(k))/dim1(k);
        if Strain(k) < 0
          n1(k,1)=n1(k,1)+1;
        else
          n2(k,1)=n2(k,1)+1;
        end
      end
    end

    %% Check if exceeding maximimum allowed dimensions
    Area=max([dim1(1),dim2(1)]) * max([dim1(3),dim2(3)]);
    if Area > MaxArea
      propflag=0;
      savedata=0;
    end

    A=normr(Lat1.vectors)';
    B=normr(Lat2.vectors)';

    if propflag==1;
      %% Compute the disorientation
      k=0;
      [~,~,iSym1]=size(Lat1.Sym);
      [~,~,iSym2]=size(Lat2.Sym);

      Angle_i=10000;
      Angle=10000;
      RotMat_i=zeros(3);
      RotMat=RotMat_i;
      PlaneList=zeros(iSym1*iSym2,6);

      for ii=1:iSym1
        A1=(Lat1.Sym(:,:,ii)*A);
        Plane1=(A(2,:)*Lat1.Sym(:,:,ii));
        for jj=1:iSym2
          k=k+1;          
          B1=Lat2.Sym(:,:,jj)*B;
          Plane2=(B(2,:)*Lat2.Sym(:,:,jj));
          PlaneList(k,:)=[Plane1 Plane2];
          RotMat_i=A1*(B1)';
          Angle_i=acos((RotMat_i(1,1)+RotMat_i(2,2)+RotMat_i(3,3)-1)/2);
          if Angle_i<Angle;
            Angle=Angle_i;
            RotMat=RotMat_i;
          end
        end
      end

     PlaneList=[PlaneList;PlaneList(:,4:6) PlaneList(:,1:3)];
      
      %% Check for uniqueness
      if count>1
        CompareAngle=abs(GBAngles-Angle)<1e-7;
        % Check disorientation angle
        if sum(CompareAngle)>0
          % Check grain boundary plane
          CompareID=find(CompareAngle==1);
          if savedata==1
            for k=1:length(CompareID);
              Compare=GBorientations(CompareID(k)).PlaneList;
              [tmp,~]=size(Compare);
              if  sum(sum((repmat([A(2,:) B(2,:)],tmp,1)-Compare).^2')<1e-10) > 1
                savedata=0;
                break
              end
            end
          end
        end
      end
      
      %% Save data to GB structure
      if savedata==1;  
        GBorientations(count).CommonAxis=[Lat1.vectors(RotAxis1,:)/Lattice1...
                              Lat2.vectors(RotAxis2,:)/Lattice2];
        GBorientations(count).Planes=[Lat1.vectors(2,:)/Lattice1 Lat2.vectors(2,:)/Lattice2];   
        GBorientations(count).PlaneList=PlaneList;
        GBorientations(count).Lat(1)=Lat1;
        GBorientations(count).Lat(2)=Lat2;
        GBorientations(count).Lat(1).n=n1;
        GBorientations(count).Lat(2).n=n2;  
        GBorientations(count).Lat(1).dim=dim1;
        GBorientations(count).Lat(2).dim=dim2;
        GBorientations(count).Strain=Strain;
        GBorientations(count).Area=Area;
        GBorientations(count).RotMat=RotMat;
        GBorientations(count).Angle=Angle*180/pi;
        GBAngles(count)=Angle;
        count=count+1;
      end
    end      

  end
end 

nGBs=count-1;

if nGBs==0 
  error(['No GB orientations found. '  ...
     'Try increasing dis_tol or deg_tol, or decrease symmetry.'])
end

if Verbose==true
  fprintf('\nFound %d grain boundary orientations.',nGBs)     
end

clear count i j Orientation GBaxis GBp Val
clear propflag A B Angle C_Angle ID k nGBs nOrientation

%% Save general info

% Sorted Index
[~, id_Area]=sort([GBorientations.Area]);
[~, id_Angle]=sort([GBorientations.Angle]);

GBorientations(1).Info.GBType=GBType;
GBorientations(1).Info.Directions=Directions;
GBorientations(1).Info.Strain_Tol=Strain_Tol;
GBorientations(1).Info.MaxArea=MaxArea;
GBorientations(1).Info.Lattice1=Lattice1;
GBorientations(1).Info.Basis1=Basis1;
GBorientations(1).Info.Species1=GrainRot1(1).Info.Species;
GBorientations(1).Info.Masses1=GrainRot1(1).Info.Masses;
GBorientations(1).Info.Axis1=Axis1;
GBorientations(1).Info.Direction1=Direction1;
GBorientations(1).Info.MaxMiller1=MaxMiller1;
GBorientations(1).Info.Dis_Tol1=Dis_Tol1;
GBorientations(1).Info.Deg_Tol1=Deg_Tol1;
GBorientations(1).Info.Symmetry1=Symmetry1;
GBorientations(1).Info.Sym_Tol1=Sym_Tol1;
GBorientations(1).Info.Lattice2=Lattice2;
GBorientations(1).Info.Basis2=Basis2;
GBorientations(1).Info.Species2=GrainRot2(1).Info.Species;
GBorientations(1).Info.Masses2=GrainRot2(1).Info.Masses;
GBorientations(1).Info.Axis2=Axis2;
GBorientations(1).Info.Direction2=Direction2;
GBorientations(1).Info.MaxMiller2=MaxMiller2;
GBorientations(1).Info.Dis_Tol2=Dis_Tol2;
GBorientations(1).Info.Deg_Tol2=Deg_Tol2;
GBorientations(1).Info.Symmetry2=Symmetry2;
GBorientations(1).Info.Sym_Tol2=Sym_Tol2;


save(['GBorient_' Species1{:} '_' num2str(Axis1(1)) '_' ...
      num2str(Axis1(2)) '_' num2str(Axis1(3)) '_' Direction1 ...
      '_' num2str(MaxMiller1) '_' Species2{:} '_' ...
      num2str(Axis2(1)) '_' num2str(Axis2(2)) '_' ...
      num2str(Axis2(3)) '_' Direction2 '_' ...
      num2str(MaxMiller2) '.mat'],...
      'GBorientations')


else
  fprintf('Reusing previously computed grain boundary orientations\n')
end

end

