function GB_Construct(GBorientations,Const_Options)


% Creating grain boundaries using grain boundary orientaitons 
% and options to fine tune the structure (i.e. translations, 
% atom deletions, terminating plane adjustments)


%% Setup variables

% Read information stored in GBorientation structure
Info=fieldnames(GBorientations(1).Info);
for  i=1:length(Info)
	eval([Info{i} '=' 'GBorientations(1).Info.(Info{i});']);
end

% Read optional paramaters stored in Options structure
if exist('Const_Options')==1
  Opts=fieldnames(Const_Options);
  for  i=1:length(Opts)
    eval([Opts{i} '=' 'Const_Options.(Opts{i});']);
  end
end

w = whos;
for a = 1:length(w)
  vars.(w(a).name) = eval(w(a).name);
end
 
vars_pass=GB_Vars(vars,'GB_Construct');

varNames=fieldnames(vars_pass);
for  i=1:length(varNames)
	eval([varNames{i} '=' 'vars_pass.(varNames{i});']);
end


%% Selection of GBs to study

% Create GBs within a disorientaiton angle range
if Angle_Range(1)>=0
  Selection=[GBorientations.Angle]>=Angle_Range(1) & [GBorientations.Angle]<=Angle_Range(2);
  nGBmax=sum(Selection);
  if nGBmax<1
    error('No GBs found within Angle_Range')
  end

  GBSort=find(Selection);
  [~,GBAreaSort]=sort([GBorientations(GBSort).Area]);
  GBSort=GBSort(GBAreaSort);

  if nGBs>nGBmax
    nGBs=nGBmax;
  end
end

GBid=GBSort(1:nGBs);


%% Naming conventions:
Dir_Base=[GBType '__'];

if length(Species1)==length(Species2)
  if strcmp(Species1,Species2)==1;
    Dir_Base=[Dir_Base Species1{:} '__'];
    if strcmp(GBType,'twist') || strcmp(GBType,'tilt');
      if norm(Axis1-Axis2)<10e-9 % tilt/twist same axis
        Dir_Base=[Dir_Base  num2str(Axis1(1)) '_' ...
                  num2str(Axis1(2)) '_' num2str(Axis1(3))];
      else
       Dir_Base=[Dir_Base num2str(Axis1(1)) '_' ...
                 num2str(Axis1(2)) '_' num2str(Axis1(3)) '-' ...
                 num2str(Axis2(1)) '_' num2str(Axis2(2)) '_' ...
                 num2str(Axis2(3))];
      end
    else % mixed
      Dir_Base=[Dir_Base Direction1 '_' num2str(Axis1(1)) ...
                '_' num2str(Axis1(2)) '_' num2str(Axis1(3)) ...
                '-' Direction2 '_' num2str(Axis2(1)) '_' ...
                num2str(Axis2(2)) '_' num2str(Axis2(3))];
    end
  else
    Dir_Base=[Dir_Base Species1{:} '-' Species2{:} '__'];
    if strcmp(GBType,'twist') || strcmp(GBType,'tilt');
      if norm(Axis1-Axis2)<10e-9 % tilt/twist same axis
        Dir_Base=[Dir_Base num2str(Axis1(1)) '_' ...
                  num2str(Axis1(2)) '_' num2str(Axis1(3))];
      else
       Dir_Base=[Dir_Base num2str(Axis1(1)) '_' ...
                 num2str(Axis1(2)) '_' num2str(Axis1(3)) '-' ...
                 num2str(Axis2(1)) '_' num2str(Axis2(2)) '_' ...
                 num2str(Axis2(3))];
      end
    else % mixed
      Dir_Base=[Dir_Base Direction1 '_' num2str(Axis1(1)) ...
                '_' num2str(Axis1(2)) '_' num2str(Axis1(3)) ...
                '-' Direction2 '_' num2str(Axis2(1)) '_' ...
                num2str(Axis2(2)) '_' num2str(Axis2(3))];
    end
  end
else 
  Dir_Base=[Dir_Base Species1{:} '-' Species2{:} '__'];
  if strcmp(GBType,'twist') || strcmp(GBType,'tilt');
    if norm(Axis1-Axis2)<10e-9 
      Dir_Base=[Dir_Base num2str(Axis1(1)) '_' ...
                num2str(Axis1(2)) '_' num2str(Axis1(3))];
    else
      Dir_Base=[Dir_Base num2str(Axis1(1)) '_'... 
                num2str(Axis1(2)) '_' num2str(Axis1(3)) '-' ...
                num2str(Axis2(1)) '_' num2str(Axis2(2)) '_' ...
                num2str(Axis2(3))];
   end
  else % mixed
    Dir_Base=[Dir_Base Direction1 '_' num2str(Axis1(1)) '_' ...
              num2str(Axis1(2)) '_' num2str(Axis1(3)) '-' ...
              Direction2 '_' num2str(Axis2(1)) '_' ...
              num2str(Axis2(2)) '_' num2str(Axis2(3))];
  end
end

File_Base=@(id)[Dir_Base '__' ...
              num2str(GBorientations(id).Angle,'%.2f') '_deg' ];

%% Setup stoichiometry calculations 
Lat1Types=[unique(Basis1(:,1))
           unique(Basis1(:,1))+max(Basis1(:,1))];
for i=1:length(unique(Basis1(:,1)))
  Stoich1_ideal(i,1)=sum(ismember(Basis1(:,1), ...
                         Lat1Types(i)))/length(Basis1);
end

Lat2Types=[unique(Basis2(:,1))
           unique(Basis2(:,1))+max(Basis2(:,1))];
for i=1:length(unique(Basis2(:,1)))
  Stoich2_ideal(i,1)=sum(ismember(Basis2(:,1), ...
                         Lat2Types(i)))/length(Basis2);
end
Lat2Types=Lat2Types+2*max(Basis1(:,1));
      
GBTypes=[unique(Basis1(:,1))+max(Basis1(:,1))
        unique(Basis2(:,1))+2*max(Basis1(:,1))+max(Basis2(:,1))];
Species=[repmat(Species1,[1,2]) repmat(Species2,[1,2]) ]';

Masses=[repmat(reshape(Masses1,prod(size(Masses1)),1),[2 1])
        repmat(reshape(Masses2,prod(size(Masses2)),1),[2 1])];
Types=[Lat1Types;Lat2Types];

nTypes1=length(Lat1Types);
nTypes2=length(Lat2Types);
nTypes=length(Types);

Ideal_Stoich=[Stoich1_ideal;Stoich1_ideal;...
              Stoich2_ideal;Stoich2_ideal];



%% Prepare For Transfer
if Archive==true
  if ismember(1,Write_Style)
      if exist([Dir_Base '_Data'],'dir')~=7
        mkdir([Dir_Base '_Data'])
      end
  end
  if ismember(2,Write_Style)
      if exist([Dir_Base '_Dump'],'dir')~=7
        mkdir([Dir_Base '_Dump'])
      end
  end
  if ismember(3,Write_Style)
      if exist([Dir_Base '_Car'],'dir')~=7
        mkdir([Dir_Base '_Car'])
      end
  end
end

%% Looping over GB orientaitons
GBcount=0;
File_List={};
for id=GBid                    % Looping over GB set
  %% Adjust naming convention to not write over files
  GBcount=GBcount+1;
  File_Base_id=[File_Base(id) Suffix];
  [Nfiles,~]=size(File_List);
  File_Base_id0=File_Base_id;
  j=1;
  for i = 1:Nfiles
      while strcmp(File_Base_id,File_List(i,:)) ==1
        File_Base_id=[File_Base_id0 '__' num2str(j)];
        j=j+1;
    end
  end
  File_List(GBcount,:)={File_Base_id};

  %% Construct Grain Boundaries
  Planes=GBorientations(id).Planes;
  CommonAxis=GBorientations(id).CommonAxis;
  Angle=GBorientations(id).Angle;
  Strain=GBorientations(id).Strain;
  GBLat1=round(GBorientations(id).Lat(1).vectors/Lattice1);
  GBLat2=round(GBorientations(id).Lat(2).vectors/Lattice2);
  nLat1=GBorientations(id).Lat(1).n;
  nLat2=GBorientations(id).Lat(2).n;
  dimLat1=GBorientations(id).Lat(1).dim;
  dimLat2=GBorientations(id).Lat(2).dim;
  
  
  fprintf('\n\n# Working On: [%d %d %d]',CommonAxis(1:3));
  fprintf('%s/[%d %d %d]%s',Direction1,CommonAxis(4:6),...
                            Direction2);
  fprintf(', (%d %d %d)/(%d %d %d), %.4f %s ',Planes,Angle,GBType);
  fprintf('grain boundary.\n');
  
  fid=[];
  Head1={};
  if Verbose==true
    fid=[1];
  end  
  if ismember(1,Write_Style)
    fid_data=2;
    fid=[fid,fid_data];
    Head1={};
  end

  if Verbose==true || ismember(1,Write_Style)
    Head1=mfprintf(fid,Head1,'#\n#      GB Index Number: %d\n',id);
    Head1=mfprintf(fid,Head1,'#     Type of boundary: [%d %d %d]%s/[%d %d %d]%s, (%d %d %d)/(%d %d %d) %s\n', CommonAxis(1:3),Direction1,CommonAxis(4:6),Direction2,Planes,GBType);
    Head1=mfprintf(fid,Head1,'# Misorientation angle: %6.4f degrees\n',Angle);
    Head1=mfprintf(fid,Head1,'#    GB overlap cutoff: %7.4f angstroms\n',Overlap_Tol);
    Head1=mfprintf(fid,Head1,'#                 -Top  Lat1-   -Bottom Lat2-\n');
    Head1=mfprintf(fid,Head1,'# Orientation: %4d %4d %4d %4d %4d %4d\n',GBLat1(1,:),GBLat2(1,:));
    Head1=mfprintf(fid,Head1,'#              %4d %4d %4d %4d %4d %4d\n',GBLat1(2,:),GBLat2(2,:));
    Head1=mfprintf(fid,Head1,'#              %4d %4d %4d %4d %4d %4d\n',GBLat1(3,:),GBLat2(3,:)); 
    Head1=mfprintf(fid,Head1,'#  Dimension:     N  Distance        N  Distance\n'); 
    Head1=mfprintf(fid,Head1,'#          x: %5d %9.4f    %5d %9.4f\n',nLat1(1),dimLat1(1),nLat2(1),dimLat2(1));
    Head1=mfprintf(fid,Head1,'#          y: %5d %9.4f    %5d %9.4f\n',nLat1(2),dimLat1(2),nLat2(2),dimLat2(2));
    Head1=mfprintf(fid,Head1,'#          z: %5d %9.4f    %5d %9.4f\n#\n',nLat1(3),dimLat1(3),nLat2(3),dimLat2(3));
    Head1=mfprintf(fid,Head1,'#    Strain  (Lat1-Lat2)/Lat1 :\n'); 
    Head1=mfprintf(fid,Head1,'#          x:  %12.11f\n',Strain(1));  
    Head1=mfprintf(fid,Head1,'#          y:  %12.11f\n',Strain(2));  
    Head1=mfprintf(fid,Head1,'#          z:  %12.11f\n#\n',Strain(3));
    
    Head1=mfprintf(fid,Head1,'# Min gamma Surface Dimension :\n') ; 
    Head1=mfprintf(fid,Head1,'#         gx:  %12.11f\n',max(GBorientations(id).Lat(1).norms(1),GBorientations(id).Lat(2).norms(1)));  
    Head1=mfprintf(fid,Head1,'#         gy:  %12.11f\n',max(GBorientations(id).Lat(1).norms(2),GBorientations(id).Lat(2).norms(2)));   
    Head1=mfprintf(fid,Head1,'#         gz:  %12.11f\n#\n', max(GBorientations(id).Lat(1).norms(3),GBorientations(id).Lat(2).norms(3)));     
    
    Head1=mfprintf(fid,Head1,'#      Types:');
    for s=1:length(Types)
      Head1=mfprintf(fid,Head1,'%10d  ' ,Types(s));
    end
    Head1=mfprintf(fid,Head1,'\n#    Species:');
    for s=1:length(Species)
      Head1=mfprintf(fid,Head1,'%10s  ' ,Species{s,:});
    end
    Head1=mfprintf(fid,Head1,'\n#     Masses:');
    for s=1:length(Masses)
      Head1=mfprintf(fid,Head1,'%10.6f  ' ,Masses(s));
    end
    Head1=mfprintf(fid,Head1,'\n#\n#  Lat1Types:');
    for s=1:prod(size(Lat1Types))
      Head1=mfprintf(fid,Head1,'%4d ',Lat1Types(s));
    end
    Head1=mfprintf(fid,Head1,'\n#  Lat2Types:');
    for s=1:prod(size(Lat2Types))
      Head1=mfprintf(fid,Head1,'%4d ',Lat2Types(s));
    end
    Head1=mfprintf(fid,Head1,'\n#    GBTypes:');
    for s=1:length(GBTypes)
      Head1=mfprintf(fid,Head1,'%4d ',GBTypes(s));
    end
    Head1=mfprintf(fid,Head1,'\n');
  end
    
  %% Create Lattices
  % Determine the min and maximum translations of the Basis 
  %  required to fill the simulation space

  Lx=max([GBorientations(id).Lat(1).dim(1) ...
          GBorientations(id).Lat(2).dim(1)]);
  Lz=max([GBorientations(id).Lat(1).dim(3) ...
          GBorientations(id).Lat(2).dim(3)]);
  
  % Use periodic distance in y-direction if desired
  Head1=mfprintf(fid,Head1,'#\n#  FullyPeriodic: %d \n',FullyPeriodic);
  if FullyPeriodic==true
    Tmp1=ceil(NormSlabDim(1)/GBorientations(id).Lat(1).norms(2));
    Tmp2=ceil(NormSlabDim(2)/GBorientations(id).Lat(2).norms(2));
    NormSlabDim=[Tmp1*GBorientations(id).Lat(1).norms(2) ...
                 Tmp2*GBorientations(id).Lat(2).norms(2)];
    Vacuum=0;
    clear Tmp1 Tmp2
  end
  
  % Create Lattice 1
  Ly=-NormSlabDim(1)-Plane_Shift(1);

  for k=[1 3]
      if Strain(k)<0
          Strain1(k)=-Strain(k);
      else 
          Strain1(k)=0;
      end
  end

  AtomDat1=GB_FillRegion(Lattice1,Basis1,...
                         GBorientations(id).Lat(1),Lx,Ly,Lz,0,Strain1);

  % Shift grain inward to expose different terminating plane
  AtomDat1(:,3)=AtomDat1(:,3)+Plane_Shift(1);                       
  AtomDat1(AtomDat1(:,3)>0,:)=[];
  
  % Create Lattice 2    
  Ly=NormSlabDim(2)+Plane_Shift(2);
  
  for k=[1 3]
      if Strain(k)>0
          Strain2(k)=Strain(k);
      else 
          Strain2(k)=0;
      end
  end
  
  AtomDat2=GB_FillRegion(Lattice2,Basis2,...
                         GBorientations(id).Lat(2),Lx,Ly,Lz,...
                         2*max(Basis1(:,1)),Strain2);
                       
  % Shift grain inward to expose different terminating plane                       
  AtomDat2(:,3)=AtomDat2(:,3)-Plane_Shift(2); 
  AtomDat2(AtomDat2(:,3)<0,:)=[];
                       
  AtomData=[AtomDat1;AtomDat2];
  Initial_AtomTypes=AtomData(:,1);
         
  %% Assigning atom types to the various regions
  
  AtomData(:,1)=Initial_AtomTypes;
  GBSelect1=abs(AtomData(:,3))<=GBregion & ismember(AtomData(:,1),Lat1Types);
  AtomData(GBSelect1,1)=Initial_AtomTypes(GBSelect1,1)+max(Basis1(:,1));  
  GBSelect2=abs(AtomData(:,3))<=GBregion & ismember(AtomData(:,1),Lat2Types);
  AtomData(GBSelect2,1)=Initial_AtomTypes(GBSelect2,1)+max(Basis2(:,1));  
 
  if Verbose==true
    Head1=mfprintf(fid,Head1,'#\n#  Initial Simulation Conents (GB region):\n');
    Head1=VerboseDetail(Stoich,AtomData,Lat1Types,Lat2Types,GBTypes,...
                  Species1,Species2,Stoich1_ideal,...
                  Stoich2_ideal,fid,Head1);
  end
  
  clear x y z buff nLat1 nLat2 Rot1 Rot2 BoundsLat1 BoundsLat2
  clear nAtoms nAtomsGB nLat1 nLat2 nLat1GB nLat2GB s
  
  
  %% Adjusting for overlapping periodic boundary conditions
  PBC_Overlap_Tol=1;
  [AtomData,Overlap_Id]=GB_WrapPBC(AtomData,[Lx 0 Lz],PBC_Overlap_Tol);
  nPBC_overlap=length(Overlap_Id);
  Initial_AtomTypes(Overlap_Id)=[];

  if Verbose==true || ismember(1,Write_Style)
    Head1=mfprintf(fid,Head1,'#\n#  Removing %d overlapping atoms at PBC :  \n',...
            nPBC_overlap);
   Head1=VerboseDetail(Stoich,AtomData,Lat1Types,Lat2Types,...
                  GBTypes,Species1,Species2,Stoich1_ideal,...
                  Stoich2_ideal,fid,Head1);
  end
  
  clear check delta m n1 n2 overlap_tol overlap_id buff Edge1 
  clear Edge2 Lat1Select Lat2Select k nPBC_overlap nLat1 nLat2  
  clear overlap_id nAtoms nAtomsGB nLat1 nLat2 nLat1GB nLat2GB s 
  clear nGB_Overlap

  %% Translations 
  
  Lat1=ismember(AtomData(:,1),Lat1Types);
  Lat2=ismember(AtomData(:,1),Lat2Types);
  
  AtomData(Lat1,2)=AtomData(Lat1,2)+Translate(1);
  AtomData(Lat1,3)=AtomData(Lat1,3)+Translate(2);
  AtomData(Lat1,4)=AtomData(Lat1,4)+Translate(3);
  AtomData(Lat2,2)=AtomData(Lat2,2)+Translate(4);
  AtomData(Lat2,3)=AtomData(Lat2,3)+Translate(5);
  AtomData(Lat2,4)=AtomData(Lat2,4)+Translate(6);
  
  Wrap=1;
  while Wrap>0
    Wrap1=AtomData(:,2)<0;
    AtomData(Wrap1,2)=AtomData(Wrap1,2)+Lx;
    Wrap2=AtomData(:,2)>Lx;
    AtomData(Wrap2,2)=AtomData(Wrap2,2)-Lx;
    Wrap3=AtomData(:,4)<0;
    AtomData(Wrap3,4)=AtomData(Wrap3,4)+Lz;
    Wrap4=AtomData(:,4)>Lz;
    AtomData(Wrap4,4)=AtomData(Wrap4,4)-Lz;
    Wrap=sum(Wrap1+Wrap2+Wrap3+Wrap4);
  end
  
  
  %% Re-assigning atom types to the various regions
  
  AtomData(:,1)=Initial_AtomTypes;
  GBSelect1=abs(AtomData(:,3))<=GBregion & ismember(AtomData(:,1),Lat1Types);
  AtomData(GBSelect1,1)=Initial_AtomTypes(GBSelect1,1)+max(Basis1(:,1));  
  GBSelect2=abs(AtomData(:,3))<=GBregion & ismember(AtomData(:,1),Lat2Types);
  AtomData(GBSelect2,1)=Initial_AtomTypes(GBSelect2,1)+max(Basis2(:,1));  
  
  if sum(GBSelect1)==0 
    error('Grain 1 y-translation must be less than %f (GBregion)',GBregion)
  end
  if sum(GBSelect1)==0
    error('Grain 2 y-translation must be less than %f (GBregion)',GBregion)
  end
  
  %% Adjusting for overlapping atoms at the GB
  overlap_id=[];
  buff=Overlap_Tol;  % only consider atoms within buffer near boundary

  Edge1=find(ismember(AtomData(:,1),Lat1Types) & AtomData(:,3)>=-buff);
  n1=length(Edge1);
  Edge2=find(ismember(AtomData(:,1),Lat2Types) & AtomData(:,3)<=buff);
  n2=length(Edge2);
  for m=Edge1'
    delta=(ones(n2,1)*AtomData(m,2:4)-AtomData(Edge2,2:4));
    check=min(sqrt(sum(delta.^2,2)));
    if check<=Overlap_Tol
      overlap_id=[overlap_id; m];
    end
  end
  
  if FullyPeriodic==true
    % Removing atoms from both GBs in fully periodic simulation
    Ly=NormSlabDim(2)+NormSlabDim(1);
    Edge1=find(AtomData(:,3)<=min(AtomData(:,2))+buff);
    n1=length(Edge1);
    Edge2=find(AtomData(:,3)>=max(AtomData(:,2))-buff);
    n2=length(Edge2);
    for m=Edge1'
      delta=(ones(n2,1)*(AtomData(m,2:4)+[0 Ly 0])-AtomData(Edge2,2:4));
      check=min(sqrt(sum(delta.^2,2)));
      if check<=PBC_Overlap_Tol
        overlap_id=[overlap_id; m];
      end
    end
  end
  
  overlap_id=unique(overlap_id);
  nGB_overlap=length(overlap_id);
  AtomData(overlap_id,:)=[];

  if Verbose==true || ismember(1,Write_Style)
    Head1=mfprintf(fid,Head1,'#\n#  Removing %d overlapping atoms at GB :  \n',nGB_overlap); 
    Head1=VerboseDetail(Stoich,AtomData,Lat1Types,Lat2Types,GBTypes,...
                  Species1,Species2,Stoich1_ideal,...
                  Stoich2_ideal,fid,Head1);
  end

  clear buff Edge1 Edge2 n1 n2 delta check
  clear GBSelect1 GBSelect2  m nLat1 nLat2
  
  
  %% Forcing stoicheometric regions by removing atoms on surfaces
  
  if ForceStoich==true
    nAtoms=length(AtomData);
    nLat1=sum(ismember(AtomData(:,1),Lat1Types));
    nLat2=sum(ismember(AtomData(:,1),Lat2Types));

    [N1,D1]=rat(Stoich1_ideal);
    [N2,D2]=rat(Stoich2_ideal);
    if sum(D1(1)==D1)~=length(D1) || sum(D2(1)==D2)~=length(D2)
      error('Cannot determine correct stoich')
    end
    D1=D1(1);
    D2=D2(1);
    
    AtomData=sortrows(AtomData,3);
    Ids=[1:nAtoms]';
    Remove=[];

    U1=floor(nLat1/D1);
    for i=1:length(Lat1Types)/2
      target=U1*N1(i);
      TypeI=ismember(AtomData(:,1), [Lat1Types(i),Lat1Types(i+length(Lat1Types)/2)]);
      actual=sum(TypeI);
      Remove_i=Ids(TypeI,:);
      if actual>target
        Remove_i=Remove_i(1:actual-target);
        Remove=[Remove;Remove_i];
      elseif actual<target
        Remove_i=Remove_i(1:target-actual);
        Remove_i=Remove_i(end)
      end
    end
    
    U2=floor(nLat2/D2);
    for i=1:length(Lat2Types)/2
      target=U2*N2(i);
      TypeI=ismember(AtomData(:,1), [Lat2Types(i),Lat2Types(i+length(Lat2Types)/2)]);
      actual=sum(TypeI);
      Remove_i=Ids(TypeI,:);
      if actual>target
        Remove_i=Remove_i(end-(actual-target-1):end);
        Remove=[Remove;Remove_i];
      elseif actual<target
        Remove_i=Remove_i(end-(target-actual-1):end);
        Remove=[Remove;Remove_i];
      end
    end    
    nRemove=length(Remove);
    AtomData(Remove,:)=[];
    
    if Verbose==true || ismember(1,Write_Style)
      Head1=mfprintf(fid,Head1,'#\n#  Removing %d atoms for Stoich :  \n',nRemove);
      Head1=VerboseDetail(Stoich,AtomData,Lat1Types,Lat2Types,GBTypes,...
        Species1,Species2,Stoich1_ideal,...
        Stoich2_ideal,fid,Head1);
    end
  end
  

         
  %% Write Files
  if ismember(1,Write_Style)
    Header.Head1=Head1;
  end
  if ismember(2,Write_Style)
    Head2={''};
    Header.Head2=Head2;
  end
  if ismember(3,Write_Style)
    Head3={sprintf('# %s GB: [%d %d %d]%s/[%d %d %d]%s, (%d %d %d)/(%d %d %d), %.4f deg\n',GBType,CommonAxis(1:3),Direction1,CommonAxis(4:6),Direction2,Planes,Angle)};
    Header.Head3=Head3;
  end
  Corners=[0, Lx, -NormSlabDim(1)-Vacuum/2, ...
           NormSlabDim(2)+Vacuum/2, 0 Lz]; 
  GB_WriteFiles(File_Base_id,Write_Style,AtomData,Species,Masses,Corners,Header,AtomStyle,Archive,Dir_Base,Verbose)
  clear Corners
  

  
  if Summary==true;
    if GBcount==1;
      fid_report=fopen([Summary_Name],'w');
      fprintf(fid_report,'%40s | ','FileName');
      fprintf(fid_report,'%27s | ','[Common Axis] direction');
      fprintf(fid_report,'%23s | ','(GB planes)');
      fprintf(fid_report,'%9s | ','Angle');
      fprintf(fid_report,'%11s | ','Overlap Tol');
      fprintf(fid_report,'%42s | ','Translations');
      fprintf(fid_report,'%15s | ','Planar Shifts'); 
      fprintf(fid_report,'%11s | ','Area');
      fprintf(fid_report,'%10s \n','Natoms');
    end
      
    fprintf(fid_report,'%40s | ',File_Base_id);
    fprintf(fid_report,'[%3d%3d%3d]%2s/[%3d%3d%3d]%2s | ', CommonAxis(1:3),Direction1,CommonAxis(4:6),Direction2);
    fprintf(fid_report,'(%3d%3d%3d)/(%3d%3d%3d) | ',Planes);
    fprintf(fid_report,'%9.4f | ',Angle);
    fprintf(fid_report,'%11.5f | ',Overlap_Tol);
    fprintf(fid_report,'%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f |', Translate);
    fprintf(fid_report,'%8.4f%8.4f | ',Plane_Shift);
    fprintf(fid_report,'%11.5f | ',Lx*Lz);
    fprintf(fid_report,'%10d\n',length(AtomData(:,1)));
    

  end
  
    
end

if Summary==true;
  fclose(fid_report);
end
  
clear DataFile DumpFile GBvarFile GBID
clear nAtoms nAtomsGB nLat1 nLat2 nLat1GB nLat2GB s


%% Bulk
GB_Bulk(GBorientations,Write_Style,AtomStyle,Archive,Dir_Base,Verbose)

%% Prepare Output For Transfer
if Archive==true
  if ismember(1,Write_Style)
    tar([Dir_Base '_Data.tgz'],[Dir_Base '_Data'])
  end
  if ismember(2,Write_Style)
    tar([Dir_Base '_Dump.tgz'],[Dir_Base '_Dump'])
  end
  if ismember(3,Write_Style)
    tar([Dir_Base '_Car.tgz'],[Dir_Base '_Car'])
  end
end

fprintf('\n\n# Finished:  Wrote %d GB Files\n',GBcount);

end


% Write verbose to screen and store header for data files
function Head1=mfprintf(fid,Head1,varargin)
  for i=1:length(fid)
    if fid(i)==1
      fprintf(fid(i), varargin{:});
    else
      Head1(length(Head1)+1)={sprintf(varargin{:})};
    end
  end
end

% Write detailed information in verbose mode
function Head1=VerboseDetail(Stoich,AtomData,Lat1Types,Lat2Types,GBTypes,Species1,Species2,Stoich1_ideal,Stoich2_ideal,fid,Head1)
  nAtoms=length(AtomData);
  nAtomsGB=sum(ismember(AtomData(:,1),GBTypes));
  nLat1=sum(ismember(AtomData(:,1),Lat1Types));
  nLat2=sum(ismember(AtomData(:,1),Lat2Types));
  nLat1GB=sum(ismember(AtomData(:,1),intersect(Lat1Types,GBTypes)));
  nLat2GB=sum(ismember(AtomData(:,1),intersect(Lat2Types,GBTypes)));

  Head1=mfprintf(fid,Head1,'#\n#       Total number of atoms: %d  (%d) \n',nAtoms,nAtomsGB);
  Head1=mfprintf(fid,Head1,'#                  Lat1 atoms: %d  (%d) \n',nLat1,nLat1GB);
  Head1=mfprintf(fid,Head1,'#                  Lat2 atoms: %d  (%d) \n',nLat2,nLat2GB);
  
  if Stoich==1
    Head1=mfprintf(fid,Head1,'#\n# Stoichiometry Lat1:  Type    Total       (GB)    Ideal   Species  \n');
    for i=1:length(Lat1Types)/2
      a=sum(ismember(AtomData(:,1), [Lat1Types(i),Lat1Types(i+length(Lat1Types)/2)]))/nLat1;
      b=sum(ismember(AtomData(:,1), Lat1Types(i+length(Lat1Types)/2)))/nLat1GB;
      Head1=mfprintf(fid,Head1,'#                       %3d   %4.4f   (%4.4f)   %4.4f   %s\n',Lat1Types(i),a,b,Stoich1_ideal(i),Species1{i});
      Stoich1(i,1)=a;
      Stoich1(i+length(Lat1Types)/2,1)=b;
    end

    Head1=mfprintf(fid,Head1,'#\n# Stoichiometry Lat2:  Type    Total       (GB)    Ideal   Species  \n');
    for i=1:length(Lat2Types)/2
      a=sum(ismember(AtomData(:,1), [Lat2Types(i),Lat2Types(i+length(Lat2Types)/2)]))/nLat2;
      b=sum(ismember(AtomData(:,1), Lat2Types(i+length(Lat2Types)/2)))/nLat2GB;
      Head1=mfprintf(fid,Head1,'#                       %3d   %4.4f   (%4.4f)   %4.4f   %s\n',Lat2Types(i),a,b,Stoich2_ideal(i),Species2{i});
      Stoich1(length(Lat1Types)+i,1)=a;
      Stoich1(length(Lat1Types)+i+length(Lat2Types)/2,1)=b;
    end  
  end
end