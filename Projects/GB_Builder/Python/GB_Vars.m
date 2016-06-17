function vars_pass = GB_Vars(vars,fun)
  
  varNames=fieldnames(vars);
  
  for  i=1:length(varNames)
    Variable=varNames{i};
    eval([(Variable) '=' 'vars.(varNames{i});'])
  end
  
  %% Ensure variables are used correctly
  for  i=1:length(varNames)
    % Adjust for GB variables(1/2)
    Variable=varNames{i};
    if strcmpi(Variable(end),'1') || strcmpi(Variable(end),'2')
      suffix=Variable(end);
    else
      suffix='';
    end
    
    switch varNames{i}
      
      case ['Lattice' suffix]
        [Lr,Lc] = size(eval(['Lattice' suffix]));
        if Lr~=3 || Lc ~=3
          error(['Lattice' suffix ' must be a 3x3 matrix.'])
        end

      case ['Basis' suffix]
        [~,Bc]=size(eval(['Basis' suffix]));
        if Bc~=4
          error(['Basis' suffix ' must contain 4 columns.'])
        end

      case ['Axis' suffix]
        [Ar,Ac]=size(eval(['Axis' suffix]));
        if prod([Ar,Ac])~=3
          error(['Axis' suffix ' must be a 3x1 or 1x3 vector.'])
        end
        if Ar==3
          tmp=reshape(eval(['Axis' suffix]),1,3);
          eval([['Axis' suffix] '=' 'tmp;'])
        end

      case ['Direction' suffix]
        if [strcmpi(eval(['Direction' suffix]),'x')  || ...
            strcmpi(eval(['Direction' suffix]),'tilt')]
          eval([['DirIdex' suffix] '=' '1;'])                  
          eval([['RH' suffix] '=' '[1 0 0;0 0 1;0 1 0];'])
          eval([['Direction' suffix] '=' '''x'';']) 
        elseif [strcmpi(eval(['Direction' suffix]),'y') || ...
                strcmpi(eval(['Direction' suffix]),'twist')]
          eval([['DirIdex' suffix] '=' '2;'])                   
          eval([['RH' suffix] '=' '[0 0 1;0 1 0;1 0 0];'])
          eval([['Direction' suffix] '=' '''y'';']) 
        elseif strcmpi(eval(['Direction' suffix]),'z') 
          eval([['DirIdex' suffix] '=' '3;'])                   
          eval([['RH' suffix] '=' '[0 1 0; 1 0 0; 0 0 1];'])
        else
          error(['Direction' suffix ...
                 ' (x, y, or z) must be set.'])
        end
        tmp=setdiff([1:3],eval(['DirIdex' suffix]));
        eval([['DirIdex' suffix '(2:3)'] '=' 'tmp;'])

      case ['MaxMiller' suffix]
        if rem(eval(['MaxMiller' suffix]),1)~=0
          error(['MaxMiller' suffix ' must be an integer.'])
        end

      case ['Basis' suffix]
        a=length(eval(['Species' suffix]));
        b=length(eval(['Masses' suffix]));
        c=length(unique(eval(['Basis' suffix '(:,1)'])));
        if a~=b || b~=c
          error(['Inconsistent number of atom Types (' ...
                  num2str(c) '), Species (' num2str(a)...
                  '), and Masses (' num2str(b) ').' ])
        end
        
      case 'NormSlabDim'
        if prod(size(NormSlabDim))==1
          NormSlabDim=[NormSlabDim NormSlabDim];
        end
        
      case 'GBSort'
        if strcmpi(GBSort,'Area')
          [~,GBSort]=sort([GBorientations.Area]);
        elseif strcmpi(GBSort,'Angle')
          [~,GBSort]=sort([GBorientations.Angle]);
        elseif strcmpi(GBSort,'Index')
          GBSort=1:length(GBorientations);
        else
          error('Invalid GB sorting technique')
        end
        
      case 'nGBs'
        if nGBs>length(GBorientations)
          nGBs=length(GBorientations);
        end

      case 'WriteStyle'
        for i=1:length(WriteStyle)
          if strcmpi(WriteStyle{i},'data')
            tmp(i)=1;
          elseif strcmpi(WriteStyle{i},'dump')
            tmp(i)=2;
          elseif strcmpi(WriteStyle{i},'car')
            tmp(i)=3;
          else
            error('Invalid WriteStyle')
          end
        end
        Write_Style=tmp;

        
      case 'Plane_Shift'
        if prod(size(Plane_Shift))~=2
          error('Plane_Shift must be 2x1 or 1x2 vector')
        end
        if Plane_Shift(1)<0 || Plane_Shift(2)<0
          error('Plane_Shift values must be greater than 0')
        end
        if exist('FullyPeriodic')==1
          if FullyPeriodic==true;
            FullyPeriodic=false;
            warning(['Cannot shift planes with fully '...
                     'periodic boundary conditions'])
          end
        end
           
      case 'Translate'
        if length(Translate)~=6 & prod(size(Translate))~=6
          error('Translate must be 6x1 or 1x6 vector')
        end  
        
      case 'Angle_Range'
        if length(Angle_Range)~=2 & prod(size(Angle_Range))~=2
          error('Translate must be 2x1 or 1x2 vector')
        end 
        
      case 'Corners'
        if length(Corners)~=6 & prod(size(Corners))~=6
          error('Corners must be 6x1 or 1x6 vector')
        end        

      case 'AtomData'
        [~,tmp]=size(AtomData);
        if tmp~=4
          error('AtomData must be a nAtoms x 4 matrix')
        end  
        
      case 'MaxEdge'
        if length(MaxEdge)~=3 & prod(size(MaxEdge))~=3
          error('MaxEdge must be a 3x1 or 3x1 vector')
        end
        
      case 'Verbose'
        if  Verbose~=1 & Verbose~=0
          error('Verbose must be true or false')
        end 
        
      case 'Summary'
        if  Summary~=1 & Summary~=0
          error('Summary must be true or false')
        end 
                    
            
      case 'Symmetry'
        if  Symmetry~=1 & Symmetry~=0
          error('Symmetry must be true or false')
        end
    end
  end

  
  %% Define required and optional variables for each function
  switch fun
    case 'GB_GrainRot'
      ReqVars={'Lattice';'Basis';'Species';'Masses';...
               'Axis';'Direction'};
      OptVars={'MaxMiller';'Dis_Tol';'Deg_Tol';'Symmetry';...
              'Sym_Tol';'Verbose'};
             
      
    case  'GB_SymOp'
      ReqVars={'Lattice';'Basis'};
      OptVars={'Symmetry','Orientation','Sym_Tol'};
      
      
    case  'GB_FillRegion'
      ReqVars={'Lattice';'Basis';'Lx';'Ly';'Lz'};
      OptVars={'TypeOffset'};
      
      
    case 'GB_Orientations'
      ReqVars={'Lattice1';'Basis1';'Species1';'Masses1';...
               'Axis1';'Direction1';'MaxMiller1';
               'Dis_Tol1';'Deg_Tol1';'Symmetry1';'Sym_Tol1';...
               'Lattice2';'Basis2';'Species2';'Masses2';...
               'Axis2';'Direction1';'MaxMiller2';
               'Dis_Tol2';'Deg_Tol2';'Symmetry2';'Sym_Tol2'};
      OptVars={'Strain_Tol';'MaxArea';'Verbose'};
      
      
    case 'GB_Construct'
      ReqVars={'Lattice1';'Basis1';'Species1';'Masses1';...
               'Axis1';'Direction1';'MaxMiller1';
               'Dis_Tol1';'Deg_Tol1';'Symmetry1';'Sym_Tol1';...
               'Lattice2';'Basis2';'Species2';'Masses2';...
               'Axis2';'Direction1';'MaxMiller2';
               'Dis_Tol2';'Deg_Tol2';'Symmetry2';'Sym_Tol2'};
      OptVars={'GBSort';'nGBs';'Overlap_Tol';'WriteStyle';...
               'NormSlabDim';'Vacuum';'GBregion';...
               'FullyPeriodic';'Stoich';'Archive';'Suffix';...
               'Verbose';'PBC_Overlap_Tol'; 'Plane_Shift';...
               'Translate';'Angle_Range';'Summary';...
               'Summary_Name';'ForceStoich';'AtomStyle'};

             
    case 'GB_WriteFiles'             
      ReqVars={'File_Base_id';'Write_Style';'AtomData';...
               'Species';'Masses';'Corners'};
      OptVars={'Header';'Archive';'Dir_Base';'Verbose'};  
        
        
    case 'GB_WrapPBC'
      ReqVars={'AtomData';'MaxEdge'};
      OptVars={'PBC_Overlap_Tol'};
            
  end

  
  %% Checking for required variables
  for i=1:length(ReqVars)
    if exist(ReqVars{i})==0
      error('%s must contain %s',fun,ReqVars{i})
    end
  end

  
  %% Set optional variable defaults
  for i=1:length(OptVars)
    switch OptVars{i}
      
      case 'MaxMiller'
        if exist('MaxMiller')==0
          MaxMiller=5;
        end  
        
      case 'Sym_Tol'
        if exist('Sym_Tol')==0
          Sym_Tol=0.001;
        end    
        
      case 'Dis_Tol'
        if exist('Dis_Tol')==0
          Dis_Tol=0.001;
        end

      case 'Deg_Tol'
        if exist('Deg_Tol')==0
          Deg_Tol=0.001;
        end                

      case 'Orientation'
        if exist('Orientation')==0
          Orientation=eye(3);
        end
            
      case 'Strain_Tol'
        if exist('Strain_Tol')==0
          Strain_Tol=1e-3;
        end
            
      case 'MaxArea'
        if exist('MaxArea')==0
          MaxArea=250*250;
        end
            
      case 'Symmetry'
        if exist('Symmetry')==0
          Symmetry=true;
        end
            
      case 'Verbose'
        if exist('Verbose')==0
          Verbose=true;
        end
        
      case 'GBSort'  
        if exist('GBSort')==0
          [~,GBSort]=sort([GBorientations.Area]);
        end
        
      case 'nGBs'
        if exist('nGBs')==0;
          nGBs=length(GBorientations);
        end
        
      case 'Overlap_Tol'
        if exist('Overlap_Tol')==0
          Overlap_Tol=0;
        end
        
      case 'WriteStyle'
        if exist('WriteStyle')==0
          Write_Style=[1 2 3];
        end
        
      case 'NormSlabDim'
        if exist('NormSlabDim')==0
          NormSlabDim=[50 50];
        end
        
      case 'Vacuum'
        if exist('Vacuum')==0
          Vacuum=20;
        end
      
      case 'GBregion'
        if exist('GBregion')==0
          GBregion=10;
        end
        
      case 'Suffix'
        if exist('Suffix')==0
          Suffix='';
        end
        
      case 'FullyPeriodic'
        if exist('FullyPeriodic')==0
          FullyPeriodic=false;
        end

      case 'PBC_Overlap_Tol'
        if exist('PBC_Overlap_Tol')==0
          PBC_Overlap_Tol=0.5;
        end
        
      case 'Stoich'
        if exist('Stoich')==0
          Stoich=true;
        end
        
      case 'ForceStoich'
        if exist('ForceStoich')==0
          ForceStoich=false;
        end
        

      case 'Plane_Shift'
        if exist('Plane_Shift')==0
          Plane_Shift=[0 0];
          if exist('FullyPeriodic')==1
            if FullyPeriodic==true;
              FullyPeriodic=false;
              warning(['Cannot shift planes with fully '...
                       'periodic boundary conditions'])
            end
          end
        end

      case 'Translate'
        if exist('Translate')==0
          Translate=[0 0 0 0 0 0]; 
        end

      case 'Summary'
        if exist('Summary')==0
          Summary=true;
        end
        
      case 'Summary_Name'
        if exist('Summary_Name')==0
          Summary_Name='Files.info';
        end
        
      case 'Angle_Range'
        if exist('Angle_Range')==0
          Angle_Range=[-1000 1000];
        end

      case 'Archive'
        if exist('Archive')==0
          Archive=true;
        end
        
      case 'Dir_Base'
        if exist('Dir_Base')==0;
          Dir_Base=pwd;
        end
        
      case 'AtomStyle'
        if exist('AtomStyle')==0;
          AtomStyle='atomic';
        end
        
      case 'TypeOffset'
        if exist('TypeOffset')==0
          TypeOffset=0;
        end
        
      case 'Header'
        if exist('Header')==0
          Header.Head1={'# Matlab Builder\n'};
          Header.Head2={''};
          Header.Head3={'# Matlab Builder\n'};
        end          
          
        if isfield(Header,'Head1')==0
          Header.Head1={'# Matlab Builder\n'};
        end        
        if isfield(Header,'Head2')==0
          Header.Head2={''};
        end   
        if isfield(Header,'Head3')==0
          Header.Head3={'# Matlab Builder\n'};
        end           
    end
  end

  clear varNames Variable suffix vars fun i Lc Lr Ac Ar Bc 
  clear tmp a b c ReqVars OptVars

  %% Pass valid variables back to functions
  w = whos;
  for a = 1:length(w)
    vars_pass.(w(a).name) = eval(w(a).name);
  end

  
  

end