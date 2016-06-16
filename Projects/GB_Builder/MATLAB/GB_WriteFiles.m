function GB_WriteFiles(File_Base_id,Write_Style,AtomData,Species,Masses,Corners,Header,AtomStyle,Archive,Dir_Base,Verbose,Bulk)

%% Setup variables

  w = whos;
  for a = 1:length(w)
    vars.(w(a).name) = eval(w(a).name);
  end
  
  if exist('Bulk')==0
    Bulk=0;
  end
  
  vars_pass=GB_Vars(vars,'GB_WriteFiles');

  varNames=fieldnames(vars_pass);
  for  i=1:length(varNames)
    eval([varNames{i} '=' 'vars_pass.(varNames{i});']);
  end

  xlo=Corners(1);
  xhi=Corners(2);
  ylo=Corners(3);
  yhi=Corners(4);
  zlo=Corners(5);
  zhi=Corners(6);

  nAtoms=length(AtomData(:,1));
  Types=unique(AtomData(:,1));
  nTypes=length(Types);
  AtomIDs=[1:nAtoms]';


  %% Writing to LAMMPS data file
  if ismember(1,Write_Style)
    
    if Bulk==0
      fid=fopen([File_Base_id '.data'],'w');
    else 
      fid=fopen([File_Base_id '.Data'],'w');
    end

    for i=1:length(Header.Head1)
      fprintf(fid,Header.Head1{:,i});
    end
    
    if strcmpi(AtomStyle,'atomic')
    fprintf(fid, '\n%d atoms\n',nAtoms);
    fprintf(fid, '%d atom types\n\n',nTypes);
    fprintf(fid,'%f %f xlo xhi\n',xlo,xhi);
    fprintf(fid,'%f %f ylo yhi\n',ylo, yhi);
    fprintf(fid,'%f %f zlo zhi\n',zlo,zhi);
    fprintf(fid,'\nAtoms\n\n');
    fprintf(fid,'%d %d %f %f %f\n',[AtomIDs AtomData]');
    fprintf(fid,'\nMasses\n\n');
    for s=1:length(Types)
      fprintf(fid,'%d %f\n',Types(s),Masses(s));
    end
    fclose(fid);

    elseif strcmpi(AtomStyle,'charge')
    fprintf(fid, '\n%d atoms\n',nAtoms);
    fprintf(fid, '%d atom types\n\n',nTypes);
    fprintf(fid,'%f %f xlo xhi\n',xlo,xhi);
    fprintf(fid,'%f %f ylo yhi\n',ylo, yhi);
    fprintf(fid,'%f %f zlo zhi\n',zlo,zhi);
    fprintf(fid,'\nAtoms\n\n');
    fprintf(fid,'%d %d %f %f %f %f\n',[AtomIDs AtomData(:,1) zeros(nAtoms,1) AtomData(:,2:4)]');
    fprintf(fid,'\nMasses\n\n');
    for s=1:length(Types)
      fprintf(fid,'%d %f\n',Types(s),Masses(s));
    end
    fclose(fid);
    end
    
      
    
    if Verbose==true
      fprintf(['#\n# Data File Written: ' File_Base_id '.data\n']);
    end
  end


  %% Writing to LAMMPS dump file
  if ismember(2,Write_Style)
    fid=fopen([File_Base_id '.dump'],'w');
    for i=1:length(Header.Head2)
      fprintf(fid,Header.Head2{:,i});
    end    
    fprintf(fid,'ITEM: TIMESTEP\n');
    fprintf(fid,'0\n');
    fprintf(fid,'ITEM: NUMBER OF ATOMS\n');
    fprintf(fid,'%d\n',nAtoms);
    fprintf(fid,'ITEM: BOX BOUNDS pp pp pp\n');
    fprintf(fid,'%f %f\n',xlo,xhi);
    fprintf(fid,'%f %f\n',ylo,yhi);
    fprintf(fid,'%f %f\n',zlo,zhi);
    fprintf(fid,'ITEM: ATOMS id type x y z\n');
    fprintf(fid,'%d %d %f %f %f\n',[AtomIDs AtomData]');
    fclose(fid);
    if Verbose==true
      fprintf(['#\n# Dump File Written: ' File_Base_id '.dump\n'])
    end
  end
 
  %% Writing to car file for Materials Studio
  if ismember(3,Write_Style)
    for s=1:nTypes
      Label(ismember(AtomData(:,1),Types(s)),1)=Species(s);
      clear Tmp1
    end
    
    fid=fopen([File_Base_id '.car'],'w');  
    fprintf(fid,'!BIOSYM archive 3\n');
    fprintf(fid,'PBC=ON\n');
    for i=1:length(Header.Head3)
      fprintf(fid,Header.Head3{:,i});
    end 
    fprintf(fid,'!DATE\n');
    fprintf(fid,'PBC %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f (P1)\n',...
            xhi-xlo, yhi-ylo, zhi-zlo,...
            90.000, 90.000,90.00);
    for i=1:nAtoms
      tmp=length(char(Label(i))');
      format1=['%s%d % ' num2str(18-tmp) '.9f %14.9f %14.9f XXXX 1'];
      format2=['      xx      %s %' num2str(8-tmp) '.3f\n'];
      format=[format1 format2];
      fprintf(fid,format,char(Label(i))',AtomData(i,1),AtomData(i,2),...
                  AtomData(i,3), AtomData(i,4),char(Label(i))',0);
    end
    fprintf(fid,'end\nend\n');
    if Verbose==true
      fprintf(['#\n# Car File Written: ' File_Base_id '.car\n']  )    
    end
    fclose(fid);
    clear Label
  end
  

  %% Prepare Output For Transfer
  if Archive==true
    if ismember(1,Write_Style)
      movefile([File_Base_id '.data'],[Dir_Base '_Data'])
    end
    if ismember(2,Write_Style)
      movefile([File_Base_id '.dump'],[Dir_Base '_Dump'])
    end
    if ismember(3,Write_Style)
      movefile([File_Base_id '.car'],[Dir_Base '_Car'])
    end
  end

  
end