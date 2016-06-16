 function [AtomData,Overlap_Id]=GB_WrapPBC(AtomData,MaxEdge,PBC_Overlap_Tol)

  %% Setup varialbes
  w = whos;
  for a = 1:length(w)
    vars.(w(a).name) = eval(w(a).name);
  end

  vars_pass=GB_Vars(vars,'GB_WrapPBC');

  varNames=fieldnames(vars_pass);
  for  i=1:length(varNames)
    eval([varNames{i} '=' 'vars_pass.(varNames{i});']);
  end

  Overlap_Id=[];
  buff=PBC_Overlap_Tol;
  
  Edge1=[];
  Edge2=[];
  if MaxEdge(1)>0  % x-axis
    Edge1=find(AtomData(:,2)<=min(AtomData(:,2))+buff);
    n1=length(Edge1);
    Edge2=find(AtomData(:,2)>=max(AtomData(:,2))-buff);
    n2=length(Edge2);
    for m=Edge1'
      delta=(ones(n2,1)*(AtomData(m,2:4)+[MaxEdge(1) 0 0])-AtomData(Edge2,2:4));
      check=min(sqrt(sum(delta.^2,2)));
      if check<=PBC_Overlap_Tol
        Overlap_Id=[Overlap_Id; m];
      end
    end
  end
  
  Edge1=[];
  Edge2=[];  
  if MaxEdge(2)>0 % y-axis
    Edge1=find(AtomData(:,3)<=min(AtomData(:,3))+buff);
    n1=length(Edge1);
    Edge2=find(AtomData(:,3)>=max(AtomData(:,3))-buff);
    n2=length(Edge2);
    for m=Edge1'
      delta=(ones(n2,1)*(AtomData(m,2:4)+[0 MaxEdge(2) 0 ])-AtomData(Edge2,2:4));
      check=min(sqrt(sum(delta.^2,2)));
      if check<=PBC_Overlap_Tol
        Overlap_Id=[Overlap_Id; m];
      end
    end  
  end
  
  Edge1=[];
  Edge2=[];  
  if MaxEdge(3)>0 % z-axis
    Edge1=find(AtomData(:,4)<=min(AtomData(:,4))+buff);
    n1=length(Edge1);
    Edge2=find(AtomData(:,4)>=max(AtomData(:,4))-buff);
    n2=length(Edge2);
    for m=Edge1'
      delta=(ones(n2,1)*(AtomData(m,2:4)+[0 0 MaxEdge(3)])-AtomData(Edge2,2:4));
      check=min(sqrt(sum(delta.^2,2)));
      if check<=PBC_Overlap_Tol
        Overlap_Id=[Overlap_Id; m];
      end
    end
  end
  
  Overlap_Id=unique(Overlap_Id);
  AtomData(Overlap_Id,:)=[];
  
end