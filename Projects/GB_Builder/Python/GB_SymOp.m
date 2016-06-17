function [ActiveSymmetry,NSymmetry] = GB_SymOp(Lattice,Basis,Symmetry,Orientation,Sym_Tol)

%% Setup variables
w = whos;
for a = 1:length(w) %1
vars.(w(a).name) = eval(w(a).name);
end %end 1
vars_pass=GB_Vars(vars,'GB_SymOp');

varNames=fieldnames(vars_pass);
for  i=1:length(varNames) %2
    eval([varNames{i} '=' 'vars_pass.(varNames{i});']);
end%end2

if Symmetry==true  %%3
  
  Rot=normr(Orientation)';

  %% Define possible symmetry operators - first time / save
  if exist('SymOpSearch.mat', 'file') ~= 2  %%4
  Reflection(:,:,1)=[1 0 0; 0 1 0; 0 0 1];
  Reflection(:,:,2)=[-1 0 0; 0 1 0; 0 0 1];
  Reflection(:,:,3)=[1 0 0; 0 -1 0; 0 0 1];
  Reflection(:,:,4)=[1 0 0; 0 1 0; 0 0 -1];
  Reflection(:,:,5)=[1 0 0; 0 -1 0; 0 0 -1];
  Reflection(:,:,6)=[-1 0 0; 0 1 0; 0 0 -1];
  Reflection(:,:,7)=[-1 0 0; 0 -1 0; 0 0 1];
  Reflection(:,:,8)=-[1 0 0; 0 1 0; 0 0 1];

  Rotation_x=@(theta)[ 1           0           0
               0          cos(theta) -sin(theta)
               0          sin(theta)  cos(theta)];
  Rotation_y=@(theta)[ cos(theta) 0           sin(theta)
               0          1           0
              -sin(theta) 0           cos(theta)];
  Rotation_z=@(theta)[ cos(theta) -sin(theta)  0
               sin(theta)  cos(theta)  0
               0           0           1];

  Theta=2*pi./[-1 -2 -3 -4 -6 1 2 3 4 6];
  [~,~,nReflection]=size(Reflection);
  nTheta=length(Theta);
  Sym=eye(3);
  for i=1:nReflection%5
    for j=1:nTheta%6
      for k=1:nTheta%7
        for m=1:nTheta%8
          SymTest=Reflection(:,:,i)*Rotation_x(Theta(j))*...
                  Rotation_y(Theta(k))*Rotation_z(Theta(m));
          nSym=length(Sym)/3;
          Compare=zeros(nSym*3,3);
          Compare(1:3:end,:)=ones(nSym,1)*SymTest(1,:);
          Compare(2:3:end,:)=ones(nSym,1)*SymTest(2,:);
          Compare(3:3:end,:)=ones(nSym,1)*SymTest(3,:);
          if sum(sum(reshape(sum((Sym-Compare).^2'),3,[]))<=1e-2)==0%9
                Sym=[Sym;SymTest];
          end%end 9
        end%end 8
      end %end 7
    end %end 6
  end%end 5
  save('SymOpSearch.mat','Sym');
  else  %% corallary 4
  load('SymOpSearch.mat')
  end  %% end 4
  nSym=length(Sym)/3;   
  
  %% Create expanded crytal for reference
  MaxN=3;
  n = 2*MaxN + 1;
  R = fullfact([n n n])-MaxN-1;
  nR = length(R);
  [nB,~]=size(Basis);   
  Expand=zeros(nR*nB,4);
  % Replicating each basis atom
  Tmp1=ones(nR,1);
  for m = 1:nB, %10
    Expand((m-1)*nR+1:m*nR,:)=[Tmp1*Basis(m,1) Tmp1*Basis(m,2:4)+R];
  end %end 10
  Ref=[Expand(:,1) Expand(:,2:end)*Lattice*Rot];
  
  [nRef,~]=size(Ref);
  Ref2=zeros(nRe f,nSym*4);

  %% Create test by applying symmetry opperation to crystal
  for n=1:nSym % 11
    Test(:,(n-1)*4+1:4*n)=[Basis(:,1) Basis(:,2:4)*Sym((n-1)*3+1:n*3,:)*Lattice*Rot];
  end %end 11

  Ref2(:,1:4:nSym*4)=Ref(:,1)*ones(1,nSym);
  Ref2(:,2:4:nSym*4)=Ref(:,2)*ones(1,nSym);
  Ref2(:,3:4:nSym*4)=Ref(:,3)*ones(1,nSym);
  Ref2(:,4:4:nSym*4)=Ref(:,4)*ones(1,nSy3
  parfor i=1:nB %12
     SymBasis(i,:)=(sum(reshape(sum(reshape((Ref2-ones(nRef,1)*Test(i,:)).^2'<=Sym_Tol,4,[])),nSym,[])'==4))>=1;
  end % end 12
  ActSym=find(sum(SymBasis)==nB);
  
  count=1;
  for n=ActSym %13
  ActiveSymmetry(:,:,count)=Sym((n-1)*3+1:n*3,:);
  count=count+1;
  end %end 13
  [~,~,NSymmetry]=size(ActiveSymmetry);

else %% corallary 3
  ActiveSymmetry=eye(3);
  NSymmetry=1;
end  %% end 3

end
