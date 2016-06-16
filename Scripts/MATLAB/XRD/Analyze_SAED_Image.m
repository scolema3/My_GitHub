% cd /home/scoleman/Research/Polymer_Branch/Cellulose/Publish/slices
% cd /data/Research/NC_Deformation/Foley/Analyze/10nm/10nm_e_comp_y_notext_notwin
close all; clear; clc

files = dir('*opt1*.png');


Rmin=0.1;   %% Inner Circle Radius (percentage of image size) 
Rmax=1.85;  %% Outter Circle Radius

if length(files)>0
for i=1:length(files)
  % Load Image
  Image=imread(files(i).name);
  
  % Mask middle region with circle
  a=size(Image,1);b=size(Image,2);c=(a+b)/2;
  Mask=zeros(a,b);
  [X,Y]=meshgrid(( (1:a)-a/2).^2 , ((1:a)-a/2).^2);
  Mask(X+Y<(Rmax*c)^2 & X+Y>(Rmin*c)^2)=1;
  Mask=uint8(Mask);
  Image=Image.*repmat(Mask,[1,1,3]);

  
  % Filtering -- creating BW image
  BW_option1=im2bw(Image(:,:,1));             % Isolate red regions
  BW_option2=im2bw(Image,graythresh(Image));  % Using whole image
  BW_option3=Image(:,:,1)-Image(:,:,2)-Image(:,:,3)>=10;         % Isolate red & fine_tune
  BW=BW_option3;
  
  
  % Find boundaries of BW image & computing properties
  [Edges,Label]=bwboundaries(BW);
  Stats = regionprops('table',Label,'all');
  Centers1 = Stats.Centroid;
  
  Scene(i).Image=Image;
  Scene(i).BW=BW;
  Scene(i).Filename=files(i).name;
  Scene(i).Edges=Edges;
  Scene(i).Centers=Centers1;
  Scene(i).Stats=Stats;
  
end
%%

clc; close all
PlotScene(Scene)

end


%%
h1=figure('color','w','position',[2561 546 1920 948]); hold on

Dmax=100;
Dmin=0;

cmap=colormap(lines(length(Scene)));
figure(h1)
for i=1:length(Scene)
  % Draw edges
  for k = 1:length(Scene(i).Edges)
    boundary = Scene(i).Edges{k};
    plot(boundary(:,2), boundary(:,1), 'color',cmap(i,:), 'LineWidth', 2)
  end
  
  % Plot centroids
  Centers = Scene(i).Centers;
  text(Centers(:,1),Centers(:,2),num2str(i),'color',cmap(i,:))
end



for i=1:length(Scene)-1
  
  j=i+1;
  X = KDTreeSearcher(Scene(i).Centers);
  Y = Scene(j).Centers;
  [IDX D]=knnsearch(X,Y);

  Centers1 = Scene(i).Centers;
  Centers2 = Scene(j).Centers;  

  Origin1=[size(Scene(i).Image,1)/2,size(Scene(i).Image,2)/2];
  Origin2=[size(Scene(j).Image,1)/2,size(Scene(j).Image,2)/2];
  
  Vec1=Centers1-ones(size(Centers1,1),1)*Origin1;
  Vec2=Centers2-ones(size(Centers2,1),1)*Origin2;
  Dcut=find(D>=Dmin & D<=Dmax);
  Pairs(i).Id=[IDX, [1:length(Y)]'];
  Pairs(i).Vec =[ X.X(IDX(Dcut),1) X.X(IDX(Dcut),2) Y(Dcut,1) Y(Dcut,2)];
  
  for k = 1:length(Dcut)
    arrow([ X.X(IDX(Dcut(k)),1) X.X(IDX(Dcut(k)),2)], [Y(Dcut(k),1) Y(Dcut(k),2)])
    Ang(i,k)=acos(dot(Vec1(IDX(Dcut(k)),:),Vec2(Dcut(k),:))/(norm(Vec1(IDX(Dcut(k)),:))*norm(Vec2(Dcut(k),:))))*180/pi;
  end
  
end


set(gca, 'YDir', 'reverse');
[a,b,~]=size(Scene(i).Image);
axis([0 a 0 b]);
axis square


%%
clc; close all


h1=figure('color','w','position',[2561 546 1920 948]); hold on

Dmax=100;
Dmin=0;

cmap=colormap(lines(length(Scene)));
figure(h1)
for i=1:length(Scene)
  % Draw edges
  for k = 1:length(Scene(i).Edges)
    boundary = Scene(i).Edges{k};
    plot(boundary(:,2), boundary(:,1), 'color',cmap(i,:), 'LineWidth', 2)
  end
  
  % Plot centroids
  Centers = Scene(i).Centers;
  text(Centers(:,1),Centers(:,2),num2str(i),'color',cmap(i,:))
end

set(gca, 'YDir', 'reverse');
[a,b,~]=size(Scene(i).Image);
axis([0 a 0 b]);
axis square

pause

[x y]=ginput(2);

xmin=min(x);
xmax=max(x);
ymin=min(y);
ymax=max(y);

h1=figure('color','w','position',[2561 546 1920 948]); hold on

for i=1:length(Scene)
  select=Scene(i).Centers(:,1) > xmin & Scene(i).Centers(:,1) < xmax & Scene(i).Centers(:,2) > ymin  & Scene(i).Centers(:,2) <ymax;
  SelCenter(i,1)=mean(Scene(i).Centers(select,1))
  SelCenter(i,2)=mean(Scene(i).Centers(select,2))
  
  SelCenterNorm=normr(SelCenter(i,:));
  atan(SelCenterNorm(2)/SelCenterNorm(1))*180/pi;
  AngleDegMat(i)=acos(SelCenter(i,1)/norm(SelCenter(i,:)))*180/pi;
  AngleDegMat2(i)=atan(SelCenterNorm(2)/SelCenterNorm(1))*180/pi;
end

plot(AngleDegMat,'o')

AngleDegMat
AngleDegMat2

%%  
