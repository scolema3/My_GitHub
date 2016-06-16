close all;
clear all;
clc

% Intensity threshold for keeping data.
Threshold=0;

% VTK data file of intensities.
% cd /home/scoleman/Research/Polymer/PACM/2D
fid=fopen('gs5_compression_All_0_saed.0.vtk');

% Read in data and dimensions.
head=textscan(fid,'%s %*[^\n]',4);
dimensions=textscan(fid,'%*s %f %f %f',1);
dimensions=cell2mat(dimensions);
aspect_ratio=textscan(fid,'%*s %f %f %f',1);
aspect_ratio=cell2mat(aspect_ratio);
origin=textscan(fid,'%*s %f %f %f',1);
origin=cell2mat(origin);
point_data=textscan(fid,'%*s %f',1);
point_data=cell2mat(point_data);
head=textscan(fid,'%s %*[^\n]',2);

% Create Index
Index=fullfact([dimensions(1) dimensions(2) dimensions(3)]);
Index=Index-ones(length(Index),1)*[(dimensions(1)+1)/2 (dimensions(2)+1)/2 (dimensions(3)+1)/2];
Index(:,1)=Index(:,1)*aspect_ratio(1);
Index(:,2)=Index(:,2)*aspect_ratio(2);
Index(:,3)=Index(:,3)*aspect_ratio(3);

% Read in data by chunks
formatSpec = '%f';
BlockSize=100000;
Intensity=[];  % matrix for: [x y z intensity]
count=0;
while ~feof(fid)
  	tmp = textscan(fid,formatSpec,BlockSize);
    index = [(count*BlockSize)+1:(count)*BlockSize+length(tmp{1})]';
    rawdata = [Index((count*BlockSize)+1:(count)*BlockSize+length(tmp{1}) ,:) tmp{:}];
    keep = rawdata(:,4)>=Threshold;
    Intensity = [Intensity;rawdata(keep,:)];
    count = count+1;
end


fclose(fid);


%% testing visualization

plot3(Intensity(:,1),Intensity(:,2),Intensity(:,3),'.')

slice=find(Intensity(:,1)==0);

scatter3(Intensity(slice,2),Intensity(slice,3),log10(Intensity(slice,4)),[],log10(Intensity(slice,4)))
axis square