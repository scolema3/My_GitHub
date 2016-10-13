%% Initialization
clc;clear;close all

% column 1= 2 theta (degrees)
% column 2= intentisty (normalized Imax=100)
ExpData=[   
   10.9233  100.0000
   14.2830   94.9431
   10.7528   92.6100
   16.2784   89.0624
   16.5512   88.3934
   15.1528   78.9204
   15.5962   78.4505
   13.1063   73.8211
    9.4908   58.0397
    8.9450   43.7658
    8.4163   29.2142
    6.5574   24.1237
    3.1977   14.1050
    2.3279   13.6405];
ExpData=[  
  10.9233  100.0000
   10.7528   81.4938
   14.2830   42.1657
    9.4908   35.4058
    8.9280   19.6745
    6.5400   26.1632
   16.2784   22.5120
   16.5512   19.5285
   13.1063   11.6837
    8.4163    7.4066
   13.8900    0.3964
   15.1528    3.3799
   15.5962    2.0238
   15.0200    8.3664
   15.2900    0.7720
    9.6780   12.7686
   10.1600    9.1383
   10.3800    8.3664
   11.2300    1.8151
   12.1700    0.8763
   12.5800    1.3770
   13.4100    1.0015
   16.8900    1.7317
   17.0100    1.0432]    
    
% XRD files 
files = dir('*.xrd');

% Data Filter Options

scale=1;               % Scale results to maximum int value in all XRD files
theta_window=0.06;     % Window for filter/smooth (2theta)

% Data Fitting Options

Npeaks=25;              % Number of peaks to fit
min_peak_dist=0.15;       % Minimum value for placement of consecutive peaks (2theta)

fit=0;                 % 0/1  : Perform fitting or not
whole=0;               % 0/1  : Perform whole pattern fitting or not
manual=0;              % 0/1  : Turn off manual peak placement
peakshape=13;          % 1-33 : See below
background=0;          % 0-3  : See below 

NumTrials_ind=1;       % Number of starting  for each individual fit
NumTrials_whole=1;     % Trials to fit for whole peak

addPeaks=0;            % 0/1  : Add peak to improve goodness of fit
rmPeaks=0;             % 0/1  : Remove peaks to improve goodness of fit 
plot_on=0;             % 0/1  : Plot off/on

% Display Results

display_result=1;      % 0/1  : Display all plots
save_result=1;         % 0/1  : Save plot


% Save location / intensity of peaks 

save_peak_info=1;
save_peak_file='XRD_peak_data.txt';


% Analysis Options 

lambda=1.541838;       % Radiation Wavelength (angstroms)
K=1;                   % Williamson-Hall factor


% "peakshape" specifies the peak shape of the model: (1=Gaussian (default),
% 2=Lorentzian, 3=logistic distribution, 4=Pearson, 5=exponentionally
% broadened Gaussian; 6=equal-width Gaussians; 7=Equal-width Lorentzians;
% 8=exponentionally broadened equal-width Gaussian, 9=exponential pulse,
% 10=up-sigmoid (logistic function), 11=Fixed-width Gaussian,
% 12=Fixed-width Lorentzian; 13=Gaussian/ Lorentzian blend; 14=Bifurcated
% Gaussian, 15=Breit-Wigner-Fano, 16=Fixed-position Gaussians;
% 17=Fixed-position Lorentzians; 18=exponentionally broadened Lorentzian;
% 19=alpha function; 20=Voigt profile; 21=triangular; 22=multiple shapes;
% 23=down-sigmoid; 25=lognormal; 26=slope; 27=Gaussian first derivative; 
% 28=polynomial; 29=piecewise linear; 30=variable-alpha Voigt; 31=variable
% time constant ExpGaussian; 32=variable Pearson; 33=independently-variable
% Gaussian/Lorentzian blend


% 'background' sets the baseline correction mode:
% background=0 (default) does not subtract baseline from data segment; 
% background=1 interpolates a linear baseline from the edges of the data 
%              segment and subtracts it from the signal (assumes that the 
%              peak returns to the baseline at the edges of the signal); 
% background=2 is like mode 1 except that it computes a quadratic curved baseline; 
% background=3 compensates for a flat baseline without reference to the 
%              signal itself (best if the peak does not return to the
%              baseline at the edges of the signal).


% Collect information
clear XRD
maxY=0;
for i = 1:length(files)
  fid = fopen(files(i).name,'r');
  tmp = textscan(fid,'%*f %f %f %*f', 'headerlines',4);
  fclose(fid);
  XRD{i}.filename = files(i).name(1:end-4);
  XRD{i}.rawtheta = tmp{1};
  XRD{i}.rawintensity = tmp{2};
  x = tmp{1}; y = tmp{2};
  XRD{i}.theta = x;
  XRD{i}.intensity = y;
  
  % Filter & Smooth Data
  windowsize=round(theta_window/((max(x)-min(x))/length(x)));
  b=(1/windowsize)*ones(1,windowsize);
  a=1;
  y1=filter(b,a,y);
  y2=smooth(x,y1,windowsize,'lowess');

  if exist('ThetaRange')
    select=[x>=ThetaRange(1) & x<=ThetaRange(2)];
    x=x(select);
    y=y(select);
    y1=y1(select);
    y2=y2(select);
    XRD{i}.theta = x;
    XRD{i}.intensity = y;
  else 
    ThetaRange(1)=min(XRD{i}.rawtheta);
    ThetaRange(2)=max(XRD{i}.rawtheta);
  end
  
  XRD{i}.windowsize = windowsize;
  XRD{i}.theta_f = x;
  XRD{i}.intensity_f = y1;
  XRD{i}.theta_fs = x;
  XRD{i}.intensity_fs = y2;
  
  maxY=max(max(y2),maxY);
end

if scale==1
  Scale=100/maxY;
else
  Scale=1;
end

clear x y tmp fid files i

if background==0
  Fbaseline=@(coef,xx)0;
elseif background==1
  Fbaseline=@(coef,xx)polyval(coef,xx);
elseif background==2
  Fbaseline=@(coef,xx)polyval(coef,xx);
elseif background==3
  Fbaseline=@(coef,xx)coef;
end

% Filter and smooth data

if plot_on==1;
  f1=figure('Name','Bin,Smooth,Fit','Color','w','Position',[2560 650 1250 750]);
end

if save_peak_info==1;
  fid=fopen(save_peak_file,'w');
end

for i =1:length(XRD)
  x = XRD{i}.theta;
  y = XRD{i}.rawintensity*Scale;
  y1 = XRD{i}.intensity_f*Scale;
  y2 = XRD{i}.intensity_fs*Scale;
  
  
  if plot_on==1;
    clf    
    plot(x,y,'.k'); hold on
    plot(x,y1,'.g'); hold on
    plot(x,y2,'-b');
    hl=legend(XRD{i}.filename,'filtered','smoothed');
    set(hl,'interpret','none','box','off')
    xlabel('2\theta (degree)')
    ylabel('Intensity (%)')
    title('Bin,Smooth,Fit')
    if scale==1
      axis([ThetaRange(1) ThetaRange(2) -5 105])
    end    
    drawnow
  end

  %% Find peak maxima
  [pks,locs] = findpeaks(y2,...
               'MinPeakHeight',max(y2)*0.015,...
               'MinPeakDistance',round(min_peak_dist/((max(x)-min(x))/length(x))),...
               'npeaks',Npeaks,'sortstr','descend');
             
   XRD{i}.PeakLocInt=[x(locs) y2(locs)];
   
   
   if save_peak_info==1;
     fprintf(fid,'\n%s\n',XRD{i}.filename);
     fprintf(fid,'%f %f\n',XRD{i}.PeakLocInt');
   end

  if plot_on==1;
    plot(x(locs),pks,'*r')
    drawnow
  end

if fit~=0;    
    
    if manual==1;
      
      f1=figure('Name','Bin,Smooth,Fit','Color','w','Position',[2560 650 1250 750])
      plot(x,y,'.k'); hold on
      plot(x,y1,'.g'); hold on
      plot(x,y2,'-b');
      hl=legend(XRD{i}.filename,'filtered','smoothed');
      set(hl,'interpret','none','box','off')
      xlabel('2\theta (degree)')
      ylabel('Intensity (%)')
      title('Bin,Smooth,Fit')
      if scale==1
        axis([ThetaRange(1) ThetaRange(2) -5 105])
      end
      drawnow
      
      title({'LEFT CLICK TO SELECT PEAKS '; 'RIGHT CLICK TO DE-SELECT PEAKS';' PRESS ENTER WHEN FINSIHED'});
      [x,y,B] = ginput;
      
      PeakAdd=[];
      PeakRemove=[];
      % Move autmoated peaks
      for j=1:length(x)
        % Remove peaks
        tol=1.0^2;
        select=(x(j)*ones(size(locs))-x(locs)).^2<=tol;
        PeakRemove=[PeakRemove;find(select==1)];
        
        % Add back peaks with right click
        if B(j)==1
          [tmpint tmploc]=max((((XRD{i}.theta_fs)-x(j)).^2<=tol).*XRD{i}.intensity_fs);
          PeakAdd=[PeakAdd;tmploc tmpint];
        end
      end
      
      PeakRemove=unique(PeakRemove);
      pks(PeakRemove)=[];
      locs(PeakRemove)=[];
      
      if prod(size(PeakAdd))>1
        locs=[locs;PeakAdd(:,1)];
        pks=[pks; PeakAdd(:,2)];
      end
      
      clf
      plot(x,y,'.k'); hold on
      plot(x,y2,'-b');
      plot(x(locs),pks,'*r')
      hl=legend(XRD{i}.filename);
      set(hl,'interpret','none','box','off')
      xlabel('2\theta  (degree)')
      ylabel('Intensity  (%)')
      title('Bin,Smooth,Fit')
      drawnow
    end
    
    
    XRD{i}.pks=pks;
    XRD{i}.locs=locs;
    
    % Find peak bounds
    Npeaksi=length(XRD{i}.pks);
    y=XRD{i}.intensity_fs*Scale;

    [locsort,n1] = sort(locs(1:Npeaksi));
    pksort = pks(n1);
    clear xbound
    for j = 1:Npeaksi+1
        if j == 1
            n=1:locsort(j);
            xbound(j) = round(mean(n));
        elseif j == Npeaksi+1
            n = locsort(j-1):length(x);
            xbound(j) = round(mean(n));
        else
            n = locsort(j-1):locsort(j);
            xbound(j) = round(mean(n));
        end
    end
    if plot_on==1;
      plot(x(xbound),y(xbound),'*b')
      drawnow
    end
    
    % Fit individual peaks to smoothed data
    ExtraPeak=0; J=0; Draw=1;
    for j = 1:Npeaksi
      J=J+1;
      n = xbound(j):xbound(j+1);
      signal = [x(n),y(n)];
      center = x(locsort(j));
      NumPeaks = 1;
      StartGuess=0;
      extra=0;
      [extra,fval] = fminsearch(@(j) peakfitfun(j,x(n),y(n),NumPeaks,peakshape,NumTrials_ind,background,StartGuess),50);
      [FitResults,GOF,baseline,coeff,BestStart,xi,yi]=peakfit([x(n),y(n)],0,0,NumPeaks,peakshape,extra,NumTrials_ind,StartGuess,background,0,0);
      XRD{i}.numpeaks(J) = NumPeaks;
      XRD{i}.peakshape(J) = peakshape;
      XRD{i}.extra_fs(J) = extra;
      XRD{i}.baseline_fs(J,:)= baseline;
      XRD{i}.FitResults_fs_i(J,:) = FitResults;
      XRD{i}.rsquare_fs_i(J) = GOF(2);
      XRD{i}.Fit_fs_i_xi(J,:)=xi;
      XRD{i}.Fit_fs_i_yi(J,:)=yi;      

      % Add Peaks to improve GOF
      GOFmin=0.95;
      if GOF(2)<=GOFmin & addPeaks
        NumPeaks = 2;
        StartGuess=0;
        [extra,fval] = fminsearch(@(j) peakfitfun(j(1:2),x(n),y(n),NumPeaks,[peakshape peakshape],NumTrials_ind,background,StartGuess),[50 50]);
        [FitResults,GOF,baseline,coeff,BestStart,xi,yi]=peakfit([x(n),y(n)],0,0,NumPeaks,[peakshape peakshape],extra,NumTrials_ind,StartGuess,background,0,0);
        XRD{i}.numpeaks(J:J+1) = [NumPeaks NumPeaks];
        XRD{i}.peakshape(J:J+1) = [peakshape peakshape];
        XRD{i}.extra_fs(J:J+1) = extra;
        XRD{i}.baseline_fs(J:J+1,:)= [baseline;baseline];
        XRD{i}.FitResults_fs_i(J:J+1,:) = FitResults;
        XRD{i}.rsquare_fs_i(J:J+1) = [GOF(2) GOF(2)];
        XRD{i}.Fit_fs_i_xi(J:J+1,:)=[xi;xi];
        XRD{i}.Fit_fs_i_yi(J:J+1,:)=yi;
        ExtraPeak=ExtraPeak+1;
        J=J+1;
      end
      
      % Remove peaks that are fitting to background or individual
      % data point
      Widtol=[0.1 15]; % Peak Width tolerence
      BG=(FitResults(:,4)<=Widtol(1) | FitResults(:,4)>=Widtol(2) | FitResults(:,2)<=ThetaRange(1) | FitResults(:,2)>=ThetaRange(2)) & rmPeaks;
      if sum(BG)>0 & Npeaksi>2
      XRD{i}.numpeaks(J) = [];        
        XRD{i}.peakshape(J) = [];
        XRD{i}.extra_fs(J) = [];
        XRD{i}.baseline_fs(J,:)= baseline;
        XRD{i}.FitResults_fs_i(J,:) = [];
        XRD{i}.rsquare_fs_i(J) = [];
        XRD{i}.Fit_fs_i_xi(J,:)=[];
        XRD{i}.Fit_fs_i_yi(J,:)=[];
        J=J-1;
        Draw=0;
      end
      
      if plot_on==1;
        if Draw~=0
          if NumPeaks>1
            plot(xi,sum(yi)+Fbaseline(baseline,xi),'--m')
          else
            plot(xi,yi+Fbaseline(baseline,xi),'--m')
          end
          drawnow
        end
      end
    end
  
    Draw=1;
    % Fit individual peaks to orig data
    x=XRD{i}.theta_f;
    y=XRD{i}.intensity_f*Scale;
    J=0;
    for j = 1:Npeaksi
      J=J+1;
      
      if XRD{i}.numpeaks(J)==2
        J2=J+1;
      else
        J2=J;
      end
      
      n = xbound(j):xbound(j+1);
      [~,I1]=min(abs(XRD{i}.theta_f-min(x(n))));
      [~,I2]=min(abs(XRD{i}.theta_f-max(x(n))));
      n=[I1:I2];

      signal = [x(n),y(n)];
      center = x(locsort(j));
      
      NumPeaks = XRD{i}.numpeaks(J);
      PeakShape= XRD{i}.peakshape(J:J2);
      
      StartGuess=reshape([XRD{i}.FitResults_fs_i(J:J2,2),XRD{i}.FitResults_fs_i(J:J2,4)]',1,[]);
      [extra,fval] = fminsearch(@(j) peakfitfun(j,x(n),y(n),NumPeaks,PeakShape,NumTrials_ind,background,StartGuess),XRD{i}.extra_fs(J:J2));      
      [FitResults,GOF,baseline,coeff,BestStart,xi,yi]=peakfit([x(n),y(n)],0,0,NumPeaks,PeakShape,extra,NumTrials_ind,StartGuess,background,0,0);
      XRD{i}.peakshape_f(J:J2) = PeakShape;
      XRD{i}.extra_f(J:J2) = extra;
      XRD{i}.baseline_f(J,:)= baseline;
      XRD{i}.baseline_f(J2,:)= baseline;
      XRD{i}.FitResults_f_i(J:J2,:) = FitResults;
      XRD{i}.rsquare_f_i(J:J2) = GOF(2);
      XRD{i}.Fit_f_i_xi(J,:)=xi;
      XRD{i}.Fit_f_i_xi(J2,:)=xi;      
      XRD{i}.Fit_f_i_yi(J:J2,:)=yi;      

      J=J2;
      if plot_on==1;
        if NumPeaks>1
          plot(xi,sum(yi)+Fbaseline(baseline,xi),'--g')
        else
          plot(xi,yi+Fbaseline(baseline,xi),'--g')
        end
        drawnow
      end
    end
    
    if whole==1
    Npeaksi=Npeaksi+ExtraPeak;
    
    % Refine Fit Whole Pattern
    peakshapes=ones(1,Npeaksi)*peakshape;
    extras=XRD{i}.extra_fs;
    StartGuess = reshape([XRD{i}.FitResults_fs_i(:,2),XRD{i}.FitResults_fs_i(:,4)]',1,[]);        
    
    % Fit to binned and smoothed data
    x=XRD{i}.theta_fs;
    y=XRD{i}.intensity_fs*Scale;
    [FitResults_fs,GOF_fs,baseline_fs,coeff,BestStart,xi_fs,yi_fs]=peakfit([x,y],0,0,Npeaksi,peakshapes,extras,NumTrials_whole,StartGuess,background,0,0);
    if plot_on==1;
      plot(xi_fs,sum(yi_fs),'-y')
      drawnow
    end  

    peakshapes=ones(1,Npeaksi)*peakshape;
    extras=XRD{i}.extra_f;
    StartGuess = reshape([XRD{i}.FitResults_f_i(:,2),XRD{i}.FitResults_f_i(:,4)]',1,[]);        
    
    x=XRD{i}.theta_f;
    y=XRD{i}.intensity_f*Scale;
    [FitResults_f,GOF_f,baseline_f,coeff,BestStart,xi_f,yi_f]=peakfit([x,y],0,0,Npeaksi,peakshapes,extras,NumTrials_whole,StartGuess,background,0,0);
    
    if plot_on==1;
      plot(xi_f,sum(yi_f)+Fbaseline(baseline,xi_f),'-c')
      drawnow
    end

    XRD{i}.Fit_fs_peakshape = peakshapes;
    XRD{i}.Fit_fs_extra = extras;
    XRD{i}.Fit_fs_baseline = baseline_fs;
    XRD{i}.Fit_fs_FitResults=FitResults_fs;
    XRD{i}.Fit_fs_rsquare=GOF_fs(2);
    XRD{i}.Fit_fs_xi=xi_fs;
    XRD{i}.Fit_fs_yi=yi_fs;

    XRD{i}.Fit_f_peakshape = peakshapes;
    XRD{i}.Fit_f_extra = extras;
    XRD{i}.Fit_f_baseline = baseline_f;
    XRD{i}.Fit_f_FitResults=FitResults_f;
    XRD{i}.Fit_f_rsquare=GOF_f(2);
    XRD{i}.Fit_f_xi=xi_f;
    XRD{i}.Fit_f_yi=yi_f;
  
     clear  h hl x1 y1 y2
    end
  end
end

if save_peak_info==1;
  fclose(fid);
end

% Display Results - NON-Fitted
%%

if display_result==1
  close all; clc
  figure('name','XRD','Color','w', 'position',[1000 50 1500 1500])
  cmap=colormap(lines(length(XRD)));
  
  linw=3;
  fs0=18;
  fs1=24;
  fs2=24;
  fname='Arial';
  
  doffset=.95;
  crop=0.025;
  xloc=.60;
  yloc=0.7*doffset;
  
  
  % Adjust order of line profiles using ID
  ID=[1:length(XRD)];
  LabelID=[ID];
  
  
  for i=1:length(ID);
    hold on
    if i==1
      offset=0;
      plot(XRD{ID(i)}.theta_f,XRD{ID(i)}.intensity_fs/max(XRD{ID(i)}.intensity_fs)+offset,'color',cmap(i,:),'linewidth',linw);
    else
      offset=offset+doffset;
      plot(XRD{ID(i)}.theta_f,XRD{ID(i)}.intensity_fs/max(XRD{ID(i)}.intensity_fs)+offset,'color',cmap(i,:),'linewidth',linw);
    end
    text((max(XRD{ID(i)}.theta_f)-min(XRD{ID(i)}.theta_f))*xloc+min(XRD{ID(i)}.theta_f),offset+yloc,XRD{ID(i)}.filename(1:end),'interpret','none','color',cmap(i,:),'fontsize',fs1)
  end
  
  axis([min(XRD{1}.theta_f) max(XRD{1}.theta_f) -crop offset+1+crop])
  xlabel(['2',char(952),char(176)],'Interpreter','tex','fontsize',fs2,'Fontname',fname);
  ylabel(['Intensity'],'Interpreter','tex','fontsize',fs2,'Fontname',fname)
  set(gca,'fontsize',fs0,'ytick',[],'xminortick','on','tickdir','out','linewidth',2,'Fontname',fname)
  
  if save_result==1
    Image_Name=XRD{1}.filename(1:end);
    expfig=['export_fig ' Image_Name ' -eps -pdf -png']
    eval(expfig)
  end
end



%% Display Results - Fitted
% 
% if display_result==1 && fit==1
%   % XRD profiles
%   if background==0
%     Fbaseline=@(coef,xx)0;
%   elseif background==1
%     Fbaseline=@(coef,xx)polyval(coef,xx);
%   elseif background==2
%     Fbaseline=@(coef,xx)polyval(coef,xx);
%   elseif background==3
%     Fbaseline=@(coef,xx)coef;
%   end
%   
% 
%   figure('name','XRD','Color','w', 'position',[1000 50 1500 1500])
%   cmap=colormap(lines(length(XRD)));
%   
%   linw=3;
%   fs0=18;
%   fs1=24;
%   fs2=24;
%   fname='Arial';
%   
%   doffset=.95;
%   crop=0.025;
%   xloc=.60;
%   yloc=0.7*doffset;
%   
%   
%   % Adjust order of line profiles using ID
%   ID=[1:length(XRD)];
%   LabelID=[ID];
%   
%   
%   for i=1:length(ID);
%     hold on
%     if i==1
%       offset=0;
%       if size(XRD{1}.Fit_fs_i_yi,1)==1
%         y=(XRD{i}.Fit_fs_i_yi+Fbaseline(XRD{i}.baseline_fs,XRD{i}.Fit_fs_i_xi));
%       else
%         y=sum(XRD{i}.Fit_fs_yi+Fbaseline(XRD{i}.Fit_fs_baseline,XRD{i}.Fit_fs_i_xi));
%       end
%       plot(XRD{i}.Fit_fs_i_xi,(y-min(y))/max(y)+offset,'color',cmap(i,:),'linewidth',linw);
%     else
%       offset=offset+doffset;
%       if size(XRD{1}.Fit_fs_i_yi,1)==1
%         y=(XRD{i}.Fit_fs_i_yi+Fbaseline(XRD{i}.baseline_fs,XRD{i}.Fit_fs_i_xi));
%       else
%         y=sum(XRD{i}.Fit_fs_yi+Fbaseline(XRD{i}.Fit_fs_baseline,XRD{i}.Fit_fs_i_xi));  
%       end      
%       plot(XRD{i}.Fit_fs_i_xi,(y-min(y))/max(y)+offset,'color',cmap(i,:),'linewidth',linw);
%     end
%     text((max(XRD{ID(i)}.theta_f)-min(XRD{ID(i)}.theta_f))*xloc+min(XRD{ID(i)}.theta_f),offset+yloc,XRD{ID(i)}.filename(1:end),'interpret','none','color',cmap(i,:),'fontsize',fs1)
%   end
%   
%   axis([min(XRD{1}.theta_f) max(XRD{1}.theta_f) -crop offset+1+crop])
%   xlabel(['2',char(952),char(176)],'Interpreter','tex','fontsize',fs2,'Fontname',fname);
%   ylabel(['Intensity'],'Interpreter','tex','fontsize',fs2,'Fontname',fname)
%   set(gca,'fontsize',fs0,'ytick',[],'xminortick','on','tickdir','out','linewidth',2,'Fontname',fname)
%   
%   if save_result==1
%     Image_Name=XRD{1}.filename(1:end);
%     expfig=['export_fig ' Image_Name ' -eps -pdf -png']
%     eval(expfig)
%   end
% end
% 



%%  Compare simulated data to experimental data

% Here all peaks have the same weight, but this can be adjusted
% later if desired
Weight=ones(length(ExpData),1);

Weight=ExpData(:,2)/max(ExpData(:,2));

Obj=zeros(length(XRD),1);
for i=1:length(XRD)
  for j=1:length(ExpData)
    MinDist=min(abs(XRD{i}.PeakLocInt(:,1)-ExpData(j,1)));
    Obj(i)=Obj(i)+MinDist*Weight(j);
  end
end

[~,IDobj]=sortrows(Obj);


if save_peak_info==1;
  fid=fopen(['LIST_' save_peak_file],'w');
  fprintf(fid,'FileName   | Objective Function');
  for i=1:length(IDobj)
    fprintf(fid,'\n%s %3.2f',XRD{IDobj(i)}.filename,Obj(IDobj(i)));
  end
  fid=fclose(fid);
end




if display_result==1
  close all; clc
  figure('name','XRD','Color','w', 'position',[-1566 204 1508 1490])
  cmap=colormap(lines(length(XRD)));
  
  linw=3;
  fs0=18;
  fs1=24;
  fs2=24;
  fname='Arial';
  
  doffset=.95;
  crop=0.025;
  xloc=.60;
  yloc=0.7*doffset;
  
  
  % Adjust order of line profiles using ID
  N=5
  ID=IDobj([1:N]);
  LabelID=[ID];
  
    
  for i=1:length(ExpData)
      plot([ExpData(i,1) ExpData(i,1)],[-100 100],'linewidth',ExpData(i,2)/ExpData(1,2)*10); hold on
  end
  
  for i=1:length(ID);
    hold on
    if i==1
      offset=0;
      plot(XRD{ID(i)}.theta_f,XRD{ID(i)}.intensity_fs/max(XRD{ID(i)}.intensity_fs)+offset,'color',cmap(i,:),'linewidth',linw);
    else
      offset=offset+doffset;
      plot(XRD{ID(i)}.theta_f,XRD{ID(i)}.intensity_fs/max(XRD{ID(i)}.intensity_fs)+offset,'color',cmap(i,:),'linewidth',linw);
    end
    text((max(XRD{ID(i)}.theta_f)-min(XRD{ID(i)}.theta_f))*xloc+min(XRD{ID(i)}.theta_f),offset+yloc,XRD{ID(i)}.filename(1:end),'interpret','none','color',cmap(i,:),'fontsize',fs1)
  end
  

  
  axis([min(XRD{1}.theta_f) max(XRD{1}.theta_f) -crop offset+1+crop])
  xlabel(['2',char(952),char(176)],'Interpreter','tex','fontsize',fs2,'Fontname',fname);
  ylabel(['Intensity'],'Interpreter','tex','fontsize',fs2,'Fontname',fname)
  set(gca,'fontsize',fs0,'ytick',[],'xminortick','on','tickdir','out','linewidth',2,'Fontname',fname)
  
  if save_result==1
    Image_Name=XRD{1}.filename(1:end);
    expfig=['export_fig ' Image_Name ' -eps -pdf -png']
    eval(expfig)
  end
end

