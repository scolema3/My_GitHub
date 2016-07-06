% script to plot data from LAMMPS minimization routne that performs a
% hydrostatic loading and unloading on the strucutre in attempt to find a
% global energy minima 

close all
clc

prefix='B11Cp_CBC_tri_4x4x4_';
natoms=14580;

file=[prefix 'Min1.txt'];
D1=dlmread(file,' ',1,0);

file=[prefix 'Load1.txt'];
D2=dlmread(file,' ',1,0);

file=[prefix 'Min2.txt'];
D3=dlmread(file,' ',1,0);

file=[prefix 'Load2.txt'];
D4=dlmread(file,' ',1,0);

file=[prefix 'Min3.txt'];
D5=dlmread(file,' ',1,0);

file=[prefix 'Load3.txt'];
D6=dlmread(file,' ',1,0);

% file=[prefix 'Min4.txt'];
D7=dlmread(file,' ',1,0);



           %1      2      3    4    5                       6       7         8     9     10    11    12    13    14   15   16   17       18     19    20    21      22     23
axislabel={'step' 'temp' 'ke' 'pe' 'Cohessive Energy (eV)' 'press' 'P (GPa)' 'pxx' 'pyy' 'pzz' 'pxy' 'pxz' 'pyz' 'lx' 'ly' 'lz' 'Volume'  'a'   'b'   'c'   'alpha' 'beta' 'gamma'};
xid=17;        % index for x-data
yid=5;         % index for y-data
NormX=natoms;
NormY=1;

figure('Position',[1000 569 1010 913],'Color','w')

cmap=colormap(lines(7));
LnW=2;
plot(D1(1,xid)/NormX,D1(1,yid)/NormY,'go', 'HandleVisibility','off','LineWidth',LnW,'color',cmap(1,:)); hold on
plot(D1(:,xid)/NormX,D1(:,yid)/NormY,'-.', 'DisplayName','Minimization 1','LineWidth',LnW,'color',cmap(1,:)); hold on
plot(D2(:,xid)/NormX,D2(:,yid)/NormY,'--', 'DisplayName','Load Cycle 1','LineWidth',LnW,'color',cmap(2,:)); hold on
plot(D3(:,xid)/NormX,D3(:,yid)/NormY,'--', 'DisplayName','Minimization 2','LineWidth',LnW,'color',cmap(3,:)); hold on
plot(D4(:,xid)/NormX,D4(:,yid)/NormY,'--', 'DisplayName','Load Cycle 2','LineWidth',LnW,'color',cmap(4,:)); hold on
plot(D5(:,xid)/NormX,D5(:,yid)/NormY,'--', 'DisplayName','Minimization 3','LineWidth',LnW,'color',cmap(5,:)); hold on
plot(D6(:,xid)/NormX,D6(:,yid)/NormY,'--', 'DisplayName','Load Cycle 3','LineWidth',LnW,'color',cmap(6,:)); hold on
plot(D7(:,xid)/NormX,D7(:,yid)/NormY,'--', 'DisplayName','Minimization 4','LineWidth',LnW,'color',cmap(7,:)); hold on
plot(D7(end,xid)/NormX,D7(end,yid)/NormY,'ro', 'HandleVisibility','off','LineWidth',LnW); hold on
h=legend('show');
set(h,'box','off','location','northwest','fontsize',14);
xlabel([axislabel{xid} ' per Atom' ],'fontsize',14)
ylabel([axislabel{yid} ],'fontsize',14) 

expfig=['export_fig ' prefix 'alldata' ' -eps -pdf -png'];
eval(expfig)


clf
plot(D2(:,xid)/NormX,D2(:,yid)/NormY,'--', 'DisplayName','Load Cycle 1','LineWidth',LnW,'color',cmap(2,:)); hold on
plot(D4(:,xid)/NormX,D4(:,yid)/NormY,'--', 'DisplayName','Load Cycle 2','LineWidth',LnW,'color',cmap(4,:)); hold on
plot(D6(:,xid)/NormX,D6(:,yid)/NormY,'--', 'DisplayName','Load Cycle 3','LineWidth',LnW,'color',cmap(6,:)); hold on
plot(D7(end,xid)/NormX,D7(end,yid)/NormY,'ro', 'HandleVisibility','off','LineWidth',LnW); hold on
h=legend('show');
set(h,'box','off','location','northwest','fontsize',14);
xlabel([axislabel{xid} ' per Atom' ],'fontsize',14)
ylabel([axislabel{yid} ],'fontsize',14) 


x=[D1(:,xid); D2(:,xid) ;D3(:,xid) ;D4(:,xid) ;D5(:,xid); D6(:,xid); D7(:,xid)]/NormX;
y=[D2(:,yid) ;D3(:,yid) ;D4(:,yid) ;D5(:,yid); D6(:,yid); D7(:,yid)]/NormY;
axis([min(x),max(x),min(y),max(y)])


arrow([D6(1:end-1,xid)/NormX,D6(1:end-1,yid)/NormY],[D6(2:end,xid)/NormX,D6(2:end,yid)/NormY],'Length',15,'Width',0,'tipangle',25,'baseangle',30); hold on

plot(D2(:,xid)/NormX,D2(:,yid)/NormY,'--', 'DisplayName','Load Cycle 1','LineWidth',LnW,'color',cmap(2,:)); hold on
plot(D4(:,xid)/NormX,D4(:,yid)/NormY,'--', 'DisplayName','Load Cycle 2','LineWidth',LnW,'color',cmap(4,:)); hold on
plot(D6(:,xid)/NormX,D6(:,yid)/NormY,'--', 'DisplayName','Load Cycle 3','LineWidth',LnW,'color',cmap(6,:)); hold on
plot(D7(end,xid)/NormX,D7(end,yid)/NormY,'ro', 'HandleVisibility','off','LineWidth',LnW); hold on
expfig=['export_fig ' prefix 'load_data' ' -eps -pdf -png'];
eval(expfig)

