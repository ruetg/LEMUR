clear mex;
clear all;
Z=zeros(500,500);
Z=Z+rand(size(Z));
%esh(Z)
Z1=Z;

path = pwd;
cd ../matlab_setup


addpaths
cd(path)

UL.k=zeros(size(Z1))+2e-6;%stream power k 
UL.Zi=Z+100;%initial grid - set this to whatever you wish or load in a diferent grid 
[mm,nn] = size(Z1);%Size of grid
UL.dx = 1000;%Grid distance (m)
UL.dy = 1000;
UL.t=30e6;%total time (yr)
UL.dt=1e5;%time step (yr)
UL.display=1;%display topo during run?
UL.Udt = 1e3;%Time to switch between uplift maps, if using cell array (yr)
UL.flex =20e3;%effective elastic thickness (m)
UL.kd = .00;%Diffusion coefficient (m2/yr)
UL.m=.5;% Drainage area exponent in stream power 
UL.n=1;%stream power n

BC=zeros(size(Z1)); %Boundarycondition 
BC(:,end)=0;
BC(end,:)=1;
BC(:,1)=0;
BC(1,:)=1;
UL.BC=find(BC==1);%BC must be a list of outlet linear indices
% 
% %Alternatively, use BC=-1 to make all nodes below sea level (transiently)
%outlets
%UL.BC=-1;
UL.evaprate = 3; %(m/yr)
UL.wdt=UL.dt;%Record data every x years
UL.firstcall=0;
UL.U=zeros(mm,nn);
UL.U(150:end-150,:)=.00025;
%UL.Zi=Z2;
UL.srho=2400; %Sediment density (kg/m3)
UL.deposit=0; %Use marine deposition?
UL.drawdt=1;% Plot topography every x timesteps
UL.undercapacity=ones(size(UL.Zi));% Where to apply undercapacity model (1=undercapacity, 0=stream power).
UL.l=1; %Length scale in undercapacity model (m)
UL.ks =1; %Sediment transport coefficient in undercapacity model;
UL.sinkfill = true; %Do not pair false with marine deposition - massive slowdown
%A value of 2 indicates precip-based filling, a value of 1 indicates sediment-based
%filling
UL.precip=1;
UL.maxareasinkfill = 0;
%%
load('../data/Z1_andes.mat')
Z2=Z1;
%imagesc(Z2);shading interp;demcmap(Z2)
%imagesc(SS(1).ero);colormap(custmap(2,1));
[caxis([-2 2])];
c = colorbar;title(c,'erosion net (m')
UL.Zi = Z2*.75+100;

%%
clear mex;
close all
UL.drawdt=1;% Plot topography every x timesteps
k = 5e-6;


UL.l = 1;
UL.precip = .2;
UL.t = 1e6;
UL.dt=1e4;
c=1;
UL.evaprate = 2.5; %(m/yr)
UL.kd=0;
UL.massconservativesinkfill = 3; %mass conservative sinkfill? - this is not recommended to be paired with marine deposition
maxsed=[];
disp('here')
U=UL.U;
UL.U(100:200,:)=8*U(250,250);
UL.k=k*UL.precip^.5;
UL.k_sed=1.0001*UL.k;
UL.ks=1;
[Z2,FD,SS,f] = lemur(UL);
c=c+1;
sed=getacc(SS(end).ero,FD);
