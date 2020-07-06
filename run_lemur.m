clear mex;
clear all;
load('lemlat');
%load('ll');
Z1=Z1;%double(calbuco_dem(1:end,1:end));
 Z1=zeros(500,502);
 Z1(1:250,:)=-5000;
 Z1(251:end,:)=500;
Z1=Z1+rand(size(Z1))*.01;% DEM has large flat areas?
%Z1(Z1<0)=-1000;
UL.k=zeros(size(Z1))+5e-5;%stream power k 
UL.Zi=Z1;%initial grid - set this to whatever you wish or load in a diferent grid 
[mm,nn] = size(Z1);%Size of grid
UL.dx = 3000;%Grid distance (m)
UL.dy = 3000;
UL.t=20e4;%total time (yr)
UL.dt=1e4;%time step (yr)
UL.display=1;%display topo during run?
UL.Udt = 1e3;%Time to switch between uplift maps, if using cell array (yr)
UL.flex =.0000001e3;%effective elastic thickness (m)
UL.kd = .001;%Diffusion coefficient (m2/yr)
UL.m=.5;% Drainage area exponent in stream power 
UL.n=1;%stream power n

% 
% BC=zeros(size(Z1)); %Boundarycondition 
% BC(:,end)=1;
% BC(end,:)=1;
% BC(:,1)=1;
% BC(1,:)=1;
% UL.BC=find(BC==1);%BC must be a list of outlet linear indices

%Alternatively, use BC=-1 to make all nodes below sea level (transiently)
%outlets
UL.BC=-1;
UL.wdt=UL.dt;%Record data every x years
UL.firstcall=0;
UL.U=zeros(mm,nn)+.00000;
%UL.U(150:end,:)=.0;
%UL.Zi=Z2;
UL.srho=2400; %Sediment density (kg/m3)
UL.deposit=true; %Use marine deposition?
UL.drawdt=10;% Plot topography every x timesteps
UL.undercapacity=zeros(size(UL.Zi))+0;% Where to apply undercapacity model (1=undercapacity, 0=stream power).
UL.l=0; %Length scale in undercapacity model (m)
UL.ks =0; %Sediment transport coefficient in undercapacity model;
UL.sinkfill = true; %Do not pair false with marine deposition - massive slowdown
UL.massconservativesinkfill = false; %mass conservative sinkfill? - this is not recommended to be paired with marine deposition
UL.maxareasinkfill = 1000000;
[Z2,FD,SS] = lemur_wrapper(UL); 
