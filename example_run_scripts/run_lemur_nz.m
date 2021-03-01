clear mex;
clear all;

path = pwd;
cd ../matlab_setup
addpaths
cd(path)

load('../data/nz_uplift');
load('../data/nz_topo');
upm=upm2;
clear upm2;
upm2=cell(1,122);
c=1;
for i = 5:10:length(upm)
    upm2{c} = upm{i}/2+upm{i+1}/2;
    c=c+1;
end


    
%load('lemlat');
%load('ll');
Z1=z;
Z1(:,end+1) = Z1(:,end);

%Z1=Z1+rand(size(Z1))*1;% DEM has large flat areas?

%Z1(Z1<0)=-1000;
UL.k=zeros(size(Z1))+5e-5;%stream power k 

UL.Zi=Z1;%initial grid - set this to whatever you wish or load in a diferent grid 
[mm,nn] = size(Z1);%Size of grid
UL.dx = 900;%Grid distance (m)
UL.dy = 900;
UL.t=122e3;%total time (yr)
UL.dt=1e3;%time step (yr)
UL.display=0;%display topo during run?
UL.Udt = 1e3;%Time to switch between uplift maps, if using cell array (yr)
UL.flex = 30e3;%effective elastic thickness (m)
UL.kd = .01;%Diffusion coefficient (m2/yr)
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
UL.k_sed = UL.k;
UL.BC=-1;
UL.wdt=1000;
UL.firstcall=0;
UL.U=upm2;
%UL.U(150:end,:)=.0;
%UL.Zi=Z2;
UL.precip = 1;
UL.srho=2400; %Sediment density (kg/m3)
UL.deposit=true; %Use marine deposition?
UL.drawdt=100;% Plot topography every x timesteps
UL.undercapacity=ones(size(UL.Zi))+1;% Where to apply undercapacity model (1=undercapacity, 0=stream power).
UL.l=1; %Length scale in undercapacity model (m)
UL.ks =.001; %Sediment transport coefficient in undercapacity model;
UL.sinkfill = true; %Do not pair false with marine deposition - massive slowdown
UL.massconservativesinkfill = false; %mass conservative sinkfill? - this is not recommended to be paired with marine deposition
%A value of 1 means fill with sediment and a value of 2 means fill with
%runoff
UL.maxareasinkfill = 1000000;
UL.evaprate = 0;
[Z2,FD,SS] = lemur(UL); 
