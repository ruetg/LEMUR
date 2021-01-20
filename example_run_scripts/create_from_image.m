    clear mex;

    path = pwd;
      
    cd ../matlab_setup
    addpaths
    cd(path)
    img = flipud(double(rgb2gray(imread('../data/thumbnail_image.png'))));
    img(img<=200)=1;
    img(img>200)=0;
    img=imresize(img,290*2,884*2);
   % img=img(:,1:end-1);
    %load('lemlat');
    clear UL
    Z1=rand(size(img));
    
    UL.k=1e-6;%stream power k 
    UL.Zi=Z1;%initial grid - set this to whatever you wish or load in a diferent grid 
    [mm,nn] = size(Z1);%Size of grid

    UL.dx = 1500;%Grid distance (m)
    UL.dy = 1500;
    UL.t=50e6;%total time (yr)
    UL.dt=1e6;%time step (yr)
    UL.display=0;%display topo during run?
    UL.Udt = 1e3;%Time to switch between uplift maps, if using cell array (yr)
    UL.flex =0e3;%effective elastic thickness (m)
    UL.kd = .0001;%Diffusion coefficient (m2/yr)
    UL.m=.5;% Drainage area exponent in stream power 
    UL.n=1;%stream power n

    % 
    BC=zeros(size(Z1)); %Boundarycondition 
    BC(:,end)=1;
    BC(end,:)=1;
    BC(:,1)=1;
    BC(1,:)=1;
    UL.BC=find(BC==1);%BC must be a list of outlet linear indices
    %Alternatively, use BC=-1 to make all nodes below sea level (transiently)
    %outlets
   % UL.BC=0;
    UL.wdt=UL.dt;%Record data every x years
    UL.firstcall=0;
    UL.U=zeros(mm,nn)+.00000;
    %UL.U(150:end,:)=.0;
    %UL.Zi=Z2;
    UL.srho=2400; %Sediment density (kg/m3)
    UL.deposit=false; %Use marine deposition?
    UL.drawdt=10;% Plot topography every x timesteps
    UL.undercapacity=zeros(size(UL.Zi))+0;% Where to apply undercapacity model (1=undercapacity, 0=stream power).
    UL.l=1; %Length scale in undercapacity model (m)
    UL.ks =1; %Sediment transport coefficient in undercapacity model;
    UL.sinkfill = true; %Do not pair false with marine deposition - massive slowdown
    UL.massconservativesinkfill = 0; %mass conservative sinkfill? - this is not recommended to be paired with marine deposition
    UL.maxareasinkfill = 1000000;
    UL.precip=1;
    UL.evaprate=4;
    UL.k_sed=UL.k+.3e-8;
    UL.U = ones(size(Z1));
    UL.Zi=UL.Zi+700;
   
    UL.k=zeros(mm,nn)+1e-6;
    UL.k(img==1)=.3e-6;
    UL.U = UL.U*3e-5;
    UL.Zi(UL.BC)=0;
    UL.k_sed=UL.k;
    disp('here')
    [Z2,FD,SS,f] = lemur(UL); 
