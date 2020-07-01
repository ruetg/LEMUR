

function [Z,FD,OUT] = lemur_wrapper(UL)

clear mex;
%Relatively constant constants
wrho = 1029;%water density
srho =UL.srho;%
sinkfilling=UL.sinkfill;
maxareasinkfill=UL.maxareasinkfill;
landsed=UL.massconservativesinkfill;
%Gather parameters - see input file for descriptions
[m,n] = size(UL.Zi);
Z=UL.Zi;
U=UL.U;
BC=UL.BC;
k=UL.k;
tt = UL.t;
dt = UL.dt;
tisoero=zeros(size(Z));
dx = UL.dx;
dy = UL.dy;
fluvialm = UL.m;
kd = UL.kd;
wdt = UL.wdt;
flex = UL.flex;
ks = UL.ks;
nval = UL.n;
drawdt=UL.drawdt;
L=UL.undercapacity;
deposit=UL.deposit;
l1=UL.l;
tsinkfill = zeros(m,n);
FD.ix = ones(1,m*n);
FD.ixc = ones(1,m*n);
display=UL.display;
firstcall=1;%First run through the time loop is the often the slowest, establish memory needed

if firstcall
    vsave=0;%Saved radius indices used in deposition function
end

tic;
%If erodibility is not a map, make it a map.
if numel(k)==1
    k = ones(m,n)*k;
    
end

%Is uplift a set of maps or just one map?
dynamicU = 0;

if isa(U,'cell')
    dynamicU = 1;
    iU = U;
    U = iU{1};
    Udt = params(10);
else
    Udt = 3;
    if length(U) ==1
        
        U = ones(m,n).*U;
    end
end

tero = zeros(m,n);%Keep track of total advective erosion
tdepo = zeros(m,n);%Keep track of total deposition
tdiffu = zeros(m,n);%Keep track of total diffusive erosion
depo=zeros(m,n);

%Is boundary condition dynamic <0?
dBC = 0;
if BC == -1
    dBC = 1;
end

ero = zeros(m,n);

%Use same format for flow direction as topotoolbox, in case topotoolbox is used to analyze
%results.  This has not been tested to see if it actually works with
%topotoolbox though
FD.I=zeros(m*n,1);
FD.R=zeros(m*n,1);
size(Z)

lemur_mex(m,n);

%The ulem mex functions must be fed all necessary parameters
lemur_mex('set','dx',dx,'dy',dy,'dt',dt);
lemur_mex('set','l',l1,'ks',ks,'kd',kd);
lemur_mex('set','m',fluvialm,'n',nval,'landsed',double(landsed),'maxareasinkfill',maxareasinkfill);

% main time loop
for t = 0:dt:tt-dt
    lemur_mex('set','u',U(:));
    
    lemur_mex('set','z',Z(:));
    
    lemur_mex('set','k',k(:));
    lemur_mex('set','undercapacity',L(:))
    Zi = Z;
    %If uplift is a set of maps, use the a new map every <Udt> years
    if and(dynamicU==1,(mod(t,Udt) == 0)) && t~=tt
        if t/Udt<length(iU)
            
            U = iU{floor(t/Udt)+1};
            if length(U) ==1
                
                U = ones(m,n).*U;
            end
            disp('using U#')
            disp(num2str(floor(t/Udt)+1));
        else
            U = zeros(m,n);
        end
    end
    %If boundary condition is dynamic, set it to all points below zero
    %elevation
    if dBC
        BC = find(Z<=0);
        % BC = [BC];
    else
        U(BC) = 0;
        Z(BC) = 0;
    end
    lemur_mex('set','bc',BC(:));
    disp(strcat(num2str(t),' years gone by'));
    
    
    sinkfill=Z;
    %Sink filling algorithm - if use "landsed" option for mass conservative
    %sinkfill.
    if sinkfilling||landsed
        lemur_mex('run','lakefill')
        Z=lemur_mex('get','z');
        Z=reshape(Z,m,n);
    end
    %keep track of sinkfill
    sinkfill = Z-sinkfill;
    sinkfill(BC)=0;
    
    %Fluvial erosion
    %%
    DEM=Z;
    m
    tic;lemur_mex('run','erode_fluvial');toc
    
    Z=lemur_mex('get','z');
    Z=reshape(Z,m,n);
    BCi = zeros(m,n);
    BCi(BC) = 1;
    Z(1,~BCi(1,:)) =  Z(2,~BCi(1,:));
    Z(end,~BCi(end,:)) =  Z(end-1,~BCi(end,:));
    Z(~BCi(:,end),end) =  Z(~BCi(:,end),end-1);
    Z(~BCi(:,1),1) =  Z(~BCi(:,1),2);
    
    ero = DEM-Z;
    
    ero(:,1)=0;
    ero(:,end)=0;
    ero(end,:)=0;
    ero(1,:)=0;
    ero(Z<0)=0;
    %%
    lemur_mex('set','ero',ero(:));
    lemur_mex('set','z',Z(:));
    tero = ero + tero;
    
    %Diffusion presently only supports square cells
    if  kd~=0
        diffu =Z;
        lemur_mex('set','seadiffusion',0);
        lemur_mex('set','kd',kd);
        
        lemur_mex('run','erode_diffusion');
        Z=lemur_mex('get','z');
        
        Z=reshape(Z,m,n);
        
        diffu = diffu-Z;
        tdiffu = tdiffu+diffu;
        
    else
        U(BC) = 0;
        diffu = 0;
    end
    
    %Deposition presently only supports square cells
    if deposit==true
        depo=Z;
        
        lemur_mex('run','deposition');
        Z=lemur_mex('get','z');
        Z=reshape(Z,m,n);
        depo=Z-depo;
        
        
        
    else
    end
    
    
    if 1
        Zc = Z;
        Zc(Zc>0)=0;
        Zi(Zi>0)=0;
        %Total loading due to sediment
        alls=-depo+ero+diffu-sinkfill;
        
        Des=zeros(size(Z))+srho;
        Dew=zeros(size(Z))+wrho;
        
        %Presently flexure due to water loading is no longer supported but will
        %be added back soon
        if flex>0
            isoall = IsoFlex(alls,flex,dx,dy,srho,100);
        else
            isoall=0;
        end
        
        
        %Calculate total deformation due to sediment
        tisoero=tisoero+isoall;
        isoall(:,1)=0;
        isoall(:,end)=0;
        isoall(end,:)=0;
        isoall(1,:)=0;
        Z = Z+isoall;
        
    end
    if sinkfilling
        sinkareas=lemur_mex('get','sinkareas');
    end
    Z = Z+U*dt;
    if sinkfilling
        depo(sinkareas>maxareasinkfill)=sinkfill(sinkareas>maxareasinkfill);
        depo(Z<0)=0;
    end
    %Display topography every drawdt timesteps
    if display && mod(t/dt,drawdt)==0
        figure(2);
        h=imagesc(Z);shading interp;demcmap(Z,2500);
        %[LON,LAT]=meshgrid(1:n,1:m);hold off;
        %dem((1:n).*.1,(1:m).*.075,Z,'latlon','zlim',[-3e4,3e4]);
        drawnow;
        hold on;
        
        %%
        
        % figure(2);drawnow;
        %               frame = getframe(gcf);
        %       im = frame2im(frame);
        %       [imind,cm] = rgb2ind(im,256);
        %       % Write to the GIF File
        %       if t == 0
        %           imwrite(imind,cm,'AA','gif', 'Loopcount',inf,'delaytime',.1);
        %       else
        %           imwrite(imind,cm,'AA','gif','WriteMode','append','delaytime',.1);
        %       end
        
    end
    
    
    nd=0;
    tsinkfill = tsinkfill+sinkfill+nd;
            tdepo = tdepo+depo;

    %Output
    %Write output every wdt years
    FD.I=lemur_mex('get','stack');
    FD.R=lemur_mex('get','rec');
    if mod(t,wdt)==0
        
        wt = floor(t/(wdt+1e-9))+1;
        OUT(wt).topo = Z;
        OUT(wt).ero = tero;
        OUT(wt).FD = FD;
        OUT(wt).depo = tdepo;
        OUT(wt).isostasy = tisoero;
        OUT(wt).diffusion = tdiffu;
        OUT(wt).sinkfill = tsinkfill;
    end
    
    firstcall=0;
    
    [~,rtrt(t/dt+1)]=max(Z(:,250));
end
%save('rtrt','rtrt');

