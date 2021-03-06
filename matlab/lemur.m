

function [Z,FD,OUT,f] = lemur(UL)

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
k_sed=UL.k_sed;
BC=UL.BC;
precip=UL.precip;
evaprate =UL.evaprate;
tt = UL.t;
dt = UL.dt;
tisoero=zeros(size(Z));
dx = UL.dx;
dy = UL.dy;
fluvialm = UL.m;
k=UL.k;
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
    ki=k(1);

%Is uplift a set of maps or just one map?
dynamicU = 0;

if isa(U,'cell')
    dynamicU = 1;
    iU = U;
    U = iU{1};
    Udt = UL.Udt;
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
[X,Y] = meshgrid(1:n,1:m);

%The ulem mex functions must be fed all necessary parameters
lemur_mex('set','dx',dx,'dy',dy,'dt',dt);
lemur_mex('set','l',l1,'ks',ks,'kd',kd);
lemur_mex('set','m',fluvialm,'n',nval,'precip',precip,'maxareasinkfill',maxareasinkfill);
lemur_mex('set','uselandsed',landsed);
lemur_mex('set','evaprate',evaprate);

% main time loop
  lemur_mex('set','firstcall',1);
  if display
   hs=tight_subplot(1,2,.05,.05,.05);
  end
km=ki;
ksedm = k_sed;
tsed=0;
for t = 0:dt:tt-dt
    tic()
    if t>tt/2
        lemur_mex('set','precip',precip*2);
                lemur_mex('set','evaprate',evaprate/2);

         km=ki*2^.5;
         k_sed = 2^.5*ksedm;

    end
%    lemur_mex('set','u',U(:));
    
    lemur_mex('set','z',Z(:));
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
        BC = find(Z < 0);
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
    tsinkfill = tsinkfill+sinkfill;
wt=1;
    %Fluvial erosion
    %%
    DEM=Z;
    %k(:)=km;
    tsed = tsed-ero+sinkfill;
    %tsed = -ero+sinkfill;
    tsed = tsed/2 +abs(tsed)/2;
    %k(tsed>0)=k_sed;
    lemur_mex('set','k',k(:));
    lemur_mex('run','erode_fluvial');
    
    disp('herex')
    
    Z=lemur_mex('get','z');
    Z=reshape(Z,m,n);
    Z=Z+sinkfill;

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
    toc()
    %Diffusion presently only supports square cells
    if numel(U)>1
        U = imresize(U,size(ero));
    end

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
        size(U)
        if numel(U) > 1
            U(BC) = 0;
            diffu = 0;
        end
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
    if sinkfilling&&landsed
        sinkareas=lemur_mex('get','sinkareas');
    end
     %modified for nz data, must interpolate because grid is too large
     size(U)
    Z = Z+U*dt;

    if sinkfilling&&landsed
        depo(sinkareas>maxareasinkfill)=sinkfill(sinkareas>maxareasinkfill);
        depo(Z<0)=0;
    end
    water=lemur_mex('get','watertot');
    %Display topography every drawdt timesteps
    FD.I=lemur_mex('get','stack');
    FD.R=lemur_mex('get','rec');
    sed=getacc(ero,FD);
    smax(t/dt+1)=max(max(sed(95:105,75:135)));
    if 1% temporary
        imagesc(flipud(Z));shading interp;axis off;
        demcmap([-1000 1000]);
        drawnow;
    end
    if display ==1 && mod(t/dt,drawdt)==0&&0
        axes(hs(2));
        hold off;

        pcolor(Z);shading interp;demcmap([-5000 4000]);
        hold on;

        x=X(k==k_sed);
        y=Y(k==k_sed);
        
        if length(x)>0
            rands = randi(length(x),[1,length(x)]);
            scatter(x(rands),y(rands),50,[1,1,0],'filled','s','markerfacealpha',.5,'markeredgealpha',.5)
%            p.Color(4)=.4;
            
        end
        
        x=X(water>0);
        y=Y(water>0);
        if length(x)>0
            rands = randi(length(x),[1,length(x)]);
            if 0
                scatter(x(rands),y(rands),50,[0,0,1],'filled','s','markerfacealpha',.5,'markeredgealpha',.5)
%                        p.Color(4)=.4;
            end

        end
        plotstrm(FD,Z,X,Y);
        set(gca,'fontsize',14);
        c=colorbar;
        title(c,'elevation (m)');

                    axes(hs(1));
        hold off
        plot(0:dt/1000:t/1000,smax*dx*dy/dt,'linewidth',2)
        xlabel('time (ka)')
        set(gca,'fontsize',14);
        ylabel('Sed. flux (m^3/yr)')
        er=ero(100:350,:);
        mean(mean(er))
        set(gcf,'units','normal','position',[0 0 .9 .9]);

        title([num2str(t),' yrs gone by']);
        drawnow;
        f(t/(drawdt*dt)+1)=getframe(gcf);
        
    end
    
    set(gcf,'units','normal','position',[0 0 .9 .9]);
    hold off;
    tdepo = tdepo+depo;
    
    %Output
    %Write output every wdt years
mean(mean(ero(Z>0)))
    if (length(wdt)>1&&ismember(t,wdt)) || ...
            (length(wdt) == 1 && mod(t,wdt)==0)
        if length(wdt) == 1
            wt = floor(t/(wdt+1e-9))+1;
        else
            wt = wt+1;
        end

        OUT(wt).topo = Z;
        OUT(wt).ero = tero;
        OUT(wt).FD = FD;
        OUT(wt).depo = tdepo;
        OUT(wt).isostasy = tisoero;
        OUT(wt).diffusion = tdiffu;
        OUT(wt).sinkfill = tsinkfill;
    end
    firstcall=0;
    set(gcf,'Units','normal')
    f(t/dt+1)=getframe(gcf);
    [~,rtrt(t/dt+1)]=max(Z(:,250));
    lemur_mex('set','firstcall',0);
end
% v=VideoWriter('~/climate_w_ero.mp4','mpeg-4');
% v.FrameRate=5;
% v.open()
% 
% v.writeVideo(f);
% v.close();
%save('rtrt','rtrt');
end
function plotstrm(FD,Z1,LON,LAT)
thres=500;
[mm,nn]=size(Z1);
acc = zeros(mm,nn)+1;
    FD.ix=(FD.I);
    FD.ixc = FD.R(FD.I);
    for i = length(FD.I):-1:1
        acc(FD.ixc(i)) = acc(FD.ixc(i)) + acc(FD.ix(i));
    end

plot(LON(acc>thres),LAT(acc>thres),'.k','markersize',2);
end

