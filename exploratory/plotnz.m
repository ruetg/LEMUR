%load ../data/
hs = tight_subplot(1,3,.05,.05,.05);
axes(hs(1));
t1 = 20;
t2 = t1-1;
er1 = SS(t1).ero(:,1:end-1)*0.1-SS(t1).depo(:,1:end-1)*0.1;
er2 = SS(t2).ero(:,1:end-1)*0.1-SS(t2).depo(:,1:end-1)*0.1;
er = er1-er2;
pcolor(LON,LAT,er);shading interp;
colormap(gca,bluewhitered(50))
hold on
plot(coastlon,coastlat,'k','linewidth',2)
set(gca,'fontsize',14)
ylabel('Lat \circ')
xlabel('Lon \circ')
title('100 Ka');


axes(hs(2));
t1 = 90;
t2 = t1-1;
er1 = SS(t1).ero(:,1:end-1)*0.1-SS(t1).depo(:,1:end-1)*0.1;
er2 = SS(t2).ero(:,1:end-1)*0.1-SS(t2).depo(:,1:end-1)*0.1;
er = er1-er2;
pcolor(LON,LAT,er/1000);shading interp;
colormap(gca,bluewhitered(50))
hold on
plot(coastlon,coastlat,'k','linewidth',2)
set(gca,'fontsize',14)
ylabel('Lat \circ')
xlabel('Lon \circ')
title('120 Ka');
title('30 Ka');


axes(hs(3));
t1 = 121;
t2 = t1-1;
er1 = SS(t1).ero(:,1:end-1)*0.1-SS(t1).depo(:,1:end-1)*0.1;
er2 = SS(t2).ero(:,1:end-1)*0.1-SS(t2).depo(:,1:end-1)*0.1;
er = er1-er2;
pcolor(LON,LAT,er);shading interp;
colormap(gca,bluewhitered(50))
hold on
plot(coastlon,coastlat,'k','linewidth',2)
set(gca,'fontsize',14)
ylabel('Lat \circ')
xlabel('Lon \circ')
title('1 Ka');

