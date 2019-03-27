%% Plot Survey in map view and depth
f1 = figure(1); clf;
set(gcf,'position',[53 292 1100 513]);
clr = lines(4);

% limits etc. 
latlim = [min(lats_ship) max(lats_ship)]+[-1 1]*.2/111; % 200m outside ship circle
lonlim = [min(lons_ship) max(lons_ship)]+[-1 1]*.2/111; % 200m outside ship circle
deplim = [mean(z_sta) 0]+[-100 20];
londist = vdist(mean(lats_ship),lonlim(1),mean(lats_ship),lonlim(2));

% make axes
ax1 = axes('pos',[0.10 0.1100 0.38 0.8150]); hold(ax1,'on')
ax2 = axes('pos',[0.57 0.1100 0.38 0.8150]); hold(ax2,'on')


plot(ax1,lons_ship,lats_ship,'ok','markerfacecolor',clr(2,:),'markersize',13,'linewidth',1);
plot(ax1,data_bad.lons,data_bad.lats,'xk','markersize',13,'linewidth',2);
plot(ax1,lon_drop,lat_drop,'sk','markerfacecolor',[0.5 0.5 0.5],'markersize',15,'linewidth',1);
plot(ax1,mean(lon_sta),mean(lat_sta),'pk','markerfacecolor',[1 1 0],'markersize',25,'linewidth',1)
grid(ax1,'on'); 
set(ax1,'fontsize',16,'linewidth',2,'box','on',...
    'xlim',lonlim,'ylim',latlim);
title(ax1,['\textbf{Drift: ',num2str(mean(drift),'%.1f'),' m ~~~  Azi: ',num2str(mean(azi),'%.1f'),'$^{\circ}$}'],'fontsize',18,'interpreter','latex');
xlabel(ax1,'Longitude','fontsize',18,'interpreter','latex');
ylabel(ax1,'Latitude','fontsize',18,'interpreter','latex');
axis(ax1,'equal');

plot(ax2,[lons_ship, ones(size(lons_ship))*mean(lon_sta)]',[z_ship, ones(size(z_ship))*mean(z_sta)]','-k','linewidth',0.5,'color',[0.5 0.5 0.5]);
h1(1) = plot(ax2,lons_ship,z_ship,'ok','markerfacecolor',clr(2,:),'markersize',13,'linewidth',1);
plot(ax2,data_bad.lons,zeros(size(data_bad.lats)),'xk','markersize',13,'linewidth',2);
h1(2) = plot(ax2,lon_drop,z_drop,'sk','markerfacecolor',[0.5 0.5 0.5],'markersize',15,'linewidth',1);
h1(3) = plot(ax2,mean(lon_sta),mean(z_sta),'pk','markerfacecolor',[1 1 0],'markersize',25,'linewidth',1);
set(ax2,'fontsize',16,'linewidth',2,'box','on',...
    'xlim',lonlim,'ylim',deplim);
ht = title(ax2,['\textbf{Z$_{sta} -$ Z$_{drop}$ = ',num2str(mean(dz_sta),'%.1f'),' m}'],'fontsize',18,'interpreter','latex'); 
set(ht,'pos',get(ht,'pos')+[0 50 0]);
xlabel(ax2,'Longitude','fontsize',18,'interpreter','latex');
ylabel(ax2,'Depth','fontsize',18,'interpreter','latex');
legend(h1,{'Ship','Drop Pt.','Solution'},'location','southeast','fontsize',13);

% correct vertical exaggeration
kmpdeg = londist./diff(lonlim);
vexag_actual = vertexag( ax2 )*kmpdeg*1000;
vexag_want = (diff(deplim)./londist)/1000;
if vexag_want./vexag_actual < 1
posvec = [axpos(ax2,1),...
          axpos(ax2,2) + axpos(ax2,4)*0.5*(1- vexag_want./vexag_actual),...
          axpos(ax2,3),...
          axpos(ax2,4).*(vexag_want./vexag_actual)];
else
posvec = [axpos(ax2,1),...
          axpos(ax2,2),...
          axpos(ax2,3)./(vexag_want./vexag_actual),...
          axpos(ax2,4)];
end
set(ax2,'pos',posvec);
% vexag_actual_new = vertexag( ax2 )*kmpdeg*1000;

