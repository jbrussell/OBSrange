%% plot the histograms of model parameters
% histogram bins
Nbins = 15;

% set up figure
f100 = figure(100); clf;
set(gcf,'position',[3 1 1208 697],'color','white');

% set up axes
ax1 = axes('pos',[0.08 0.58 0.25 0.35]);
ax2 = axes('pos',[0.38 0.58 0.25 0.35]);
ax3 = axes('pos',[0.68 0.58 0.25 0.35]);
% ax4 = axes('pos',[0.08 0.12 0.25 0.35]);
ax5 = axes('pos',[0.38 0.12 0.25 0.35]);
ax6 = axes('pos',[0.68 0.12 0.25 0.35]);

% make the histograms
[ Ncount_lat, cent_lat ] = plot_hist(ax1,lat_sta,Nbins);
title(ax1,'\textbf{Latitude}','fontsize',18,'interpreter','latex');

[ Ncount_lon, cent_lon ] = plot_hist(ax2,lon_sta,Nbins);
title(ax2,'\textbf{Longitude}','fontsize',18,'interpreter','latex');

[ Ncount_z, cent_z ] = plot_hist(ax3,z_sta,Nbins);
title(ax3,'\textbf{Depth (m)}','fontsize',18,'interpreter','latex');

% [ Ncount_TAT, cent_TAT ] = plot_hist(ax4,TAT*1000,Nbins);
% title(ax4,'\textbf{TAT (ms)}','fontsize',18,'interpreter','latex');

[ Ncount_vw, cent_vw ] = plot_hist(ax5,V_w,Nbins);
title(ax5,'\textbf{Water Velocity (m/s)}','fontsize',18,'interpreter','latex');

[ Ncount_drift, cent_drift ] = plot_hist(ax6,drift,Nbins);
title(ax6,'\textbf{Drift (m)}','fontsize',18,'interpreter','latex');

%% ticks for the histograms - in m about central value
[cent_x,cent_y] = lonlat2xy_nomap(mean(cent_lon),mean(cent_lat),cent_lon,cent_lat);
% ticks every two m
tick_x = [2*floor(min(cent_x)/2):2:2*ceil(max(cent_x)/2)];
tick_y = [2*floor(min(cent_y)/2):2:2*ceil(max(cent_y)/2)];
[~,tick_lat] = xy2lonlat_nomap(mean(cent_lon),mean(cent_lat),0,tick_y);
[tick_lon,~] = xy2lonlat_nomap(mean(cent_lon),mean(cent_lat),tick_x,0);
[~,latlim] = xy2lonlat_nomap(mean(cent_lon),mean(cent_lat),0,1.2*[min(cent_y),max(cent_y)]);
[lonlim,~] = xy2lonlat_nomap(mean(cent_lon),mean(cent_lat),1.2*[min(cent_x),max(cent_x)],0);
set(ax1,'xtick',tick_lat,'xticklabel',tick_y,'xlim',latlim)
set(ax2,'xtick',tick_lon,'xticklabel',tick_x,'xlim',lonlim)
xlabel(ax1,sprintf('m from %.5f$^{\\circ}$',mean(cent_lat)),'interpreter','latex');
xlabel(ax2,sprintf('m from %.5f$^{\\circ}$',mean(cent_lon)),'interpreter','latex');
set(ax3,'xlim',mean(cent_z) + 1.3*([min(cent_z),max(cent_z)]-mean(cent_z)));
% set(ax4,'xlim',mean(cent_TAT) + 1.3*([min(cent_TAT),max(cent_TAT)]-mean(cent_TAT)));


%% drift directions
ax7 = polaraxes('pos',[0.621 0.35 0.16 0.16]); hold on
driftaz = atan2d(x_sta,y_sta);
plot_drift_polar( ax7,driftaz );

%% prettify
% remove ax6 tick underneath the polar plot
ytl=get(ax6,'YTicklabel');set(ax6,'YTickLabel',{ytl{1:end-1},''});

%% Reposition bottom plots
dx_shift = -0.16;
dx = 1;
dy_shift = -0.04;
ax5.Position = [ax5.Position(1)+dx_shift, ax5.Position(2)+dy_shift, ax5.Position(3)*dx, ax5.Position(4)];
ax6.Position = [ax6.Position(1)+dx_shift, ax6.Position(2)+dy_shift, ax6.Position(3)*dx, ax6.Position(4)];
ax7.Position = [ax7.Position(1)+dx_shift, ax7.Position(2)+dy_shift, ax7.Position(3)*dx, ax7.Position(4)];