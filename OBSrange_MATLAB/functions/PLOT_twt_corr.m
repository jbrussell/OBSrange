%% map the ship velocities and the inversion misfits (and twt corrections) 
f3 = figure(3); clf;
set(gcf,'pos',[60 105 1100 800]);

% axes
ax1 = axes('pos',[0.1 0.56 0.37 0.38]);  hold(ax1,'on');
ax2 = axes('pos',[0.55 0.56 0.37 0.38]); hold(ax2,'on');
ax3 = axes('pos',[0.1 0.08 0.37 0.38]);  hold(ax3,'on');
ax4 = axes('pos',[0.55 0.08 0.37 0.38]);  hold(ax4,'on');

% limits
latlim = [min(lats_ship) max(lats_ship)]+[-1 1]*.2/111; % 200m outside ship circle
lonlim = [min(lons_ship) max(lons_ship)]+[-1 1]*.2/111; % 200m outside ship circle
londist = vdist(mean(lats_ship),lonlim(1),mean(lats_ship),lonlim(2));
resid_max = prctile(abs(dtwt_bs*1000),97);

% Divergent colormap
% cmap = cmocean('balance');
cmap = flip(brewermap(100,'RdBu'));

% Plot Velocities
scatter(ax1,lons_ship,lats_ship,100,sqrt(sum(v_ship.^2,1)),'o','filled','markeredgecolor',[0 0 0],'linewidth',1); hold on;
plot(ax1,mean(lon_sta),mean(lat_sta),'pk','markerfacecolor',[0.5 0.5 0.5],'markersize',20,'linewidth',1);

title(ax1,'\textbf{Survey pattern}','fontsize',18,'interpreter','latex');
xlabel(ax1,'Longitude','fontsize',18,'interpreter','latex');
ylabel(ax1,'Latitude','fontsize',18,'interpreter','latex');

set(ax1,'fontsize',16,'linewidth',2,'box','on',...
        'xlim',lonlim,'ylim',latlim); grid(ax1,'on');
cb1 = colorbar('peer',ax1,'linewidth',2);
ylabel(cb1,'\textbf{Ship velocity (m/s)}','interpreter','latex');
axis(ax1,'equal');
xlim(ax1,lonlim);
ylim(ax1,latlim);

% Plot twt correction at each site
scatter(ax2,lons_ship,lats_ship,100,dtwtcorr_bs*1000,'o','filled','markeredgecolor',[0 0 0],'linewidth',1); hold on;
plot(ax2,mean(lon_sta),mean(lat_sta),'pk','markerfacecolor',[0.5 0.5 0.5],'markersize',20,'linewidth',1);

title(ax2,'\textbf{Doppler corrections}','fontsize',18,'interpreter','latex');
xlabel(ax2,'Longitude','fontsize',18,'interpreter','latex');
ylabel(ax2,'Latitude','fontsize',18,'interpreter','latex');

set(ax2,'fontsize',16,'linewidth',2,'box','on',...
        'xlim',lonlim,'ylim',latlim); grid(ax2,'on');
colormap(ax2,cmap);
cb2 = colorbar('peer',ax2,'linewidth',2);
ylabel(cb2,'\textbf{Doppler correction (ms)}','interpreter','latex');
caxis(ax2,[-max(abs(dtwtcorr_bs)*1000) max(abs(dtwtcorr_bs)*1000)])
% cb2.Ticks = linspace(-max(cb2.Ticks),max(cb2.Ticks),11);
axis(ax2,'equal');
xlim(ax2,lonlim);
ylim(ax2,latlim);


% Plot residuals at each site
scatter(ax3,lons_ship,lats_ship,100,dtwt_bs*1000,'o','filled','markeredgecolor',[0 0 0],'linewidth',1); hold on;
plot(ax3,mean(lon_sta),mean(lat_sta),'pk','markerfacecolor',[0.5 0.5 0.5],'markersize',20,'linewidth',1);

title(ax3,'\textbf{Travel time residuals}','fontsize',18,'interpreter','latex');
xlabel(ax3,'Longitude','fontsize',18,'interpreter','latex');
ylabel(ax3,'Latitude','fontsize',18,'interpreter','latex');

set(ax3,'fontsize',16,'linewidth',2,'box','on',...
        'xlim',lonlim,'ylim',latlim); grid(ax3,'on');
colormap(ax3,cmap);
cb3 = colorbar('peer',ax3,'linewidth',2);
ylabel(cb3,'\textbf{Travel time residual (ms)}','interpreter','latex');
caxis(ax3,[-resid_max resid_max])
% cb3.Ticks = linspace(-max(cb3.Ticks),max(cb3.Ticks),max(cb3.Ticks)+1);
axis(ax3,'equal');
xlim(ax3,lonlim);
ylim(ax3,latlim);


% Plot residuals as a function of azimuth
obs_num = 1:length(dtwt_bs);
scatter(ax4,azi_locs{1},dtwt_bs*1000,100,dtwtcorr_bs*1000,'o','filled','markeredgecolor',[0 0 0],'linewidth',1); hold on

title(ax4,'\textbf{Travel time residuals}','fontsize',18,'interpreter','latex');
xlabel(ax4,'Ship Azimuth (deg)','fontsize',18,'interpreter','latex');
ylabel(ax4,'Travel time residuals (ms)','fontsize',18,'interpreter','latex');
set(ax4,'fontsize',16,'linewidth',2,'box','on'); grid(ax4,'on');
colormap(ax4,cmap)
cb = colorbar('peer',ax4,'linewidth',2);
ylabel(cb,'\textbf{Doppler correction (ms)}','interpreter','latex');
caxis(ax4,[-max(abs(dtwtcorr_bs)*1000) max(abs(dtwtcorr_bs)*1000)])
% cb.Ticks = linspace(-max(cb.Ticks),max(cb.Ticks),11);
% axis(ax4,'equal');
