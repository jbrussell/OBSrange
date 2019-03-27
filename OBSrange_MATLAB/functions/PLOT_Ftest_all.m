%% F-plot plots
f101 = figure(101); clf;
set(gcf,'position',[3 402 1208 396]);
ax1 = axes('pos',[0.06 0.15 0.25 0.76]); hold(ax1,'on');
ax2 = axes('pos',[0.37 0.15 0.25 0.76]); hold(ax2,'on');
ax3 = axes('pos',[0.68 0.15 0.30 0.76]); hold(ax3,'on');

% bounds and ticks
x_bds = [min(x_grid),max(x_grid)];
y_bds = [min(y_grid),max(y_grid)];
z_bds = [min(z_grid),max(z_grid)];
% x_tix = unique(round_level([x_bds(1)-10:5:x_bds(2)+10],5));
% y_tix = unique(round_level([y_bds(1)-10:5:y_bds(2)+10],5));
% z_tix = unique(round_level([z_bds(1)-10:5:z_bds(2)+10],5));
dtick = [1,2,4,5,10,20];
dtickx = dtick(mindex(abs((x_grid(end)-x_grid(1))/8 - dtick)));
dticky = dtick(mindex(abs((y_grid(end)-y_grid(1))/8 - dtick)));
dtickz = dtick(mindex(abs((z_grid(end)-z_grid(1))/8 - dtick)));
mimax = round_level(diff(x_bds)/2,dtickx)*[-1 1];
mimay = round_level(diff(y_bds)/2,dticky)*[-1 1];
mimaz = round_level(diff(z_bds)/2,dtickz)*[-1 1];
x_tix = [mimax(1):dtickx:mimax(2)] + median(x_grid);  x_tix_lab = [mimax(1):dtickx:mimax(2)];
y_tix = [mimay(1):dticky:mimay(2)] + median(y_grid);  y_tix_lab = [mimay(1):dticky:mimay(2)];
z_tix = [mimaz(1):dtickz:mimaz(2)] + median(z_grid);  z_tix_lab = [mimaz(1):dtickz:mimaz(2)];

col = colormap(parula);
fill(ax1,[-2000,2000,2000,-2000,-2000],[-2000,-2000,2000,2000,-2000],col(1,:))
fill(ax2,[-2000,2000,2000,-2000,-2000],[-8000,-8000,8000,8000,-8000],col(1,:))
fill(ax3,[-2000,2000,2000,-2000,-2000],[-8000,-8000,8000,8000,-8000],col(1,:))

ang = [0:0.1:2*pi];
%X-Y
contourf(ax1,Xgrd(:,:,Iz_max),Ygrd(:,:,Iz_max),P(:,:,Iz_max),'linestyle','none');
plot(ax1,x_sta_bs,y_sta_bs,'.k','MarkerSize',5); hold on;
% shading(ax1,'flat');
contour(ax1,Xgrd(:,:,Iz_max),Ygrd(:,:,Iz_max),P(:,:,Iz_max),[[0.05 0.05],[0.32 0.32]],'-w','linewidth',2); hold on;
plot(ax1,Xgrd(Ix_max,Iy_max,Iz_max),Ygrd(Ix_max,Iy_max,Iz_max),'ok','markerfacecolor',[1 0 0],'markersize',12,'linewidth',1)
set(ax1,'fontsize',16,'linewidth',2,'box','on','layer','top',...
    'xlim',x_bds,'ylim',y_bds,...
    'xtick',x_tix,'xticklabel',x_tix_lab,'ytick',y_tix,'yticklabel',y_tix_lab,...
    'TickDir','out');
xlabel(ax1,'$\delta$X (m)','fontsize',18,'interpreter','latex');
ylabel(ax1,'$\delta$Y (m)','fontsize',18,'interpreter','latex');
title(ax1,'\textbf{X-Y}','fontsize',18,'interpreter','latex');
% average values
htx1 = text(ax1,x_bds(1)+0.04*diff(axlim(ax1,1:2)),y_bds(1)+0.11*diff(axlim(ax1,3:4)),...
    sprintf('$\\mathbf{\\bar{x}}$\\textbf{ = %.1f m}',median(x_grid)),'color','white','interpreter','latex','fontsize',17);
hty1 = text(ax1,x_bds(1)+0.04*diff(axlim(ax1,1:2)),y_bds(1)+0.05*diff(axlim(ax1,3:4)),...
    sprintf('$\\mathbf{\\bar{y}}$\\textbf{ = %.1f m}',median(y_grid)),'color','white','interpreter','latex','fontsize',17);
% axis equal;

%X-Z
contourf(ax2,squeeze(Xgrd(Iy_max,:,:)),squeeze(Zgrd(Iy_max,:,:)),squeeze(P(Iy_max,:,:)),'linestyle','none');
plot(ax2,x_sta_bs,z_sta_bs,'.k','MarkerSize',5); hold on;
% shading(ax1,'flat');
contour(ax2,squeeze(Xgrd(Iy_max,:,:)),squeeze(Zgrd(Iy_max,:,:)),squeeze(P(Iy_max,:,:)),[[0.05 0.05],[0.32 0.32]],'-w','linewidth',2); hold on;
plot(ax2,Xgrd(Ix_max,Iy_max,Iz_max),Zgrd(Ix_max,Iy_max,Iz_max),'ok','markerfacecolor',[1 0 0],'markersize',12,'linewidth',1)
set(ax2,'fontsize',16,'linewidth',2,'box','on','layer','top',...
    'xlim',x_bds,'ylim',z_bds,...
    'xtick',x_tix,'xticklabel',x_tix_lab,'ytick',z_tix,'yticklabel',z_tix_lab,...
    'TickDir','out');
xlabel(ax2,'$\delta$X (m)','fontsize',18,'interpreter','latex');
ylabel(ax2,'$\delta$Z (m)','fontsize',18,'interpreter','latex');
title(ax2,'\textbf{X-Z}','fontsize',18,'interpreter','latex');
% average values
htx2 = text(ax2,x_bds(1)+0.04*diff(axlim(ax2,1:2)),z_bds(1)+0.11*diff(axlim(ax2,3:4)),...
    sprintf('$\\mathbf{\\bar{x}}$\\textbf{ = %.1f m}',median(x_grid)),'color','white','interpreter','latex','fontsize',17);
htz2 = text(ax2,x_bds(1)+0.04*diff(axlim(ax2,1:2)),z_bds(1)+0.05*diff(axlim(ax2,3:4)),...
    sprintf('$\\mathbf{\\bar{z}}$\\textbf{ = %.1f m}',median(z_grid)),'color','white','interpreter','latex','fontsize',17);
% axis equal;

%Y-Z
contourf(ax3,squeeze(Ygrd(:,Ix_max,:)),squeeze(Zgrd(:,Ix_max,:)),squeeze(P(:,Ix_max,:)),'linestyle','none');
plot(ax3,y_sta_bs,z_sta_bs,'.k','MarkerSize',5); hold on;
% shading(ax1,'flat');
contour(ax3,squeeze(Ygrd(:,Ix_max,:)),squeeze(Zgrd(:,Ix_max,:)),squeeze(P(:,Ix_max,:)),[[0.05 0.05],[0.32 0.32]],'-w','linewidth',2); hold on;
plot(ax3,Ygrd(Ix_max,Iy_max,Iz_max),Zgrd(Ix_max,Iy_max,Iz_max),'ok','markerfacecolor',[1 0 0],'markersize',12,'linewidth',1)
set(ax3,'fontsize',16,'linewidth',2,'box','on','layer','top','xdir','reverse',...
    'xlim',y_bds,'ylim',z_bds,...
    'xtick',y_tix,'xticklabel',y_tix_lab,'ytick',z_tix,'yticklabel',z_tix_lab,...
    'TickDir','out');
xlabel(ax3,'$\longleftarrow$ $\delta$Y (m)','fontsize',18,'interpreter','latex');
ylabel(ax3,'$\delta$Z (m)','fontsize',18,'interpreter','latex');
title(ax3,'\textbf{Y-Z}','fontsize',18,'interpreter','latex');
% average values
hty3 = text(ax3,y_bds(2)-0.04*diff(axlim(ax3,1:2)),z_bds(1)+0.11*diff(axlim(ax3,3:4)),...
    sprintf('$\\mathbf{\\bar{y}}$\\textbf{ = %.1f m}',median(y_grid)),'color','white','interpreter','latex','fontsize',17);
htz3 = text(ax3,y_bds(2)-0.04*diff(axlim(ax3,1:2)),z_bds(1)+0.05*diff(axlim(ax3,3:4)),...
    sprintf('$\\mathbf{\\bar{z}}$\\textbf{ = %.1f m}',median(z_grid)),'color','white','interpreter','latex','fontsize',17);
% axis equal;
colorbar('peer',ax3,'linewidth',2)
