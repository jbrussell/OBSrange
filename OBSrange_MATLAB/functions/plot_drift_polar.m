function [ mean_drift_az ] = plot_drift_polar( ax,driftaz )
%[ mean_drift_az ] = plot_drift_polar( ax,driftaz )
% 
%   Plotting script to plot polar histogram of drift direction
% 
% INPUTS:
%   ax:      handle to axes in which to plot
%   driftaz: azimuth of drift (in degrees)
% 
% OUTPUTS:
%   mean_drift_az: average drift azimuth, in degrees
% 
% Z. Eilon, 2018

% polar histogram
hpol = polarhistogram(ax,d2r(driftaz),36,'Normalization','pdf','FaceColor',[0 0.1 0.8],'linewidth',0.1);

% mean drift direction
mean_drift_az = r2d(mean_ang(d2r(driftaz)));
% create little arrow...
polarplot(ax,d2r(mean_drift_az+[0 0 7 0 -7]),max(hpol.Values)*[0 0.97 0.85 0.97 0.85],'color','red','linewidth',2)

% bounding circle
polarplot(ax,d2r([0:1:360]'),ones(361,1)*max(get(ax,'rlim')),'k','linewidth',3);

% station
polarplot(ax,0,0,'sk','markerfacecolor',[0.5 0.5 0.5],'markersize',15,'linewidth',1);

set(ax,'ThetaZeroLocation','top','ThetaDir','clockwise',...
    'thetatick',[0:90:270],'thetaticklabel',{'N','E','S','W'},...
    'rticklabel',{},'fontweight','bold','fontsize',14);




end

