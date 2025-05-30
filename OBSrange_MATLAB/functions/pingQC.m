function [ datao, datao_bad ] = pingQC( data, par )
% [ datao, datao_bad ] = pingQC( data, vp_w )
%
% Ping quality control function - flags bad data based on difference
% between predicted travel time for some preliminary model and the actual
% travel time. If the difference exceeds "res_thresh" (which is given in
% ms, not seconds!), then the data points are flagged as bad readings. This
% is to wean out bad readings from the deckbox (spurious reflections from
% the transducers that do happen to come in between the range gates but are
% not truly from the instruments). 
%
% 11/2019:
% In addition to applying a simple residual threshold, this version inverts
% for a preliminary nominal location and removes pings above a specified 
% RMS threshold.
% 
% INPUTS:
%   data:  structure with ship location and travel time data
%   vp_w:  default water velocity (usually 1500m/s)
%   res_thresh: residual threshold, in ms, above which trave time data wil
%               be rejeted - set this high (~500ms) so that only obvious
%               outliers are thrown out, and you don't lose lots of data
%               from a station that has drifted even a significant distance
% 
% J. Russell, Z. Eilon, S. Mosher, 2018

lat_drop = data.lat_drop;
lon_drop = data.lon_drop;
z_drop = data.z_drop;
lats_ship = data.lats;
lons_ship = data.lons;
twt = data.twt;
vp_w = par.vp_w;
res_thresh = par.res_thresh;
rms_thresh = par.rms_thresh;
TAT = par.TAT;
% dtwt_thresh = par.dtwt_thresh/1e3;
alt = data.alt;

[ x_ship, y_ship ] = lonlat2xy_nomap( lon_drop, lat_drop, lons_ship, lats_ship );
% [ x_drop, y_drop ] = lonlat2xy_nomap( lon_drop, lat_drop, lon_drop, lat_drop );
x_drop=0;
y_drop=0;

twt_pre = calcTWT(x_drop, y_drop, z_drop, 0, TAT, x_ship, y_ship, 0, vp_w);
dtwt = twt - twt_pre;

I_bad = abs(dtwt)*1000 > res_thresh | alt' < 0;
if_fix_oneway_tt = 0;
if (sum(I_bad)/length(I_bad))>0.9
    % MUST BE 1-way travel time??!!
    twt = 2*twt;
    dtwt = twt - twt_pre;
    I_bad = abs(dtwt)*1000 > res_thresh;
    if_fix_oneway_tt = 1;
end
if (sum(I_bad)/length(I_bad))>0.9
    error('Something very wrong with travel times...');
end
if if_fix_oneway_tt
    warning('Travel times were apparently one-way and have been multiplied by 2!');
end

% Do quick inversion and remove bad pings based on RMS of2 residuals for best location
m0_strt(1,1) = x_drop; %x0;
m0_strt(2,1) = y_drop; %y0;
m0_strt(3,1) = z_drop; %z0;
m0_strt(4,1) = 0; %dvp;
z_ship = zeros(size(x_ship));
azi_ship = atan2d(x_ship,y_ship);
azi_ship(azi_ship<0)=azi_ship(azi_ship<0)+360;
M = length(m0_strt);
v_ship = zeros(3,sum(~I_bad));
H = eye(M, M) .* diag([par.dampx, par.dampy, par.dampz, par.dampdvp]);
[ m_final,models,v,N,R,Cm ] = ...
        inv_newtons( par,m0_strt,twt(~I_bad),...
                    x_ship(~I_bad),y_ship(~I_bad),z_ship(~I_bad),...
                    v_ship,H,[]);
dtwt_solv = models(end).dtwt;                
I_bad_rms = abs(dtwt_solv - median(dtwt_solv)) >  rms_thresh*rms(dtwt_solv);
% I_bad_rms = I_bad_rms | (abs(dtwt_solv) >  dtwt_thresh);

% Index all bad pings
ping_nums = [1:length(dtwt)];
ping_nums_bad = ping_nums(I_bad);
ping_nums_solv = ping_nums(~I_bad);
ping_nums_rms_bad = ping_nums_solv(I_bad_rms);
I_bad_all = sort([ping_nums_bad, ping_nums_rms_bad]);
I_good_all = setdiff(ping_nums,I_bad_all);

% Plot
if 0
    figure(50); clf;
    plot(azi_ship(~I_bad),dtwt(~I_bad)*1000,'ok'); hold on;
    plot(azi_ship(I_bad),dtwt(I_bad)*1000,'o','color',[0.8 0.8 0.8]); hold on;
    plot(azi_ship(~I_bad),dtwt_solv*1000,'om');
    plot(azi_ship(I_bad_all),dtwt(I_bad_all)*1000,'.','color',[1 0 0]);
    xlabel('azimuth (deg)');
    ylabel('residual (ms)');
end

datao.lat_drop = lat_drop;
datao.lon_drop = lon_drop;
datao.z_drop = z_drop;
datao.lats = lats_ship(I_good_all);
datao.lons = lons_ship(I_good_all);
datao.t_ship = data.t_ship(I_good_all);
datao.twt = twt(I_good_all);
datao.sta = data.sta;
datao.alt = data.alt(I_good_all);

datao_bad.lats = lats_ship(I_bad_all);
datao_bad.lons = lons_ship(I_bad_all);
datao_bad.t_ship = data.t_ship(I_bad_all);
datao_bad.twt = twt(I_bad_all);

end

