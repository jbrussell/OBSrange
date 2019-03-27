% MAIN script to run location inversion for an OBS deployment
% 
% Requires: 
%   Acoustic survey files in some directory 
%   Inputs below
% 
% Uses two-way travel time information and ship coordinates from an OBS 
% survey to invert for station location on the seafloor (Lat, Lon, Depth),
% static correction to the sound velocity through the
% water column (dvp), and the velocity of the ship in the radial
% direction of the survey circle (vr0).
%
% Josh Russell & Zach Eilon 4/16/18

clear; close all;

%% INPUTS
% path to project

% JOSH
projpath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA_v2/'; % DATA

% ZACH
% projpath = '~/Work/OBSrange/projects/PacificORCA/';

% STEVE
% projpath = '~/Seismo/projects/OBSrange/projects/PacificORCA/';

% TESTING
% projpath = '../Tests';

% path to survey data from the project directory
datapath = './'; 
% path to output directory from project directory (will be created if it does not yet exist)
outdir = './OUT_OBSrange/'; 
% Put a string station name here to only consider that station. 
% Otherwise, to locate all stations, put ''
onesta = 'EC03'; %''; %'EC03';

%% Parameters
ifsave = 1; % Save results to *.mat?
ifplot = 1; % Plot results?

par = struct([]);
par(1).vp_w = 1500; % Assumed water velocity (m/s)
par.N_bs = 1000; % Number of bootstrap iterations
par.E_thresh = 1e-5; % RMS reduction threshold for inversion

% Traveltime correction parameters
% ==>  +1 if location is RECEIVE, -1 if location is SEND, 0 if no correction
par.if_twtcorr = 0; % Apply a traveltime correction to account for ship velocity?
par.npts_movingav = 1; %5; % number of points to include in moving average smoothing of ship velocity (1 = no smoothing);

% Raybending correction parameters
par.if_raycorrect = 1; % Apply a traveltime correction to account for rays bending?
%  NOTE - if you choose to do this you can either input your own
%  depth-soundspeed profile (make sure "SSP_stationname.txt" files are
%  sitting in the directory defined below, where "stationname" is the same
%  as "sta" in the code below) OR our code will calculate one for you based
%  on the location and the date of the survey (using monthly means and a
%  decadal SSP projection)
par.sspfiledir = './OUT_OBSrange/station_soundspeed_profiles';

% GPS-transponder offset
% if this is not known, set both of these to zero
par.dforward = 0; % units of meters, positive means transponder is further forward than GPS 
par.dstarboard = 0; % units of meters, positive means transponder is further starboard than GPS

% Ping QC -- Remove pings > ping_thresh ms away from neighbor
ifQC_ping = 1; % Do quality control on pings? (We strongly recommend this!)
res_thresh = 500; % (ms) Will filter out pings with residuals > specified magnitude

% TAT - Define turnaround which remains fixed (transponder-specific value)
par.TAT = 0.013; % (s) [13 ms for Edgetech instruments]

% Norm damping for each model parameter (damping towards starting model)
% Larger values imply more damping towards the starting model.
par.dampx = 0;
par.dampy = 0;
par.dampz = 0; %1e3; %0
par.dampdvp = 5e-8; %1e3; %5e-8

% Global norm damping for stabilization
par.epsilon = 1e-10;

%% ===================================================================== %%
%% ================ NOT ADVISED TO EDIT BELOW THIS LINE ================ %%
%% ===================================================================== %%

% prepend functions directory to MATLAB path
functionspath = [pwd,'/functions'];
addpath(functionspath);
functionsSSPpath = [pwd,'/functions/ocean_profiles_obsrange'];
addpath(functionsSSPpath);

%% Load 2-way Travel Time Data
wd = pwd;
if ~exist(projpath)
    mkdir(projpath);
end
cd(projpath);
files = dir([datapath,'/*.txt']);
stas = unique(strtok({files.name},{'_','.txt'}));
Nstas = length(stas);
for is = 1:Nstas
sta = stas{is};
if ~isempty(onesta), if ~strcmp(sta,onesta),continue; end; end
fprintf('===========================\nLocating OBS %s\n\n',sta);
%% load all the data    
%     rawdatafile = sprintf('%s%s.txt',datapath,sta);
rawdatfile = dir([datapath,sta,'*']);data = [];
for irf = 1:length(rawdatfile)
	if any(regexp(rawdatfile(irf).name,'orrect')) 
		continue;
	end
data = load_pings([datapath,rawdatfile(irf).name]);        
end
if isempty(data) || isempty(data.lon_drop) || isempty(data.lat_drop)
	continue;
end
if ifQC_ping
	[dataQC, data_bad] = pingQC(data, par.vp_w, res_thresh);
	data = dataQC;
    N_badpings = length(data_bad.twt);
    fprintf('\nNumber of bad pings removed: %d\n',N_badpings);
end
lat_drop = data.lat_drop;
lon_drop = data.lon_drop;
z_drop = data.z_drop;
lats_ship = data.lats;
lons_ship = data.lons;
t_ship = data.t_ship;
twt = data.twt;

if par.if_raycorrect
    % Make soundspeed profiles (or load user's) and ray-line correction matrix 
    % NOTE - if user profile, make sure it extends down to z_drop depth
    dt_rvl = ray_correct_makegrid([lat_drop,lon_drop],par,sta,z_drop,t_ship(1),[functionspath,'/ocean_profiles_obsrange']);
    % note these will be two-way correction times, to be SUBTRACTED from
    % actual data in order to remove the effect of ray bending.
else
    dt_rvl = [];
end

% Set origin of coordinate system to be lat/lon of drop point
olon = lon_drop;
olat = lat_drop;

Nobs = length(twt);
z_ship = zeros(Nobs,1); % ship is always at surface

% Convert Lon/Lat to x/y
[ x_ship, y_ship ] = lonlat2xy_nomap( olon, olat, lons_ship, lats_ship );
[ x_drop, y_drop ] = lonlat2xy_nomap( olon, olat, lon_drop, lat_drop );

% Calculate velocity of ship
v_ship = pt_veloc( x_ship, y_ship, z_ship, t_ship );
v_ship = [moving_average(v_ship(1,:),par.npts_movingav)'; moving_average(v_ship(2,:),par.npts_movingav)'; moving_average(v_ship(3,:),par.npts_movingav)'];

% Account for GPS-transponder offset
survcog = atan2d(v_ship(1,:),v_ship(2,:));
[dx,dy] = GPS_transp_correction(par.dforward,par.dstarboard,survcog');
x_ship = x_ship + dx;
y_ship = y_ship + dy;


%% Set up initial model
m0_strt(1,1) = x_drop; %x0;
m0_strt(2,1) = y_drop; %y0;
m0_strt(3,1) = z_drop; %z0;
m0_strt(4,1) = 0; %dvp;

%% Do Bootstrap Inversion
Nobs = length(twt);
M = length(m0_strt);
% damping matrix
H = eye(M, M) .* diag([par.dampx, par.dampy, par.dampz, par.dampdvp]);

[xmat_ship_bs, ymat_ship_bs, zmat_ship_bs, vmat_ship_bs, twtmat_bs, indxs] = bootstrap(x_ship, y_ship, z_ship, v_ship, twt, par.N_bs-1);

% prepare bootstrap result vectors
dtwt_mat = zeros(size(indxs));
twtcorr_mat = zeros(size(indxs));
dtwtcorr_mat = zeros(size(indxs));
vr_mat = zeros(size(indxs));
x_sta = nan(par.N_bs,1); y_sta = nan(par.N_bs,1); z_sta = nan(par.N_bs,1);
lon_sta = nan(par.N_bs,1); lat_sta = nan(par.N_bs,1);
TAT = nan(par.N_bs,1);
dvp = nan(par.N_bs,1);
V_w = nan(par.N_bs,1);
E_rms = nan(par.N_bs,1);
v_eff = nan(par.N_bs,1);
dx_drift = nan(par.N_bs,1);
dy_drift = nan(par.N_bs,1);
dz_sta = nan(par.N_bs,1);
drift = nan(par.N_bs,1);
azi  = nan(par.N_bs,1);
dist = cell(par.N_bs,1);
azi_locs = cell(par.N_bs,1);
models = cell(par.N_bs,1);
R_mat = zeros(M,M,par.N_bs);
Cm_mat = zeros(M,M,par.N_bs);
tic
for ibs = 1:par.N_bs
    x_ship_bs = xmat_ship_bs(:,ibs);
    y_ship_bs = ymat_ship_bs(:,ibs);
    z_ship_bs = zmat_ship_bs(:,ibs);
    v_ship_bs = [vmat_ship_bs(1).amatbs(:,ibs)'; vmat_ship_bs(2).amatbs(:,ibs)'; vmat_ship_bs(3).amatbs(:,ibs)'];
    twt_bs = twtmat_bs(:,ibs);
    
    [ m_final,models,v,N,R,Cm ] = ...
        inv_newtons( par,m0_strt,twt_bs,...
                    x_ship_bs,y_ship_bs,z_ship_bs,...
                    v_ship_bs,H,dt_rvl);

    x_sta(ibs) = m_final(1);
    y_sta(ibs) = m_final(2);
    z_sta(ibs) = m_final(3);
    TAT(ibs) = par.TAT;
    dvp(ibs) = m_final(4);
    V_w(ibs) = par.vp_w + dvp(ibs);
    E_rms(ibs) = models(end).E;
    v_eff(ibs) = v;
    dtwt_mat(:,ibs) = models(end).dtwt;
    twtcorr_mat(:,ibs) = models(end).twt_corr;
    dtwtcorr_mat(:,ibs) = models(end).dtwtcorr;
    vr_mat(:,ibs) = models(end).vr;
    [lon_sta(ibs), lat_sta(ibs)] = xy2lonlat_nomap(olon, olat, x_sta(ibs), y_sta(ibs));
    R_mat(:,:,ibs) = R;
    Cm_mat(:,:,ibs) = Cm;
    
    % Calculate OBS drift distance and azimuth
    dx_drift(ibs) = m_final(1) - x_drop;
    dy_drift(ibs) = m_final(2) - y_drop;
    dz_sta(ibs) = m_final(3) - z_drop;
    drift(ibs) = sqrt(dx_drift(ibs).^2 + dy_drift(ibs).^2);
    azii = -atan2d(dy_drift(ibs),dx_drift(ibs))+90;
    azii(azii<0) = azii(azii<0) + 360;
    azi(ibs) = azii;
    
    % Calculate ship distance azimuth from OBS
    dx_ship = x_ship_bs-m_final(1);
    dy_ship = y_ship_bs-m_final(2);
    dist{ibs} = sqrt(dx_ship.^2 + dy_ship.^2);
    azi_locss = -atan2d(dy_ship , dx_ship) + 90;
    azi_locss(azi_locss<0) = azi_locss(azi_locss<0) + 360;
    azi_locs{ibs} = azi_locss;
end
toc
dtwt_bs = mean(unscramble_randmat(dtwt_mat,indxs),2);
twtcorr_bs = mean(unscramble_randmat(twtcorr_mat,indxs),2);
dtwtcorr_bs = mean(unscramble_randmat(dtwtcorr_mat,indxs),2);
vr_bs = mean(unscramble_randmat(vr_mat,indxs),2);

range = sqrt( (mean(x_sta)-x_ship).^2 + (mean(y_sta)-y_ship).^2 + (mean(z_sta)-z_ship).^2 );

%% HISTOGRAMS OF MODEL PARAMETERS
if ifplot
PLOT_histograms_all
end

%% F-test for uncertainty using grid search
tic
Ftest_res = f_test_gridsearch(par,x_ship,y_ship,x_sta,y_sta,z_sta,V_w,TAT,v_eff,twtcorr_bs,ifplot);        
toc
%% Print some results
fprintf('\nStation: %s',data.sta);
fprintf('\nlat:   %.5f deg (%f) \nlon:   %.5f deg (%f) \nx:     %.1f m (%.1f) \ny:    %.1f m (%.1f) \ndepth: %.1f m (%.1f) \nTAT:   %.1f ms \nv_H20: %.1f m/s (%.1f)',mean(lat_sta),std(lat_sta)*2,mean(lon_sta),std(lon_sta)*2,mean(x_sta),std(x_sta)*2,mean(y_sta),std(y_sta)*2,mean(z_sta),std(z_sta)*2,mean(TAT)*1000,mean(V_w),std(V_w)*2);
fprintf('\nDrift Lon: %.1f m (%.1f) \nDrift Lat: %.1f m (%.1f) \nDrift:    %.1f m (%.1f) \nDrift Azi: %.1f deg (%.1f)\ndz: %.1f m (%.1f)\n',mean(dx_drift),std(dx_drift)*2,mean(dy_drift),std(dx_drift)*2,mean(drift),std(drift)*2,mean(azi),std(azi)*2,mean(dz_sta),std(dz_sta)*2);
fprintf('\nRMS:  %.3f ms (%.3f)\n',mean(E_rms)*1000,std(E_rms)*2*1000);

%% PLOTTING
if ifplot
	%% Plot Misfit
	PLOT_misfit
	%% Geographic PLOTTING
	PLOT_survey
	PLOT_twt_corr
    %% Model Resolution & Covariance
    PLOT_resolution_covariance
    drawnow;
end


%% Save output
% output directory
if par.if_twtcorr == 1
    modified_outdir = [outdir,'OUT_wcorr_xrec/'];
elseif par.if_twtcorr == -1
    modified_outdir = [outdir,'OUT_wcorr_xsend/'];
elseif par.if_twtcorr == 0
    modified_outdir = [outdir,'OUT_nocorr/'];
end  
if ~exist(modified_outdir)
	mkdir(modified_outdir);
end

% Save textfile
fid = fopen([modified_outdir,'/',data.sta,'_location.txt'],'w');
fprintf(fid,'Bootstrap inversion results (2sigma uncertainty)');
fprintf(fid,'\nStation: %s',data.sta);
fprintf(fid,'\nLat:   %.5f deg (%f) \nLon:   %.5f deg (%f) \nX:     %.1f m (%.1f) \nY:    %.1f m (%.1f) \nDepth: %.1f m (%.1f) \nTAT:   %.1f ms (%f) \nWater Vel.: %.1f m/s (%f)',mean(lat_sta),std(lat_sta)*2,mean(lon_sta),std(lon_sta)*2,mean(x_sta),std(x_sta)*2,mean(y_sta),std(y_sta)*2,mean(z_sta),std(z_sta)*2,mean(TAT)*1000,std(TAT)*1000*2,mean(V_w),std(V_w)*2);
fprintf(fid,'\nDrift Lon: %f m (%f) \nDrift Lat: %f m (%f) \nDrift:    %f m (%f) \nDrift Azi: %f deg (%f)\ndz: %f m (%f)\n',mean(dx_drift),std(dx_drift)*2,mean(dy_drift),std(dy_drift)*2,mean(drift),std(drift)*2,mean(azi),std(azi)*2,mean(dz_sta),std(dz_sta)*2);
fprintf(fid,'\nRMS:  %.4f ms (%f)\n',mean(E_rms)*1000,std(E_rms)*2*1000);
fprintf(fid,'\nBad pings Removed: %d',N_badpings);
fprintf(fid,'\n===================================================\n');
fprintf(fid,'%10s %10s %15s %15s %15s %15s \n','Lat','Lon','Range (m)','Residual (s)','Doppler Vel. (m/s)','TWT corr. (ms)');
for ii = 1:Nobs
	fprintf(fid,'%3d: %10f %10f %10f %10f %10f %10f\n',ii,lats_ship(ii),lons_ship(ii),range(ii),dtwt_bs(ii)*1000,vr_bs(ii),dtwtcorr_bs(ii)*1000);
end
fclose(fid);

if ifsave && ifplot
% Save plots
if ~exist([modified_outdir,'/plots/'])
	mkdir([modified_outdir,'/plots/']);
end
save2pdf([modified_outdir,'/plots/',data.sta,'_1_OBSlocation.pdf'],f1,500)
save2pdf([modified_outdir,'/plots/',data.sta,'_2_misfit.pdf'],f2,500)
save2pdf([modified_outdir,'/plots/',data.sta,'_3_VelCorrs.pdf'],f3,500)
save2pdf([modified_outdir,'/plots/',data.sta,'_4_bootstrap.pdf'],f100,500)
save2pdf([modified_outdir,'/plots/',data.sta,'_5_Ftest.pdf'],101,500)
save2pdf([modified_outdir,'/plots/',data.sta,'_6_Resolution_Covariance.pdf'],f103,500)
end

if ifsave
    datamat.par = par;
    datamat.sta = data.sta;
	datamat.drop_lonlatz = [data.lon_drop,data.lat_drop,data.z_drop];
	datamat.lons_ship = lons_ship;
	datamat.lats_ship = lats_ship;
	datamat.x_ship = x_ship;
	datamat.y_ship = y_ship;
	datamat.z_ship = z_ship;
	datamat.v_ship = v_ship;
	datamat.databad = data_bad;
	datamat.dtwt_bs = dtwt_bs;
	datamat.twtcorr_bs = twtcorr_bs;
	datamat.dtwtcorr_bs = dtwtcorr_bs;
	datamat.lon_sta_bs = lon_sta;
	datamat.lat_sta_bs = lat_sta;
	datamat.x_sta_bs = x_sta;
	datamat.y_sta_bs = y_sta;
	datamat.z_sta_bs = z_sta;
	datamat.drift_bs = drift;
	datamat.azi_bs = azi;
	datamat.TAT_bs = TAT;
	datamat.V_w_bs = V_w;
	datamat.E_rms = E_rms;
    datamat.Ftest_res = Ftest_res;
	datamat.loc_xyz = [mean(x_sta),mean(y_sta),mean(z_sta)];
	datamat.loc_lolaz = [mean(lon_sta),mean(lat_sta),mean(z_sta)];
	datamat.mean_drift_az = [mean(drift) r2d(mean_ang(d2r(azi)))];
    datamat.R_mat = R_mat;
    datamat.Cm_mat = Cm_mat;
	if ~exist([modified_outdir,'/mats'])
		mkdir([modified_outdir,'/mats']);
	end
	save([modified_outdir,'/mats/',data.sta,'_data.mat'],'datamat');
end


end

% message if no success
if is==Nstas && ~exist('rawdatfile','var')
    fprintf('*********************\nStation %s does not seem to exist in that folder\n',onesta);
end
cd (wd)
