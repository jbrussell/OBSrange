% Convert Brandon's airgun picks made on SPOBS at Me GUSTA! experiment to
% equivalent one-way travel time "pings". We will multiply by two in the
%
% Josh Russell, 2025
% main code.

clear; close all;

functionspath = ['../functions'];
addpath(functionspath);

% stationcode = '1AS16';
stationcode = '1AS23';
% stationcode = '5S18';

% Setup paths
path2picks = './active_source_cruise_files/picks/';
path2shotlogs = './active_source_cruise_files/shotlogs/';
path2geometry = './active_source_cruise_files/geometry/';
path2pingfile_out = './survey_files/';

% Get line number and station name
toks = strsplit(stationcode,'S');
line = toks{1};
station = ['S',toks{2}];

%% Load tables

% Load picks
pickfile = [path2picks,'/',stationcode,'_directwave.csv'];
picks = readtable(pickfile);

% Load shotlog
shotlogfile = [path2shotlogs,'/','*_OBS',line,'.shotlog.csv'];
temp = dir(shotlogfile);
shotlog = readtable([temp.folder,'/',temp.name]);

% Load geometry
geometryfile = [path2geometry,'/','*_OBS',line,'_locations.csv'];
temp = dir(geometryfile);
geometry = readtable([temp.folder,'/',temp.name]);

%% Get drop lat,lon,depth from geometry file

ista = find(strcmp(geometry.OBS,stationcode));
lat_drop = geometry.DeployLatitude(ista);
lon_drop = geometry.DeployLongitude(ista);
z_drop = geometry.elevation(ista);


%% Associate picks with shot lat/lon

Npicks = size(picks,1);
tt = nan(Npicks,1);
sourceLat = nan(Npicks,1);
sourceLon = nan(Npicks,1);
datetime = NaT(Npicks,1);
for ii = 1:Npicks
    % Get one-way travel time
    tt(ii) = picks.TIME(ii);
    
    % Index shot using FFID
    FFID = picks.FFID(ii);
    ishot = find(shotlog.shotnumber==FFID);
    sourceLat(ii) = shotlog.sourceLat(ishot);
    sourceLon(ii) = shotlog.sourceLon(ishot);
    v_datetime(ii) = shotlog.date(ishot) + shotlog.time(ishot);
end

[~,isrt] = sort(v_datetime);
v_datetime = v_datetime(isrt);
sourceLat = sourceLat(isrt);
sourceLon = sourceLon(isrt);
tt = tt(isrt);


%% Plot travel time data and track
figure(1); clf;
set(gcf,'color','w','position',[370         550        1151         471]);

subplot(1,2,1); box on; hold on;
scatter(sourceLon,sourceLat,50,datenum(v_datetime),'o','filled');
plot(lon_drop,lat_drop,'pk','linewidth',1.5,'markerfacecolor','r','MarkerSize',20);
xlabel('Longitude');
ylabel('Latitude');
set(gca,'fontsize',15,'linewidth',1.5);

subplot(1,2,2); box on; hold on;
scatter(v_datetime,tt,50,datenum(v_datetime),'o','filled');
xlabel('Date');
ylabel('Travel time picks (ms)');
set(gca,'fontsize',15,'linewidth',1.5);

%% Write picks to ping file

pingfile_out = [path2pingfile_out,'/',stationcode,'.txt'];

build_pingfile(stationcode,lat_drop,lon_drop,z_drop,tt,sourceLat,sourceLon,v_datetime,pingfile_out)

