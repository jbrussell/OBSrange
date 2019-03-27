% levSSPs = getLevSSPs (typeSSP, lons, lats, dataBaseDir);
%
% interpolate SSPS derived from Levitus DataBase of typeSSP (as described below)
%  to get Sound Speed Profiles at (lons, lats)
%
% typeSSP = 0 indicates annual Levitus,
% typeSSP = 1:12 indicates Levitus for month 1:12
% 
% typeVAR = 1 indicates sound speed
%           2 indicates temperature
%           3 indicates salinity
%           4 indicates buoyancy
% 
% UNITS:  sound speed in m/s, depth in positive meters.
% 
% Modified from original for simpler application (only sound speed) by 
%  Zach Eilon on 01/25/19

function levSSPs = getLevSSPs_obsrange(typeSSP, lons, lats, dataBaseDir)

qnq='';

if isunix, slsh = '/'; else, slsh = '\'; end

typeStr = {'ann';'jan';'feb';'mar';'apr';'may';'jun';'jul';'aug';'sep';'oct';'nov';'dec'};

levFileName = sprintf ('%s%slev%s_%s.mat',dataBaseDir, slsh, qnq,typeStr{typeSSP+1});


eval (['load ' levFileName ])
% this loads in sound speed: (1000 + c*0.01) is sound speed
% UNITS:  sound speed in m/s, depth in positive meters.

levDataName = sprintf ('%s%slev_latlonZ.mat',dataBaseDir, slsh);
eval (['load ' levDataName])
lon=reshape(lon,360,180); %#ok<NODEF>
lat=reshape(lat,360,180); %#ok<NODEF>
lat=lat*0.1;
lon=lon*0.1;
% also loads in z, the levitus depths.
% this is a function so z stays in this subroutine.

% interpolate one depth at a time.
bb=[];

[n, m]=size(lons);
if m>n
lons=lons';
lats=lats';
end

for i=1:33
c1=c(:,i);
c1=reshape(c1,360,180);
bb=[bb interp2(lon',lat',c1',lons,lats,'*cubic')];
end

% convert bb to real sound speed values
% rows of bb are the sound speed profiles
bb=1000 + bb*0.01;

levSSPs = cell (size (lons));
for k=1:length(lons)
   levSSPs{k} = num2cell(bb(k,:));
end
