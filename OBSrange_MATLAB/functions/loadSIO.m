function [ datao ] = loadSIO( filename )
% [ datao ] = loadSIO( filename )
% 
% Load SIO OBS location file "filename" for location comparison - see
% example directory for format expected for these files.
% 
% J. Russell, 2018

fid = fopen(filename);
iping = 0;

% Skip to final model results
while 1
    row = fgetl(fid);
    if length(row) == length('Individual residuals (based at final drop point): ')
        if row == 'Individual residuals (based at final drop point): '
            break
        end
    end
end

while 1
    row = fgetl(fid);
    if isempty(row)
        break
    else
        iping = iping+1;
        scan = textscan(row,'%f)  Lat:%f%f(%f),  Lon:%f%f(%f), range: %fIndividual residual: %f');
        iobs(iping) = scan{1};
        lat(iping) = scan{4};
        lon(iping) = scan{7};
        range(iping) = scan{8};
        res(iping) = scan{9};
    end
end

row = fgetl(fid);
scan = textscan(row,'Initial Drop:     Lat: %f %f (%f),  Lon: %f %f (%f), depth: %f');
lat_drop = scan{3};
lon_drop = scan{6};
z_drop = scan{7};
row = fgetl(fid);
scan = textscan(row,'Final Drop:     Lat: %f %f (%f),  Lon: %f %f (%f), depth: %f');
lat_sta = scan{3};
lon_sta = scan{6};
z_sta = scan{7};
for i=1:4; row = fgetl(fid); end
scan = textscan(row,'                sqrt(sum(Residual^2/N)  : %f');
E_rms = scan{1};
for i=1:2; row = fgetl(fid); end
scan = textscan(row,'Offset Distance: Lat=%f meters, Lon=%f meters, (r=%f meters, angle=%f)');
dy_drift = scan{1};
dx_drift = scan{2};
drift = scan{3};
azi = scan{4}*-1 + 90;
azi(azi<0) = azi(azi<0) + 360;

fclose(fid);

% Data structure
datao.iobs = iobs;
datao.lat = lat;
datao.lon = lon;
datao.range = range;
datao.res = res/1000;
datao.lat_drop = lat_drop;
datao.lon_drop = lon_drop;
datao.z_drop = -z_drop;
datao.lat_sta = lat_sta;
datao.lon_sta = lon_sta;
datao.z_sta = -z_sta;
datao.E_rms = E_rms/1000;
datao.dy_drift = dy_drift;
datao.dx_drift = dx_drift;
datao.drift = drift;
datao.azi = azi;
end

