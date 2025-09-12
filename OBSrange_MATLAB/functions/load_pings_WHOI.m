function [ data ] = load_pings_WHOI( filename )
%[ data ] = load_pings( filename )
%
% Load ping data in file called "filename", with format based on WHOI
% acoustic survey codes. 
% See Hooft_Galapagos_Plume_T104_44_20230329_145238.txt
% for expected format. 
%
% Anant Hariharan, 2025

fid = fopen(filename);
iping = 0;

% Read header info
for i=1:6 row = fgetl(fid); end
hdr = textscan(row,'# Instrument: %s');
sta = cell2mat(hdr{1}); 
row = fgetl(fid); 
hdr = textscan(row,'# Drop-Point Latitude: %f%f%s');
if isempty(hdr{2})&&isempty(hdr{3})
    lat_drop = hdr{1};
elseif isempty(hdr{3})&&~isempty(hdr{2})
    lat_drop = sign(hdr{1})*(abs(hdr{1}) + hdr{2}/60);
else
    lat_drop = (hdr{1}+hdr{2}/60)*(-1).^(1+double(strcmp(hdr{3},'N')));
end
row = fgetl(fid);
hdr = textscan(row,'# Drop-Point Longitude: %f%f%s');
if isempty(hdr{2})&&isempty(hdr{3})
    lon_drop = hdr{1};
elseif isempty(hdr{3})&&~isempty(hdr{2})
    lon_drop = sign(hdr{1})*(abs(hdr{1}) + hdr{2}/60);
else    
    lon_drop = (hdr{1}+hdr{2}/60)*(-1).^(1+double(strcmp(hdr{3},'E')));
end

% lon_drop = hdr{1};
row = fgetl(fid);
hdr = textscan(row,'# Drop-Point Depth (m): %f');
z_drop = -hdr{1};
for i=1:8; row = fgetl(fid); end

% Read pings
while 1
    row = fgetl(fid);
    if row == -1
        break
    elseif  any(regexp(row,'#'))
        disp('Bad ping denoted by ''*''');
        continue
    elseif length(row)<10
        continue
    elseif strcmp(row(1:5),'Event')
        continue
    else
        iping = iping+1;
        scan = textscan(row,'%s %s %s %s %f %f %s %f %f %s %s %f');
        % check full line of text!
        if isempty(scan{end}), continue; end
        

        % What information do we want? 
        % Latitude, Longitude, 2-Way travel time, Ship time, Station, 
        % Drop Longitude, Drop Latitude, Drop Depth, and 'alt', which is
        % likely not used? 

        sysdate = scan{1};


        dt = datetime(sysdate, 'InputFormat', 'yyyy-MM-dd');
        % Extract Julian day (day of year)
        JDAY = day(dt, 'dayofyear');


        timeutc = strrep(scan{2},',','');
        timeutc=timeutc{1};
        HOUR = str2num(timeutc(1:2));
        MIN = str2num(timeutc(4:5));
        SEC = str2num(timeutc(7:8));
      
        rng_m = scan{3};
        rng_ms_str = strrep(scan{4},',','');
        rng_ms = str2num(rng_ms_str{1});
        msec(iping) = rng_ms;
        
        lat_deg(iping) = scan{5};
        lat_min(iping) = scan{6};
        N_S(iping) = cell2mat(strrep(scan{7},',',''));
        lon_deg(iping) = scan{8};
        lon_min(iping) = scan{9};
        E_W(iping) = cell2mat(strrep(scan{10},',',''));

        t_ship(iping,1) = JDAY*24*60*60 + HOUR*60*60 + MIN*60 + SEC;

    end
end
fclose(fid);

% Convert to degrees
lat = lat_deg+lat_min/60;
lon = lon_deg+lon_min/60;


% Convert from N-S, E-W to +/-
lat(N_S=='S') = lat(N_S=='S')*-1;
lon(E_W=='W') = lon(E_W=='W')*-1;

% data
data.lats = lat';
data.lons = lon';
data.twt = msec'/1000;
data.t_ship = t_ship;
data.sta = sta;
data.lon_drop = lon_drop;
data.lat_drop = lat_drop;
data.z_drop = z_drop;
data.alt = 99999*ones(size(data.twt))';

end

