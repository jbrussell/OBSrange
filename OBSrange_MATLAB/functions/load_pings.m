function [ data ] = load_pings( filename )
%[ data ] = load_pings( filename )
%
% Load ping data in file called "filename", with format based on SIO
% acoustic survey codes. See example directory for expected format. 
%
% J. Russell, 2018

fid = fopen(filename);
iping = 0;

% Read header info
for i=1:3; row = fgetl(fid); end
hdr = textscan(row,'Site: %s');
sta = cell2mat(hdr{1}); 
for i=1:2; row = fgetl(fid); end
hdr = textscan(row,'Drop Point (Latitude): %f%f%s');
if isempty(hdr{2})&&isempty(hdr{3})
    lat_drop = hdr{1};
elseif isempty(hdr{3})&&~isempty(hdr{2})
    lat_drop = sign(hdr{1})*(abs(hdr{1}) + hdr{2}/60);
else
    lat_drop = (hdr{1}+hdr{2}/60)*(-1).^(1+double(strcmp(hdr{3},'N')));
end
% lat_drop = hdr{1};
row = fgetl(fid);
hdr = textscan(row,'Drop Point (Longitude): %f%f%s');
if isempty(hdr{2})&&isempty(hdr{3})
    lon_drop = hdr{1};
elseif isempty(hdr{3})&&~isempty(hdr{2})
    lon_drop = sign(hdr{1})*(abs(hdr{1}) + hdr{2}/60);
else    
    lon_drop = (hdr{1}+hdr{2}/60)*(-1).^(1+double(strcmp(hdr{3},'E')));
end
% lon_drop = hdr{1};
row = fgetl(fid);
hdr = textscan(row,'Depth (meters): %f');
z_drop = -hdr{1};
for i=1:3; row = fgetl(fid); end

% Read pings
while 1
    row = fgetl(fid);
    if row == -1
        break
    elseif  any(regexp(row,'*'))
        disp('Bad ping denoted by ''*''');
        continue
    elseif length(row)<10
        continue
    elseif strcmp(row(1:5),'Event')
        continue
    else
        iping = iping+1;
        scan = textscan(row,'%f%s%s%f%f%s%s%f%f%s%s%f%s%f:%f:%f:%f:%f');
        % check full line of text!
        if isempty(scan{18}), continue; end
        msec(iping) = scan{1};
        lat_deg(iping) = scan{4};
        lat_min(iping) = scan{5};
        N_S(iping) = cell2mat(scan{6});
        lon_deg(iping) = scan{8};
        lon_min(iping) = scan{9};
        E_W(iping) = cell2mat(scan{10});
        YEAR = scan{14};
        JDAY = scan{15};
        HOUR = scan{16};
        MIN = scan{17};
        SEC = scan{18};
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

end

