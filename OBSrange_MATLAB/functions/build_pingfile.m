function [] = build_pingfile(stationcode,lat_drop,lon_drop,z_drop,tt,lat_source,lon_source,datetime_source,pingfile_out)
%[] = build_pingfile(stationcode,lat_drop,lon_drop,z_drop,tt,lat_source,lon_source,datetime_source,pingfile_out)
%
% Generate a ping file with format based on SIO acoustic survey codes. 
% Requires lat/lon/z of drop position, GPS coordinates and timing of pings,
% and travel times.
% 
% INPUTS
% stationcode:     station name
% lat_drop:        latitude of drop position
% lon_drop:        longitude of drop position
% z_drop:          water depth at drop position
% tt:              traveltime in milliseconds (two-way if acoustic transponder, one-way if active source)
% lat_source:      latitude of "source" (i.e., GPS/transponder, air gun, etc.)
% lon_source:      longitude of "source" (i.e., GPS/transponder, air gun, etc.)
% datetime_source: datetime vector containing timing of GPS coordinates
% pingfile_out:    output name of .txt file
% 
% J. Russell, 2025

fid = fopen(pingfile_out,'w');
fprintf(fid,'Ranging data taken on:  %s\n',string(datetime_source(1),'yyyy-MM-dd HH:mm:ss.SSSSSS')); %2018-04-23 23:03:18.383000
fprintf(fid,'Cruise:                 obs-cruise\n');
fprintf(fid,'Site:                   %s\n',stationcode);
fprintf(fid,'Instrument:             \n');
fprintf(fid,'Drop Point (Latitude):  %.5f\n',lat_drop); %-4.88241
fprintf(fid,'Drop Point (Longitude): %.5f\n',lon_drop); %-132.68907
fprintf(fid,'Depth (meters):         %d\n',z_drop);
fprintf(fid,'Comment:                \n');
fprintf(fid,'==================================================\n\n');

% Build ping lines
for ii = 1:length(tt)
    % Degrees to minutes
    D_lat = floor(abs(lat_source(ii))) * sign(lat_source(ii));
    M_lat = 60*(lat_source(ii)-D_lat);
    if lat_source(ii)>=0
        NS_lat = 'N';
    else
        NS_lat = 'S';
    end
    if lon_source(ii)>180
        lon_source(ii) = lon_source(ii)-360;
    end
    D_lon = floor(abs(lon_source(ii))) * sign(lon_source(ii));
    M_lon = 60*(lon_source(ii)-D_lon);
    if lon_source(ii)>=0
        EW_lon = 'E';
    else
        EW_lon = 'W';
    end
    shot_timestamp = string(datetime_source(ii),'yyyy:D:HH:mm:ss');
    fprintf(fid,'%d msec. Lat: %d %.4f %s  Lon: %d %.4f %s  Alt: %.2f Time(UTC): %s\n',tt(ii),abs(D_lat),abs(M_lat),NS_lat,abs(D_lon),abs(M_lon),EW_lon,0,shot_timestamp);
end
fclose(fid);

end

