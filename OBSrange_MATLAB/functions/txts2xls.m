function [ ] = txts2xls( path2txts,ofile )
% Build an excel spreadsheet from individual .txt files after running
% OBSrange
%
% 2019

% ofile = 'OldORCA_stalocs.xls';
% files = dir('./OldORCA/output/OUT_nocorr/txts/*.txt');
files = dir([path2txts,'/*.txt']);
for ifil = 1:length(files)
    file = [files(ifil).folder,'/',files(ifil).name];
    fid = fopen(file);
    fgetl(fid);
    row = fgetl(fid);
    temp = textscan(row,'Station: %s');
    sta{ifil,1} = temp{1}{:};
    row = fgetl(fid);
    temp = textscan(row,'Lat:   %f deg (%f) ');
    lat(ifil,1) = temp{1};
    row = fgetl(fid);
    temp = textscan(row,'Lon:   %f deg (%f) ');
    lon(ifil,1) = temp{1};
    for ii=1:2; fgetl(fid); end
    row = fgetl(fid);
    temp = textscan(row,'Depth: %f m (%f) ');
    depth(ifil,1) = temp{1};
    fgetl(fid);
    row = fgetl(fid);
    temp = textscan(row,'Water Vel.: %f m/s (%f)');
    vp(ifil,1) = temp{1};
    for ii=1:2; fgetl(fid); end
    row = fgetl(fid);
    temp = textscan(row,'Drift:    %f m (%f) ');
    drift(ifil,1) = temp{1};
    row = fgetl(fid);
    temp = textscan(row,'Drift Azi: %f deg (%f)');
    azi(ifil,1) = temp{1};
    fclose(fid);
end

vars = {'Station','Lat','Lon','Depth','Vp','Drift','Azimuth'};
units = {'','deg','deg','m','m/s','m','deg'};
T = table(sta, lat, lon, depth, vp, drift, azi,'VariableNames',vars );
T.Properties.VariableUnits = units;
writetable(T,ofile,'Sheet',1);

end