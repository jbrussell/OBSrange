function [ssp,z] = load_sspfile(sspfile)
% [ssp,z] = load_sspfile(sspfile)
% 
%  Loads a soundspeed-depth profile from file "sspfile"
%  
%  sspfile expected to have one header file, followed by n lines
%  corresponding to columns of:
%     (1) depth in meters 
%     (2) ssp in m/s
% 
% Z. Eilon, 01/2019

fid = fopen(sspfile,'r');
DAT = textscan(fid,'%f %f','headerlines',1);
fclose(fid);

z = DAT{1};
ssp = DAT{2};

end

