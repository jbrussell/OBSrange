% function [SSP,z]=getlev_obsrange(glat,glon,typessd,ifofile,dataBaseDir)
%
% A routine to extract SSPs, temperature, salinity, or buoyancy profiles 
% along a given point or points along a path (glat,glon)
%
% Uses "getLevSSPs.m", a program to load the variables and interpolate
% onto the desired location.
%
% Gets the World Ocean Atlas files from directory '/home/dushaw/mapprogram/dataBase'
% Modify this program to set the directory where files lev_ann.mat and  
% lev_latlonZ.mat are located.
% These files are preprocessed to save disk space, and have all shortened
% profiles of the original data filled in by nearest neighbor data.
% All profiles extend to 5500 m depth, and there is world wide coverage
% (yes, including continents)
%
%   See:  http://staff.washington.edu/dushaw/WOA/
%   Data originally from:  http://www.nodc.noaa.gov/OC5/indprod.html
%
% Data at points where the original WOA has data should be unchanged from the WOA.
% The annual mean World Ocean Atlas is used.
%
%
% glat,glon are a single poin
% 
% typeSSP = 0 indicates annual Levitus,
% typeSSP = 1:12 indicates Levitus for month 1:12
% typeSSP = 13:16 indicates Levitus for Winter, Spring, Summer, or Fall  
%
% Returns P, a matrix of profiles: 33X(No. of points)
% and z, the standard 33 World Ocean Atlas depths
%
% Interpolation from World Ocean Atlas grid points to 
% desired points by cubic spline interpolation horizontally.
%
% Writes out file "ofile"; if sound speed is called, 
% this can be used in the RAY or EIGENRAY programs.
%
% The format of file templev.ssp is two columns: 
% first,            "-1  range(m)', 
% followed by   "depth(n) soundspeed (n)" for n=1:33
% and repeated for points glat(m), glon(m)
%
% UNITS:  sound speed in m/s, depth in positive meters.
% 
% Modified from original for simpler application (only sound speed) by 
%  Zach Eilon on 01/25/19

function [ssp,z,typeStr]=getlev_obsrange(glat,glon,typeSSP,ofile,dataBaseDir)

if nargin < 3
    typeSSP = 0;
end

if nargin < 4
    ofile = false;
end

if nargin < 5
    dataBaseDir='./ssp_09';
else
    dataBaseDir = [dataBaseDir,'/ssp_09'];
end

if glon(1) < 0
   glon=glon+360;
end

stdDpts = [0 10 20 30 50 75 100 125 150 200 250 300 400 500 600 700 800 900 ...
           1000 1100 1200 1300 1400 1500 1750 2000 2500 3000 3500 4000 4500 ...
           5000 5500]; 

typeStr = {'Annual'; 'January';     'February'; 'March';    'April'; ...
                     'May';         'June';     'July';     'August'; ...
                      'September';  'October';  'November'; 'December'};

typeStr = typeStr{typeSSP+1};

% obtains ssps from levitus database
predSSPs = getLevSSPs_obsrange(typeSSP, glon, glat, dataBaseDir);

ssp = [predSSPs{1}{:}]';
z=stdDpts(:);

% Finally print SSPs
if ofile~=false
    disp (['Sound speed saved in file ',ofile]);
    fid = fopen (ofile, 'w');
    fprintf (fid, 'depth(m) ssp(m/s)\n');
    for l=1:33
        fprintf (fid, '%5d %10.2f\n', z(l), ssp(l));
    end
    fclose (fid);
end


