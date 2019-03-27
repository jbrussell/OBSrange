function [dE,dN] = GPS_transp_correction(dforward,dstarboard,COG)
%[dE,dN] = GPS_transp_correction(dforward,dstarboard,COG)
% 
%  This function computes the location corrections for the case where the
%  GPS and hull/side transponder are not co-located. 
%  The values dE and dN are distances in METERS that should be ADDED to the
%  locations that come from the GPS to give the location of the
%  transponder. They are positive East and North, respectively.
%  The values dforward and dstarboard are positive if the transponder is
%  further forward (towards the prow) of the ship than the GPS, and if the
%  transponder is further starboard than the GPS, respectively. COG is the
%  ship heading ("course over ground") in DEGREES.
% 
%  For instance, in the Chagall-esque artwork below, for a ship (sailing 
%  up the page) with transponder (T) and GPS (G): 
%   the value of dforward is positive (+6 ASCII units) and 
%   the value of dstarboard is negative (-5 ASCII units)
%         __
%       /    \  
%      /      \
%     /        \
%    |  T       | 
%    |          | 
%    |          | 
%    |          | 
%    |          | 
%    |          | 
%    |       G  | 
%    |__________| 
%       
% Z. Eilon, 01/2019      

% calculate azimuth from GPS to transponder in ship ref frame
theta = atan2d(dstarboard,dforward);

% calculate azimuth from GPS to transponder in geographic ref frame
phi = theta + COG;

% calculate absolute horizontal distance from GPS to transponder
dr = sqrt(dforward.^2 + dstarboard.^2);

% calculate East and North offset of transponder from GPS
dE = dr(:).*sind(phi(:));
dN = dr(:).*cosd(phi(:));




end

