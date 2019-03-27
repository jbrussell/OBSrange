function [ G ] = buildG( x, y, z, dvp, x_ship, y_ship, z_ship, vp_w, Nobs, M)
% [ G ] = buildG( x, y, z, dvp, x_ship, y_ship, z_ship, vp_w, Nobs, M)
% 
% Build the G matrix for 2-way travel time inversion, based on instrument
% position for this iteration (x,y,z), water speed (vp_w) and perturbation
% to water speed (dvp), and ship position (x_ship,y_ship,z_ship). This is
% done for Nobs observations and M (=4, by default) data.
% 
% J. Russell, 2018

G = zeros(Nobs,M);

% Shortcut for distance...
D = sqrt((x_ship-x).^2 + (y_ship-y).^2 + (z_ship-z).^2);

% Setup G Matrix  
% dti/dx
G(:,1) = -(x_ship-x) .* 2 ./ (vp_w+dvp) ./ D;
% dti/dy
G(:,2) = -(y_ship-y) .* 2 ./ (vp_w+dvp) ./ D;
% dti/dz
G(:,3) = -(z_ship-z) .* 2 ./ (vp_w+dvp) ./ D;
% dti/ddvp
G(:,4) = -2 .* D ./ (vp_w+dvp).^2;

end

