function [ TWT ] = calcTWT( x, y, z, dvp, TAT, x_ship, y_ship, z_ship, vp_w )
% [ TWT ] = calcTWT( x, y, z, dvp, TAT, x_ship, y_ship, z_ship, vp_w )
%
% Calculate two-way travel time based on simple straight-line rays
%
% J. Russell, 2018

TWT = 2 .* sqrt((x_ship-x).^2+(y_ship-y).^2+(z_ship-z).^2) ./ (vp_w+dvp) + TAT;

end

