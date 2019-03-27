function [ v ] = pt_veloc( x, y, z, t )
% [ v ] = pt_veloc( x, y, z, t )
% 
% Calculate velocity between points in 3-D - "x,y,z" are vectors of points
% visited by the ship (units of meters) at times "t" (expecting times in
% units of seconds). "v" is a column matrix with velocities in x,y, and z
% directions in columns 1,2, and 3 respectively (velocities then in m/s).
%
% J. Russell, 2018


% Check for stationary positions
dt = diff(t);
Idt0 = (dt==0);

% Calculate velocity
v_half = [diff(x)./dt, diff(y)./dt, diff(z)./dt]';
v_half(:,Idt0') = 0; % replace stationary positions with zero velocity;
v = (v_half(:,1:end-1) + v_half(:,2:end))/2;
v = [v_half(:,1), v, v_half(:,end)]; % Replace endpoints

end

