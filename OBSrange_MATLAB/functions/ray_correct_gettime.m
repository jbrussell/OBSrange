function dT_ray_v_line = ray_correct_gettime(dt_rvl,r)
% dT_ray_v_line = ray_correct_gettime(dt_rvl,r)
% 
% Function to extract ray correction times (to account for ray bending)
% from a lookup table of correction times as a function of horizontal
% offset (spaced every 5 meters) and depth (spaced every 20 m).
% 
% These times are given in units of seconds and are POSITIVE if ray is
% slower than straight line approx.
% 
% We should therefore SUBTRACT these from the observed travel times to
% correct from bent rays to straight lines. Having turned data from rays to
% lines, the code - which assumes lines - can then invert for position.
% 
% Z. Eilon, 2019


dh_ship = sqrt(sum(r(1:2,:).^2));
dz_ship = r(3,:);

% extract values from the offsets and depths that are closest to the
% solved-for values (interpolation would be more accurate, but slower)
indx=mindex(dt_rvl.Dx_grid_m(:),dh_ship);
indz=mindex(dt_rvl.dz_grid_km(:),dz_ship/1000);

dT_ray_v_line = dt_rvl.dT_grid_ms(indx',indz');
dT_ray_v_line = dT_ray_v_line(:,1)/1000; % convert from ms to s


end

