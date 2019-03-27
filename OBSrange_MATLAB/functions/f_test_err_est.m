function [f_test_err] = f_test_err_est(P,p_levels,x_grid,x_sta,y_grid,y_sta,z_grid,z_sta)

[Pz_max, Iz_max] = max(max(max(P)));
[Py_max, ~] = max(max(P(:,:,Iz_max)));

f_test_err = struct([]);
hf = figure;
for ip = 1:length(p_levels)
    aval = 1 - p_levels(ip)/100;
    xyz_str = ['xyz',num2str(p_levels(ip))];
    xy_str = ['xy',num2str(p_levels(ip))];
    
    c_3D = isosurface(x_grid-mean(x_sta),y_grid-mean(y_sta),z_grid-mean(z_sta),P./Pz_max,aval);
    c_3D_amp = sqrt(sum(c_3D.vertices.^2,2));
    c_2D = contour(x_grid-mean(x_sta),y_grid-mean(y_sta),P(:,:,Iz_max)./Py_max,aval*[1,1])';
    c_2D=c_2D(2:end,:); % not sure why first row of this contour is junk.
    c_2D_amp = sqrt(sum(c_2D.^2,2));

    f_test_err(1).(xyz_str) = mean(c_3D_amp);
    f_test_err(1).(xy_str) = mean(c_2D_amp);

end
close(hf);delete(hf); 


end

