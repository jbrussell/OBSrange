function [Ftest_res] = f_test_gridsearch(par,x_ship,y_ship,x_sta_bs,y_sta_bs,z_sta_bs,V_w_bs,TAT_bs,v_eff_bs,twtcorr_bs,ifplot)
% [Ftest_res] = f_test_gridsearch(par,x_ship,y_ship,x_sta_bs,y_sta_bs,z_sta_bs,V_w_bs,TAT_bs,v_eff_bs,twtcorr_bs,ifplot)
%  

% some parms

% set grid spacing
ngridpts = 40;
grdswidthf = 4.5; 

% work out some things
dvp_bs = V_w_bs - par.vp_w;

D = max([ std(x_sta_bs) std(y_sta_bs) std(z_sta_bs) ]*grdswidthf);
Dx = D;
Dy = D; 
Dz = D;
dx = 2*Dx/ngridpts;
dy = 2*Dy/ngridpts;
dz = 2*Dz/ngridpts;
x_grid = [mean(x_sta_bs)-Dx:dx:mean(x_sta_bs)+Dx];
y_grid = [mean(y_sta_bs)-Dy:dy:mean(y_sta_bs)+Dy];
z_grid = [mean(z_sta_bs)-Dz:dz:mean(z_sta_bs)+Dz];
Nx = length(x_grid);
Ny = length(y_grid);
Nz = length(z_grid);

Nobs = length(x_ship);
z_ship = zeros(size(x_ship));

% residual of bootstrap mean result
twt_pre_bs = calcTWT(mean(x_sta_bs), mean(y_sta_bs), mean(z_sta_bs), mean(dvp_bs), mean(TAT_bs), x_ship, y_ship, z_ship, par.vp_w);
resid_bs = twtcorr_bs-twt_pre_bs;

% Determine the eigenvectors for z_sta, V_w, and TAT
X = [z_sta_bs, V_w_bs, TAT_bs];
[V, ~] = eig(X'*X);
eigvec1 = V(:,1); % Closest to TAT axis
eigvec2 = V(:,2); % Closest to V_w axis
eigvec3 = V(:,3); % Closest to z_sta axis
eig3_z = eigvec3(1);
eig3_vw = eigvec3(2);
eig3_TAT = eigvec3(3);

% % Grid search
[Xgrd,Ygrd,Zgrd] = meshgrid(x_grid,y_grid,z_grid);

dzs =  Zgrd - mean(z_sta_bs);
dvw = (eig3_vw/eig3_z)*dzs;
dTAT = (eig3_TAT/eig3_z)*dzs;
% make 4-d grids to do all in one...
% each grid will be a [Nx,Ny,Nz,Nobs] 4-D matrix
Xgrd4 = repmat(Xgrd,[1 1 1 Nobs]); Ygrd4 = repmat(Ygrd,[1 1 1 Nobs]); Zgrd4 = repmat(Zgrd,[1 1 1 Nobs]); 
dvw4 =  repmat(dvw,[1 1 1 Nobs]); dTAT4 = repmat(dTAT,[1 1 1 Nobs]); 
x_ship4 = permute(reshape(permute(repmat(x_ship',[Nx,Ny,Nz,1]),[1 3 2]),[Nx Ny Nobs Nz]),[1 2 4 3]);
y_ship4 = permute(reshape(permute(repmat(y_ship',[Nx,Ny,Nz,1]),[1 3 2]),[Nx Ny Nobs Nz]),[1 2 4 3]);
z_ship4 = permute(reshape(permute(repmat(z_ship',[Nx,Ny,Nz,1]),[1 3 2]),[Nx Ny Nobs Nz]),[1 2 4 3]);
twt_pre_gs_grd4 = calcTWT(Xgrd4, Ygrd4, Zgrd4, mean(dvp_bs)+dvw4, mean(TAT_bs)+dTAT4, x_ship4, y_ship4, z_ship4, par.vp_w);
twtcorr_bs_grd4 = permute(reshape(permute(repmat(twtcorr_bs',[Nx,Ny,Nz,1]),[1 3 2]),[Nx Ny Nobs Nz]),[1 2 4 3]);
resid_gs4 = twtcorr_bs_grd4-twt_pre_gs_grd4;
resid_bs4 = permute(reshape(permute(repmat(resid_bs',[Nx,Ny,Nz,1]),[1 3 2]),[Nx Ny Nobs Nz]),[1 2 4 3]);
chi2_gs = sum((resid_gs4.^2),4);
chi2_bs = sum((resid_bs4.^2),4);
E_gs = sqrt(chi2_gs./Nobs);

% notes from Ftest_dof function:
%[ P ] = ftest( res1,parms1,res2,parms2 )
%
% function to perform a formal F-test on two sets of residuals of fits to
% data (res1, res2) that have respectively been fitted using parms1,parms2
% (where these are integers equal to the number of parameters used in each
% fit)
% IMPORTANTLY: 
%   This will test whether the second fit (yielding res2) is statistically
%   superior to the first fit, but used more parms. I.e.:
%          sum(abs(res2)) < sum(abs(res1)) 
%   and            parms2 > parms1
%
% the degrees of freedom for each fit are therefore:
% 1) N - parms1
% 2) N - parms2
% where N is the number of data, or length(res*)
% The residuals are just equal to dobs-dpred, so we square and sum these to
% get the chi^2 values (Ea)
% 
% Z. Eilon, 2015
%
% J. Russell : This version takes the degrees of freedom as input (instead of model parameters)
% for the case where v = N - M is not accurate.
% P > 1 : Model 1 actually fits the data better than Model 2 (model 1 smaller chi^2)
% P = 1 : Model 1 fits the data same as model 2 (same chi^2)
% P < 1 : Model 2 fits the data better than model 1 (model 2 smaller chi^2)
% P < 0.05 : Model 2 fits the data better model 1 with 95% confidence


% calculate chi^2 sums
% Ea_1 = sum(res1.^2);
% Ea_2 = sum(res2.^2);
% 
% Fobs = (Ea_1/v1)/(Ea_2/v2);
% 
% P = 1 - ( fcdf(Fobs,v1,v2) - fcdf(1/Fobs,v1,v2) );

Fobs = chi2_gs./chi2_bs;

v_eff = mean(v_eff_bs);

P = 1 - ( fcdf(Fobs,v_eff,v_eff) - fcdf(1/Fobs,v_eff,v_eff) );


[Pz_max, Iz_max] = max(max(max(P)));
[Py_max, Iy_max] = max(max(P(:,:,Iz_max)));
[Px_max, Ix_max] = max(P(:,Iy_max,Iz_max));

% statistics on F-test surface
f_test_err = f_test_err_est(P,[68 95],x_grid,x_sta_bs,y_grid,y_sta_bs,z_grid,z_sta_bs);

Ftest_res = struct('x_grid',x_grid,'y_grid',y_grid,'z_grid',z_grid,...
                  'Pstat',P,'Erms',E_gs,'uncertainties',f_test_err);

if ifplot
PLOT_Ftest_all
end 

end

