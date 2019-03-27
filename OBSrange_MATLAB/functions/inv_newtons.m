function [ m_final,models,v_eff,N,R,Covm ] = inv_newtons( par,m_start,twt,x_ship,y_ship,z_ship,v_ship,H,dt_rvl)
%[ m_final,models,v_eff ] = inv_newton( par,m_start,twt,x_ship,y_ship,z_ship,v_ship,H )
% 
%   Conduct interative inversion for OBS location and other parameters
%   using Newton's method to step down an error gradient according to
%   partial derivatives of the forward model with respect to the model
%   parameters. 
% 
% INPUTS:
%   m_start:  starting model (x,y,z location of OBS, and dvp)
%   x_ship:   x position of ship relative to drop point, in m
%   y_ship:   y position of ship relative to drop point, in m
%   z_ship:   z position of ship relative to drop point, in m (should be zero!)
%   v_ship:   absolute velocity of ship, in m/s
%   H:        norm damping matrix, with weights set in the main script
%   dt_rvl:   structure with correction times for ray-bending effects
% OUTPUTS:
%   m_final:  output model after Newton's method iterations cease 
%   models:   output structure with all models from the iterative inversion
%             (and their associated data errors)
%   v_eff:    effective number of degrees of freedom for the inverted data
%   N:        data resolution matrix
%   R:        model resolution matrix
%   Covm:     model covariance matrix
% 
% J. Russell & Z. Eilon, 2018

vp_w = par.vp_w;
E_thresh = par.E_thresh;
epsilon = par.epsilon;
if_twtcorr = par.if_twtcorr;

Nobs = length(twt);
M = length(m_start);

dE = 1000; % initialize
iiter = 0;
models = struct([]);
m1 = m_start(:);
while dE > E_thresh
    % Update from previous iteration
    m0 = m1; 
    
    iiter = iiter+1;
    x0 = m0(1);
    y0 = m0(2);
    z0 = m0(3);
    dvp0 = m0(4);

    % Apply correction to two-way travel time due to ship velocity
    r = [x_ship-x0, y_ship-y0, z_ship-z0]'; % ship distance from OBS (maybe use drop point instead of updating)
    rmag = sqrt(sum(r.^2,1));
    r_hat = r./[rmag; rmag; rmag]; % direction vectors from OBS to ship
    vr = sum(v_ship.*r_hat,1)'; % ship velocity in the radial direction
    dr = vr .* twt; % difference in distance from ship to OBS before and after ping
    dtwtcorr = dr./(vp_w + dvp0); % traveltime correction due to velocity
    
    % apply correction (+ive, -ive, or not at all)
    % (+) if logging ship location at receive time (*)
    % (-) if logging ship location at transmit time
    twt_corr = twt + if_twtcorr*dtwtcorr; % Apply correction to the data
    
    % Apply ray bending correction
    if par.if_raycorrect && ~isempty(dt_rvl)
        dT_ray_v_line = ray_correct_gettime(dt_rvl,r); % positive if ray slower than line
        twt_corr = twt_corr - dT_ray_v_line; 
    end

    % Build the G matrix
    G = buildG_noTAT( x0, y0, z0, dvp0, x_ship, y_ship, z_ship, vp_w, Nobs, M);

    % Set up norm damping for each parameter
    h = zeros(M,1);

    % Predicted twt for this iteration
    twt_pre = calcTWT(x0, y0, z0, dvp0, par.TAT, x_ship, y_ship, z_ship, vp_w);

    % Residual
    dtwt = twt_corr-twt_pre;
    E = sqrt(sum(dtwt.^2)/Nobs); % RMS error

    % Least squares solution
    J = eye(M,M)*epsilon^0.5; % global damping
    F = [G; H; J];
    f = [dtwt; h; zeros(M,1)];
    Finv = (F'*F)\F';
    m1 = m0 + Finv*f;
    
    Ginv = (G'*G + H'*H + J'*J)\G';

    % Calculate the effective degrees of freedom (for F-test)
    v_eff = length(f) - trace(F*Finv);
    % Exclude 0=0 constraint equations from v_eff
    v_eff = v_eff - sum(sum(H,2)==0) - sum(sum(J,2)==0);
    
    models(iiter).m = m0;
    models(iiter).E = E;
    models(iiter).dtwt = dtwt;
    models(iiter).dtwtcorr = dtwtcorr;
    models(iiter).twt_corr = twt_corr;
    models(iiter).vr = vr;

    if iiter > 1
        dE = models(iiter-1).E - models(iiter).E;
    end
end
m_final = m0;

% Data Resolution
N = G*Ginv;
% Model Resolution
R = Ginv*G;
% Model Covariance
sigma2_d = dtwt'*dtwt / v_eff;
Covm = sigma2_d * (Ginv*Ginv');
end

