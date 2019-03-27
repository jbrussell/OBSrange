function dT_ray_v_line_grid = ray_correct_makegrid(lalo,par,sta,z_sta,t_julday_s,ocean_profile_dir)

z_sta = abs(z_sta);
if z_sta > 6 % infer depth in m when should be in km
    z_sta = z_sta/1000;
end

%% Write (if not there) ssp file for this station
% check ssp file directory exists
if exist(par.sspfiledir,'dir')~=7
    mkdir(par.sspfiledir);
end

sspfile = [par.sspfiledir,'/SSP_',sta,'.txt'];

% if ssp depth file does not exist for this station...
if exist(sspfile,'file')~=2

    %% month of deployment
    % t_julday_s is seconds since year beginning
    % good enough to round to nearest center of month
    month_deploy = ceil(12*t_julday_s/(365.25*24*3600));
    % write ssp file
    addpath(ocean_profile_dir);
    getlev_obsrange(lalo(1),lalo(2),month_deploy+1,sspfile,ocean_profile_dir);
end

%% Load ssp file
[ssp,z] = load_sspfile(sspfile);
v_profile = struct('ssp',ssp/1000,'z',z/1000);


%% Ray shooting 

% loop through ray parameters making grid of distances and depths, with
% straight line versions of each, getting values of two-way dT delay
ps = 0.01:0.01:1.5;

zs = z_sta + [-200:20:200]/1000 ; % do shooting for depths 200m above and below starting depth

zmax = max(zs);

Np = length(ps);
Nz = length(zs);

Dr_ray = nan(Np,Nz);
Dx = nan(Np,Nz);
t_ray = nan(Np,Nz);
Dr_lin = nan(Np,Nz);
t_lin = nan(Np,Nz);
for ip = 1:Np
    p = ps(ip);   
    try
        %% compute rays
        [rx,rz,Dr,rt,rv] = shootrays(p, v_profile, zmax);
        % get intersection of this ray at each depth 
        for iz = 1:Nz
            z = zs(iz);
            % Ray vs. straight-line 
            Dr_ray(ip,iz) = interp1(rz,Dr,z); % total length of ray (in km)
            Dx(ip,iz) = interp1(rz,rx,z); % horizontal offset of ray at zmax (in km)
            t_ray(ip,iz) = interp1(rz,rt,z); % time along ray (s)
            Dr_lin(ip,iz) = sqrt(Dx(ip,iz).^2 + z.^2); % straight line start to end
            t_lin(ip,iz) = Dr_lin(ip,iz) / harmmean(rv(rz<=z)); % time along straight line
        end % loop on depths
    end % loop on try
end % loop on rayps

% difference in ray length( metres)
dDr_m = 1000*(Dr_ray-Dr_lin);

% difference in ray travel time (ms) (TWO-way)
dt_ray_ms = 2*1000*(t_ray-t_lin); 
% to be clear: this number is POSITIVE if ray is slower than straight line
% approx. We should therefore SUBTRACT this from the observed travel times
% to correct them from bent rays to straight lines. Having turned data from
% rays to lines the code, which assumes lines, can then invert for position

%% add vertical ray
Dx = [zeros(1,Nz);Dx];
dDr_m = [zeros(1,Nz);dDr_m];
dt_ray_ms = [zeros(1,Nz);dt_ray_ms];

%% interpolate tt error onto regular mesh
Dx_grid = [0:0.005:min([4,max(Dx)])]'*1000; % grid in m, spaced every 5m, up to 4 km offset
% make grid of differential travel times - columns are each depth, rows are
% different offsets
dT_grid = nan(length(Dx_grid),Nz);
for iz = 1:Nz
indx = ~isnan(dt_ray_ms(:,iz)); % ignore offsets when rays go imaginary
    dT_grid(:,iz) = interp1(Dx(indx,iz)*1000,dt_ray_ms(indx,iz),Dx_grid,'linear',nan); % times in ms
end

dT_ray_v_line_grid = struct('Dx_grid_m',Dx_grid,...
                       'dz_grid_km',zs,...
                       'dT_grid_ms',dT_grid);
                   
%% save ray correction file
corfile = [par.sspfiledir,'/cor_rayline_',sta,'.mat'];
save(corfile,'dT_ray_v_line_grid');

end

