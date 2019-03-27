'''
FUNCTION ray_correct.py

A function and subfunctions to compute travel-time differences between straight
rays and refracting rays by ray-tracing through a 1D velocity profile.

Zach Eilon and Stephen Mosher 2/19/19
'''
# Import modules and functions
import os
import time
import pickle
import numpy as np
from scipy.stats import hmean
from scipy.interpolate import interp1d
from funcs import get, shootrays

from IPython.core.debugger import Tracer

def load_ssp(ssp_fname):
  # Declare lists for velocities and depths. 
  ssp, z = [], []
  
  # Open SSP and read in velocities and depths from each line. Skip first line.
  f = open(ssp_fname, 'r')
  lines = f.readlines()
  for i, line in enumerate(lines):
    if i == 0:
      continue
    z.append(line.split()[0])
    ssp.append(line.split()[1])
  f.close()

  return ssp, z

def makegrid(lat, lon, z_sta, stn, t, ssp_dir):
  # Path for levitus files.
  levitus_dir = './ocean_profiles/ssp_09/'

  # Take absolute drop depth estimate.
  z_sta = abs(z_sta)

  # This infers whether sation depth is in [m] or [km]. Adjusts to [km].
  if z_sta > 6:
    z_sta = z_sta/1000

  # Name current station sound-speed profile (ssp) if supplied by user.
  ssp_user_fname = ssp_dir + 'SSP_' + stn + '.txt'

  # Name if previously calculated by OBSrange.
  ssp_stock_fname = ssp_dir + 'SSP_' + stn + '.txt'

  # Check if SSP for current station exists (whether it was provided by the
  # user or computed by OBSrange) create if otherwise.
  if not os.path.exists(ssp_stock_fname):

    if not os.path.exists(ssp_user_fname):
      # Compute month of the deployment from t.
      month = time.gmtime(t).tm_mon 
  
      # Compute approptiate sound-speed-depth profile from Levitus database.
      ssp, z = get.lev_based_ssp(lat, lon, month + 1, levitus_dir, 'SSP_' + stn + '.txt', ssp_dir)
    
    else:
      ssp, z = load_ssp(ssp_user_fname)

  else:
    ssp, z = load_ssp(ssp_stock_fname)
  
  cor_fpath = ssp_dir + 'cor_rayline_' + stn +'.pkl'
  if not os.path.exists(cor_fpath):
    # SSP, and Z info come as Python lists. Convert to arrays and desired units.
    ssp = np.array([float(val) for val in ssp])
    z = np.array([float(val) for val in z])
    v_profile = np.array([z, ssp])/1000
    
    ############################### Ray-Shooting ###############################
    
    # Basic ray-shooting parameters.
    ps = np.arange(0.01, 1.5 + 0.01, 0.01)       # Ray parameters
    zs = np.arange(-200, 200 + 20, 20)/1000      # Depths 200 meters above and 
    zs += z_sta                                  #    below the starting point 
    zmax = max(zs)                               # Max depth
    Np = len(ps)                                 # Number of ray parameters
    Nz = len(zs)                                 # Number of depths
    Dr_ray = np.zeros(shape=(Np, Nz))            # Refracting ray lengths
    Dx = np.zeros(shape=(Np, Nz))                # Refracting ray horizontal dists
    t_ray = np.zeros(shape=(Np, Nz))             # Refracting ray travel-times
    Dr_lin = np.zeros(shape=(Np, Nz))            # Linear ray lengths
    t_lin = np.zeros(shape=(Np, Nz))             # Linear ray travel-times
  
    # Fill distance and time arrays with nans by default.
    Dr_ray[:] = np.nan
    Dx[:] = np.nan
    t_ray[:] = np.nan
    Dr_lin[:] = np.nan
    t_lin[:] = np.nan
  
    # Loop through ray parameters.
    print('\n Ray Tracing... ')
    for i, ray in enumerate(ps):
      # Compute distances and times. (unless rays turn complex)
      rx, rz, Dr, rt, rv = shootrays.shootrays(ray, v_profile, zmax)
  
      # Situation for when rays go complex.
      if (np.array([rx, rz, Dr, rt, rv]) == 0).all():
        continue
      
      # Get intersection of this ray at each depth.
      for j, depth in enumerate(zs):
        f_Dr = interp1d(rz, Dr)
        Dr_ray[i,j] = f_Dr(depth)
      
        f_Dx = interp1d(rz, rx)
        Dx[i,j] = f_Dx(depth)
        
        f_t = interp1d(rz, rt)
        t_ray[i,j] = f_t(depth)
        
        Dr_lin[i,j] = np.sqrt(Dx[i,j]**2 + depth**2)
        t_lin[i,j] = Dr_lin[i,j] / hmean(rv[rz <= depth]) 
    
  
    # Difference in ray length [m].
    dDr_m = 1000 * (Dr_ray - Dr_lin)
    
    # Difference in ray two-way-travel-time [ms].
    dt_ray_ms = 2 * 1000 * (t_ray - t_lin) 
    # to be clear: this number is POSITIVE if ray is slower than straight line
    # approx. We should therefore SUBTRACT this from the observed travel times
    # to correct them from bent rays to straight lines. Having turned data from
    # rays to lines the code, which assumes lines, can then invert for position
    
    # Add a vertical ray.
    Dx = np.vstack([np.zeros(Nz), Dx])
    dDr_m = np.vstack([np.zeros(Nz), dDr_m])
    dt_ray_ms = np.vstack([np.zeros(Nz), dt_ray_ms])
  
    # Interpolate tt error onto regular grid spaced every 5m up to 4km offset.
    minimum = np.min(np.hstack([4, np.nanmax(Dx, axis=0)]))
    Dx_grid = np.arange(0, minimum + 0.005, 0.005) * 1000
    
    # Make a grid of differential tts, columns are depths, rows are offsets.
    dT_grid = np.zeros(shape=(len(Dx_grid), Nz))
    dT_grid[:] = np.nan
  
    for i, depth in enumerate(zs):
      indx = ~np.isnan(dt_ray_ms[:,i])
      f_dT = interp1d(Dx[indx == True, i]*1000, dt_ray_ms[indx == True, i])
      dT_grid[:,i] = f_dT(Dx_grid)
  
    # Put results in a dictionary and write to disk
    dT_ray_v_line_grid = {'Dx_grid_m': Dx_grid,
                          'dz_grid_km': zs,
                          'dT_grid_ms': dT_grid}
  
    f = open(cor_fpath, 'wb')
    pickle.dump(dT_ray_v_line_grid, f)
    f.close()
  
    return dT_ray_v_line_grid

  else:
    dT_ray_v_line_grid = pickle.load(open(cor_fpath, 'rb'))
    return dT_ray_v_line_grid
