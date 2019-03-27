'''
FUNCTION SET get.py

A pair of routines to extract SSPs, temperature, salinity, or buoyancy profiles 
along a given point or points along a path (lat, lon). Uses interp_global_SSP(),
a function to load the variables and interpolate onto the desired location.

Gets the World Ocean Atlas files from directory 'lev_db_dir'.

These files are preprocessed to save disk space and have shortened profiles from
the original data filled in by nearest neighbor data. All profiles extend to
5500 m depth, and there is world wide coverage (yes, including continents)

  See:  http://staff.washington.edu/dushaw/WOA/
  Data originally from:  http://www.nodc.noaa.gov/OC5/indprod.html

Data at points where the original WOA has data should be unchanged from the 
WOA. The annual mean World Ocean Atlas is used.
 
  type_SSP = 0 indicates annual Levitus,
  type_SSP = 1:12 indicates Levitus for month 1:12

UNITS: sound speed in [m/s], depth in positive [meters].  

Returns P, a matrix of profiles: 33 * (No. of points) and z, the standard 33
World Ocean Atlas depths.

Interpolation from World Ocean Atlas grid points to desired points by cubic
spline interpolation horizontally.

Writes out file "stn_ssp_fname".

Modified from original for simpler application (only sound speed) by: 
Zach Eilon and Stephen Mosher 01/25/19
'''
# Import modules and functions
import os
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import interp2d

# Interpolate global SSPs from Levitus database to get local station SSP.
def interp_global_SSP(type_str, lat, lon, lev_db_dir):
  # Levitus file name. 
  lev_fname = lev_db_dir + 'lev_' + type_str + '.mat'
  
  # Levitus coordinate file.
  lev_coords = lev_db_dir + 'lev_latlonZ.mat' 

  # Load in global sound speed profiles and their coordinates.
  global_ssps = loadmat(lev_fname)['c']
  coords = loadmat(lev_coords)

  # Scipy interpolation functions skip meshgrid generation, instead they simply
  # require the coordinate paramaters in monotonic order (as opposed to MATLAB).
  lats = coords['lat'].reshape((180, 360)).T * 0.1
  lons = coords['lon'].reshape((180, 360)).T * 0.1
  lats = lats[0,:]
  lons = lons[:,0]
  
  # Declare an array for the local SSP.
  local_ssp = np.ndarray(shape=(global_ssps.shape[1]))
  
  # Loop through depth slices of global SSPs and interpolate at each depth.
  for i in range(len(local_ssp)):
    depth_val = global_ssps[:,i].reshape((180, 360)).T
    f = interp2d(lats, lons, depth_val, kind='cubic')
    local_ssp[i] = f(lon, lat) * 0.01 + 1000

  return local_ssp

# Main function.
def lev_based_ssp(lat, lon, type_SSP, lev_db_dir, stn_ssp_fname, ssp_dir):
  # Longitude convention conversion for Levitus database.
  if lon < 0:
    lon = lon + 360

  # Standard World Ocean Atlans depths for SSPs.
  z = [0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500,
       600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 
       2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500]
  
  # Keys for ssp_09 files.
  type_str = ['ann', 'jan', 'feb', 'mar', 'apr', 'may', 'jun', \
              'jul', 'aug', 'sep', 'oct', 'nov', 'dec']

  # Current file type.
  type_str = type_str[type_SSP]

  # Interpolate global SSP file from Levitus database to get local SSP.
  ssp = interp_global_SSP(type_str, lon, lat, lev_db_dir)

  # Check if output folder for generated SSPs exists.
  path = ssp_dir
  if not os.path.exists(path):
    os.mkdir(path)

  # Save local SSP to disk.
  print('\n Sound speed profile saved to ' + path + stn_ssp_fname)
  txt_fle = open(path + stn_ssp_fname, 'w')
  txt_fle.write('depth(m) ssp(m/s)\n')
  for depth, velocity in zip(z, ssp):
    txt_fle.write('{:>5}    {:>4.2f}\n'.format(depth, velocity))
  txt_fle.close()

  # Return
  return ssp, z
