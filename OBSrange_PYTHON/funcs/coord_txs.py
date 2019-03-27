'''
FUNCTION SET coord_txs.py

A set of functions to transform geographic coordinates to xy-coords and vice
versa. Uses the pymap3d module to perform the transformations. Both functions
are based off the WGS84 reference ellipsoid.

Stephen M. 4/23/18
'''
# Import modules and functions.
import pymap3d as pm

def latlon2xy(lat0, lon0, lat, lon):
  e, n, u = pm.geodetic2enu(lat=lat, lon=lon, h=0, lat0=lat0, lon0=lon0, h0=0)
  return e, n

def xy2latlon(x, y, lat0, lon0):
  lat, lon, elv = pm.enu2geodetic(e=x, n=y, u=0, lat0=lat0, lon0=lon0, h0=0)
  return lat, lon