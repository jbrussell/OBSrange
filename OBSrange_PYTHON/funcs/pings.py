'''
FUNCTION SET pings.py

A set functions to load ping data from the .txt file created during the sensor
survey, format it, and perform quality control. 

First function format() formats data in the .txt file.

Output of the second function load() is a dictionary. Fields are:

  1) Drop latitude                             -- in decimal degrees
  2) Drop longitude                            -- in decimal degrees
  3) Drop depth                                -- meters
  4) Latitude of each survey point             -- converted to decimal degress
  5) Longitude of each survey point            -- converted to decimal degrees
  6) Time each survey point was received       -- converted to seconds
  7) Two-way travel-time of each survey point  -- in msecs
  8) Station id                                -- string

Outputs of the third function qc() are two dictionaries into which the good and
bad data points are segregated. Bad points are considered as those for which the
two-way travel time exceeds a time threshold (in msec). The threshold is based 
on a reasonable amount of time needed to traverse the Euclidean distance between
the survey point and the sensor drop coordinates based on an average sound speed
in water (vpw).                                 

Josh R. & Zach E. & Stephen M. 4/23/18
'''
# Import modules and functions
import numpy as np
from funcs import coord_txs, calc

def format(lines):
  for i, line in enumerate(lines):
    # Remove any leading whitespace from lines
    lines[i] = line.strip()
  return lines

def load(txt_file):
  # Create a python list of the lines in the survey .txt file.
  lines = list( open(txt_file, mode='r') )
  
  # Format lines.
  lines = format(lines)

  # Populate single-valued fields.
  sta = lines[2].split(' ')[-1]
  lat_drop = float( lines[4].split(' ')[-1] )
  lon_drop = float( lines[5].split(' ')[-1] )
  z_drop = float( lines[6].split(' ')[-1] )
  
  # Initialize python lists for lats, lons, ts, and twts.
  lats = []
  lons = []
  ts = []
  twts = []
  
  # Loop through lines of the .txt file to fill the above lists, skip bad lines.
  for line in lines[10:]:
    if line.split(' ')[0] == 'Event': # Most bad lines start 'Event skipped'
      continue
    if line.split(' ')[0][0] == '*': # Some bad lines start with '*'
      continue
    else:
      lats.append( line.split('Lat:')[1].split('Lon:')[0].strip() )
      lons.append( line.split('Lon:')[1].split('Alt:')[0].strip() ) 
      ts.append( line.split('Time(UTC): ')[-1].strip() )
      twts.append( line.split('msec')[0].strip() )

  # Format lats, lons, ts, and twts to desired units.
  for i, (lat, lon, t, twt) in enumerate( zip(lats, lons, ts, twts) ):
    # Account for N-S / E-W when converting coordinate to decimal degrees.
    if lat.split(' ')[2] == 'S':
      N = -1
    else:
      N = 1
    if lon.split(' ')[2] == 'W':
      E = -1
    else:
      E = 1
    lats[i] = (float(lat.split(' ')[0]) + float(lat.split(' ')[1])/60) * N
    lons[i] = (float(lon.split(' ')[0]) + float(lon.split(' ')[1])/60) * E
    
    jday = float(t.split(':')[1])
    hour = float(t.split(':')[2])
    mint = float(t.split(':')[3])
    secd = float(t.split(':')[4]) 
    ts[i] = jday * 24 * 60 * 60 + hour * 60 * 60 + mint * 60 + secd
    
    twts[i] =  float(twt) * 1e-3

  # Package everything into a dictionary and return. Convert lists to np arrays.
  data = {
          'lat_drop' : lat_drop,
          'lon_drop' : lon_drop,
          'z_drop' : z_drop*-1,
          'lats' : np.array(lats),
          'lons' : np.array(lons),
          'ts' : np.array(ts),
          'twts' : np.array(twts),
          'sta' : sta
          }
          
  return data

def qc(data, vpw, thresh):
  # Grab required variables to calculate time to traverse the Euclidean dist.
  lat0 = data['lat_drop']
  lon0 = data['lon_drop']
  z0   = data['z_drop']
  lats = data['lats']
  lons = data['lons']
  twts = data['twts']
  
  # Convert geographic coordinates to Cartesian.
  xs, ys = coord_txs.latlon2xy(lat0, lon0, lats, lons)
 
  # Calculate theoretical two-way travel-times and residuals.
  theor_twts = calc.twt(x0=0,y0=0,z0=z0,xs=xs,ys=ys,zs=0,vpw=vpw,dvp=0,tat=0)
  dtwts = twts - theor_twts

  # Find which survey points exceed the time threshold (units msecs).
  is_bad = abs(dtwts) * 1000 > thresh
  
  # Segregate bad and good data points into separate dictionaries then return.
  data_good = {
              'lat_drop' : lat0,
              'lon_drop' : lon0,
              'z_drop' : z0,
              'lats' : data['lats'][np.where(is_bad==False)],
              'lons' : data['lons'][np.where(is_bad==False)],
              'ts' : data['ts'][np.where(is_bad==False)],
              'twts' : data['twts'][np.where(is_bad==False)],
              'sta' : data['sta']
               }

  data_bad = {
             'lat_drop' : lat0,
             'lon_drop' : lon0,
             'z_drop' : z0,
             'lats' : data['lats'][np.where(is_bad==True)],
             'lons' : data['lons'][np.where(is_bad==True)],
             'ts' : data['ts'][np.where(is_bad==True)] ,
             'twts' : data['twts'][np.where(is_bad==True)],
             'sta' : data['sta']
             }

  return data_good, data_bad