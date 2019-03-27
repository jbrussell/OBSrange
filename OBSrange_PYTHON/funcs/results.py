'''
CLASS results.py

A class to initialize various vectors and matrices, as well as an output
dictionary, all for holding results of the bootstrap inversion.

Stephen M. 04/23/18
'''
# Import modules and functions
import numpy as np

class results:
  '''
  Initialization
  '''
  def __init__(self, M=1, N=1, P=1):  
    # Vectors
    self.lats    = np.zeros(M)   # sensor latitudes
    self.lons    = np.zeros(M)   # sensor longitudes
    self.tats    = np.zeros(M)   # turn-around times
    self.xs      = np.zeros(M)   # sensor x positions
    self.ys      = np.zeros(M)   # sensor y positions
    self.zs      = np.zeros(M)   # sensor z positions
    self.dzs     = np.zeros(M)   # changes in depth
    self.vpws    = np.zeros(M)   # sound speeds for water
    self.dvps    = np.zeros(M)   # perturbations to water sound speeds
    self.drifts  = np.zeros(M)   # sensor drift distance
    self.dxdrfts = np.zeros(M)   # changes to drift in x-direction
    self.dydrfts = np.zeros(M)   # changes to drift in y-direction
    self.azs     = np.zeros(M)   # sensor drift azimuth
    self.E_rms   = np.zeros(M)   # RMS error terms
    self.v_effs  = np.zeros(M)   # effective degrees of freedom for F-test
    
    # Matrices
    self.dists   = np.ndarray(shape=(N, M))    # ship distances from OBS
    self.az_locs = np.ndarray(shape=(N, M))    # ship azimuths from OBS
    self.twts    = np.ndarray(shape=(N, M))    # travel-times
    self.corrs   = np.ndarray(shape=(N, M))    # travel-time corrections
    self.dtwts   = np.ndarray(shape=(N, M))    # travel-time residuals
    self.vrs     = np.ndarray(shape=(N, M))    # radial velocities of ship
    self.dat_res = np.ndarray(shape=(N, N, M)) # data resolution matrices
    self.mod_res = np.ndarray(shape=(P, P, M)) # model resolution matrices
    self.mod_cov = np.ndarray(shape=(P, P, M)) # model covariance matrices
    
    # Output dictionary
    self.models  = {}                          # holds each model updates
  
  '''
  A function to update the results object.
  '''
  def update(self, i, m0, vpw, E, v, dt, tts, corr, vr, x0, y0, z0, xs, ys, \
             dat_res, mod_res, mod_cov, tat):
    self.xs[i]      = m0[0]
    self.ys[i]      = m0[1]
    self.zs[i]      = m0[2]
    self.dvps[i]    = m0[3]
    self.vpws[i]    = vpw + self.dvps[i]
    self.tats[i]    = tat
    self.E_rms[i]   = E
    self.v_effs[i]  = v
    self.dtwts[:,i] = dt
    self.twts[:,i]  = tts
    self.corrs[:,i] = corr
    self.vrs[:,i]   = vr
    self.dat_res[:,:,i] = dat_res
    self.mod_res[:,:,i] = mod_res
    self.mod_cov[:,:,i] = mod_cov
    
    # Calculate OBS drift distance and azimuth.
    self.dxdrfts[i] = m0[0] - x0
    self.dydrfts[i] = m0[1] - y0
    self.dzs[i] = m0[2] - z0
    self.drifts[i] = np.sqrt(self.dxdrfts[i]**2 + self.dydrfts[i]**2)
    az = -(np.degrees( np.arctan2(self.dydrfts[i], self.dxdrfts[i]) )) + 90
    if az < 0:
      az += 360
    self.azs[i] = az
    
    # Calculate ship distance and azimuths from OBS.
    dxs = xs - m0[0]
    dys = ys - m0[1]
    self.dists[:,i] = np.sqrt(dxs**2 + dys**2)
    az_locs = - (np.degrees(np.arctan2(dys , dxs))) + 90
    az_locs[az_locs < 0 ] += 360
    self.az_locs[:,i] = az_locs