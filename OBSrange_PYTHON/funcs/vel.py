'''
FUNCTION SET velocity.py 

A function to calculate the velocity vectors of the ship at survey points with
positions (xs, ys, zs) recorded at times ts. 

Another function smooths the ship's velocity vectors.

Zach E. & Josh R. & Stephen M. 4/23/18
'''
# Import modules and functions
import numpy as np
import scipy.linalg as LA

def vector(xs, ys, zs, ts):
  # Check for stationary positions.
  dts = np.diff(ts)
  is_dts_0 = dts == 0

  # Avoid zero division. Replace later.
  dts[is_dts_0 == True] = 1
  
  # Calculate half velocity.
  v_half = np.array([np.diff(xs)/dts, np.diff(ys)/dts, np.diff(zs)/dts]).T

  # Replace stationary points with zeros.
  v_half[is_dts_0] = [0,0,0]

  # Reconstruct a point since np.diff(X) throws away a point.
  v = (v_half[0:-1] + v_half[1:])/2
  
  # Calculate full velocity with reconstructed point. Then return.
  v = np.insert(v, 0, v_half[0], axis=0)
  v = np.insert(v, len(v), v_half[-1], axis=0)

  return v

def smooth(vel, npts):
  # Decompose velocity vectors.
  vx = vel[:,0]
  vy = vel[:,1]
  vz = vel[:,2]

  # Make sure npts is odd.
  if npts % 2 != 1:
    npts = npts + 1

  dy = npts//2
  L = vel.shape[0]

  column = np.zeros(L-npts)
  column = np.insert(column, 0, 1)
  row = np.ones(npts)
  row = np.concatenate([row, np.zeros(L-npts)])

  # mid full moving average section.
  C = LA.toeplitz(c=column, r=row)
  
  # taper section (wish I could think of a way to avoid the loop!).
  D = np.zeros([dy,L])
  f = 2 * np.arange(1, dy + 1) - 1
  
  for i in range(len(f)):
    D[i, 0:f[i]] = 1

  # complete moving average filter (not normalised).
  G = np.concatenate([D, C, np.rot90(D, 2)])

  # Multiply and normalize.
  vx = np.dot(G,vx)/np.sum(G, axis=1)
  vy = np.dot(G,vy)/np.sum(G, axis=1)
  vz = np.dot(G,vz)/np.sum(G, axis=1)
  
  return np.vstack([vx, vy, vz]).T