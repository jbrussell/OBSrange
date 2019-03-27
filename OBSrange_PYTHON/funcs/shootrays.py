'''
FUNCTION shootrays.py 

1D ray tracing function. All units in km or km/s

Z. Eilon 01/2019 (modified from code by Brandon Schmandt)
Translated to Python and further modified by S. Mosher 01/2019
'''
# Import modules and functions
import numpy as np
from numpy.lib.scimath import sqrt
from scipy.interpolate import interp1d

from IPython.core.debugger import Tracer

# Given two arrays, this function finds the closest value in the smaller array
# to each value in the larger array.
def find_ilay(r, R):
  N = len(r)
  ilay = np.zeros(N, dtype=int)

  for i, depth in enumerate(r):
    idx = np.where(R - depth > 0)[0][0]
    ilay[i] = idx

  return ilay

# Linear interpolation function.
def linterp(X, Y, XI):
  YI = np.zeros(len(XI))
  ilay = find_ilay(XI, X)
  
  # Linearly interpolate.
  YI = (XI - X[ilay-1])*(Y[ilay] - Y[ilay-1])/(X[ilay] - X[ilay-1]) + Y[ilay-1]

  # Sort out coincident elements.
  olap = np.intersect1d(X, XI)
  for i, val in enumerate(olap):
    YI[XI == olap[i]] = np.mean(Y[X == olap[i]])

  return YI

def shootrays(p, v_profile, zmax, dr=0.001, vdz=0.001):
  # Anonymous functions
  eta = lambda u, p: sqrt(u**2 - p**2)
  gradv_dist = lambda b, u1, u2, p: (eta(u1,p)/(b*u1*p)) - (eta(u2,p)/(b*u2*p))
  constv_dist = lambda u, dz, p: (p * dz) / eta(u, p)

  zz1 = np.arange(0, zmax, vdz)
  zz2 = np.arange(vdz, zmax + vdz, vdz)
  
  v1  = linterp(v_profile[0], v_profile[1], zz1)
  v2  = linterp(v_profile[0], v_profile[1], zz2)

  # Assume flat Earth (no radius terms)
  u1 = 1 / v1
  u2 = 1 / v2 
  
  dv = v2 - v1
  
  dz = zz2 - zz1
  b = dv / dz
  const_indx = (b == 0)

  X = np.zeros(len(v1), dtype=complex)
  
  X[const_indx == True] = \
                constv_dist( u1[const_indx == True], dz[const_indx == True], p)
  
  X[const_indx != True] = \
  gradv_dist(b[const_indx!=True], u1[const_indx!=True], u2[const_indx!=True], p)
  

  # When rays go complex skip out.
  if np.imag(X).any():
    rx = 0
    rz = 0
    Dr = 0
    rt = 0
    rv = 0

    return rx, rz, Dr, rt, rv

  else:
    X = np.real(X)
    X = np.hstack([0, X])
    Xd = np.cumsum(X)
    rayz = np.hstack([zz1[0], zz2]) 
    rayz = rayz[0:len(X)]
    dz = dz[0:len(X) - 1]
   
    # Calculate Dr
    Dr1 = np.sqrt(X[1:]**2 + dz**2)
    Dr1 = np.hstack([0, Dr1])
    Dr2 = np.cumsum(Dr1)
    Dr = np.unique( np.hstack([np.arange(0, max(Dr2), dr), max(Dr2)]))
    fx = interp1d(Dr2, Xd)
    fz = interp1d(Dr2, rayz)
    rx = fx(Dr)
    rz = fz(Dr)
  
    # Calculate travel-time (JBR)
    dDr = Dr[1:len(Dr)] - Dr[0:len(Dr)-1]
    dDr = np.hstack([0, dDr])
    v = np.hstack([v1[0], v2])
    fv = interp1d(Dr2, v)
    rv = fv(Dr)
    dt = dDr / rv
    rt = np.cumsum(dt)
  
    return rx, rz, Dr, rt, rv