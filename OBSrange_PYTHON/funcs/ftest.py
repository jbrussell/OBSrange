'''
FUNCTION ftest.py

Function to perform a formal F-test on two sets of residuals of fits to data
(res1, res2) that have respectively been fitted using parms1, parms2 (where 
these are integers equal to the number of parameters used in each fit).

IMPORTANTLY: 
  This will test whether the second fit (yielding res2) is statistically
  superior to the first fit, but using more parms.
    I.e.:          sum(abs(res2)) < sum(abs(res1)) 
    and            parms2 > parms1
  The degrees of freedom for each fit are therefore:
    1) N - parms1
    2) N - parms2
  where N is the number of data, or len(res*). The residuals are just equal to
  dobs - dpred, so we square and sum these to get the chi^2 values (Ea)

Z. Eilon

J. Russell:
  This version takes the degrees of freedom as input (instead of model 
  parameters) for the case where v = N - M is not accurate.
  P > 1 : Model 1 actually fits the data better than Model 2 (model 1 smaller
          chi^2)
  P = 1 : Model 1 fits the data same as model 2 (same chi^2)
  P < 1 : Model 2 fits the data better than model 1 (model 2 smaller chi^2)
  P < 0.05 : Model 2 fits the data better model 1 with 95% confidence

S. Mosher: 
  Translation from MATLAB to Python
'''
# Import modules and functions
import numpy as np
import scipy.linalg as LA
from scipy.stats import f
from funcs import coord_txs, calc

def dof(res1, v1, res2, v2):

  # Calculate chi**2 sums.
  Ea_1 = np.sum(res1**2, axis=0)
  Ea_2 = np.sum(res2**2, axis=0)
  
  Fobs = (Ea_1/v1)/(Ea_2/v2)
  P = 1 - ( f.cdf(Fobs, v1, v2) - f.cdf(1/Fobs, v1, v2) )
  
  return P

def test(R, coords, lat0, lon0, vpw):
  # Grab coords.
  xs = coords[1][0]
  ys = coords[1][1]
  zs = coords[1][2]

  # Grab intermediate variables from results.
  xM = np.mean(R.xs)
  yM = np.mean(R.ys)
  zM = np.mean(R.zs)
  tatM = np.mean(R.tats)
  dvpM = np.mean(R.dvps)
  v_effM = np.mean(R.v_effs)
  
  # Set grid size and grids.
  ngridpts = 40
  D = max(np.std(R.xs), np.std(R.ys), np.std(R.zs)) * 4
  Dx = D 
  Dy = D 
  Dz = D 
  dx = 2*Dx/ngridpts 
  dy = 2*Dy/ngridpts 
  dz = 2*Dz/ngridpts
  xg = np.linspace(xM - Dx, xM + Dx, ngridpts)
  yg = np.linspace(yM - Dy, yM + Dy, ngridpts)
  zg = np.linspace(zM - Dz, zM + Dz, ngridpts)
  assert len(xg) == len(yg) == len(zg)

  lat_grid, lon_grid = coord_txs.xy2latlon(xg, yg, lat0, lon0)
  Nx = len(xg)
  Ny = len(yg)
  Nz = len(zg)

  # Bootstrap residuals.
  twt_pre = calc.twt(xM, yM, zM, xs, ys, zs, vpw, dvpM, tatM)
  resid_bs = R.twts - twt_pre
  
  # Determine the eigenvectors for z, vpw, and tat.
  X = np.vstack([R.zs, R.vpws, R.tats]).T
  V = LA.eigh(np.dot(X.T, X))[1]
  
  eigvec1 = V[:,0] # Closest to TAT axis
  eigvec2 = V[:,1] # Closest to V_w axis
  eigvec3 = V[:,2] # Closest to z axis
  eig3_z = eigvec3[0]
  eig3_vw = eigvec3[1]
  eig3_tat = eigvec3[2]

  # Create spatial mesh and tensors (faster implementation than nested loops).
  Xg, Yg, Zg = np.meshgrid(xg, yg, zg)
  
  xtensor = np.repeat(Xg[np.newaxis,:,:,:], len(xs), axis=0)
  ytensor = np.repeat(Yg[np.newaxis,:,:,:], len(ys), axis=0)
  ztensor = np.repeat(Zg[np.newaxis,:,:,:], len(zs), axis=0)
  
  xtensor_minus_xs = np.ndarray(xtensor.shape)
  ytensor_minus_ys = np.ndarray(ytensor.shape)
  ztensor_minus_zs = np.ndarray(ztensor.shape)
  
  dz = np.ndarray(shape=xtensor.shape)
  dvw = np.ndarray(shape=xtensor.shape)
  dtat = np.ndarray(shape=xtensor.shape)
  
  # Intermediate calculations for grid search residuals and eigen-scaling.
  for i, (x,y,z) in enumerate(zip(xs,ys,zs,)):
    xtensor_minus_xs[i,:,:,:] = Xg - xs[i]
    ytensor_minus_ys[i,:,:,:] = Yg - ys[i]
    ztensor_minus_zs[i,:,:,:] = Zg - zs[i]
    
    dz[i,:,:,:] = Zg - zM

  # Eigenvector scaling.
  dvw = ((eig3_vw/eig3_z) * dz) + dvpM 
  dtat = (eig3_tat/eig3_z) * dz 

  # Compute travel times (grid search residuals).
  tt = 2*np.sqrt(xtensor_minus_xs**2 + ytensor_minus_ys**2 + ztensor_minus_zs**2)
  twt_pre_gs = (tt/(vpw+dvw)) + (dtat+tatM)
  
  resid_gs = np.ndarray(shape=xtensor.shape)
  for i, r in enumerate(R.twts):
    resid_gs[i,:,:,:] = r - twt_pre_gs[i,:,:,:]

  # P-statistic.
  P = np.ndarray(shape=(Nx, Ny, Nz))
  E_gs = np.ndarray(shape=(Nx, Ny, Nz))

  P = dof(resid_gs, v_effM, resid_bs, v_effM) 
  E = np.sqrt( np.sum(resid_gs**2) / len(resid_gs))

  xmax, ymax, zmax = np.where( P == np.amax(P) )

  return xg, yg, zg, Xg, Yg, Zg, P, xmax[0], ymax[0], zmax[0], E