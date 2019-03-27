'''
FUNCTION SET bootstrap.py

The first function gets indices for balanced resampling of data. It takes in a
numpy array (data) as well as the number of bootstrap iterations (N). The output
is an array of the same set of indices randomly permuted N times.

Second function re-indexes supplied data arrays according to previously permuted
indices. The overall reult of the first two functions is balanced random 
resampling, where each data point occurs N times in the bootstraping scheme.

The third function inverts the results of the previous two functions (i.e. it 
re-orders an array by its proper index order).

The fifth function is a simple helper function to get the correct ray bending
corrections for the current station. 

The fifth function performs the inversion.

Josh R. & Zach E. & Stephen M. 4/23/18
'''
#Import modules and functions
import numpy as np
import numpy.linalg as LA
from funcs import results, calc, coord_txs

def get_bs_indxs(data, N):
  # Create an array of indices for the number of data points.
  Ndata = len(data)
  idxs = np.arange(0, Ndata)

  # Create N-1 copies of indices. Nth copy added later.
  idxs = np.tile(idxs, N-1)

  # Randomly permute indices and reshape in a 2D array.
  rand_idxs = np.random.permutation(idxs)
  rand_idxs = np.reshape(rand_idxs, (Ndata, N-1))

  # Add a column of unpermuted indices at the beginning of the array.
  rand_idxs = np.c_[np.arange(Ndata), rand_idxs]

  # Doing it this way ensures each index occurs exactly N times.
  
  return rand_idxs

def sampling(x, y, z, v, t, N):
  # Get randomly sampled indices.
  idxs = get_bs_indxs(x, N)

  # Index parameter data by bootstrap indices then return.
  xmatbs = x[idxs]
  ymatbs = y[idxs]
  zmatbs = z[idxs]
  vmatbs = v[idxs]
  tmatbs = t[idxs]

  # Small adjustment to axes of vmatbs, given vmatbs has a different shape.
  vmatbs = np.swapaxes(vmatbs, 1, -1,)

  return xmatbs, ymatbs, zmatbs, vmatbs, tmatbs, idxs

def unscramble(data, idxs):
  # Initialze a new array for the unscrambled data.
  unscrambled_data = np.ndarray(shape=data.shape)
  Ninds = idxs.shape[0]
  
  # Re-index by proper index order.
  for index in np.arange(Ninds):
    unscrambled_data[index] = data[index == idxs]
  
  return unscrambled_data

def ray_correct_get_time(dt_rvl, r):
  # Loads in ray correction grids.
  x_grid = dt_rvl['Dx_grid_m']
  z_grid = dt_rvl['dz_grid_km']
  t_grid = dt_rvl['dT_grid_ms']

  # Reshape for later operations.
  x_grid = x_grid.reshape(len(x_grid), 1)
  z_grid = z_grid.reshape(len(z_grid), 1)

  # Break apart ship position into xy and z.
  dxy_ship = np.sqrt(np.sum(r[0:2,:]**2, axis=0)) 
  dz_ship = r[2,:]

  # Find closest indices to get appropriate ray correction times.
  indx = np.argmin(abs(x_grid - dxy_ship), axis=0)
  indz = np.argmin(abs(z_grid - dz_ship/1000), axis=0)

  return t_grid[[indx, indz]]/1000

def inv(X, Y, Z, V, T, R, parameters, m0_strt, coords, M, dt_rvl):
  # Grab required parameters.
  vpw0     = parameters[0]    
  dvp0     = parameters[1]
  tat      = parameters[2]       
  N_bs     = parameters[3]    
  E_thresh = parameters[4]
  dampx    = parameters[6]   
  dampy    = parameters[7]   
  dampz    = parameters[8]   
  dampdvp  = parameters[9]
  eps      = parameters[10] 
  twtcorr  = parameters[15]
  raycorr  = parameters[16]   
  Nobs     = X.shape[0]
  x0       = coords[0][0]
  y0       = coords[0][1]
  z0       = coords[0][2]
  lat0     = coords[0][3]
  lon0     = coords[0][4]
  xs       = coords[1][0]
  ys       = coords[1][1]
  zs       = coords[1][2]

  # Loop over a random sample of model parameters.
  for i in range(N_bs):
    # Create a dictionary for the output models from each bootstrap iteration. 
    R.models[str(i)] = {}
    
    # Grab a set of randomly sampled model parameters for current bootstrap it. 
    xbs = X[:,i]
    ybs = Y[:,i]
    zbs = Z[:,i]
    vbs = V[:,:,i]
    twtbs = T[:,i]

    # Initialize.
    m1 = m0_strt
    dE = 1000
    j = 0 
    # Iterate models until RMS stabilizes.
    while dE > E_thresh:
      # Increase j
      j +=1

      # Set current model parameters
      m0 = m1
      x = m0[0]
      y = m0[1]
      z = m0[2]
      dvp = m0[3]

      # Apply correction to two-way travel-times due to ship velocity.
      ctd, cns, vr = calc.tt_corr(x, y, z, xbs, ybs, zbs, vbs, vpw0, dvp, twtbs)
      if twtcorr:
        twts = ctd # ctd = corrected, cns = corrections
      else:
        twts = twtbs

      # Apply ray bending correction.
      if raycorr and dt_rvl:
        r = np.array([xs - x0, ys - y0, zs - z0])
        dT_ray_v_line = ray_correct_get_time(dt_rvl, r)
        twts = twts - dT_ray_v_line

      # Build the G matrix.
      G = calc.G(x, y, z, dvp, xbs, ybs, zbs, vpw0, Nobs, M)
      
      # Set up norm damping for each parameter.
      H = np.eye(M) * np.diag([dampx, dampy, dampz, dampdvp])
      h = np.zeros(M)
      
      # Calculate predicted travel-times for this iteration and residuals.
      twt_pre = calc.twt(x, y, z, xbs, ybs, zbs, vpw0, dvp, tat)
      dtwt = twts - twt_pre

      # Calculate RMS error.
      E = np.sqrt( np.sum(dtwt**2) / Nobs)
      
      # Least squares solution.
      J = np.eye(M) * np.sqrt(eps)  
      F = np.concatenate([G, H, J])
      f = np.concatenate([dtwt, h, np.zeros(M)])
      Finv = LA.solve(np.dot((F.T), F), F.T)

      # Model update. 
      m1 = m0 + np.dot(Finv, f)

      # Invert G.
      Ginv = LA.solve(np.dot(G.T, G) + np.dot(H.T, H) + np.dot(J.T, J), G.T)
      
      # Calculate the effective degrees of freedom (for F-test)
      v = len(f) - np.trace(np.dot(F, Finv))
      
      # Exclude 0=0 constraint equations from v_eff
      v = v - np.sum(np.sum(H,axis=1) == 0) - np.sum(np.sum(J, axis=1) == 0)
      
      # Record output of current iteration in m0, E, and dtwt.
      R.models[str(i)][str(j)] = {'m': m0, 'E': E, 'dtwt': dtwt}
      
      # Update inversion RMS.
      if j > 1:
        dE = R.models[str(i)][str(j-1)]['E'] - R.models[str(i)][str(j)]['E']
    
    # Model resolution and covariance
    dat_res = np.matmul(G, Ginv)
    mod_res = np.matmul(Ginv, G)
    mod_cov = np.matmul(Ginv, Ginv.T)

    # Update results object with stabilized results of current bootstrap it.
    R.update(i, m0, vpw0, E, v, dtwt, twts, cns, vr, x0, y0, z0, xs, ys, \
             dat_res, mod_res, mod_cov, tat)
    
    # Convert stabilized result coords of current iteration back to lat lon.
    R.lats[i], R.lons[i] = coord_txs.xy2latlon(R.xs[i], R.ys[i], lat0, lon0)

  # Return filled results.
  return R