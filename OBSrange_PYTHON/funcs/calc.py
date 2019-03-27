'''
FUNCTION SET calc.py

Functions for several intermediate calculations. Artwork courtesy of Zach.

Josh R. & Zach E. & Stephen M. 4/23/18
'''
# Import modules and functions
import numpy as np

'''
A function to calculate two-way travel-times between a series of points
(xs, ys, zs) and a single reference point (x0, y0, z0) given an average sound 
speed for water (vpw), a sound speed perturbation (dvp), and a constant turn-
around-time for the sensor.
'''
def twt(x0, y0, z0, xs, ys, zs, vpw, dvp, tat):
  twt = 2 * np.sqrt((xs-x0)**2 + (ys-y0)**2 + (zs-z0)**2) / (vpw+dvp) + tat
  return twt

'''
A function to correct measured travel-times by accounting for the ship's radial
velocity.
'''
def tt_corr(x0, y0, z0, xs, ys, zs, vs, vp, dvp, tts):
  # Calculate ship's position unit vector.
  r = np.array([xs - x0, ys - y0, zs - z0])
  r_hat = r / np.sqrt( np.sum(r**2, axis=0) )

  # Correct travel-times for ship's radial velocity.
  vr = np.sum(vs.T * r_hat, axis=0)
  dr = vr * tts
  corrections = dr / (vp + dvp)
  corrected = tts + corrections # (+/-) if logging ship location at
                                # (receive/transmit) time.
  # Return. 
  return corrected, corrections, vr

'''
This function computes the location corrections for the case where the GPS and
hull/side transponder are not co-located. The values dE and dN are distances in
METERS that should be ADDED to the locations that come from the GPS to give the
location of the transponder. They are positive East and North, respectively. The
values dforward and dstarboard are positive if the transponder is further 
forward (towards the prow) of the ship than the GPS, and if the transponder is
further starboard than the GPS, respectively. COG is the ship heading 
("course over ground") in DEGREES.
 
For instance, in the Chagall-esque artwork below, for a ship (sailing up the
page) with transponder (T) and GPS (G):
  the value of dforward is positive (+6 ASCII units) and 
  the value of dstarboard is negative (-5 ASCII units)
       __
     /    \  
    /      \
   /        \
  |  T       | 
  |          | 
  |          | 
  |          | 
  |          | 
  |          | 
  |       G  | 
  |__________| 
'''
def GPS_transp_correction(dforward, dstarboard, COG):
  # Calculate azimuth from GPS to transponder in ship reference frame.
  theta = np.rad2deg(np.arctan2(dstarboard, dforward))
  
  # Calculate azimuth from GPS to transponder in geographic reference frame.
  phi = theta + COG
  
  # Calculate absolute horizontal distance from GPS to transponder.
  dr = np.sqrt(dforward**2 + dstarboard**2)
  
  # Calculate East and North offset of transponder from GPS.
  dE = dr * np.rad2deg(np.sin(phi))
  dN = dr * np.rad2deg(np.cos(phi))

  return dE, dN

'''
A function to build the G matrix for the inversion.
'''
def G(x, y, z, dvp, x_ship, y_ship, z_ship, vpw, Nobs, M):
  # Initialize G matrix.
  G = np.ndarray(shape=(Nobs, M))
  
  # Anonymous function for calculating distances.
  D = lambda x0, y0, z0, x, y, z : np.sqrt((x0-x)**2 + (y0-y)**2 + (z0-z)**2)

  # First column: dti/dx
  G[:,0] = -(x_ship-x) * 2 / (vpw+dvp) / D(x_ship, y_ship, z_ship, x, y, z)
  # Second column: dti/dy
  G[:,1] = -(y_ship-y) * 2 / (vpw+dvp) / D(x_ship, y_ship, z_ship, x, y, z)
  # Third column: dti/dz
  G[:,2] = -(z_ship-z) * 2 / (vpw+dvp) / D(x_ship, y_ship, z_ship, x, y, z)
  # Fifth column: dti/ddvp
  G[:,3] = -2 * D(x_ship, y_ship, z_ship, x, y, z) / (vpw+dvp)**2

  return G