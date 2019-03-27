'''
FUNCTION txt.py

A function to format and write inversion results into a .txt file.

Stephen M. 4/23/18
'''
# Import modules and functions
from numpy import std as S
from numpy import mean as M
from numpy import sqrt as sq

def build(fle, R):
  # Fancy characters (degree sign, lower-case sigma, plus-minus)
  deg = u'\N{DEGREE SIGN}'
  sig = u'\u03C3'
  pm = u'\u00b1'

  # Anonymous functions and formatting strings to make things more tractable.
  dst = lambda rx,ry,rz,x,y,z : sq((M(rx)-x)**2 + (M(ry)-y)**2 + (M(rz)-z)**2)
  
  a = '{:>11} {: >11.5f} {}{: .5f}\n'
  fmt1 = lambda s1, f1, s2, f2 : a.format(s1, f1, s2, f2)

  b = '{:>3} {: >11.5f}  {: <11.5f} {: <11.5f} {: >8.5f}  {: >8.5f}  {: >8.5f}\n'
  fmt2 = lambda f1,f2,f3,f4,f5,f6,f7 : b.format(f1,f2,f3,f4,f5,f6,f7)

  c = '{:>14}   {:<11} {:<11} {:<10} {:<9} {}'

  # .txt header and main results.
  fle.write('Bootstrap Inversion Results ('+pm+' 2' + sig + ' uncertainty)\n')
  fle.write('\n')
  fle.write('{:>11} {} \n'.format('Station:', R['sta']))
  fle.write(fmt1('Lat:', M(R['lat_sta']), deg+' '+pm, 2 * S(R['lat_sta'])))
  fle.write(fmt1('Lon:', M(R['lon_sta']), deg+' '+pm, 2 * S(R['lon_sta'])))
  fle.write(fmt1('X:', M(R['x_sta']), 'm '+pm, 2 * S(R['x_sta'])))
  fle.write(fmt1('Y:', M(R['y_sta']), 'm '+pm, 2 * S(R['y_sta'])))
  fle.write(fmt1('Depth:', M(R['z_sta']), 'm '+pm, 2 * S(R['z_sta'])))
  fle.write(fmt1('Water Vel.:', M(R['vpws']), 'm '+pm, 2 * S(R['vpws'])))
  fle.write(fmt1('Drift:', M(R['drifts']), 'm '+pm, 2 * S(R['drifts'])))
  fle.write(fmt1('Drift Az:', M(R['azs']), 'm '+pm, 2 * S(R['azs'])))
  fle.write(fmt1('dz:', M(R['dzs']), 'm '+pm, 2 * S(R['dzs'])))
  fle.write('\n')
  fle.write(fmt1('RMS:', M(R['E_rms']*1000), 'ms '+pm, 2 * S(R['E_rms'])*1000))
  fle.write('\n')
  fle.write('{} {} \n'.format('Bad Pings Removed:', R['Nbad']))
  fle.write('{:=<80} \n'.format(''))

  # .txt survey points.
  fle.write( c.format(
                       'Lat (' + deg + ')',
                       'Lon (' + deg +')',
                       'Range (m)',
                       'Resid (s)',
                       'Vr (m/s)',
                       'TWT corr. (ms)\n' 
                       )
            )

  for i in range(len(R['svy_lats'])):
    fle.write(fmt2(
                   i+1, 
                   R['svy_lats'][i], 
                   R['svy_lons'][i], 
                   dst(R['x_sta'][i],
                       R['y_sta'][i],
                       R['z_sta'][i],
                       R['svy_xs'][i],
                       R['svy_ys'][i],
                       R['svy_zs'][i]),
                   R['dtwts'][i] * 1000,
                   R['vrs'][i],
                   R['corrs'][i] * 1000))
  fle.close()