'''
FUNCTION SET output.py

A function to save the results outputted from locate.instruments(). 

A second function which checks to see if a given survey file has already been
processed.

Stephen M. 4/23/18
'''
# Import modules and functions
import os
import pickle
from funcs import txt

def out(results, figs, pkls, plts, txts):
  # Print statement.
  print('\n Writing results ...')

  # Create and write a .pkl file containing the final results.
  pkl_path = pkls + results['sta'] + '_out.pkl'
  f = open(pkl_path, 'wb')
  pickle.dump(results, f) 
  f.close()
  
  # Create and write a .txt file containing the final results.
  txt_path = txts + results['sta'] + '_location.txt'
  txt_fle = open(txt_path, 'w')
  txt.build(txt_fle, results)

  # Save figures in the output directory.
  for i, fig in enumerate(figs):
    fig_path = plts + results['sta'] + '_fig' + str(i+1) +'.pdf'
    fig.savefig(fig_path, bbox_inches='tight', format='pdf')

def exists(fle, txt_out):
  # Use the .txt file output as a proxy for whether file has been processed.
  sta = fle.split('/')[-1].split('.txt')[0]
  txt_path = txt_out + sta + '_location.txt'
  
  if os.path.exists(txt_path):
    print('\n Survey file ' + sta + ' processed. Skipping ...')
    return True
  else:
    return False