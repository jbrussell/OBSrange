'''
Function fetch.py

Return a python list of absolute file paths to files within a given directory.
Can also specify unix style wildcards via a matchkey.

Stephen M.
'''
# Import modules and functions
import glob

def data_paths(data_dir, matchkey):
  print(' Listing files in %s' % data_dir)
  files = glob.glob(data_dir + matchkey)
  
  # Return file paths.
  return sorted(files)