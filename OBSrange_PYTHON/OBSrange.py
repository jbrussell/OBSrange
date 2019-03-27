'''
SCRIPT OBSrange.py

This script sets parameters and loops through several survey files to locate
instruments on the seafloor via a bootstrap inversion procedure. The main 
function that gets called is locate.instruments(). The guts of the code are 
contained in that function and several functions within. The results, written in
an output directory, consist of .pkl and .txt files of the inversion results, as
well as several figures. 

Stephen M. & Zach E. & Josh R. 4/23/18
'''
# Import modules and functions
import os
import pickle
from funcs import fetch, locate, output

############################# Inversion Parameters #############################

vpw = 1500.0            # Assumed velocity of sound in water [m/s]
dvp = 0                 # Assumed perturbation to vpw [m/s]
tat = 0.014             # Assumed sensor turn-around time [s]
N_bs = 1000             # Number of bootstrap iterations
E_thresh = 1e-5         # RMS reduction threshold for inversion
npts = 1                # Npts in moving avg of ship vel. (1 = no smoothing)
dampx = 0               # Norm damping for each model parameter             
dampy = 0               #             .
dampz = 0               #             .    
dampdvp = 5e-8          #             .      
eps = 1e-10             # Global norm damping for stabilization
QC = True               # Option to perform quality control on survey points
res_thresh = 500        # Threshold [ms] beyond which survey points are tossed
dforward = 0            # GPS-transp offset [m] (+ means trans. further forward)
dstarboard = 0          # GPS-transp offset [m] (+ means trans. further stboard)
twtcorr = False   # Option to apply a travel-time correction for ship velocity
raycorr = True    # Option to apply a travel-time correction due to ray bending:
                  #   If you choose to correct for ray beding you can either
                  #   input your own depth-soundspeed profile OR our code will 
                  #   calculate one for you see the README for further details.

parameters = [vpw, dvp, tat, N_bs, E_thresh, npts, dampx, dampy, dampz, dampdvp,
              eps, QC, res_thresh, dforward, dstarboard, twtcorr, raycorr]

################################ Directory Setup ###############################

survey_dir = '../new_project/survey_files/'
output_dir = '../new_project/output/'
ssp_dir = '../new_project/output/station_sound_speed_profiles'

###################### Run Inversion for Each Survey File ######################

# Grab survey files.
survey_fles = fetch.data_paths(survey_dir, matchkey='*.txt')

# Create output sub-directories.
out_pkls = output_dir + 'data_pkls/'
out_plts = output_dir + 'plots/'
out_txts = output_dir + 'data_txts/'
#out_ssps = output_dir + 'station_sound_speed_profiles/'
if not os.path.exists(output_dir):
  os.mkdir(output_dir)
  os.mkdir(out_pkls)
  os.mkdir(out_plts)
  os.mkdir(out_txts)
if not os.path.exists(ssp_dir):
    os.mkdir(ssp_dir)

# Perform locations for each survey site then write results to output.
for survey_fle in survey_fles:
  results, figs = locate.instruments(survey_fle, parameters, ssp_dir)
  output.out(results, figs, out_pkls, out_plts, out_txts)
  
##################################### FIN ######################################
