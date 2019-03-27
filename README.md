# OBSrange
### A robust, efficient, open-source OBS location code for the marine geophysics community.

OBSrange is a set of scripts written in both MATLAB and Python for precisely locating ocean bottom seismometers (OBSs) using acoustic transponder data. The code reads data from acoustic ranging surveys that measure two-way travel times from ship to instrument as the ship runs a survey pattern. Using these data, imported in survey files, the code inverts for instrument locations, and depth averaged sound speeds in water. Additionally, OBSrange generates several figures visualizing these results as well as estimates of parameter uncertainties. For a detailed description of the inversion algorithm, synthetic tests, and results from the 2018 Young Pacific ORCA experiment, please refer to the paper under the "Files" tab.

OBSrange features:
* Compute precise locations for OBSs on the seafloor from acoustic survey travel time data
* Produce and plot detailed and comprehensive information about location uncertainties
* Solve for depth and mean water sound speed
* Correct travel times for "doppler" ship movement effects
* Correct for a known horizontal offset between shipboard GPS and transponder
* Correct travel times for ray bending due to refraction through the water column using either a custom 1D sound-speed profile or a geographically appropriate profile chosen from the World Ocean Atlas database

-------
Russell, J.B., Z. Eilon, & S. Mosher (2019) OBSrange: A new tool for the precise remote location of Ocean Bottom Seismometers, SRL, xx, xx-xx.