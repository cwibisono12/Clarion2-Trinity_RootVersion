# Clarion2-Trinity_RootVersion

This version of data processing is intended for a user who wants to analyze Clarion2-Trinity using CERN Root software. The output of this program is root objects that consists of TTree and some of TH2 histograms. 

To compile user can run ./build_root.sh

To link user can run ./link_root.sh

TTree that are generated from this sort-code contains Ge data after add-back and Compton Suprression. If a user wants to generate the very raw data, user can modify some of the source codes. Also, the data from GAGG is the post-processing data after peak and tail reconstruction from the trace or QDC sum. Therefore the TTree from this sorting code is already in the form of Physics Data (No need to build or reconstuct Detectors).


#Description of TTree:

The structure of TTree from this sort-code generates the following:
Gemult;GeID;Ge-Energy;Ge-Time;Ge-Edopp;
GAGGmult;GAGGID;GAGG-Time;GAGG-peak;GAGG-tail;GAGG-traceint


#2D-Histograms:

The default sorting code from this software also generate 2D Histogram from raw channel and calibrated channel. The Y axis is ID and the X axis is channel. User can utilize object e_raw to calibrate Ge crystals. To check calibration, user can utilize object e_cal.


#Running the Sort Code:
./Clarion_root -up [time-ordered-file] [calibration file] [mapping file] [gagg_calibration_file] [banana gate file] [scaling factor for Doppler correction] [Doppler correction option] [reaction channel information] [output root file]

#Example:
./Clarion_root -up calEu152.evt.to cal_fsu.ca3 id_fsu_Jun9.map gagg_calib_proton.txt protoncutex.root 0.8 1 proton.txt calEu152.root

All of necessary parameter can be found under param directory. Note that if we only have Germanium detector without GAGG, this sorting code still requires the charged particle parameter files though, they do not affect the output of TTree. If there are no charged particle detectors, the reading of gaggmult in TTree would be zero and no other branches related to GAGG will be generated.

#Generating GAGG-PID and making particle cuts file for channel selection:
PIDs for GAGG are created separately. Gamma gate condition can also be made to give more restrictive requirement so that the separation between charged particles are more distinct.  

Sorting Code above assumes that particle cuts have been made and the cuts were stored to a root file via TCutG class (see param/protoncutex.root for example of particle cuts file). The particle cuts file only contain TCutG objects. The number of TCutG objects has to be equal to the number of reconstructed GAGG. The name of TCutG objects has to follow the convention of the name used as in param/protoncutex.root file. 

Particle cuts file can be made after the PIDs were generated. Current standard PID for the analysis is trace-integral for the y axis and tail to peak ratio for the x axis. In order to produce the PIDs, one needs to compile differently with the one previously mentioned [see Clarion2-Trinity repository on this github as well].

To compile: ./build_pid_root.sh
To link: ./link_pid_root.sh

The procedure for running the program is the same except that now the executable name is clarion_pid_root. Also, there is an additional argument next to reaction channel parameter (before root output file) that is the information about the gamma gate energy window that we want to put to create gamma-gated pid.

#Example of how to generate PID:

./clarion_pid_root -up data.evt.to cal_fsu.ca3 id_fsu_Jun9.map gagg_calib_proton.txt protoncutex.root 0.8 1 proton.txt gamgatefile.txt out.root

All of GAGG PIDs will be stored on out.root file. Cuts can be made by calling each of objects inside out.root via ROOT Graphics editor by using View->Toolbar->Cut and store and rename each cut objects for each PIDs to a file according to naming condition as mentioned.
