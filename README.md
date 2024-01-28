# Clarion2-Trinity_RootVersion

This version of data processing is intended for a user who wants to analyze Clarion2-Trinity using CERN Root software. The output of this program is root objects that consists of TTree and some of TH2 histograms. 

To compile user can run ./build_root.sh

To link user can run ./link_root.sh

TTree that are generated from this sort-code contains Ge data after add-back and Compton Suprression. If a user wants to generate the very raw data, user can modify some of the source codes. Also, the data from GAGG is the post-processing data after peak and tail reconstruction from the trace or QDC sum. Therefore the TTree from this sorting code is already in the form of Physics Data (No need to build or reconstuct the Detector).


#Description of TTree:

The structure of TTree from this sort-code generates the following:
Gemult;GeID;Ge-Energy;Ge-Time;Ge-Edopp;
GAGGmult;GAGGID;GAGG-Time;GAGG-peak;GAGG-tail;GAGG-traceint


#2D-Histograms:

The dafault sorting code from this software also generate 2D Histogram from raw channel and calibrated channel. The Y axis is ID and the X axis is channel. User can utilize object e_raw to calibrate Ge crystals. To check calibration, user can utilize object e_cal.


#Running the Sort Code:
./Clarion_root -up [time-ordered-file] [calibration file] [mapping file] [gagg_calibration_file] [banana gate file] [scaling factor for Doppler correction] [Doppler correction option] [reaction channel information] [output root file]

#Example:
./Clarion_root -up calEu152.evt.to cal_fsu.ca3 id_fsu_Jun9.map gagg_calib_proton.txt 2d_all_proton.txt 0.8 1 proton.txt calEu152.root

All of necessary parameter can be found under param directory. Note that if we only have Germanium detector without GAGG, this sorting code still requires the charged particle parameter files though, they do not affect the output of TTree. If there are no charged particle detectors, the reading of gaggmult in TTree would be zero and no other branches related to GAGG will be generated.


#Added Support for NSCLDAQ data format:
This repository also supports the data format generated from the NSCLDAQ which have been event-built. Directory is located at src/nscl/. The output would be in forms of TTree where the structure of TTree can be found at src/nscl/pxitree.h

After compiling and linking (1.build_evtnsclroot.sh 2.link_evtnsclroot.sh) the program can be run via:

./pixie16nscl_root .evtfile .rootfile 
