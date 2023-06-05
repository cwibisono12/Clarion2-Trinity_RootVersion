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

