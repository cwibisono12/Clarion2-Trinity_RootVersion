# Clarion2-Trinity_RootVersion

This version of data processing is intended for a user who wants to analyze Clarion2-Trinity using CERN Root software. The output of this program is root objects that consist TTree and some of TH2 histograms. 

To compile user can run ./build_root.sh

To link user can run ./link_root.sh

TTree that are generated from this sort-code contains Ge data after add-back and Compton Suprression. If a user want to generate the very raw data, user can modify some of the source codes. Also, the data from GAGG is the post-processing data after peak and tail reconstruction from the trace or QDC sum.

#Description of TTree:
The structure of TTree from this sort-code generates the following:
Gemult;GeID;Ge-Energy;Ge-Time;Ge-Edopp;
GAGGmult;GAGGID;GAGG-Time;GAGG-peak;GAGG-tail;GAGG-traceint
