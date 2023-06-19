g++ -o clarion_pid_root main_pid.o pxi16reader.o param.o detmap.o germanium.o gaggpid.o gaggpid_dec.o kinmat.o analysis.o spec.o gaggcut.o -lm `root-config --cflags --libs`

