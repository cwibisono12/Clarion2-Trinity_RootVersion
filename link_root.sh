g++ -o clarion_root main_ev3.o pxi16reader.o param.o detmap.o germanium.o gagg.o kinmat.o bangate.o write.o analysis.o spec.o to2root.o clariontree.o -lm `root-config --cflags --libs`

