dualBEAGLE

~/msdir/ms 500 1 -s 8 > ms.sim    # sim with ms
./ms2BEAGLE ms.sim ms.sim.bgl     # convert to BEAGLE format




# compile statements
gcc ms2BEAGLE.c -o ms2BEAGLE

