dualBEAGLE

~/msdir/ms 500 1 -s 8 > ms.sim    # sim with ms
~/msdir/ms 500 1 -t 10.04 -s 8 -r 100.0 2501 > ms.sim  # alternate sim
./ms2BEAGLE ms.sim ms.sim.bgl     # convert to BEAGLE format




# compile statements
gcc ms2BEAGLE.c -o ms2BEAGLE

