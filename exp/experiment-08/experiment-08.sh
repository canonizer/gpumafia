#! /bin/bash
# to launch from this directory only

# done with dimension-first layout
# TODO: set to 1 if the data is not generated yet
GENERATE=1

# compile for profiling
make -C ../../cppmafia

# compile clugen, if not already
make -C ../../utils/clugen

# clean cluster data and generate new data
DATA_DIR=~/try/mafia/cluster-10m-10d
DATA_FILE=$DATA_DIR/cluster-10.dat

mkdir -p $DATA_DIR
../../utils/clugen/bin/clugen -n10000000 -d10 -k10 $DATA_FILE

#export GOMP_CPU_AFFINITY="0-32"

# generate profile data for different number of OpenMP threads
for nt in {1..32}; do
		OMP_NUM_THREADS=$nt ../../cppmafia/bin/cppmafia --timing $DATA_FILE \
				>> time-$nt.log
done
