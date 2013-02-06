#! /bin/bash
# to launch from this directory only

# done with dimension-first layout
# TODO: set to 1 if the data is not generated yet
GENERATE=1
DATA_DIR=~/try/mafia/cluster-n10m-d20-m5
DEVICE=k20x

# compile for profiling
make -C ../../cppmafia

# compile clugen, if not already
make -C ../../utils/clugen

# clean cluster data and generate new data
mkdir -p $DATA_DIR
if [ $GENERATE == 1 ]; then
		rm -rf $DATA_DIR/*
		for k in {3..16}; do
				../../utils/clugen/bin/clugen -n10000000 -d20 -k$k -m3 -N \
						$DATA_DIR/cluster-$k.dat
		done
fi

for k in {3..16}; do
		OMP_NUM_THREADS=1 ../../cppmafia/bin/cppmafia --timing -d \
				$DATA_DIR/cluster-$k.dat > time-$DEVICE-$k.log
done
