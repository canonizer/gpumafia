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
DATA_DIR=~/try/mafia/cluster-n100k-d30

mkdir -p $DATA_DIR
if [ $GENERATE == 1 ]; then
		rm -rf $DATA_DIR/*
		for k in {2..17}; do
				../../utils/clugen/bin/clugen -n100000 -d30 -k$k \
						$DATA_DIR/cluster-$k.dat
		done
fi

# generate profile data; also alternate between bitmap and direct
for s in {opt,noopt}; do
		sopt="--seq"
		if [ $s == noopt ]; then
				sopt="--seq --no-bitmap --no-set-dedup --no-set-gen"
		fi
		for k in {2..17}; do
				../../cppmafia/bin/cppmafia --timing $sopt \
						$DATA_DIR/cluster-$k.dat > time-$s-$k.log
		done
done
