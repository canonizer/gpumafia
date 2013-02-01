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
DATA_DIR=~/try/mafia/cluster-n10m-d20-m5

mkdir -p $DATA_DIR
if [ $GENERATE == 1 ]; then
		rm -rf $DATA_DIR/*
		for k in {3..16}; do
				../../utils/clugen/bin/clugen -n10000000 -d20 -k$k -m3 -N \
						$DATA_DIR/cluster-$k.dat
		done
fi

# generate profile data; also alternate between bitmap and direct
for s in {seq,par6,par24,fermi}; do
#for s in kepler; do
		sopt=""
		nt=1
		if [ $s == seq ]; then
				sopt="--seq"
		elif [ $s == par6 ]; then
				nt=6
		elif [ $s == par24 ]; then
				nt=24
		elif [ $s == fermi -o $s == kepler ]; then
				sopt="-d"
				nt=1
		fi
		for k in {3..16}; do
				OMP_NUM_THREADS=$nt ../../cppmafia/bin/cppmafia --timing $sopt \
						$DATA_DIR/cluster-$k.dat > time-$s-$k.log
		done
done
