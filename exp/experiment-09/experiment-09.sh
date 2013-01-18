#! /bin/bash
# to launch from this directory only

# done with dimension-first layout
# TODO: set to 1 if the data is not generated yet
GENERATE=0

# compile for profiling
make -C ../../cppmafia

# compile clugen, if not already
make -C ../../utils/clugen

# clean cluster data and generate new data
DATA_DIR=~/try/mafia/cluster-n10m-d10

mkdir -p $DATA_DIR
if [ $GENERATE == 1 ]; then
		rm -rf $DATA_DIR/*
		for k in {2..10}; do
				../../utils/clugen/bin/clugen -n10000000 -d10 -k$k \
						$DATA_DIR/cluster-$k.dat
		done
fi

# generate profile data; also alternate between bitmap and direct
for s in {seq,par,dev}; do
		sopt=''
		if [ $s == seq ]; then
				sopt=--seq
		elif [ $s == dev ]; then
				sopt="--seq --device"
		fi
		for k in {2..10}; do
				cp cppmafia-dev.sh cppmafia-tmp-0.sh
				echo mpiexec -np '$NSLOTS' ../../cppmafia/bin/cppmafia --timing \
						$sopt	$DATA_DIR/cluster-$k.dat >> cppmafia-tmp-0.sh
				# submit and wait for completion
				TASK_ID=`msub -qdevel cppmafia-tmp-0.sh`
				while qstat $TASK_ID 2>&1>/dev/null; do 
						sleep 1
				done
				# get the time data
				mv cppmafia-dev-0.out time-$s-$k.log
				rm cppmafia-dev-0.err
		done
done
