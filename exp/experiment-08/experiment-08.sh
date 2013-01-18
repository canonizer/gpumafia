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

mkdir -p $DATA_DIR
../../utils/clugen/bin/clugen -n10000000 -d10 -k10 $DATA_DIR/cluster-10.dat

#export GOMP_CPU_AFFINITY="0-23"

# generate profile data for different number of OpenMP threads
for nt in {1..24}; do
		cp cppmafia-omp.sh cppmafia-tmp-0.sh
		echo export OMP_NUM_THREADS=$nt >> cppmafia-tmp-0.sh
		echo mpiexec -np '$NSLOTS' ../../cppmafia/bin/cppmafia --timing \
				$DATA_DIR/cluster-10.dat >> cppmafia-tmp-0.sh
		# submit and wait for completion
		TASK_ID=`msub -qdevel cppmafia-tmp-0.sh`
		while qstat $TASK_ID 2>&1>/dev/null; do 
				sleep 1
		done
		# get the time data
		mv cppmafia-omp-0.out time-$nt.log
		rm cppmafia-omp-0.err
done
