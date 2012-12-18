#!/bin/bash -x
#MSUB -l nodes=1:ppn=24:gpus=1
#MSUB -l walltime=00:01:30
#MSUB -e ./cppmafia-dev-0.err
# if keyword omitted : default is submitting directory
#MSUB -o ./cppmafia-dev-0.out
# if keyword omitted : default is submitting directory
#MSUB -v tpt=1
# for OpenMP/hybrid jobs only

### start of jobscript
export OMP_NUM_THREADS=1
# for OpenMP jobs only
cd $PBS_O_WORKDIR
echo "workdir: $PBS_O_WORKDIR"
NSLOTS=1
echo "running on $NSLOTS cpus ..."
