#!/bin/bash -x
#MSUB -l nodes=1:ppn=24:gpus=1:performance
#MSUB -l walltime=00:02:00
#MSUB -e ./cppmafia-omp-0.err
# if keyword omitted : default is submitting directory
#MSUB -o ./cppmafia-omp-0.out
# if keyword omitted : default is submitting directory
#MSUB -v tpt=1
# for OpenMP/hybrid jobs only
# for OpenMP jobs only
cd $PBS_O_WORKDIR
echo "workdir: $PBS_O_WORKDIR"
NSLOTS=1
echo "running on $NSLOTS cpus ..."
