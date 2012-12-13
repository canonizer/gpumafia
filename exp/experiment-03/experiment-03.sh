#! /bin/bash
# to launch from this directory only

# done with point-first layout

# compile for profiling
make -C ../../cppmafia profil

# compile clugen, if not already
make -C ../../utils/clugen

# clean cluster data and generate new data
DATA_DIR=/private/adinetz/cluster-n1m-d10
mkdir -p $DATA_DIR
rm -rf $DATA_DIR/*
for k in {2..10}; do
    ../../utils/clugen/bin/clugen -n1000000 -d10 -k$k \
				$DATA_DIR/cluster-$k.dat
done

# generate profile data; also alternate between bitmap and direct
for s in {direct,bitmap}; do
		sopt=''
		if [ $s == direct ]; then
				sopt=--no-bitmap
		fi
		for k in {2..10}; do
				../../cppmafia/bin/cppmafia -n100 $sopt	\
						$DATA_DIR/cluster-$k.dat
				gprof ../../cppmafia/bin/cppmafia -q gmon.out | head -100 > prof-$s-$k.log
				rm gmon.out
		done
done
