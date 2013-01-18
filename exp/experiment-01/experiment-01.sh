#! /bin/bash
# to launch from this directory only

# done in point-first layout

# compile for profiling
make -C ../../cppmafia profil

# compile clugen, if not already
make -C ../../utils/clugen

# clean cluster data and generate new data
rm -rf /private/adinetz/cluster-n20k-20d/*
for k in {8..17}; do
    ../../utils/clugen/bin/clugen -n20000 -d20 -k$k \
				/private/adinetz/cluster-n20k-20d/cluster-$k.dat
done

# generate profile data
for k in {8..17}; do
		../../cppmafia/bin/cppmafia -n100	--no-set-dedup --no-bitmap --seq\
				/private/adinetz/cluster-n20k-20d/cluster-$k.dat
		gprof ../../cppmafia/bin/cppmafia -q gmon.out | head -50 > prof-$k.log
		rm gmon.out
done
