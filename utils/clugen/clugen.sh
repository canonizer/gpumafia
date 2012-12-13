#! /bin/bash
# generates a number of .dat files with specific parameters

n=20000
d=20

for k in {2..20}; do
		bin/clugen -n$n -d$d -k$k /private/adinetz/cluster-n20k-20d/cluster-$k.dat
done
