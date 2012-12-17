include makefile.def
# set OpenMP option
OMP=-fopenmp
ifeq ($(CXX), icc)
	OMP=-openmp
else ifeq ($(CXX), nvcc)
	OMP=-Xcompiler -fopenmp
endif

NAME=cppmafia
SRC_C=src/*.cpp
SRC_CU=src/*.cu
SRC_H=src/*.h
SRC=$(SRC_C) $(SRC_CU) $(SRC_H)
TGT=bin/$(NAME)
TMP=*~ src/*~ $(TGT) *.o tst/data/cluster-*.idx tst/data/cluster-*.dat *.out *.err

build: $(TGT) 
$(TGT):	$(SRC) makefile.def
	$(CXX) -arch=sm_20 -O3 $(OMP) $(SRC_C) $(SRC_CU) -o $(TGT)

profil:
	$(CXX) -arch=sm_20 -O3 $(OMP) -pg $(SRC_C) -o $(TGT)

#test run on a local node
run: $(TGT)
	$(TGT) --verbose tst/data/cluster.dat

#test run in the cluster
batch: $(TGT)
	msub cppmafia.sh

clean:
	rm -f $(TMP)
