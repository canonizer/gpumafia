NAME=cppmafia
SRC_C=src/*.cpp
SRC_H=src/*.h
SRC=$(SRC_C) $(SRC_H)
TGT=bin/$(NAME)
TMP=*~ src/*~ $(TGT) *.o tst/data/cluster-*.idx tst/data/cluster-*.dat

build: $(TGT)
$(TGT):	$(SRC)
	g++ -O3 -fopenmp $(SRC_C) -o $(TGT)

profil:
	g++ -O3 -fopenmp -pg $(SRC_C) -o $(TGT)

run: $(TGT)
	$(TGT) --verbose tst/data/cluster.dat

clean:
	rm -f $(TMP)
