NAME=clugen
SRC_C=src/*.cpp
SRC=$(SRC_C) src/*.h
TGT=bin/$(NAME)
TMP=$(TGT) *~ src/*~ 

build:	$(TGT)
$(TGT):	$(SRC)
	g++ $(SRC_C) -o $(TGT)

debug: 
	g++ $(SRC_C) -g -o $(TGT)

clean:
	rm -f $(TMP)
