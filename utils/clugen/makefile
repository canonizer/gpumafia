NAME=clugen
SRC_C=src/*.cpp
SRC=$(SRC_C)
TGT=bin/$(NAME)
TMP=$(TGT) *~ src/*~ 

build:	$(TGT)
$(TGT):	$(SRC)
	g++ $(SRC_C) -o $(TGT)

clean:
	rm -f $(TMP)
