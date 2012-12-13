NAME=clugen
SRC_C=src/*.c
SRC=$(SRC_C)
TGT=bin/$(NAME)
TMP=$(TGT) *~ src/*~ 

build:	$(TGT)
$(TGT):	$(SRC)
	gcc -std=gnu99 $(SRC_C) -o $(TGT)

clean:
	rm -f $(TMP)
