include makefile.def

# set OpenMP option
OMP=-fopenmp
ifeq ($(CXX), icc)
	OMP=-openmp
else ifeq ($(CXX), nvcc)
	OMP=-Xcompiler -fopenmp
endif

# set the flags
FLAGS := 
FLAGS += $(OMP)
FLAGS += -O3

NAME=cppmafia
SRC_C=src/*.cpp
SRC_CU=src/*.cu
SRC_H=src/*.h
SRC=$(SRC_C) $(SRC_H)
TGT=bin/$(NAME)
TMP=*~ src/*~ $(TGT) *.o tst/data/cluster-*.idx *.out \
 *.err

# set device support
ifeq ($(CXX), nvcc)
	GENCODE_SM10    := -gencode arch=compute_10,code=sm_10
	GENCODE_SM20    := -gencode arch=compute_20,code=sm_20
	GENCODE_SM30    := -gencode arch=compute_30,code=sm_30
	GENCODE_SM35    := -gencode arch=compute_35,code=sm_35
	GENCODE_FLAGS   := $(GENCODE_SM20) $(GENCODE_SM35)
#	FLAGS += -DKEPLER_OPTIMIZATIONS $(GENCODE_FLAGS) -lineinfo
	FLAGS += $(GENCODE_FLAGS) -lineinfo
	FLAGS += -DMAFIA_USE_DEVICE
	SRC += $(SRC_CU)
	SRC_C += $(SRC_CU)
endif

all: build

build: $(TGT) 
$(TGT):	$(SRC) makefile makefile.def
	$(CXX) $(FLAGS) $(SRC_C) -o $(TGT)

profil:
	$(CXX) $(FLAGS) -pg $(SRC_C) -o $(TGT)

debug:
	$(CXX) $(FLAGS) -g $(SRC_C) -o $(TGT)

#test run on a local node
run: $(TGT)
	$(TGT) --verbose tst/data/cluster.dat

clean:
	rm -f $(TMP)

#installs the program
install:
	cp $(TGT) $(PREFIX)/bin/$(NAME)

uninstall:
	rm -f $(PREFIX)/bin/$(NAME)
