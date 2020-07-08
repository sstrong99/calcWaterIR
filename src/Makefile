#On Midway, first load modules:
#module load gcc/4.8 cuda/8.0 mkl/2017

exe     = calcIR
CXX     = g++
NVCC    = nvcc

HOST=$(firstword $(subst -, ,$(shell hostname)))
ifeq ($(HOST),midway)
	SOFTDIR=/home/strong1/install
else
	SOFTDIR=/usr/local
endif
## $(info $$SOFTDIR is [${SOFTDIR}])

CUDAROOT?=$(CUDADIR)

INC     = -I$(SOFTDIR)/include/xdrfile -I$(MKLROOT)/include
FLAGS	= -O2 -march=native -Wall -std=c++11 -MMD -MP -fopenmp -DUSEOMP
DEBUGFLAGS = -g -std=c++11 -Wall -fopenmp -MMD -MP #-DUSEOMP #-DDEBUG -pg
LIBS    = -lxdrfile -lfftw3f -lm -llapack -lblas -ldl -lgomp -lpthread
LIBDIRS = -L$(SOFTDIR)/lib -L$(SOFTDIR)/magma/lib
GPUINC  = -I$(SOFTDIR)/magma/include -I$(CUDAROOT)/include $(INC)
GPULIBS = -lmagma -Xcompiler "-fopenmp" $(LIBS)

#MMD and MP generate .d files
OBJ_DIR = OBJ
SRCS  = main.cpp traj.cpp xtcTraj.cpp calcIR.cpp integrateF.cpp adamsBashforth.cpp exactDiag.cpp timer.cpp histogram.cpp calcDists.cpp calcW.cpp initTraj.cpp input.cpp calcIR_TAA.cpp calculation.cpp

OBJS := $(SRCS:%.cpp=$(OBJ_DIR)/%.o)
DEPS := $(OBJS:%.o=%.d)
OBJS_GPU := $(SRCS:%.cpp=$(OBJ_DIR)/%_g.o)
DEPS_GPU := $(OBJS:%.o=%_g.d)
OBJS_DGPU := $(SRCS:%.cpp=$(OBJ_DIR)/%_dg.o)
DEPS_DGPU := $(OBJS:%.o=%_dg.d)
OBJS_DEBUG := $(SRCS:%.cpp=$(OBJ_DIR)/%_d.o)
DEPS_DEBUG := $(OBJS:%.o=%_d.d)

exeg :=  $(exe)
exedg :=  "$(exe)_debugGPU"
execpu := "$(exe)_cpu"
exed := "$(exe)_debug"

gpu: $(exeg)
dgpu: $(exedg)
cpu: $(execpu)
debug: $(exed)

all: gpu dgpu cpu debug

#link
$(execpu): $(OBJS)
	$(CXX) $(FLAGS) -o $(execpu) $^ $(LIBDIRS) $(LIBS)
$(exeg): $(OBJS_GPU)
	$(NVCC) -Xcompiler "$(FLAGS)" -o $(exeg) $^ $(LIBDIRS) $(GPULIBS)
$(exedg): $(OBJS_DGPU)
	$(NVCC) -G -Xcompiler "$(DEBUGFLAGS)" -o $(exedg) $^ $(LIBDIRS) $(GPULIBS)
$(exed): $(OBJS_DEBUG)
	$(CXX) $(DEBUGFLAGS) -o $(exed) $^ $(LIBDIRS) $(LIBS)

#compile
$(OBJ_DIR)/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(INC) $(FLAGS) -c -o $@ $<
$(OBJ_DIR)/%_g.o: %.cpp
	mkdir -p $(dir $@)
	$(NVCC) $(GPUINC) -Xcompiler "$(FLAGS)" -DUSEGPU -c -o $@ $<
$(OBJ_DIR)/%_dg.o: %.cpp
	mkdir -p $(dir $@)
	$(NVCC) $(GPUINC) -Xcompiler "$(DEBUGFLAGS)" -DUSEGPU -c -o $@ $<
$(OBJ_DIR)/%_d.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(INC) $(DEBUGFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm -rf "$(OBJ_DIR)"

-include $(DEPS)

#with help from https://spin.atomicobject.com/2016/08/26/makefile-c-projects/
