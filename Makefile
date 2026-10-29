# Multi-GPU build of Q-POP-Diffraction.
#
# Compiles src/heffte/qpop_diffraction_mgpu.cu, with the headers beside it, using
# nvcc into the MPI program bin/qpop_diffraction_mgpu, linked against HeFFTe (distributed FFTs,
# cuFFT backend), cuFFT, parallel HDF5 and MPI.  `qpop-diffraction --mgpu` looks for the program
# on PATH, then in bin/.  See docs/installation.md.
#
#   make HEFFTE_DIR=/path/to/heffte HDF5_DIR=/path/to/parallel/hdf5 MPI_DIR=/path/to/mpi
#
# Heffte_ROOT, HDF5_DIR and MPICH_DIR from the environment are used when set.
# nvcc does not compile through the MPI compiler wrapper, so MPI is passed explicitly.

HEFFTE_DIR ?= $(Heffte_ROOT)
HDF5_DIR   ?=
MPI_DIR    ?= $(MPICH_DIR)
CUDA_ARCH  ?= sm_80
CUDA_HOME  ?= $(patsubst %/bin/nvcc,%,$(shell command -v nvcc))

HEFFTE_LIB := $(if $(HEFFTE_DIR),$(firstword $(wildcard $(HEFFTE_DIR)/lib64 $(HEFFTE_DIR)/lib)))
# cuFFT sits with the CUDA libraries, or in math_libs/<version> of the NVIDIA HPC SDK.
CUFFT_DIR  ?= $(patsubst %/,%,$(dir $(firstword $(wildcard $(CUDA_HOME)/lib64/libcufft.so \
                  $(CUDA_HOME)/../../math_libs/$(notdir $(CUDA_HOME))/lib64/libcufft.so))))

# Cray MPICH aborts CUDA programs that do not link its GPU transport layer (GTL).
ifdef CRAY_MPICH_ROOTDIR
GTL_LIBS ?= -L$(CRAY_MPICH_ROOTDIR)/gtl/lib -lmpi_gtl_cuda -Xlinker -rpath=$(CRAY_MPICH_ROOTDIR)/gtl/lib
endif

TARGET := bin/qpop_diffraction_mgpu
SRC    := src/heffte/qpop_diffraction_mgpu.cu
DEPS   := $(wildcard src/heffte/*.hpp src/heffte/*.cuh)

FLAGS := -O3 -std=c++17 -arch=$(CUDA_ARCH) -cudart shared
INC   := -I$(HEFFTE_DIR)/include -I$(HDF5_DIR)/include -I$(MPI_DIR)/include
LIBS  := -L$(HEFFTE_LIB) -lheffte -L$(CUFFT_DIR) -lcufft \
         -L$(HDF5_DIR)/lib -lhdf5 \
         -L$(MPI_DIR)/lib -lmpi $(GTL_LIBS) \
         -Xlinker -rpath=$(HEFFTE_LIB):$(CUFFT_DIR):$(HDF5_DIR)/lib:$(MPI_DIR)/lib

all: $(TARGET)

$(TARGET): $(SRC) $(DEPS)
	@test -n "$(HEFFTE_LIB)" || (echo "HeFFTe not found: pass HEFFTE_DIR=/path/to/heffte"; exit 1)
	@test -n "$(HDF5_DIR)" || (echo "Parallel HDF5 not found: pass HDF5_DIR=/path/to/parallel/hdf5"; exit 1)
	@test -n "$(MPI_DIR)" || (echo "MPI not found: pass MPI_DIR=/path/to/mpi"; exit 1)
	@test -n "$(CUFFT_DIR)" || (echo "cuFFT not found: pass CUFFT_DIR=/path/to/dir/with/libcufft.so"; exit 1)
	@mkdir -p bin
	nvcc $(FLAGS) $(INC) $< -o $@ $(LIBS)

clean:
	rm -f $(TARGET)

.PHONY: all clean
