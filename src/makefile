#===============================================================================
# User Options
#===============================================================================

COMPILER  = gnu
OPTIMIZE  = yes
DEBUG     = no
PROFILE   = no
STATUS    = yes
PAPI      = no

#===============================================================================
# Program name & source code list
#===============================================================================

program = rsbench

source = \
main.c \
xs_kernel.c \
io.c \
init.c \
CUDA_init.c\
material.c \
utils.c \
CUDA_Driver.cu

OBJS = \
main.o \
xs_kernel.o\
io.o \
init.o \
CUDA_init.o\
material.o\
utils.o\
CUDA_Driver.o

HEADERS = \
rsbench.h #\
#cudaHeader.h 

obj = $(source:.c=.o)

#===============================================================================
# Sets Flags
#===============================================================================

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  CC = gcc
endif

# BG/Q gcc Cross-Compiler
ifeq ($(MACHINE),bluegene)
  CC = /bgsys/drivers/toolchain/V1R1M2/gnu-linux/bin/powerpc64-bgq-linux-gcc
endif

# Standard Flags
CFLAGS := -fopenmp -std=gnu99 -I/usr/local/cuda/include/ -ffast-math -I.
CPPFLAGS := -fopenmp -I/usr/local/cuda/include/ -ffast-math -I.
GCC_O_FLAGS =
NVCC_FLAGS = -arch sm_13 -O3 # -Xptxas=-v -DDERIVED_ATOMIC_FUNCTIONS_H
# Linker Flags
LDFLAGS = -L/usr/local/cuda/lib64/ -lm -lcudart 

# Debug Flags
ifeq ($(DEBUG),yes)
  GCC_O_FLAGS += -g
  LDFLAGS += -g
  NVCC_FLAGS += -g -G -lineinfo 
endif

# Profiling Flags
ifeq ($(PROFILE),yes)
  GCC_O_FLAGS += -pg
  LDFLAGS  += -pg
endif

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
  GCC_O_FLAGS += -O3
endif

# Calculation Status % Complete Indicator
ifeq ($(STATUS),yes)
  GCC_O_FLAGS += -DSTATUS
endif

# PAPI Performance Counters
ifeq ($(PAPI),yes)
  GCC_O_FLAGS += -DPAPI
  source += papi.c
  LDFLAGS += -lpapi
endif

#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(OBJS) rsbench.h
	nvcc -Xcompiler "$(CFLAGS)" $(GCC_O_FLAGS)  $(OBJS) -o $@ $(LDFLAGS)

CUDA_Driver.o: CUDA_Driver.cu $(HEADERS)
	nvcc  $(NVCC_FLAGS) -Xcompiler="$(CPPFLAGS)" -c CUDA_Driver.cu -o CUDA_Driver.o -I.

%.o: %.c
	$(CC) $(CFLAGS) $(GCC_O_FLAGS) -c $< -o $@

clean:
	rm -rf rsbench rsbench.dSYM $(OBJS)

edit:
	vim -p $(source) papi.c rsbench.h

run:
	./rsbench
bgqrun:
	qsub -t 10 -n 1 -O test rsbench
