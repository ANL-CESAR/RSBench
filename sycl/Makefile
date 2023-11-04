#===============================================================================
# User Options
#===============================================================================

CC        = llvm
OPTIMIZE  = yes
DEBUG     = no
PROFILE   = no
TARGET ?=
SM_VERSION ?= 
OFFLOAD_ARCH ?=

#===============================================================================
# Program name & source code list
#===============================================================================

program = RSBench

source = \
main.cpp \
simulation.cpp\
io.cpp \
init.cpp \
material.cpp \
utils.cpp

obj = $(source:.cpp=.o)

#===============================================================================
# Sets Flags
#===============================================================================

# Standard Flags
CFLAGS := -std=c++17 -Wall 

# Linker Flags
LDFLAGS = -lm

# LLVM Compiler
ifeq ($(CC),llvm)
  CC = clang++
  CFLAGS += -fsycl
endif

# IntelOneAPI Compiler
ifeq ($(CC),oneapi)
  CC = icpx
  CFLAGS += -fsycl
endif

# AdaptiveCpp
ifeq ($(CC),adaptivecpp)
  CC = acpp
endif

ifeq ($(TARGET),CUDA)
  ifeq ($(CC),acpp)
	CFLAGS += --acpp-targets="cuda:sm_${SM_VERSION}"
  else
	CFLAGS += -fsycl-targets=nvptx64-nvidia-cuda \
			  -Xsycl-target-backend --cuda-gpu-arch=sm_${SM_VERSION}
  endif
endif

ifeq ($(TARGET),HIP)
  ifeq ($(CC),acpp)
	CFLAGS += --acpp-targets="hip:${OFFLOAD_ARCH}"
  else
	CFLAGS += -fsycl-targets=amdgcn-amd-amdhsa \
			  -Xsycl-target-backend --offload-arch=${OFFLOAD_ARCH}
  endif
endif

# Debug Flags
ifeq ($(DEBUG),yes)
  CFLAGS += -g
  LDFLAGS  += -g
endif

# Profiling Flags
ifeq ($(PROFILE),yes)
  CFLAGS += -pg
  LDFLAGS  += -pg
endif

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
  CFLAGS += -O3
endif

#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(obj) rsbench.h Makefile
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.cpp rsbench.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf rsbench $(obj)

edit:
	vim -p $(source) rsbench.h

run:
	./rsbench
