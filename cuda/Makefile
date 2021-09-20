#===============================================================================
# User Options
#===============================================================================

COMPILER    = nvidia
OPTIMIZE    = yes
DEBUG       = no
PROFILE     = no
SM_VERSION  = 80

#===============================================================================
# Program name & source code list
#===============================================================================

program = rsbench

source = \
main.cu \
simulation.cu\
io.cu \
init.cu \
material.cu \
utils.cu

obj = $(source:.cu=.o)

#===============================================================================
# Sets Flags
#===============================================================================

# Standard Flags
CFLAGS := 

# NVIDIA Compiler
ifeq ($(COMPILER),nvidia)
  CC = nvcc
  CFLAGS += -Xcompiler -Wall -Xcompiler -O3 -arch=sm_$(SM_VERSION) -std=c++14
endif

# Linker Flags
LDFLAGS = -lm

# Debug Flags
ifeq ($(DEBUG),yes)
  CFLAGS += -g -G
  LDFLAGS += -g -G
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

$(program): $(obj) rsbench.cuh Makefile
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.cu rsbench.cuh Makefile
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf rsbench $(obj)

edit:
	vim -p $(source) rsbench.cuh

run:
	./rsbench
