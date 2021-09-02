#===============================================================================
# User Options
#===============================================================================

COMPILER    = amd
OPTIMIZE    = yes
DEBUG       = no
PROFILE     = no
SM_VERSION  = 70

#===============================================================================
# Program name & source code list
#===============================================================================

program = rsbench

source = \
main.hip \
simulation.hip\
io.hip \
init.hip \
material.hip \
utils.hip

obj = $(source:.hip=.o)

#===============================================================================
# Sets Flags
#===============================================================================

# Standard Flags
CFLAGS := 

# Regular gcc Compiler
ifeq ($(COMPILER),nvidia)
  CC = nvcc
  CFLAGS += -Xcompiler -Wall -Xcompiler -O3 -arch=sm_$(SM_VERSION) -std=c++11
endif

# AMD
ifeq ($(COMPILER),amd)
  CC = hipcc
  CFLAGS += -std=c++14
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

$(program): $(obj) rsbench.h Makefile
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.hip rsbench.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf rsbench $(obj)

edit:
	vim -p $(source) rsbench.h

run:
	./rsbench
