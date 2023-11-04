#===============================================================================
# User Options
#===============================================================================

COMPILER    = amd
OPTIMIZE    = yes
DEBUG       = no
PROFILE     = no
OFFLOAD_ARCH ?= gfx90a

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
CFLAGS := 

# AMD
ifeq ($(COMPILER),amd)
  CC = hipcc
  CFLAGS += -std=c++14 --offload-arch=${OFFLOAD_ARCH}
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

%.o: %.cpp rsbench.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj)

edit:
	vim -p $(source) rsbench.h

run:
	./rsbench
