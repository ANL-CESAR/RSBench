#===============================================================================
# User Options
#===============================================================================

COMPILER  = gnu
OPTIMIZE  = yes
DEBUG     = no
PROFILE   = no
STATUS    = no
PAPI      = yes

#===============================================================================
# Program name & source code list
#===============================================================================

program = resobench

source = \
main.c \
xs_kernel.c \
io.c \
init.c \
material.c \
utils.c

obj = $(source:.c=.o)

#===============================================================================
# Sets Flags
#===============================================================================

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  CC = gcc
endif

# Standard Flags
CFLAGS := -fopenmp -std=gnu99

# Linker Flags
LDFLAGS = -lm

# Debug Flags
ifeq ($(DEBUG),yes)
  CFLAGS += -g
  LDFLAGS += -g
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

# Calculation Status % Complete Indicator
ifeq ($(STATUS),yes)
  CFLAGS += -DSTATUS
endif

# PAPI Performance Counters
ifeq ($(PAPI),yes)
  CFLAGS += -DPAPI
  source += papi.c
  LDFLAGS += -lpapi
endif

#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(obj) resobench.h
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf resobench $(obj)

edit:
	vim -p $(source) resobench.h

run:
	./resobench
