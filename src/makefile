#===============================================================================
# User Options
#===============================================================================

COMPILER  = gnu
OPTIMIZE  = yes
DEBUG     = no
PROFILE   = no

#===============================================================================
# Program name & source code list
#===============================================================================

program = multibench

source = \
main.c \
io.c

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

#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(obj) multibench.h
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf multibench $(obj)

edit:
	vim -p $(source) multibench.h

run:
	./multibench
