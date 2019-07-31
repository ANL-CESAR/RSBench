#===============================================================================
# User Options
#===============================================================================

# Compiler can be set below, or via environment variable
ifeq ($(CC),)
CC        = gcc
endif
OPTIMIZE  = yes
DEBUG     = no
PROFILE   = no

#===============================================================================
# Program name & source code list
#===============================================================================

program = rsbench

source = \
main.c \
simulation.c\
io.c \
init.c \
material.c \
utils.c

obj = $(source:.c=.o)

#===============================================================================
# Sets Flags
#===============================================================================

# Make sure we didn't pick up lower case cc
ifeq (cc,$(CC))
CC    = gcc
endif

# Standard Flags
CFLAGS := -std=gnu99 -Wall

# gcc Compiler
ifneq (,$(findstring gcc,$(CC)))
  CFLAGS += -fopenmp -flto
endif

# LLVM Compiler
ifneq (,$(findstring clang,$(CC)))
  CFLAGS += -fopenmp -flto
endif

# Intel Compiler
ifneq (,$(findstring intel,$(CC)))
  CFLAGS += -qopenmp -ipo
endif

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

$(program): $(obj) rsbench.h Makefile
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.c rsbench.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf rsbench $(obj)

edit:
	vim -p $(source) rsbench.h

run:
	./rsbench
