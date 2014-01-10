==============================================================================

           __  __         _  _    _  ____                      _     
          |  \/  | _   _ | || |_ (_)| __ )   ___  _ __    ___ | |__  
          | |\/| || | | || || __|| ||  _ \  / _ \| '_ \  / __|| '_ \ 
          | |  | || |_| || || |_ | || |_) ||  __/| | | || (__ | | | |
          |_|  |_| \__,_||_| \__||_||____/  \___||_| |_| \___||_| |_|
                                   
                         
                                   Version 0

==============================================================================
Contact Information
==============================================================================

Organization:     Center for Exascale Simulation of Advanced Reactors (CESAR)
                  Argonne National Laboratory

Development Lead: John Tramm <jtramm@mcs.anl.gov>

==============================================================================
What is MultiBench?
==============================================================================

A mini-app to represent the multipole resonance representation lookup
cross section algorithm.

==============================================================================
Quick Start Guide
==============================================================================

Download----------------------------------------------------------------------

	For the most up-to-date version of MultiBench, we recommend that you
	download from our git repository. This can be accomplished via
	cloning the repository from the command line, or by downloading a zip
	from our github page. Alternatively, you can download a tar file from
	the CESAR website directly.

	Git Repository Clone:
		
		Use the following command to clone MultiBench to your machine:

		>$ git clone https://github.com/jtramm/MultiBench.git

		Once cloned, you can update the code to the newest version
		using the following command (when in the MultiBench directory):

		>$ git pull
	
	Git Zip Download:

		Simply use the "zip download" option on our webpage at:

		https://github.com/jtramm/MultiBench

Compilation-------------------------------------------------------------------

	To compile MultiBench with default settings, use the following
	command:

	>$ make

Running MultiBench---------------------------------------------------------------

	To run MultiBench with default settings, use the following command:

	>$ ./MultiBench

	For non-default settings, MultiBench supports the following command line
	options:	

	Usage: ./MultiBench <options>
	Options include:
	  -t <threads>     Number of OpenMP threads to run
	  -s <size>        Size of H-M Benchmark to run (small, large, XL, XXL)
	  -l <lookups>     Number of Cross-section (XS) lookups
	Default (no arguments given) is equivalent to: -s large -l 15000000

	-t <threads>

		Sets the number of OpenMP threads to run. By default, MultiBench
		will run with 1 thread per hardware core. If the architecture
		supports hyperthreading, multiple threads will be run per
		core.

		If running in MPI mode, this will be the number of threads
		per MPI rank.

	-s <size>
		
		Sets the size of the Hoogenboom-Martin reactor model. There
		are four options: 'small', 'large', 'XL', and 'XXL'. By default,
		the 'large' option is selected. 

		The H-M size corresponds to the number of nuclides present
		in the fuel region.  The small version has 34 fuel nuclides,
		whereas the large version has 321 fuel nuclides. This
		significantly slows down the runtime of the program as the
		data structures are much larger, and more lookups are required
		whenever a lookup occurs in a fuel material.  Note that the
		program defaults to "Large" if no specification is made.

		The additional size options, "XL" and "XXL", do not directly correspond
		to any particular physical model. They are similar to the H-M
		"large" option, except the number of gridpoints per nuclide
		has been increased greatly. This creates an extremely
		large energy grid data structure (XL: 120GB, XXL: 252GB), which is
		unlikely to fit on a single node, but is useful for experimentation
		purposes on novel architectures.

	-l <lookups>
		
		Sets the number of cross-section (XS) lookups to perform. By
		default, this value is set to 15,000,000. Users may want to
		increase this value if they wish to extend the runtime of
		MultiBench, perhaps to produce more reliable performance counter
		data - as extending the run will decrease the percentage of
		runtime spent on initialization.

==============================================================================
Debugging, Optimization & Profiling
==============================================================================

There are also a number of switches that can be set in the makefile.

Here is a sample of the control panel at the top of the makefile:

COMPILER  = gnu
OPTIMIZE  = yes
DEBUG     = no
PROFILE   = no

-> Optimization enables the -O3 optimization flag.

-> Debugging enables the -g flag.

-> Profiling enables the -pg flag.

==============================================================================
