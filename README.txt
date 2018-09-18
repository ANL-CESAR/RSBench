==============================================================================
                    _____   _____ ____                  _     
                   |  __ \ / ____|  _ \                | |    
                   | |__) | (___ | |_) | ___ _ __   ___| |__  
                   |  _  / \___ \|  _ < / _ \ '_ \ / __| '_ \ 
                   | | \ \ ____) | |_) |  __/ | | | (__| | | |
                   |_|  \_\_____/|____/ \___|_| |_|\___|_| |_|
                         
                                   Version 11

==============================================================================
Contact Information
==============================================================================

Organization:     Center for Exascale Simulation of Advanced Reactors (CESAR)
                  Argonne National Laboratory

Development Lead: John Tramm   <jtramm@anl.gov>
                  Amanda Lund  <alund@anl.gov>
                  Ron Rahaman  <rahaman@anl.gov>

==============================================================================
What is RSBench?
==============================================================================

RSBench is mini-app to represent the multipole resonance representation lookup
cross section algorithm for Monte Carlo neutron transport.

More information can be found in the following two publications:

https://doi.org/10.1007/978-3-319-15976-8_3

https://doi.org/10.1109/Co-HPC.2014.9

==============================================================================
Quick Start Guide
==============================================================================

Download----------------------------------------------------------------------

	For the most up-to-date version of RSBench, we recommend that you
	download from our git repository. This can be accomplished via
	cloning the repository from the command line, or by downloading a zip
	from our github page. Alternatively, you can download a tar file from
	the CESAR website directly.

	Git Repository Clone:
		
		Use the following command to clone RSBench to your machine:

		>$ git clone https://github.com/jtramm/RSBench.git

		Once cloned, you can update the code to the newest version
		using the following command (when in the RSBench directory):

		>$ git pull
	
	Git Zip Download:

		Simply use the "zip download" option on our webpage at:

		https://github.com/jtramm/RSBench

Compilation-------------------------------------------------------------------

	To compile RSBench with default settings, use the following
	command:

	>$ make

Running RSBench---------------------------------------------------------------

	To run RSBench with default settings, use the following command:

	>$ ./RSBench

	For non-default settings, RSBench supports the following command line
	options:	

	Usage: ./RSBench <options>
	Options include:
	  -m <simulation method>   Simulation method (history, event)
	  -t <threads>     Number of OpenMP threads to run
	  -s <size>        Size of H-M Benchmark to run (small, large)
	  -l <lookups>     Number of Cross-section (XS) lookups per particle
	  -p <particles>   Number of particle histories to simulate
	  -P <poles>       Average Number of Poles per Nuclide
	  -w <windows>     Average Number of Windows per Nuclide
	  -d               Disables Temperature Dependence (Doppler Broadening)
	Default is equivalent to: -s large -p 300000 -l 34 -P 1000 -w 100

	-m <simulation method>

		Sets the simulation method, either "history" or "event". These
		options represent the history based or event based algorithms
		respectively. The default is the history based method. These two
		methods represent different methods of parallelizing the Monte
		Carlo transport method. In the history based method, the central
		mode of parallelism is expressed over particles, which each require
		some number of macroscopic cross sections to be executed in series
		and in a dependent order. The event based method expresses its
		parallelism over a large pool of independent macroscopic cross
		section lookups that can be executed in any order without dependence.
		They key difference between the two methods is the dependence/independence
		of the macroscopic cross section loop.

	-t <threads>

		Sets the number of OpenMP threads to run. By default, RSBench
		will run with 1 thread per hardware core. If the architecture
		supports hyperthreading, multiple threads will be run per
		core.

		If running in MPI mode, this will be the number of threads
		per MPI rank.

	-s <size>
		
		Sets the size of the Hoogenboom-Martin reactor model. There
		are two options: 'small' and 'large'. By default,
		the 'large' option is selected. 

		The H-M size corresponds to the number of nuclides present
		in the fuel region.  The small version has 34 fuel nuclides,
		whereas the large version has 321 fuel nuclides. This
		significantly slows down the runtime of the program as the
		data structures are much larger, and more lookups are required
		whenever a lookup occurs in a fuel material.  Note that the
		program defaults to "Large" if no specification is made.

	-p <particles>

		Sets the number of particle histories to simulate.
		By default, this value is set to 300,000. Users may want to
		increase this value if they wish to extend the runtime of
		XSBench, perhaps to produce more reliable performance counter
		data - as extending the run will decrease the percentage of
		runtime spent on initialization. Real MC simulations in a full
		application may use up to several billion particles per generation,
		so there is great flexibility in this variable.

	-l <lookups>

		Sets the number of cross-section (XS) lookups to perform per particle.
		By default, this value is set to 34, which represents the average
		number of XS lookups per particle over the course of its lifetime in
		a light water reactor problem. Users should only alter this value if
		they are trying to capture the behavior of a different type of reactor
		(e.g., one with a fast spectrum), where the number of lookups per
		history may be different.

	-P <poles>
		
		This is the average number of poles per nuclide. Default is
		set to 1000.

	-w <windows>

		This is the number of windows per nuclide. Default is set to
		100.

	-d
		This flag disables Doppler broadening (temperature dependence) in the
	    calculation. Doppler broadening represents a calculation where the
		temperature of the materials in the reactor are considered and
		Doppler broadening is accomplished via the evaluation of the Faddeeva
		function for each pole within a window, which is accomplished by using
	    a fast approximation method for the Faddeeva function. Disabling
	    this feature represents a neutronics calculation done at 0K.

==============================================================================
Debugging, Optimization & Profiling
==============================================================================

There are also a number of switches that can be set in the makefile.

Here is a sample of the control panel at the top of the makefile:

COMPILER  = gnu
OPTIMIZE  = yes
DEBUG     = no
PROFILE   = no
STATUS    = yes
VERIFY    = no

-> Optimization enables the -O3 optimization flag.

-> Debugging enables the -g flag.

-> Profiling enables the -pg flag.

-> STATUS enables calculation completion % printout status text.
   You may want to disable this if doing batch or scripted runs.
   Does not affect performance.

-> VERIFY enables the verification mode of RSBench (as described below).

==============================================================================
Verification Support
==============================================================================

RSBench has the ability to verify that consistent and correct results are
achieved. This mode is enabled by altering the "VERIFY" setting to 'yes' in
the makefile, i.e.:

VERIFY = yes

Once enabled, the code will generate a hash of the results and display it
with the other data once the code has completed executing. This hash can
then be verified against hashes that other versions or configurations of
the code generate. For instance, running RSBench with 4 threads vs 8 threads
(on a machine that supports that configuration) should generate the
same hash number. Changing the model / run parameters should NOT generate
the same hash number (i.e., increasing the number of particles, number
of gridpoints, etc, will result in different hashes).

Note that the verification mode runs a little slower, due to need to hash
each macroscopic cross section result. Therefore, performance measurements
should generally not be made when verification mode is on. Rather,
verification mode should be used to ensure any changes to the code have not
broken it, and then be disabled before performance metrics are recorded.

==============================================================================
Citing RSBench
==============================================================================

Papers citing the RSBench program in general should refer to:

Tramm J.R., Siegel A.R., Forget B., Josey C, "Performance Analysis of
a Reduced Data Movement Algorithm for Neutron Cross Section Data in
Monte Carlo Simulations," Prsented at EASC 2014 - Solving Software
Challenges for Exascale, Stockholm. https://doi.org/10.1007/978-3-319-15976-8_3

Bibtex Entry:

@inproceedings{Tramm:rs,
author="Tramm, John R.
and Siegel, Andrew R.
and Forget, Benoit
and Josey, Colin",
title="Performance Analysis of a Reduced Data Movement Algorithm for Neutron Cross Section Data in Monte Carlo Simulations",
booktitle = {{EASC} 2014 - Solving Software Challenges for Exascale},
address = {Stockholm},
year = "2014",
url = "https://doi.org/10.1007/978-3-319-15976-8_3"
}

==============================================================================
