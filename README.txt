==============================================================================
                    _____   _____ ____                  _     
                   |  __ \ / ____|  _ \                | |    
                   | |__) | (___ | |_) | ___ _ __   ___| |__  
                   |  _  / \___ \|  _ < / _ \ '_ \ / __| '_ \ 
                   | | \ \ ____) | |_) |  __/ | | | (__| | | |
                   |_|  \_\_____/|____/ \___|_| |_|\___|_| |_|
                         
                                   Version 12

==============================================================================
Contact Information
==============================================================================

Organization:             Computational Science Division
                          Argonne National Laboratory

Development Lead:         John Tramm   <jtramm@anl.gov>
                  
Contributing Developers:  Amanda Lund  <alund@anl.gov>
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

Selecting A Source Version----------------------------------------------------

	RSBench has been implemented in multiple different languages to target
	a variety of computational architectures and accelerators. The available
	implementations can be found in the "RSBench/src" directory:

	1. RSBench/src/openmp-threading

	This is the "default" version of RSBench that is appropriate for serial
	and multicore CPU architectures. The method of parallelism is via the
	OpenMP threading model.

	2. RSBench/src/openmp-offload

	This method of parallelism uses OpenMP 4.5 (or newer) to map program
	data to a remote accelerator memory space and run targeted kernels on
	the accelerator. This method of parallelism could be used for a wide
	variety of architectures (besides CPUs) that support OpenMP 4.5 targeting.

	NOTE: The Makefile will likely not work by default and will need to be
	adjusted to utilize your OpenMP accelerator compiler.

	3. RSBench/src/cuda

	This version of RSBench is written in CUDA for use with NVIDIA GPU
	architectures.

	NOTE: You will likely want to specify in the makefile the SM version
	for the card you are running on.

	4. RSBench/src/opencl

	This version of RSBench is written in OpenCL, and can be used for CPU,
	GPU, FPGA, or other architectures that support OpenCL. It was written
	with GPUs in mind, so if running on other architectures you may need to
	heavily re-optimize the code. You will also likely need to edit the
	makefile to supply the path to your OpenCL compiler.
	
	4. RSBench/src/sycl

	This version of RSBench is written in SYCL, and can be used for CPU,
	GPU, FPGA, or other architectures that support OpenCL and SYCL.
	It was written with GPUs in mind, so if running on other architectures
	you may need to heavily re-optimize the code. You will also likely need
	to edit the makefile to supply the path to your SYCL compiler.

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
	  -k <kernel>      Optimized kernel id selector. The baseline kernel "0" is default.
	Default is equivalent to: -s large -p 300000 -l 34 -P 1000 -w 100 -k 0

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
		core. This command is only available with the openmp-threading
		programming model.

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

	-k <kernel>

		For some of the XSBench code-bases (e.g., openmp-threading and cuda)
		there are several optimized variants of the main kernel. All source
		bases run basically the same "baseline" kernel as default. Optimized
		kernels can be selected at runtime with this argument. Default is
		"0" for the baseline, other variants are numbered 1, 2, ... etc. 
		People interested in implementing their own optimized variants are
		encouraged to use this interface for convenience rather than writing
		over the main kernel. The baseline kernel is defined at the top of
		the "simulation.c" source file, with the other variants being defined
		towards the end of the file after a large comment block delineation.
		The optimized variants are related to different ways of sorting
		the sampled values such that there is less thread divergence and
		much better cache re-usage when executing the lookup kernel on
		contiguous sorted elements.

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
Verification Support
==============================================================================

Legacy versions of RSBench had a special "Verification" compiler flag option
to enable verification of the results. However, a much more performant and
portable verification scheme was developed and is now used for all
configurations -- therefore, it is not necessary to compile with or without
the verification mode as it is always enabled by default.

RSBench generates a hash of the results at the end of the simulation and displays
it with the other data once the code has completed executing. This hash can
then be verified against hashes that other versions or configurations of
the code generate. For instance, running RSBench with 4 threads vs 8 threads
(on a machine that supports that configuration) should generate the
same hash number. Running on GPU vs CPU should not change the hash number.
However, changing the model / run parameters is expected to generate a totally
different hash number (i.e., increasing the number of particles, number
of gridpoints, etc, will result in different hashes). However, changing
the type of lookup performed (e.g., nuclide, unionized, or hash) should result
in the same hash being generated. Changing the simulation mode (history or
event) will generate different hashes. 

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
