![RSBench](docs/img/logo.png)

[![Latest Github release](https://img.shields.io/github/release/ANL-CESAR/RSBench.svg)](https://github.com/ANL-CESAR/RSBench/releases/latest)
[![Build Status](https://travis-ci.com/ANL-CESAR/RSBench.svg?branch=master)](https://travis-ci.com/ANL-CESAR/RSBench)
[![Published in Proceedings of EASC 2014](https://img.shields.io/badge/Published%20in-Proceedings%20of%20EASC%202014-167DA4.svg)](https://doi.org/10.1007/978-3-319-15976-8_3)

RSBench is a mini-app representing a key computational kernel of the Monte Carlo neutron transport algorithm. Specifically, RSBench represents the multipole method of performing continuous energy macroscopic neutron cross section lookups. The mulitpole method is a recently developed strategy for building microscopic cross section data "on-the-fly" that requires orders of magnitude less memory storage as compared to traditional methods (e.g., those represented in [XSBench](https://github.com/ANL-CESAR/XSBench)). RSBench serves as a useful performance stand-in for full neutron transport applications like [OpenMC](https://github.com/openmc-dev/openmc) that support multipole cross section representations.

## Table of Contents

1. [Selecting a Programming Model](#selecting-a-programming-model)
2. [Compilation](#Compilation)
3. [Running RSBench / Command Line Interface](#Running-RSBench)
4. [Verification Support](#Verification-Support)
5. [Theory & Algorithms](#Theory-&-Algorithms)
	* [Transport Simulation Styles](#Transport-Simulation-Styles)
		- [History-Based Transport](#History-Based-Transport)
		- [Event-Based Transport](#Event-Based-Transport)
	* [The Multiple Cross Section Lookup Method](#The-Multipole-Cross-Section-Lookup-Method)
		* [Faddeeva Function Evaluation](#Faddeeva-Function-Evaluation)
6. [Optimized Kernels](#Optimized-Kernels)
7. [Citing RSBench](#Citing-RSBench)
8. [Development Team](#Development-Team)

## Selecting A Programming Model

RSBench has been implemented in multiple different languages to target a variety of computational architectures and accelerators. The available implementations can be found in their own directories:

1. **RSBench/openmp-threading**
This is the "default" version of RSBench that is appropriate for serial and multicore CPU architectures. The method of parallelism is via the OpenMP threading model.

2. **RSBench/openmp-offload**
This method of parallelism uses OpenMP 4.5 (or newer) to map program data to a remote accelerator memory space and run targeted kernels on the accelerator. This method of parallelism could be used for a wide variety of architectures (besides CPUs) that support OpenMP 4.5 targeting. NOTE: The Makefile will likely not work by default and will need to be adjusted to utilize your OpenMP accelerator compiler.

3. **RSBench/cuda**
This version of RSBench is written in CUDA for use with NVIDIA GPU architectures. NOTE: You will likely want to specify in the makefile the SM version for the card you are running on.

4. **RSBench/opencl**
This version of RSBench is written in OpenCL, and can be used for CPU, GPU, FPGA, or other architectures that support OpenCL. It was written with GPUs in mind, so if running on other architectures you may need to heavily re-optimize the code. You will also likely need to edit the makefile to supply the path to your OpenCL compiler.

4. **RSBench/sycl**
This version of RSBench is written in SYCL, and can be used for CPU, GPU, FPGA, or other architectures that support OpenCL and SYCL. It was written with GPUs in mind, so if running on other architectures you may need to heavily re-optimize the code. You will also likely need to edit the makefile to supply the path to your SYCL compiler.

5. **RSBench/hip**
This version of RSBench is written in HIP for use with GPU architectures. This version is derived from CUDA using an automatic conversion tool with only a few small manual changes.

## Compilation

To compile RSBench with default settings, navigate to your selected source directory and use the following command:

```bash
make
```

 You can alter compiler settings in the included Makefile. Alternatively, for the OpenMP threading version of RSBench you may specify a compiler via the CC environment variable and then making as normal, e.g.:
```bash
export CC=clang
make
```

### Debugging, Optimization & Profiling

There are also a number of switches that can be set in the makefile. Here is a sample of the control panel at the top of the makefile:

```make
OPTIMIZE = yes
DEBUG    = no
PROFILE  = no
```
- Optimization enables the -O3 optimization flag.
- Debugging enables the -g flag.
- Profiling enables the -pg flag. When profiling the code, you may
wish to significantly increase the number of lookups (with the -l
flag) in order to wash out the initialization phase of the code.

## Running RSBench

To run RSBench with default settings, use the following command:
```bash
./RSBench
```
For non-default settings, RSBench supports the following command line options:

| Argument    |Description | Options     | Default
|-------------|------------|---------------|------------|
|-t           |# of OpenMP threads| integer value | System Default|
|-m           |Simulation method| history, event| history|
|-s | Problem Size | small, large | large|
-p | # of particle histories (if running using "history" method)| integer value | 500,000 |
-l | # of Cross-section (XS) lookups. If using using history based method, this is lookups per particle history. If using event-based method, this is total lookups. | integer value | (History: 34) (Event: 17,000,000) |
-p | # of avg poles per nuclide | integer value | 1,000 |
-w | # of windows per nuclide | integer value | 100 |
-d | Flag to disable Doppler broadening | ||
-k | Optimized kernel ID | integer value | 0

- **-t [threads]**
Sets the number of OpenMP threads to run. By default, RSBench will run with 1 thread per hardware core. If the architecture supports hyperthreading, multiple threads will be run per core. If running in MPI mode, this will be the number of threads per MPI rank. This argument is only used by the OpenMP threading version of RSBench.

- **-m [simulation method]**
Sets the simulation method, either "history" or "event". These options represent the history based or event based algorithms respectively. The default is the history based method. These two methods represent different methods of parallelizing the Monte Carlo transport method. In the history based method, the central mode of parallelism is expressed over particles, which each require some number of macroscopic cross sections to be executed in series and in a dependent order. The event based method expresses its parallelism over a large pool of independent macroscopic cross section lookups that can be executed in any order without dependence. They key difference between the two methods is the dependence/independence of the macroscopic cross section loop. See the [Transport Simulation Styles](#Transport-Simulation-Styles) section for more information.

- **-s [size]**
Sets the size of the Hoogenboom-Martin reactor model. There are four options: 'small', 'large', 'XL', and 'XXL'. By default, the 'large' option is selected. The H-M size corresponds to the number of nuclides present in the fuel region. The small version has 34 fuel nuclides, whereas the large version has 321 fuel nuclides. This significantly slows down the runtime of the program as the data structures are much larger, and more lookups are required whenever a lookup occurs in a fuel material. Note that the program defaults to "Large" if no specification is made. The additional size options, "XL" and "XXL", do not directly correspond to any particular physical model. They are similar to the H-M "large" option, except the number of gridpoints per nuclide has been increased greatly. This creates an extremely large energy grid data structure (XL: 120GB, XXL: 252GB), which is unlikely to fit on a single node, but is useful for experimentation purposes on novel architectures.

- **-p [particles]**
Sets the number of particle histories to simulate. By default, this value is set to 500,000. Users may want to increase this value if they wish to extend the runtime of RSBench, perhaps to produce more reliable performance counter data - as extending the run will decrease the percentage of runtime spent on initialization. Real MC simulations in a full application may use up to several billion particles per generation, so there is great flexibility in this variable.

- **-l [lookups]**
Sets the number of cross-section (XS) lookups to perform per particle. By default, this value is set to 34, which represents the average number of XS lookups per particle over the course of its lifetime in a light water reactor problem. Users should only alter this value if they are trying to capture the behavior of a different type of reactor (e.g., one with a fast spectrum), where the number of lookups per history may be different.

- **-p [poles]**
Sets the average number of poles per nuclide.

- **-w [windows]**
Sets the number of windows per nuclide. 

- **-d**
This flag disables Doppler broadening, in effect making it single temperature 0K simulation. Doppler broadening is enabled by default. While Doppler broadening is likely to be used in most cases in a full application like OpenMC, in may be useful to disable in some cases so as to compare to  traditional single temperature lookup algorithms.

- **-k [kernel]**
For some of the RSBench code-bases (e.g., openmp-threading and cuda) there are several optimized variants of the main kernel. All source bases run basically the same "baseline" kernel as default. Optimized kernels can be selected at runtime with this argument. Default is "0" for the baseline, other variants are numbered 1, 2, ... etc. People interested in implementing their own optimized variants are encouraged to use this interface for convenience rather than writing over the main kernel. The baseline kernel is defined at the top of the "Simulation.c" source file, with the other variants being defined towards the end of the file after a large comment block delineation. The optimized variants are related to different ways of sorting the sampled values such that there is less thread divergence and much better cache re-usage when executing the lookup kernel on contiguous sorted elements. More details can be found in the [Optimized Kernels](#Optimized-Kernels) section.

## Verification Support

Legacy versions of RSBench had a special "Verification" compiler flag option to enable verification of the results. However, a much more performant and portable verification scheme was developed and is now used for all configurations -- therefore, it is not necessary to compile with or without the verification mode as it is always enabled by default. RSBench generates a hash of the results at the end of the simulation and displays it with the other data once the code has completed executing. This hash can then be verified against hashes that other versions or configurations of the code generate. For instance, running RSBench with 4 threads vs 8 threads (on a machine that supports that configuration) should generate the same hash number. Running on GPU vs CPU should not change the hash number. However, changing the model / run parameters is expected to generate a totally different hash number (i.e., increasing the number of particles, number of gridpoints, etc, will result in different hashes). However, changing the type of lookup performed (e.g., nuclide, unionized, or hash) should result in the same hash being generated. Changing the simulation mode (history or event) will generate different hashes.

## Theory & Algorithms

### Transport Simulation Styles

#### History-Based Transport

The default simulation model used in RSBench is the "history-based" model. In this model, parallelism is expressed over independent particle histories, with each particle being simulated in a serial fashion from birth to death:

```c
for each particle do		   // Independent
	while particle is alive do // Dependent
		Move particle to collision site
		Process particle collision
```

This method of parallelism is very memory efficient, as the total number of particles that must be kept in memory at once is equivalent to the total number of active threads being run in the simulation. However, as there are many different types of collision events, the history-based model means that there is no natural SIMD style parallelism available for work happening between different threads.

#### Event-Based Transport

An alternative simulation model is the "event-based" model. In this model, parallelism is instead expressed over different collision (or "event") types. To facilitate this, all particles in the simulation are stored in memory at once. Each event kernel is executed in parallel on vectors of particles that currently require that event to be executed:

```c
Get vector of source particles
while any particles are alive do	     // Dependent
	for each living particle do          // Independent
		Move particle to collision site
	for each living particle do          // Independent
		Process particle collision
	Sort/consolidate surviving particles
```
This method of parallelism is requires more memory and requires an extra stream compaction kernel to sort and organize the particles periodically to ready them for the different event kernels. The benefit of this model is that kernels can potentially be execute in a SIMD manner and with higher cache efficiency due to the potential to sort particles by material and energy. On CPU architectures, the costs of sorting and buffering particles typically outweigh the benefits of the event-based model, but on accelerator architectures the tradeoff has been found to usually be more favorable.

### The Multipole Cross Section (XS) Lookup Method

RSBench represents the multipole macroscopic cross section lookup kernel. This kernel is responsible for adding together microscopic cross section data from all nuclides present in the material the neutron is travelling through, given a certain energy:

<p align="center"> <img src="docs/img/XS_equation.svg" alt="XS_Lookup_EQ" width="500"/> </p>

Macroscopic cross section data is typically required for multiple reaction channels "c", such as the total cross section, fission cross section, etc.

Historically, cross section data has been stored in pointwise format, sometimes requiring in excess of 100,000 energy level data points be stored for a single nuclide. There are a variety of methods for performing cross section lookups on traditional pointwise data, as represented in the mini-app  [XSBench](https://github.com/ANL-CESAR/XSBench). However, a more memory and bandwidth efficient cross section representation method has recently been developed known as the "multipole" format that models the quantum mechanical resonances (or "poles") that underly the pointwise data. By this method, the resonances can be modeled mathematically and assembled on-the-fly while storing only a fraction of the data that is required for the traditional pointwise format. The tradeoff is that a greatly increased amount of floating point work must be performed when expanding the quantum mechanical "residues" into useable cross section data.

More information regarding the mathematics and equations used in the multipole method can be found in:

>C. Josey, P. Ducru, B. Forget, K. Smith, Windowed multipole for cross section Doppler broadening, Journal of Computational Physics, Volume 307, 2016, Pages 715-727. https://doi.org/10.1016/j.jcp.2015.08.013

#### Faddeeva Function Evaluation

Doppler broadening of resonance data requires evaluation of the complex error function, also known as the Faddeeva Functon. This function is not typically available as a language supplied or standard library intrinsic, so use of the multipole method requires implementing our own or adding a library dependency. Typical libraries that implement the Faddeeva function (e.g., the [MIT Faddeeva Package](http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package) written by Steven Johnson), break up the phase space of the function into many different areas, with each area using its own evaluation technique.  This minimizes the number of floating point operations that must be performed, but can create a lot of branching which often precludes high SIMD efficiency. An alternative lightweight formulation is available, known as the [Fast Nuclear Faddeeva (FNF) Function](https://github.com/jtramm/FNF). FNF is a much simpler source implementation that only has one possible branch, allowing it to be sepcialized for use in complex phase space commonly seen in light water reactor simulations. FNF is therefore used in RSBench to minimize code complexity while providing acceptable accuracy for the use case of neutron transport.

## Optimized Kernels

If using the event-based model, we will be executing the lookup kernel in RSBench across all particles at once. While SIMD execution is possible using this method, typically issues can arrise that greatly reduce SIMD efficiency. In particular, different materials in the simulation have very different numbers of nuclides in them. For instance, spent fuel has 300+ nuclides, while moderator regions only have 10 or so nuclides. This creates a significant load imbalance across lanes in a SIMD vector, as some particles may only need a few iterations to complete all nuclides while others would need hundreds. Therefore, efficient SIMD execution of the event-based model is not possible without some optimizations.

One promising optimization for the event-based model is to perform a key-value sort of particles: first by material, and then by energy within each material. The first sort by material allows for adjacent particles in the vector to typically reside in the same type of material -- meaning that they will require the same number of nuclide lookup iterations. Then, the energy sort means that adjacent particles in the vector will be located close in energy space -- potentially allowing for many adjacent particles to access the same energy indices in each nuclide and therefore perform many or all of the same branching operations and read the same cache lines into memory at the same time. Once sorted, separate event kernels are then called for each material in the simulation. These two sorts can potentially boost both SIMD efficieny and cache efficiency, with effects being amplified as more particles are simulated at each event stage. The downside to this optimization is the introduction of the key-value particle sorting operations, which can be costly and potentially outweight any gains due to improved SIMD efficiency and cache performance.

We have implemented this optimization in both the OpenMP threading and CUDA models. They are not enabled by default, but must be enabled with the "-k 1" or "-k 6" flags if running with OpenMP and CUDA respectively. These optimizations have not yet been implemented in the other programming models due to the lack of an efficient parallel sorting function being easily available without having to create an external library dependency.

## Citing RSBench

Papers citing the RSBench program in general should refer to:

>Tramm J.R., Siegel A.R., Forget B., Josey C, "Performance Analysis of a Reduced Data Movement Algorithm for Neutron Cross Section Data in Monte Carlo Simulations," Presented at EASC 2014 - Solving Software Challenges for Exascale, Stockholm. [https://doi.org/10.1007/978-3-319-15976-8_3](https://doi.org/10.1007/978-3-319-15976-8_3)

Bibtex Entry:

```bibtex
@inproceedings{Tramm:rs,
author="Tramm, John R. and Siegel, Andrew R. and Forget, Benoit and Josey, Colin",
title="Performance Analysis of a Reduced Data Movement Algorithm for Neutron Cross Section Data in Monte Carlo Simulations",
booktitle = {{EASC} 2014 - Solving Software Challenges for Exascale},
address = {Stockholm},
year = "2014",
url = "https://doi.org/10.1007/978-3-319-15976-8_3"
}
```

## Development Team
Authored and maintained by John Tramm ([@jtramm](https://github.com/jtramm)) with help from Ron Rahaman, Amanda Lund, and other [contributors](https://github.com/ANL-CESAR/RSBench/graphs/contributors).
