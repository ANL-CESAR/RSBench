#include "rsbench.h"

////////////////////////////////////////////////////////////////////////////////////
// BASELINE FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////
// All "baseline" code is at the top of this file. The baseline code is a simple
// implementation of the algorithm, with only minor CPU optimizations in place.
// Following these functions are a number of optimized variants,
// which each deploy a different combination of optimizations strategies. By
// default, RSBench will only run the baseline implementation. Optimized variants
// must be specifically selected using the "-k <optimized variant ID>" command
// line argument.
////////////////////////////////////////////////////////////////////////////////////

void run_event_based_simulation(Input input, SimulationData data, unsigned long * vhash_result )
{
	int i = 0;
	printf("Beginning baseline event based simulation...\n");
	unsigned long verification = 0;

	// Main simulation loop over macroscopic cross section lookups
	#pragma omp parallel for schedule(dynamic, 1000) default(none) shared(input, data) reduction(+:verification)
	for( i = 0; i < input.lookups; i++ )
	{
		// Set the initial seed value
		uint64_t seed = STARTING_SEED;	

		// Forward seed to lookup index (we need 2 samples per lookup)
		seed = fast_forward_LCG(seed, 2*i);

		// Randomly pick an energy and material for the particle
		double E = LCG_random_double(&seed);
		int mat  = pick_mat(&seed);

		double macro_xs[4] = {0};

		calculate_macro_xs( macro_xs, mat, E, input, data ); 

		// For verification, and to prevent the compiler from optimizing
		// all work out, we interrogate the returned macro_xs_vector array
		// to find its maximum value index, then increment the verification
		// value by that index. In this implementation, we prevent thread
		// contention by using an OMP reduction on it. For other accelerators,
		// a different approach might be required (e.g., atomics, reduction
		// of thread-specific values in large array via CUDA thrust, etc)
		double max = -DBL_MAX;
		int max_idx = 0;
		for(int x = 0; x < 4; x++ )
		{
			if( macro_xs[x] > max )
			{
				max = macro_xs[x];
				max_idx = x;
			}
		}
		verification += max_idx+1;
	}

	*vhash_result = verification;
}

void run_history_based_simulation(Input input, SimulationData data, unsigned long * vhash_result )
{
	int p = 0;
	printf("Beginning history based simulation...\n");
	unsigned long verification = 0;

	// Main simulation loop over particle histories
	#pragma omp parallel for schedule(dynamic, 1000) default(none) shared(input, data) reduction(+:verification)
	for( p = 0; p < input.particles; p++ )
	{
		// Set the initial seed value
		uint64_t seed = STARTING_SEED;	

		// Forward seed to lookup index (we need 2 samples per lookup)
		seed = fast_forward_LCG(seed, p * input.lookups * 2 * 4);

		// Randomly pick an energy and material for the particle
		double E = LCG_random_double(&seed);
		int mat  = pick_mat(&seed);

		// Loop over macroscopic cross section events. This loop is dependent!
		// I.e., This loop must be executed sequentially,
		// as each lookup depends on results from the previous lookup.
		for( int i = 0; i < input.lookups; i++ )
		{
			double macro_xs[4] = {0};

			calculate_macro_xs( macro_xs, mat, E, input, data ); 

			// For verification, and to prevent the compiler from optimizing
			// all work out, we interrogate the returned macro_xs_vector array
			// to find its maximum value index, then increment the verification
			// value by that index. In this implementation, we prevent thread
			// contention by using an OMP reduction on it. For other accelerators,
			// a different approach might be required (e.g., atomics, reduction
			// of thread-specific values in large array via CUDA thrust, etc)
			double max = -DBL_MAX;
			int max_idx = 0;
			for(int x = 0; x < 4; x++ )
			{
				if( macro_xs[x] > max )
				{
					max = macro_xs[x];
					max_idx = x;
				}
			}
			verification += max_idx+1;

			// Randomly pick next energy and material for the particle
			// Also incorporates results from macro_xs lookup to
			// enforce loop dependency.
			// In a real MC app, this dependency is expressed in terms
			// of branching physics sampling, whereas here we are just
			// artificially enforcing this dependence based on altering
			// the seed
			uint64_t n_forward = 0;
			for( int x = 0; x < 4; x++ )
				if( macro_xs[x] > 1.0 )
					n_forward++;
			if( n_forward > 0 )
				seed = fast_forward_LCG(seed, n_forward);

			E   = LCG_random_double(&seed);
			mat = pick_mat(&seed);
		}
	}
	*vhash_result = verification;
}


void calculate_macro_xs( double * macro_xs, int mat, double E, Input input, SimulationData data ) 
{
	// zero out macro vector
	for( int i = 0; i < 4; i++ )
		macro_xs[i] = 0;

	// for nuclide in mat
	for( int i = 0; i < data.num_nucs[mat]; i++ )
	{
		double micro_xs[4];
		int nuc = data.mats[mat * data.max_num_nucs + i];

		if( input.doppler == 1 )
			calculate_micro_xs_doppler( micro_xs, nuc, E, input, data );
		else
			calculate_micro_xs( micro_xs, nuc, E, input, data);

		for( int j = 0; j < 4; j++ )
		{
			macro_xs[j] += micro_xs[j] * data.concs[mat * data.max_num_nucs + i];
		}
		// Debug
		/*
		printf("E = %.2lf, mat = %d, macro_xs[0] = %.2lf, macro_xs[1] = %.2lf, macro_xs[2] = %.2lf, macro_xs[3] = %.2lf\n",
		E, mat, macro_xs[0], macro_xs[1], macro_xs[2], macro_xs[3] );
		*/
	}

	// Debug
	/*
	printf("E = %.2lf, mat = %d, macro_xs[0] = %.2lf, macro_xs[1] = %.2lf, macro_xs[2] = %.2lf, macro_xs[3] = %.2lf\n",
	E, mat, macro_xs[0], macro_xs[1], macro_xs[2], macro_xs[3] );
	*/
}

// No Temperature dependence (i.e., 0K evaluation)
void calculate_micro_xs( double * micro_xs, int nuc, double E, Input input, SimulationData data)
{
	// MicroScopic XS's to Calculate
	double sigT;
	double sigA;
	double sigF;
	double sigE;

	// Calculate Window Index
	double spacing = 1.0 / data.n_windows[nuc];
	int window = (int) ( E / spacing );
	if( window == data.n_windows[nuc] )
		window--;

	// Calculate sigTfactors
	RSComplex sigTfactors[4]; // Of length input.numL, which is always 4
	calculate_sig_T(nuc, E, input, data, sigTfactors );

	// Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
	Window w = data.windows[nuc * data.max_num_windows + window];
	sigT = E * w.T;
	sigA = E * w.A;
	sigF = E * w.F;

	// Loop over Poles within window, add contributions
	for( int i = w.start; i < w.end; i++ )
	{
		RSComplex PSIIKI;
		RSComplex CDUM;
		Pole pole = data.poles[nuc * data.max_num_poles + i];
		RSComplex t1 = {0, 1};
		RSComplex t2 = {sqrt(E), 0 };
		PSIIKI = c_div( t1 , c_sub(pole.MP_EA,t2) );
		RSComplex E_c = {E, 0};
		CDUM = c_div(PSIIKI, E_c);
		sigT += (c_mul(pole.MP_RT, c_mul(CDUM, sigTfactors[pole.l_value])) ).r;
		sigA += (c_mul( pole.MP_RA, CDUM)).r;
		sigF += (c_mul(pole.MP_RF, CDUM)).r;
	}

	sigE = sigT - sigA;

	micro_xs[0] = sigT;
	micro_xs[1] = sigA;
	micro_xs[2] = sigF;
	micro_xs[3] = sigE;
}

// Temperature Dependent Variation of Kernel
// (This involves using the Complex Faddeeva function to
// Doppler broaden the poles within the window)
void calculate_micro_xs_doppler( double * micro_xs, int nuc, double E, Input input, SimulationData data )
{
	// MicroScopic XS's to Calculate
	double sigT;
	double sigA;
	double sigF;
	double sigE;

	// Calculate Window Index
	double spacing = 1.0 / data.n_windows[nuc];
	int window = (int) ( E / spacing );
	if( window == data.n_windows[nuc] )
		window--;

	// Calculate sigTfactors
	RSComplex sigTfactors[4]; // Of length input.numL, which is always 4
	calculate_sig_T(nuc, E, input, data, sigTfactors );

	// Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
	Window w = data.windows[nuc * data.max_num_windows + window];
	sigT = E * w.T;
	sigA = E * w.A;
	sigF = E * w.F;

	double dopp = 0.5;

	// Loop over Poles within window, add contributions
	for( int i = w.start; i < w.end; i++ )
	{
		Pole pole = data.poles[nuc * data.max_num_poles + i];

		// Prep Z
		RSComplex E_c = {E, 0};
		RSComplex dopp_c = {dopp, 0};
		RSComplex Z = c_mul(c_sub(E_c, pole.MP_EA), dopp_c);

		// Evaluate Fadeeva Function
		RSComplex faddeeva = fast_nuclear_W( Z );

		// Update W
		sigT += (c_mul( pole.MP_RT, c_mul(faddeeva, sigTfactors[pole.l_value]) )).r;
		sigA += (c_mul( pole.MP_RA , faddeeva)).r;
		sigF += (c_mul( pole.MP_RF , faddeeva)).r;
	}

	sigE = sigT - sigA;

	micro_xs[0] = sigT;
	micro_xs[1] = sigA;
	micro_xs[2] = sigF;
	micro_xs[3] = sigE;
}

void calculate_sig_T( int nuc, double E, Input input, SimulationData data, RSComplex * sigTfactors )
{
	double phi;

	for( int i = 0; i < 4; i++ )
	{
		phi = data.pseudo_K0RS[nuc * input.numL + i] * sqrt(E);

		if( i == 1 )
			phi -= - atan( phi );
		else if( i == 2 )
			phi -= atan( 3.0 * phi / (3.0 - phi*phi));
		else if( i == 3 )
			phi -= atan(phi*(15.0-phi*phi)/(15.0-6.0*phi*phi));

		phi *= 2.0;

		sigTfactors[i].r = cos(phi);
		sigTfactors[i].i = -sin(phi);
	}
}

// This function uses a combination of the Abrarov Approximation
// and the QUICK_W three term asymptotic expansion.
// Only expected to use Abrarov ~0.5% of the time.
RSComplex fast_nuclear_W( RSComplex Z )
{
	// Abrarov 
	if( c_abs(Z) < 6.0 )
	{
		// Precomputed parts for speeding things up
		// (N = 10, Tm = 12.0)
		RSComplex prefactor = {0, 8.124330e+01};
		double an[10] = {
			2.758402e-01,
			2.245740e-01,
			1.594149e-01,
			9.866577e-02,
			5.324414e-02,
			2.505215e-02,
			1.027747e-02,
			3.676164e-03,
			1.146494e-03,
			3.117570e-04
		};
		double neg_1n[10] = {
			-1.0,
			1.0,
			-1.0,
			1.0,
			-1.0,
			1.0,
			-1.0,
			1.0,
			-1.0,
			1.0
		};

		double denominator_left[10] = {
			9.869604e+00,
			3.947842e+01,
			8.882644e+01,
			1.579137e+02,
			2.467401e+02,
			3.553058e+02,
			4.836106e+02,
			6.316547e+02,
			7.994380e+02,
			9.869604e+02
		};

		RSComplex t1 = {0, 12};
		RSComplex t2 = {12, 0};
		RSComplex i = {0,1};
		RSComplex one = {1, 0};
		RSComplex W = c_div(c_mul(i, ( c_sub(one, fast_cexp(c_mul(t1, Z))) )) , c_mul(t2, Z));
		RSComplex sum = {0,0};
		for( int n = 0; n < 10; n++ )
		{
			RSComplex t3 = {neg_1n[n], 0};
			RSComplex top = c_sub(c_mul(t3, fast_cexp(c_mul(t1, Z))), one);
			RSComplex t4 = {denominator_left[n], 0};
			RSComplex t5 = {144, 0};
			RSComplex bot = c_sub(t4, c_mul(t5,c_mul(Z,Z)));
			RSComplex t6 = {an[n], 0};
			sum = c_add(sum, c_mul(t6, c_div(top,bot)));
		}
		W = c_add(W, c_mul(prefactor, c_mul(Z, sum)));
		return W;
	}
	else
	{
		// QUICK_2 3 Term Asymptotic Expansion (Accurate to O(1e-6)).
		// Pre-computed parameters
		RSComplex a = {0.512424224754768462984202823134979415014943561548661637413182,0};
		RSComplex b = {0.275255128608410950901357962647054304017026259671664935783653, 0};
		RSComplex c = {0.051765358792987823963876628425793170829107067780337219430904, 0};
		RSComplex d = {2.724744871391589049098642037352945695982973740328335064216346, 0};

		RSComplex i = {0,1};
		RSComplex Z2 = c_mul(Z, Z);
		// Three Term Asymptotic Expansion
		RSComplex W = c_mul(c_mul(Z,i), (c_add(c_div(a,(c_sub(Z2, b))) , c_div(c,(c_sub(Z2, d))))));

		return W;
	}
}

double LCG_random_double(uint64_t * seed)
{
	const uint64_t m = 9223372036854775808ULL; // 2^63
	const uint64_t a = 2806196910506780709ULL;
	const uint64_t c = 1ULL;
	*seed = (a * (*seed) + c) % m;
	return (double) (*seed) / (double) m;
}	

uint64_t LCG_random_int(uint64_t * seed)
{
	const uint64_t m = 9223372036854775808ULL; // 2^63
	const uint64_t a = 2806196910506780709ULL;
	const uint64_t c = 1ULL;
	*seed = (a * (*seed) + c) % m;
	return *seed;
}	

uint64_t fast_forward_LCG(uint64_t seed, uint64_t n)
{
	const uint64_t m = 9223372036854775808ULL; // 2^63
	uint64_t a = 2806196910506780709ULL;
	uint64_t c = 1ULL;

	n = n % m;

	uint64_t a_new = 1;
	uint64_t c_new = 0;

	while(n > 0) 
	{
		if(n & 1)
		{
			a_new *= a;
			c_new = c_new * a + c;
		}
		c *= (a + 1);
		a *= a;

		n >>= 1;
	}

	return (a_new * seed + c_new) % m;
}

// Complex arithmetic functions

RSComplex c_add( RSComplex A, RSComplex B)
{
	RSComplex C;
	C.r = A.r + B.r;
	C.i = A.i + B.i;
	return C;
}

RSComplex c_sub( RSComplex A, RSComplex B)
{
	RSComplex C;
	C.r = A.r - B.r;
	C.i = A.i - B.i;
	return C;
}

RSComplex c_mul( RSComplex A, RSComplex B)
{
	double a = A.r;
	double b = A.i;
	double c = B.r;
	double d = B.i;
	RSComplex C;
	C.r = (a*c) - (b*d);
	C.i = (a*d) + (b*c);
	return C;
}

RSComplex c_div( RSComplex A, RSComplex B)
{
	double a = A.r;
	double b = A.i;
	double c = B.r;
	double d = B.i;
	RSComplex C;
	double denom = c*c + d*d;
	C.r = ( (a*c) + (b*d) ) / denom;
	C.i = ( (b*c) - (a*d) ) / denom;
	return C;
}

double c_abs( RSComplex A)
{
	return sqrt(A.r*A.r + A.i * A.i);
}

// Fast (but inaccurate) exponential function
// Written By "ACMer":
// https://codingforspeed.com/using-faster-exponential-approximation/
// We use our own to avoid small differences in compiler specific
// exp() intrinsic implementations that make it difficult to verify
// if the code is working correctly or not.
double fast_exp(double x)
{
  x = 1.0 + x * 0.000244140625;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  return x;
}

// Implementation based on:
// z = x + iy
// cexp(z) = e^x * (cos(y) + i * sin(y))
RSComplex fast_cexp( RSComplex z )
{
	double x = z.r;
	double y = z.i;

	// For consistency across architectures, we
	// will use our own exponetial implementation
	//double t1 = exp(x);
	double t1 = fast_exp(x);
	double t2 = cos(y);
	double t3 = sin(y);
	RSComplex t4 = {t2, t3};
	RSComplex t5 = {t1, 0};
	RSComplex result = c_mul(t5, (t4));
	return result;
}	

////////////////////////////////////////////////////////////////////////////////////
// Parallel Quicksort Key-Value Sorting Algorithms
////////////////////////////////////////////////////////////////////////////////////
//
// These algorithms are based on the parallel quicksort implementation by
// Eduard Lopez published at https://github.com/eduardlopez/quicksort-parallel
//
// Eduard's original version was for an integer type quicksort, but I have modified
// it to form two different versions that can sort key-value pairs together without
// having to bundle them into a separate object. Additionally, I have modified the
// optimal chunk sizes and restricted the number of threads for the array sizing
// that XSBench will be using by default.
//
// Eduard's original implementation carries the following license, which applies to
// the following functions only:
//
//	void quickSort_parallel_internal_i_d(int* key,double * value, int left, int right, int cutoff) 
//  void quickSort_parallel_i_d(int* key,double * value, int lenArray, int numThreads)
//  void quickSort_parallel_internal_d_i(double* key,int * value, int left, int right, int cutoff)
//  void quickSort_parallel_d_i(double* key,int * value, int lenArray, int numThreads)
//
// The MIT License (MIT)
//
// Copyright (c) 2016 Eduard López
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////////
void quickSort_parallel_internal_i_d(int* key,double * value, int left, int right, int cutoff) 
{
	int i = left, j = right;
	int tmp;
	int pivot = key[(left + right) / 2];
	
	{
		while (i <= j) {
			while (key[i] < pivot)
				i++;
			while (key[j] > pivot)
				j--;
			if (i <= j) {
				tmp = key[i];
				key[i] = key[j];
				key[j] = tmp;
				double tmp_v = value[i];
				value[i] = value[j];
				value[j] = tmp_v;
				i++;
				j--;
			}
		}

	}

	if ( ((right-left)<cutoff) ){
		if (left < j){ quickSort_parallel_internal_i_d(key, value, left, j, cutoff); }			
		if (i < right){ quickSort_parallel_internal_i_d(key, value, i, right, cutoff); }

	}else{
		#pragma omp task 	
		{ quickSort_parallel_internal_i_d(key, value, left, j, cutoff); }
		#pragma omp task 	
		{ quickSort_parallel_internal_i_d(key, value, i, right, cutoff); }		
	}

}

void quickSort_parallel_i_d(int* key,double * value, int lenArray, int numThreads){

	// Set minumum problem size to still spawn threads for
	int cutoff = 10000;

	// For this problem size, more than 16 threads on CPU is not helpful
	if( numThreads > 16 )
		numThreads = 16;

	#pragma omp parallel num_threads(numThreads)
	{	
		#pragma omp single nowait
		{
			quickSort_parallel_internal_i_d(key,value, 0, lenArray-1, cutoff);	
		}
	}	

}

void quickSort_parallel_internal_d_i(double* key,int * value, int left, int right, int cutoff) 
{
	int i = left, j = right;
	double tmp;
	double pivot = key[(left + right) / 2];
	
	{
		while (i <= j) {
			while (key[i] < pivot)
				i++;
			while (key[j] > pivot)
				j--;
			if (i <= j) {
				tmp = key[i];
				key[i] = key[j];
				key[j] = tmp;
				int tmp_v = value[i];
				value[i] = value[j];
				value[j] = tmp_v;
				i++;
				j--;
			}
		}

	}

	if ( ((right-left)<cutoff) ){
		if (left < j){ quickSort_parallel_internal_d_i(key, value, left, j, cutoff); }			
		if (i < right){ quickSort_parallel_internal_d_i(key, value, i, right, cutoff); }

	}else{
		#pragma omp task 	
		{ quickSort_parallel_internal_d_i(key, value, left, j, cutoff); }
		#pragma omp task 	
		{ quickSort_parallel_internal_d_i(key, value, i, right, cutoff); }		
	}

}

void quickSort_parallel_d_i(double* key,int * value, int lenArray, int numThreads){

	// Set minumum problem size to still spawn threads for
	int cutoff = 10000;

	// For this problem size, more than 16 threads on CPU is not helpful
	if( numThreads > 16 )
		numThreads = 16;

	#pragma omp parallel num_threads(numThreads)
	{	
		#pragma omp single nowait
		{
			quickSort_parallel_internal_d_i(key,value, 0, lenArray-1, cutoff);	
		}
	}	

}

////////////////////////////////////////////////////////////////////////////////////
// Optimization 1 -- Event-based Sample/XS Lookup kernel splitting + Sorting
//                   lookups by material and energy
////////////////////////////////////////////////////////////////////////////////////
// This kernel separates out the sampling and lookup regions of the event-based
// model, and then sorts the lookups by material type and energy. The goal of this
// optimization is to allow for greatly improved cache locality, and XS indices
// loaded from memory may be re-used for multiple lookups.
//
// As efficienct sorting is key for performance, we also must implement an
// efficient key-value parallel sorting algorithm. We also experimented with using
// the C++ version of thrust for these purposes, but found that our own implemtation
// was slightly faster than the thrust library version, so for speed and
// simplicity we will do not add the thrust dependency.
////////////////////////////////////////////////////////////////////////////////////

void run_event_based_simulation_optimization_1(Input in, SimulationData SD, unsigned long * vhash_result )
{
	char * optimization_name = "Optimization 1 - Kernel splitting + full material & energy sort";
	
	printf("Simulation Kernel:\"%s\"\n", optimization_name);
	
	////////////////////////////////////////////////////////////////////////////////
	// Allocate Additional Data Structures Needed by Optimized Kernel
	////////////////////////////////////////////////////////////////////////////////
	printf("Allocating additional data required by optimized kernel...\n");
	size_t sz;
	size_t total_sz = 0;
	double start, stop;
	int i, m;

	sz = in.lookups * sizeof(double);
	SD.p_energy_samples = (double *) malloc(sz);
	total_sz += sz;
	SD.length_p_energy_samples = in.lookups;

	sz = in.lookups * sizeof(int);
	SD.mat_samples = (int *) malloc(sz);
	total_sz += sz;
	SD.length_mat_samples = in.lookups;
	
	printf("Allocated an additional %.0lf MB of data on CPU.\n", total_sz/1024.0/1024.0);
	
	////////////////////////////////////////////////////////////////////////////////
	// Begin Actual Simulation 
	////////////////////////////////////////////////////////////////////////////////
	
	////////////////////////////////////////////////////////////////////////////////
	// Sample Materials and Energies
	////////////////////////////////////////////////////////////////////////////////
	printf("Sampling event data...\n");
	#pragma omp parallel for schedule(dynamic, 1000)
	for( i = 0; i < in.lookups; i++ )
	{
		// Set the initial seed value
		uint64_t seed = STARTING_SEED;	

		// Forward seed to lookup index (we need 2 samples per lookup)
		seed = fast_forward_LCG(seed, 2*i);

		// Randomly pick an energy and material for the particle
		double p_energy = LCG_random_double(&seed);
		int mat         = pick_mat(&seed); 

		SD.p_energy_samples[i] = p_energy;
		SD.mat_samples[i] = mat;
	}
	printf("Finished sampling.\n");
	
	////////////////////////////////////////////////////////////////////////////////
	// Sort by Material
	////////////////////////////////////////////////////////////////////////////////
	
	start = get_time();

	quickSort_parallel_i_d(SD.mat_samples, SD.p_energy_samples, in.lookups, in.nthreads);

	stop = get_time();

	printf("Material sort took %.3lf seconds\n", stop-start);
	
	////////////////////////////////////////////////////////////////////////////////
	// Sort by Energy
	////////////////////////////////////////////////////////////////////////////////
	
	start = get_time();
	
	// Count up number of each type of sample. 
	int num_samples_per_mat[12] = {0};
	for( int l = 0; l < in.lookups; l++ )
		num_samples_per_mat[ SD.mat_samples[l] ]++;

	// Determine offsets
	int offsets[12] = {0};
	for( int m = 1; m < 12; m++ )
		offsets[m] = offsets[m-1] + num_samples_per_mat[m-1];
	
	stop = get_time();
	printf("Counting samples and offsets took %.3lf seconds\n", stop-start);
	start = stop;

	// Sort each material type by energy level
	int offset = 0;
	for( int m = 0; m < 12; m++ )
		quickSort_parallel_d_i(SD.p_energy_samples + offsets[m],SD.mat_samples + offsets[m], num_samples_per_mat[m], in.nthreads);

	stop = get_time();
	printf("Energy Sorts took %.3lf seconds\n", stop-start);
	
	////////////////////////////////////////////////////////////////////////////////
	// Perform lookups for each material separately
	////////////////////////////////////////////////////////////////////////////////
	start = get_time();

	unsigned long long verification = 0;

	// Individual Materials
	offset = 0;
	for( m = 0; m < 12; m++ )
	{
		#pragma omp parallel for schedule(dynamic,100) reduction(+:verification)
		for( i = offset; i < offset + num_samples_per_mat[m]; i++)
		{
			// load pre-sampled energy and material for the particle
			double E = SD.p_energy_samples[i];
			int mat  = SD.mat_samples[i]; 

			double macro_xs_vector[4] = {0};

			// Perform macroscopic Cross Section Lookup
			calculate_macro_xs( macro_xs_vector, mat, E, in, SD ); 

			// For verification, and to prevent the compiler from optimizing
			// all work out, we interrogate the returned macro_xs_vector array
			// to find its maximum value index, then increment the verification
			// value by that index. In this implementation, we prevent thread
			// contention by using an OMP reduction on the verification value.
			// For accelerators, a different approach might be required
			// (e.g., atomics, reduction of thread-specific values in large
			// array via CUDA thrust, etc).
			double max = -DBL_MAX;
			int max_idx = 0;
			for(int j = 0; j < 4; j++ )
			{
				if( macro_xs_vector[j] > max )
				{
					max = macro_xs_vector[j];
					max_idx = j;
				}
			}
			verification += max_idx+1;
		}
		offset += num_samples_per_mat[m];
	}
	
	stop = get_time();
	printf("XS Lookups took %.3lf seconds\n", stop-start);
	*vhash_result = verification;
}
