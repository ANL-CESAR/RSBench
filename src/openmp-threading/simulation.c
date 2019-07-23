#include "rsbench.h"
void run_event_based_simulation(Input input, SimulationData data, unsigned long * vhash_result )
{
	printf("Beginning event based simulation...\n");
	unsigned long verification = 0;

	// Main simulation loop over macroscopic cross section lookups
	#pragma omp parallel for schedule(dynamic, 1000) default(none) shared(input, data) reduction(+:verification)
	for( int i = 0; i < input.lookups; i++ )
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
	printf("Beginning history based simulation...\n");
	unsigned long verification = 0;

	// Main simulation loop over particle histories
	#pragma omp parallel for schedule(dynamic, 1000) default(none) shared(input, data) reduction(+:verification)
	for( int p = 0; p < input.particles; p++ )
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
