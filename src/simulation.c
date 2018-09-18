#include "rsbench.h"
void run_event_based_simulation(Input input, CalcDataPtrs data, long * abrarov_result, long * alls_result, unsigned long * vhash_result )
{
	printf("Beginning event based simulation...\n");
	long g_abrarov = 0;
	long g_alls = 0;
	unsigned long vhash = 0;
	#pragma omp parallel default(none) \
	shared(input, data) \
	reduction(+:g_abrarov, g_alls, vhash)
	{
		double * xs = (double *) calloc(4, sizeof(double));
		int thread = omp_get_thread_num();
		long abrarov = 0; 
		long alls = 0;

		#ifdef PAPI
		int eventset = PAPI_NULL; 
		int num_papi_events;
		#pragma omp critical
		{
			counter_init(&eventset, &num_papi_events);
		}
		#endif
		complex double * sigTfactors =
			(complex double *) malloc( input.numL * sizeof(complex double) );

		// This loop is independent!
		// I.e., macroscopic cross section lookups in the event based simulation
		// can be executed in any order.
		#pragma omp for schedule(guided)
		for( int i = 0; i < input.lookups; i++ )
		{
			// Particles are seeded by their particle ID
			unsigned long seed = ((unsigned long) i+ (unsigned long)1)* (unsigned long) 13371337;

			// Randomly pick an energy and material for the particle
			double E = rn(&seed);
			int mat  = pick_mat(&seed);

			#ifdef STATUS
			if( thread == 0 && i % 2000 == 0 )
				printf("\rCalculating XS's... (%.0lf%% completed)",
						(i / ( (double)input.lookups /
							   (double) input.nthreads )) /
						(double) input.nthreads * 100.0);
			#endif

			double macro_xs[4] = {0};

			calculate_macro_xs( macro_xs, mat, E, input, data, sigTfactors, &abrarov, &alls ); 

			// Results are copied onto heap to avoid some compiler
			// flags (-flto) from optimizing out function call
			memcpy(xs, macro_xs, 4*sizeof(double));

			// Verification hash calculation
			// This method provides a consistent hash accross
			// architectures and compilers.
			#ifdef VERIFICATION
			char line[256];
			sprintf(line, "%.2le %d %.2le %.2le %.2le %.2le",
				   E, mat,
				   macro_xs[0],
				   macro_xs[1],
				   macro_xs[2],
				   macro_xs[3]);
			unsigned long long vhash_local = hash(line, 10000);

			vhash += vhash_local;
			#endif
		}

		free(sigTfactors);

		// Accumulate global counters
		g_abrarov = abrarov; 
		g_alls = alls;

		#ifdef PAPI
		if( thread == 0 )
		{
			printf("\n");
			border_print();
			center_print("PAPI COUNTER RESULTS", 79);
			border_print();
			printf("Count          \tSmybol      \tDescription\n");
		}
		{
			#pragma omp barrier
		}
		counter_stop(&eventset, num_papi_events);
		#endif
	}
	*abrarov_result = g_abrarov;
	*alls_result = g_alls;
	*vhash_result = vhash;
}

void run_history_based_simulation(Input input, CalcDataPtrs data, long * abrarov_result, long * alls_result, unsigned long * vhash_result )
{
	printf("Beginning history based simulation...\n");
	long g_abrarov = 0;
	long g_alls = 0;
	unsigned long vhash = 0;
	#pragma omp parallel default(none) \
	shared(input, data) \
	reduction(+:g_abrarov, g_alls, vhash)
	{
		double * xs = (double *) calloc(4, sizeof(double));
		int thread = omp_get_thread_num();
		long abrarov = 0; 
		long alls = 0;

		#ifdef PAPI
		int eventset = PAPI_NULL; 
		int num_papi_events;
		#pragma omp critical
		{
			counter_init(&eventset, &num_papi_events);
		}
		#endif
		complex double * sigTfactors =
			(complex double *) malloc( input.numL * sizeof(complex double) );

		// This loop is independent!
		// I.e., particle histories can be executed in any order
		#pragma omp for schedule(guided)
		for( int p = 0; p < input.particles; p++ )
		{
			// Particles are seeded by their particle ID
            unsigned long seed = ((unsigned long) p+ (unsigned long)1)* (unsigned long) 13371337;

			// Randomly pick an energy and material for the particle
            double E = rn(&seed);
            int mat  = pick_mat(&seed);

			#ifdef STATUS
			if( thread == 0 && p % 35 == 0 )
				printf("\rCalculating XS's... (%.0lf%% completed)",
						(p / ( (double)input.particles /
							   (double) input.nthreads )) /
						(double) input.nthreads * 100.0);
			#endif

			// This loop is dependent!
			// I.e., This loop must be executed sequentially,
			// as each lookup depends on results from the previous lookup.
			for( int i = 0; i < input.lookups; i++ )
			{
				double macro_xs[4] = {0};

				calculate_macro_xs( macro_xs, mat, E, input, data, sigTfactors, &abrarov, &alls ); 

				// Results are copied onto heap to avoid some compiler
				// flags (-flto) from optimizing out function call
				memcpy(xs, macro_xs, 4*sizeof(double));

				// Verification hash calculation
                // This method provides a consistent hash accross
                // architectures and compilers.
                #ifdef VERIFICATION
                char line[256];
                sprintf(line, "%.2le %d %.2le %.2le %.2le %.2le",
                       E, mat,
                       macro_xs[0],
                       macro_xs[1],
                       macro_xs[2],
                       macro_xs[3]);
                unsigned long long vhash_local = hash(line, 10000);

                vhash += vhash_local;
                #endif

                // Randomly pick next energy and material for the particle
                // Also incorporates results from macro_xs lookup to
                // enforce loop dependency.
                // In a real MC app, this dependency is expressed in terms
                // of branching physics sampling, whereas here we are just
                // artificially enforcing this dependence based on altering
                // the seed
                for( int x = 0; x < 4; x++ )
				{
					if( macro_xs[x] > 0 )
                    	seed += 1337*p;
					else
						seed += 42;
				}

                E   = rn(&seed);
                mat = pick_mat(&seed);
			}
		}

		free(sigTfactors);

		// Accumulate global counters
		g_abrarov = abrarov; 
		g_alls = alls;

		#ifdef PAPI
		if( thread == 0 )
		{
			printf("\n");
			border_print();
			center_print("PAPI COUNTER RESULTS", 79);
			border_print();
			printf("Count          \tSmybol      \tDescription\n");
		}
		{
			#pragma omp barrier
		}
		counter_stop(&eventset, num_papi_events);
		#endif
	}
	*abrarov_result = g_abrarov;
	*alls_result = g_alls;
	*vhash_result = vhash;
}
