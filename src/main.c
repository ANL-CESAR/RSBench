#include "rsbench.h"

int main(int argc, char * argv[])
{
	// =====================================================================
	// Initialization & Command Line Read-In
	// =====================================================================

	int version = 10;
	double start, stop;
	unsigned long iseed = 0;
	#ifdef VERIFICATION
	iseed = 42;
	#else
	iseed = time(NULL);
	#endif
	
	// Process CLI Fields
	Input input = read_CLI( argc, argv );

	// Set number of OpenMP Threads
	omp_set_num_threads(input.nthreads); 
	
	// =====================================================================
	// Print-out of Input Summary
	// =====================================================================
	logo(version);
	center_print("INPUT SUMMARY", 79);
	border_print();
	print_input_summary(input);

	// =====================================================================
	// Prepare Pole Paremeter Grids
	// =====================================================================
	border_print();
	center_print("INITIALIZATION", 79);
	border_print();
	
	start = omp_get_wtime();
	
	// Allocate & fill energy grids
	printf("Generating resonance distributions...\n");
	int * n_poles = generate_n_poles( input, &iseed );

	// Allocate & fill Window grids
	printf("Generating window distributions...\n");
	int * n_windows = generate_n_windows( input, &iseed );

	// Get material data
	printf("Loading Hoogenboom-Martin material data...\n");
	Materials materials = get_materials( input, &iseed ); 

	// Prepare full resonance grid
	printf("Generating resonance parameter grid...\n");
	Pole ** poles = generate_poles( input, n_poles, &iseed );

	// Prepare full Window grid
	printf("Generating window parameter grid...\n");
	Window ** windows = generate_window_params( input, n_windows, n_poles, &iseed);

	// Prepare 0K Resonances
	printf("Generating 0K l_value data...\n");
	double ** pseudo_K0RS = generate_pseudo_K0RS( input, &iseed );

	CalcDataPtrs data;
	data.n_poles = n_poles;
	data.n_windows = n_windows;
	data.materials = materials;
	data.poles = poles;
	data.windows = windows;
	data.pseudo_K0RS = pseudo_K0RS;

	stop = omp_get_wtime();
	printf("Initialization Complete. (%.2lf seconds)\n", stop-start);
	
	// =====================================================================
	// Cross Section (XS) Parallel Lookup Simulation Begins
	// =====================================================================
	border_print();
	center_print("SIMULATION", 79);
	border_print();
	
	printf("Beginning Simulation.\n");
	#ifndef STATUS
	printf("Calculating XS's...\n");
	#endif
	
	#ifdef PAPI
	/* initialize papi with one thread here  */
	if ( PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT){
		fprintf(stderr, "PAPI library init error!\n");
		exit(1);
	}
	#endif	

	start = omp_get_wtime();

	unsigned long vhash = 0;

	long g_abrarov = 0; 
	long g_alls = 0;
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

		#pragma omp for schedule(dynamic)
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

	// Final hash step
	vhash = vhash % 1000000;

	stop = omp_get_wtime();
	#ifndef PAPI
	printf("\nSimulation Complete.\n");
	#endif
	// =====================================================================
	// Print / Save Results and Exit
	// =====================================================================
	border_print();
	center_print("RESULTS", 79);
	border_print();

	printf("Threads:               %d\n", input.nthreads);
	if( input.doppler)
	printf("Slow Faddeeva:         %.2lf%%\n", (double) g_abrarov/g_alls * 100.f);
	printf("Runtime:               %.3lf seconds\n", stop-start);
	printf("Lookups:               "); fancy_int(input.lookups*input.particles);
	printf("Lookups/s:             "); fancy_int((double) input.lookups*input.particles / (stop-start));
	#ifdef VERIFICATION
	unsigned long long large = 425149;
	unsigned long long small = 830419;
	if( input.HM  == LARGE )
	{
		if( vhash == large )
			printf("Verification checksum: %lu (Valid)\n", vhash);
		else
			printf("Verification checksum: %lu (WARNING - INAVALID CHECKSUM!)\n", vhash);
	}
	else if( input.HM  == SMALL )
	{
		if( vhash == small )
			printf("Verification checksum: %lu (Valid)\n", vhash);
		else
			printf("Verification checksum: %lu (WARNING - INAVALID CHECKSUM!)\n", vhash);
	}
	#endif

	border_print();

	return 0;
}
