#include "rsbench.h"

int main(int argc, char * argv[])
{
	// =====================================================================
	// Initialization & Command Line Read-In
	// =====================================================================

	int version = 3;
	int max_procs = omp_get_num_procs();
	double start, stop;
	
	srand(time(NULL));
	
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
	int * n_poles = generate_n_poles( input );

	// Allocate & fill Window grids
	printf("Generating window distributions...\n");
	int * n_windows = generate_n_windows( input );

	// Get material data
	printf("Loading Hoogenboom-Martin material data...\n");
	Materials materials = get_materials( input ); 

	// Prepare full resonance grid
	printf("Generating resonance parameter grid...\n");
	Pole ** poles = generate_poles( input, n_poles );

	// Prepare full Window grid
	printf("Generating window parameter grid...\n");
	Window ** windows = generate_window_params( input, n_windows, n_poles);

	// Prepare 0K Resonances
	printf("Generating 0K l_value data...\n");
	double ** pseudo_K0RS = generate_pseudo_K0RS( input );

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

	unsigned long seed = rand();
	int mat;
	double E;
	int i;

	#pragma omp parallel default(none) \
	private(seed, mat, E, i) \
	shared(input, data) 
	{
		double macro_xs[4];
		int thread = omp_get_thread_num();
		seed += thread;
		
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
		for( i = 0; i < input.lookups; i++ )
		{
			#ifdef STATUS
			if( thread == 0 && i % 1000 == 0 )
				printf("\rCalculating XS's... (%.0lf%% completed)",
						(i / ( (double)input.lookups /
						(double) input.nthreads )) /
						(double) input.nthreads * 100.0);
			#endif
			mat = pick_mat( &seed );
			E = rn( &seed );
			calculate_macro_xs( macro_xs, mat, E, input, data, sigTfactors ); 
		}

		free(sigTfactors);
		
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

	printf("Threads:     %d\n", input.nthreads);
	printf("Runtime:     %.3lf seconds\n", stop-start);
	printf("Lookups:     "); fancy_int(input.lookups);
	printf("Lookups/s:   "); fancy_int((double) input.lookups / (stop-start));

	border_print();

	return 0;
}
