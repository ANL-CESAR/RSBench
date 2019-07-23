#include "rsbench.h"

int main(int argc, char * argv[])
{
	// =====================================================================
	// Initialization & Command Line Read-In
	// =====================================================================

	int version = 10;
	double start, stop;
	uint64_t seed = INITIALIZATION_SEED;
	
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
	int * n_poles = generate_n_poles( input, &seed );

	// Allocate & fill Window grids
	printf("Generating window distributions...\n");
	int * n_windows = generate_n_windows( input, &seed );

	// Get material data
	printf("Loading Hoogenboom-Martin material data...\n");
	Materials materials = get_materials( input, &seed ); 

	// Prepare full resonance grid
	printf("Generating resonance parameter grid...\n");
	int max_num_poles;
	Pole * poles = generate_poles( input, n_poles, &seed, &max_num_poles );

	// Prepare full Window grid
	printf("Generating window parameter grid...\n");
	int max_num_windows;
	Window * windows = generate_window_params( input, n_windows, n_poles, &seed, &max_num_windows);

	// Prepare 0K Resonances
	printf("Generating 0K l_value data...\n");
	double * pseudo_K0RS = generate_pseudo_K0RS( input, &seed );

	CalcDataPtrs data;
	data.n_poles = n_poles;
	data.n_windows = n_windows;
	data.materials = materials;
	data.poles = poles;
	data.windows = windows;
	data.pseudo_K0RS = pseudo_K0RS;
	data.max_num_poles = max_num_poles;
	data.max_num_windows = max_num_windows;

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

	// Run Simulation
	if(input.simulation_method == HISTORY_BASED )
		run_history_based_simulation(input, data, &g_abrarov, &g_alls, &vhash );
	else if( input.simulation_method == EVENT_BASED )
		run_event_based_simulation(input, data, &g_abrarov, &g_alls, &vhash );

	// Final hash step
	vhash = vhash % 999983;

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
	int lookups = 0;
	if( input.simulation_method == HISTORY_BASED )
		lookups = input.lookups*input.particles;
	else
		lookups = input.lookups;
	printf("Lookups:               "); fancy_int(lookups);
	printf("Lookups/s:             "); fancy_int((double) lookups / (stop-start));
	#ifdef VERIFICATION
	unsigned long long large = 0;
	unsigned long long small = 0;
	if(input.simulation_method == HISTORY_BASED )
	{
		large = 841550;
		small = 252526;
	}
	else if( input.simulation_method == EVENT_BASED )
	{
		large = 846231;
		small = 251196;
	}

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
