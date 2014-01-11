#include "multibench.h"

int main(int argc, char * argv[])
{
	// =====================================================================
	// Initialization & Command Line Read-In
	// =====================================================================

	int version = 0;
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
	// Prepare Resonance Paremeter Grids
	// =====================================================================
	border_print();
	center_print("INITIALIZATION", 79);
	border_print();
	
	start = omp_get_wtime();
	
	// Allocate & fill energy grids
	printf("Generating resonance distributions...\n");
	int * n_resonances = generate_n_resonances( input );
	
	// Get material data
	printf("Loading Hoogenboom-Martin material data...\n");
	Materials materials = get_materials( input ); 

	// Prepare material resonance grid
	printf("Generating nuclide data...\n");
	double * nuclide_radii = generate_nuclide_radii( input );

	// Prepare full resonance grid
	printf("Generating resonance parameter grid...\n");
	Resonance ** resonance_params = generate_resonance_params( input, n_resonances );
	
	CalcDataPtrs data;
	data.n_resonances = n_resonances;
	data.materials = materials;
	data.nuclide_radii = nuclide_radii;
	data.resonance_params = resonance_params;

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
			calculate_macro_xs( macro_xs, mat, E, input, data ); 
		}
	}

	stop = omp_get_wtime();
	printf("\nSimulation Complete.\n");

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
