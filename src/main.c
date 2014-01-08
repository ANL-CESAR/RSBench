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
	input.n_resonances = 6000;
	input.width = 0.00001;

	// Set number of OpenMP Threads
	omp_set_num_threads(input.nthreads); 
	
	// =====================================================================
	// Calculate Estimate of Memory Usage
	// =====================================================================
	
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
	printf("\n");
	border_print();
	center_print("INITIALIZATION", 79);
	border_print();
	
	start = omp_get_wtime();
	
	// Allocate & fill energy grids
	double ** egrid = generate_egrid( input );
	
	// Get material data
	Materials materials = get_materials( input ); 

	stop = omp_get_wtime();
	printf("Time taken for initialization: %lf seconds\n", stop-start);
	
	
	// =====================================================================
	// Cross Section (XS) Parallel Lookup Simulation Begins
	// =====================================================================
	printf("\n");
	border_print();
	center_print("SIMULATION", 79);
	border_print();
	
	// OpenMP compiler directives - declaring variables as shared or private

	/*
	#pragma omp parallel default(none) \
	private(i, thread, p_energy, mat, seed) \
	shared( max_procs, n_isotopes, n_gridpoints, \
	energy_grid, nuclide_grids, lookups, nthreads, \
	mats, concs, num_nucs, mype, vhash) 
	{	
		#pragma omp critical
		{
			counter_init(&eventset, &num_papi_events);
		}
		#endif
		
		double macro_xs_vector[5];
		thread = omp_get_thread_num();
		seed   = (thread+1)*19+17;

		#pragma omp for schedule(dynamic)
		for( i = 0; i < lookups; i++ )
		{
			// Status text
			if( INFO && mype == 0 && thread == 0 && i % 1000 == 0 )
				printf("\rCalculating XS's... (%.0lf%% completed)",
						(i / ( (double)lookups / (double) nthreads )) / (double) nthreads * 100.0);
		
			// Randomly pick an energy and material for the particle
			p_energy = rn(&seed);
			mat      = pick_mat(&seed); 
			
			// debugging
			//printf("E = %lf mat = %d\n", p_energy, mat);
				
			// This returns the macro_xs_vector, but we're not going
			// to do anything with it in this program, so return value
			// is written over.
			calculate_macro_xs( p_energy, mat, n_isotopes,
			                    n_gridpoints, num_nucs, concs,
			                    energy_grid, nuclide_grids, mats,
                                macro_xs_vector );
		}

		#ifdef PAPI
		if( mype == 0 && thread == 0 )
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
	*/
	
	// =====================================================================
	// Print / Save Results and Exit
	// =====================================================================
	printf("\n");
	border_print();
	center_print("RESULTS", 79);
	border_print();

	return 0;
}
