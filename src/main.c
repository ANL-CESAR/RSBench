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
	printf("\n");
	border_print();
	center_print("INITIALIZATION", 79);
	border_print();
	
	start = omp_get_wtime();
	
	// Allocate & fill energy grids
	int * n_resonances = generate_n_resonances( input );
	
	// Get material data
	Materials materials = get_materials( input ); 

	// Prepare material resonance grid
	double * nuclide_masses = generate_nuclide_masses( input );

	// Prepare full resonance grid
	Resonance ** resonance_params = generate_resonance_params( input, n_resonances );

	stop = omp_get_wtime();
	printf("Time taken for initialization: %lf seconds\n", stop-start);
	
	// =====================================================================
	// Cross Section (XS) Parallel Lookup Simulation Begins
	// =====================================================================
	printf("\n");
	border_print();
	center_print("SIMULATION", 79);
	border_print();
	
	start = omp_get_wtime();

	unsigned long seed = time(NULL);
	for( int i = 0; i < input.lookups; i++ )
	{
		int mat = pick_mat( &seed );
		double E = rn( &seed );
		
		double macro_xs[4];
		
		//calculate_macro_xs( macro_xs, mat, E, input ); 
	}
	
	stop = omp_get_wtime();
	
	// =====================================================================
	// Print / Save Results and Exit
	// =====================================================================
	printf("\n");
	border_print();
	center_print("RESULTS", 79);
	border_print();
	
	printf("Time taken for simulation: %lf seconds\n", stop-start);

	return 0;
}
