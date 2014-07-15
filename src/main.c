#include "rsbench.h"

int main(int argc, char * argv[]) {
	dotp_driver(8);
	// =====================================================================
	// Initialization & Command Line Read-In
	// =====================================================================

	int version = 3;
	cudaEvent_t begin, end;
	float milliseconds = 0;
	cudaEventCreate(&begin);
	cudaEventCreate(&end);

	srand(time(NULL));

	// Process CLI Fields
	Input input = read_CLI( argc, argv );

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

	cudaEventRecord(begin, 0);

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
	CalcDataPtrs_d* data_d = init_data ( input, &data );
	//	free_CalcDataPtrs_d ( data_d );
	cudaEventRecord(end, 0);
	cudaEventSynchronize(end);
	cudaEventElapsedTime(&milliseconds, begin, end);
	printf("Initialization Complete 2. (%.2lf seconds)\n", milliseconds/1000);

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
	
	cudaEventRecord(begin, 0);
	top_calc_driver ( data_d, input, 200);
	cudaEventRecord(end, 0);
	cudaEventSynchronize(end);
	cudaEventElapsedTime(&milliseconds, begin, end);
	printf("\nSimulation 1 Complete.\n");
	// =====================================================================
	// Print / Save Results and Exit
	// =====================================================================
	border_print();
	center_print("RESULTS 1", 79);
	border_print();

	printf("NTPB:        %i\n", 200);
	printf ("Time for the kernel: %f ms\n", milliseconds);
	printf("Lookups:     "); fancy_int(input.lookups);
	printf("Lookups/s:   "); fancy_int((double) input.lookups / milliseconds * 1000 );

	cudaEventRecord(begin, 0);
	top_calc_driver ( data_d, input, 250);
	cudaEventRecord(end, 0);
	cudaEventSynchronize(end);
	cudaEventElapsedTime(&milliseconds, begin, end);
	printf("\nSimulation 2 Complete.\n");
	border_print();
	center_print("RESULTS 2", 79);
	border_print();

	printf("NTPB:        %i\n", 250);
	printf ("Time for the kernel: %f ms\n", milliseconds);
	printf("Lookups:     "); fancy_int(input.lookups);
	printf("Lookups/s:   "); fancy_int((double) input.lookups / milliseconds * 1000 );
	// free device memeory
	//	free_CalcDataPtrs_d ( data_d );

	cudaEventRecord(begin, 0);
	top_calc_driver ( data_d, input, 400);
	cudaEventRecord(end, 0);
	cudaEventSynchronize(end);
	cudaEventElapsedTime(&milliseconds, begin, end);
	printf("\nSimulation 3 Complete.\n");
	border_print();
	center_print("RESULTS 3", 79);
	border_print();

	printf("NTPB:        %i\n", 400);
	printf ("Time for the kernel: %f ms\n", milliseconds);
	printf("Lookups:     "); fancy_int(input.lookups);
	printf("Lookups/s:   "); fancy_int((double) input.lookups / milliseconds * 1000 );

	cudaEventRecord(begin, 0);
	top_calc_driver ( data_d, input, 500);
	cudaEventRecord(end, 0);
	cudaEventSynchronize(end);
	cudaEventElapsedTime(&milliseconds, begin, end);
	printf("\nSimulation 4 Complete.\n");
	border_print();
	center_print("RESULTS 4", 79);
	border_print();

	printf("NTPB:        %i\n", 500);
	printf ("Time for the kernel: %f ms\n", milliseconds);
	printf("Lookups:     "); fancy_int(input.lookups);
	printf("Lookups/s:   "); fancy_int((double) input.lookups / milliseconds * 1000 );

	cudaEventRecord(begin, 0);
	top_calc_driver ( data_d, input, 1000);
	cudaEventRecord(end, 0);
	cudaEventSynchronize(end);
	cudaEventElapsedTime(&milliseconds, begin, end);
	printf("\nSimulation 5 Complete.\n");
	border_print();
	center_print("RESULTS 5", 79);
	border_print();

	printf("NTPB:        %i\n", 1000);
	printf ("Time for the kernel: %f ms\n", milliseconds);
	printf("Lookups:     "); fancy_int(input.lookups);
	printf("Lookups/s:   "); fancy_int((double) input.lookups / milliseconds * 1000 );
	// =====================================================================
	// Print / Save Results and Exit
	// =====================================================================
	border_print();

	return 0;
}
