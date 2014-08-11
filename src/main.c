#include "rsbench.h"

void run_test (CalcDataPtrs_d* data_d, Input input, cudaEvent_t* begin, cudaEvent_t* end, int ntpb, int idx, int dist_type) {
	float milliseconds = 0;
	milliseconds = top_calc_driver ( data_d, ntpb, input, dist_type, begin, end);
	printf("\nSimulation %i Complete.\n", idx);
	border_print();
	char str[20];
	sprintf(str, "RESULTS\t%i", idx);
	center_print(str, 79);
	border_print();

	printf("NTPB:        %i\n", ntpb);
	printf ("Time for the kernel: %f ms\n", milliseconds);
	printf("Lookups:     "); fancy_int(input.lookups);
	printf("Lookups/s:   "); fancy_int((double) input.lookups / milliseconds * 1000 );
}

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

	size_t free, total;
	cudaMemGetInfo(&free,  &total); 
	printf ("%f %f\n", free/1024.0/1024.0, total/1024.0/1024.0);
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
	int ntpbs []= {32, 64, 96, 128, 192, 256, 512, 1024};	
	for ( int i = 3; i < 4/*sizeof(ntpbs)/sizeof(int)*/; i ++) 
		run_test (data_d, input, &begin, &end, ntpbs[i], i, 0);
	border_print();
//	for ( int i = 0; i < sizeof(ntpbs)/sizeof(int); i ++) 
//		run_test (data_d, input, &begin, &end, ntpbs[i], i, 1);
//	border_print();

	return 0;
}
