#include "rsbench.h"

int main(int argc, char * argv[])
{
	// =====================================================================
	// Initialization & Command Line Read-In
	// =====================================================================

	int version = 12;
	double start, stop;
	
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
	
	SimulationData SD = initialize_simulation( input );

	stop = omp_get_wtime();
	printf("Initialization Complete. (%.2lf seconds)\n", stop-start);
	
	// =====================================================================
	// Cross Section (XS) Parallel Lookup Simulation Begins
	// =====================================================================
	border_print();
	center_print("SIMULATION", 79);
	border_print();
	
	printf("Beginning Simulation.\n");
	printf("Calculating XS's...\n");

	start = omp_get_wtime();

	unsigned long vhash = 0;

	long g_abrarov = 0; 
	long g_alls = 0;

	// Run Simulation
	if(input.simulation_method == HISTORY_BASED )
		run_history_based_simulation(input, SD, &g_abrarov, &g_alls, &vhash );
	else if( input.simulation_method == EVENT_BASED )
		run_event_based_simulation(input, SD, &g_abrarov, &g_alls, &vhash );

	// Final hash step
	vhash = vhash % 999983;

	stop = omp_get_wtime();
	printf("\nSimulation Complete.\n");
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

	unsigned long long large = 0;
	unsigned long long small = 0;
	if(input.simulation_method == HISTORY_BASED )
	{
		large = 351414;
		small = 879238;
	}
	else if( input.simulation_method == EVENT_BASED )
	{
		large = 358421;
		small = 879382;
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

	border_print();

	return 0;
}
