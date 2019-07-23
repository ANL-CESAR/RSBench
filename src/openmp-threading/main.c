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
	// Intialize Simulation Data Structures
	// =====================================================================
	border_print();
	center_print("INITIALIZATION", 79);
	border_print();
	
	start = get_time();
	
	SimulationData SD = initialize_simulation( input );

	stop = get_time();
	printf("Initialization Complete. (%.2lf seconds)\n", stop-start);
	
	// =====================================================================
	// Cross Section (XS) Parallel Lookup Simulation Begins
	// =====================================================================
	border_print();
	center_print("SIMULATION", 79);
	border_print();

	start = get_time();

	unsigned long vhash = 0;

	// Run Simulation
	if(input.simulation_method == HISTORY_BASED )
		run_history_based_simulation(input, SD, &vhash );
	else if( input.simulation_method == EVENT_BASED )
		run_event_based_simulation(input, SD, &vhash );

	// Final hash step
	vhash = vhash % 999983;

	stop = get_time();
	printf("Simulation Complete.\n");

	// =====================================================================
	// Print / Save Results and Exit
	// =====================================================================
	border_print();
	center_print("RESULTS", 79);
	border_print();

	int is_invalid = validate_and_print_results(input, stop-start, vhash);

	border_print();

	return is_invalid;
}
