#include "multibench.h"

int main(int argc, char * argv[])
{
	// =====================================================================
	// Initialization & Command Line Read-In
	// =====================================================================

	int version = 0;
	int max_procs = omp_get_num_procs();
	int n_nuclides;
	int nthreads;
	int HM_size;
	double p_rrr;
	double p_UEG;
	
	srand(time(NULL));
	
	// Process CLI Fields
	Inputs input = read_CLI( argc, argv );
	
	// Set CLI variables
	nthreads =     input.nthreads;
	n_nuclides =   input.n_nuclides;
	lookups =      input.lookups;
	HM_size =      input.HM_size;

	// Set number of OpenMP Threads
	omp_set_num_threads(nthreads); 
	
	// =====================================================================
	// Calculate Estimate of Memory Usage
	// =====================================================================
	
	// =====================================================================
	// Print-out of Input Summary
	// =====================================================================
		
	logo(version);
	center_print("INPUT SUMMARY", 79);
	border_print();


	return 0;
}
