#include "rsbench.h"

// Prints program logo
void logo(int version)
{
	border_print();
	printf(
"                    _____   _____ ____                  _     \n"
"                   |  __ \\ / ____|  _ \\                | |    \n"
"                   | |__) | (___ | |_) | ___ _ __   ___| |__  \n"
"                   |  _  / \\___ \\|  _ < / _ \\ '_ \\ / __| '_ \\ \n"
"                   | | \\ \\ ____) | |_) |  __/ | | | (__| | | |\n"
"                   |_|  \\_\\_____/|____/ \\___|_| |_|\\___|_| |_|\n\n"
	);
	border_print();
	center_print("Developed at Argonne National Laboratory", 79);
	char v[100];
	sprintf(v, "Version: %d", version);
	center_print(v, 79);
	border_print();
}

// Prints Section titles in center of 80 char terminal
void center_print(const char *s, int width)
{
	int length = strlen(s);
	int i;
	for (i=0; i<=(width-length)/2; i++) {
		fputs(" ", stdout);
	}
	fputs(s, stdout);
	fputs("\n", stdout);
}

void border_print(void)
{
	printf(
	"==================================================================="
	"=============\n");
}

// Prints comma separated integers - for ease of reading
void fancy_int( int a )
{
    if( a < 1000 )
        printf("%d\n",a);

    else if( a >= 1000 && a < 1000000 )
        printf("%d,%03d\n", a / 1000, a % 1000);

    else if( a >= 1000000 && a < 1000000000 )
        printf("%d,%03d,%03d\n", a / 1000000, (a % 1000000) / 1000, a % 1000 );

    else if( a >= 1000000000 )
        printf("%d,%03d,%03d,%03d\n",
               a / 1000000000,
               (a % 1000000000) / 1000000,
               (a % 1000000) / 1000,
               a % 1000 );
    else
        printf("%d\n",a);
}

Input read_CLI( int argc, char * argv[] )
{
	Input input;

	// defaults to the history based simulation method
	input.simulation_method = HISTORY_BASED;
	// defaults to max threads on the system	
	input.nthreads = omp_get_num_procs();
	// defaults to 355 (corresponding to H-M Large benchmark)
	input.n_nuclides = 355;
	// defaults to 300,000
	input.particles = 300000;
	// defaults to 34
	input.lookups = 34;
	// defaults to H-M Large benchmark
	input.HM = LARGE;
	// defaults to 3000 resonancs (avg) per nuclide
	input.avg_n_poles = 1000;
	// defaults to 100
	input.avg_n_windows = 100;
	// defaults to 4;
	input.numL = 4;
	// defaults to no temperature dependence (Doppler broadening)
	input.doppler = 1;
	// defaults to baseline simulation kernel
	input.kernel_id = 0;
        // default to no binary read/write
	input.binary_mode = NONE;

	int default_lookups = 1;
	int default_particles = 1;

	// Collect Raw Input
	for( int i = 1; i < argc; i++ )
	{
		char * arg = argv[i];

		// nthreads (-t)
		if( strcmp(arg, "-t") == 0 )
		{
			if( ++i < argc )
				input.nthreads = atoi(argv[i]);
			else
				print_CLI_error();
		}
		// Simulation Method (-m)
		else if( strcmp(arg, "-m") == 0 )
		{
			char * sim_type;
			if( ++i < argc )
				sim_type = argv[i];
			else
				print_CLI_error();

			if( strcmp(sim_type, "history") == 0 )
				input.simulation_method = HISTORY_BASED;
			else if( strcmp(sim_type, "event") == 0 )
			{
				input.simulation_method = EVENT_BASED;
				// Also resets default # of lookups
				if( default_lookups && default_particles )
				{
					input.lookups =  input.lookups * input.particles;
					input.particles = 0;
				}
			}
			else
				print_CLI_error();
		}
		// lookups (-l)
		else if( strcmp(arg, "-l") == 0 )
		{
			if( ++i < argc )
			{
				input.lookups = atoi(argv[i]);
				default_lookups = 0;
			}
			else
				print_CLI_error();
		}
		// particles (-p)
		else if( strcmp(arg, "-p") == 0 )
		{
			if( ++i < argc )
			{
				input.particles = atoi(argv[i]);
				default_particles = 0;
			}
			else
				print_CLI_error();
		}
		// nuclides (-n)
		else if( strcmp(arg, "-n") == 0 )
		{
			if( ++i < argc )
				input.n_nuclides = atoi(argv[i]);
			else
				print_CLI_error();
		}
		// HM (-s)
		else if( strcmp(arg, "-s") == 0 )
		{	
			if( ++i < argc )
			{
				if( strcmp(argv[i], "small") == 0 )
					input.HM = SMALL;
				else if ( strcmp(argv[i], "large") == 0 )
					input.HM = LARGE;
				else
					print_CLI_error();
			}
			else
				print_CLI_error();
		}
		// Doppler Broadening (Temperature Dependence)
		else if( strcmp(arg, "-d") == 0 )
		{	
			input.doppler = 0;
		}
		// Avg number of windows per nuclide (-w)
		else if( strcmp(arg, "-W") == 0 )
		{
			if( ++i < argc )
				input.avg_n_windows = atoi(argv[i]);
			else
				print_CLI_error();
		}
		// Avg number of poles per nuclide (-p)
		else if( strcmp(arg, "-P") == 0 )
		{
			if( ++i < argc )
				input.avg_n_poles = atoi(argv[i]);
			else
				print_CLI_error();
		}
                // binary mode (-b)
		else if( strcmp(arg, "-b") == 0 )
		{
			char * binary_mode;
			if( ++i < argc )
				binary_mode = argv[i];
			else
				print_CLI_error();

			if( strcmp(binary_mode, "read") == 0 )
				input.binary_mode = READ;
			else if( strcmp(binary_mode, "write") == 0 )
				input.binary_mode = WRITE;
			else
				print_CLI_error();
		}
		// Kernel ID (-k)
		else if( strcmp(arg, "-k") == 0 )
		{
			if( ++i < argc )
				input.kernel_id = atoi(argv[i]);
			else
				print_CLI_error();
		}
		else
			print_CLI_error();
	}

	// Validate Input

	// Validate nthreads
	if( input.nthreads < 1 )
		print_CLI_error();
	
	// Validate n_isotopes
	if( input.n_nuclides < 1 )
		print_CLI_error();
	
	// Validate lookups
	if( input.lookups < 1 )
		print_CLI_error();
	
	// Validate lookups
	if( input.avg_n_poles < 1 )
		print_CLI_error();
	
	// Validate lookups
	if( input.avg_n_windows < 1 )
		print_CLI_error();
	
	// Set HM size specific parameters
	// (defaults to large)
	if( input.HM == SMALL )
		input.n_nuclides = 68;

	// Return input struct
	return input;
}

void print_CLI_error(void)
{
	printf("Usage: ./multibench <options>\n");
	printf("Options include:\n");
	printf("  -t <threads>            Number of OpenMP threads to run\n");
	printf("  -m <simulation method>  Simulation method (history, event)\n");
	printf("  -s <size>               Size of H-M Benchmark to run (small, large)\n");
	printf("  -l <lookups>            Number of Cross-section (XS) lookups per particle history\n");
	printf("  -p <particles>          Number of particle histories\n");
	printf("  -P <poles>              Average Number of Poles per Nuclide\n");
	printf("  -W <poles>              Average Number of Windows per Nuclide\n");
	printf("  -d                      Disables Temperature Dependence (Doppler Broadening)\n");
	printf("Default is equivalent to: -s large -m history -l 34 -p 300000 -P 1000 -W 100\n");
	printf("See readme for full description of default run values\n");
	exit(4);
}

void print_input_summary(Input input)
{
	// Calculate Estimate of Memory Usage
	size_t mem = get_mem_estimate(input);

	if( input.simulation_method == EVENT_BASED )
		printf("Simulation Method:           Event Based\n");
	else
		printf("Simulation Method:           History Based\n");
	printf("Materials:                   12\n");
	printf("H-M Benchmark Size:          ");
	if( input.HM == 0 )
		printf("Small\n");
	else
		printf("Large\n");
	if( input.doppler == 1 )
		printf("Temperature Dependence:      ON\n");
	else
		printf("Temperature Dependence:      OFF\n");
	printf("Total Nuclides:              %d\n", input.n_nuclides);
	printf("Avg Poles per Nuclide:       "); fancy_int(input.avg_n_poles);
	printf("Avg Windows per Nuclide:     "); fancy_int(input.avg_n_windows);

	int lookups = input.lookups;
	if( input.simulation_method == HISTORY_BASED )
	{
		printf("Particles:                   "); fancy_int(input.particles);
		printf("XS Lookups per Particle:     "); fancy_int(input.lookups);
		lookups *= input.particles;
	}
	printf("Total XS Lookups:            "); fancy_int(lookups);
	printf("Threads:                     %d\n", input.nthreads);
	printf("Est. Memory Usage (MB):      %.1lf\n", mem / 1024.0 / 1024.0);
}

int validate_and_print_results(Input input, double runtime, unsigned long vhash)
{
	printf("Threads:               %d\n", input.nthreads);
	printf("Runtime:               %.3lf seconds\n", runtime);
	int lookups = 0;
	if( input.simulation_method == HISTORY_BASED )
		lookups = input.lookups*input.particles;
	else
		lookups = input.lookups;
	printf("Lookups:               "); fancy_int(lookups);
	printf("Lookups/s:             "); fancy_int((double) lookups / (runtime));

	int is_invalid = 1;

	unsigned long long large = 0;
	unsigned long long small = 0;
	if(input.simulation_method == HISTORY_BASED )
	{
		large = 351485;
		small = 879693;
	}
	else if( input.simulation_method == EVENT_BASED )
	{
		large = 358389;
		small = 880018;
	}

	if( input.HM  == LARGE )
	{
		if( vhash == large )
		{
			printf("Verification checksum: %lu (Valid)\n", vhash);
			is_invalid = 0;
		}
		else
			printf("Verification checksum: %lu (WARNING - INAVALID CHECKSUM!)\n", vhash);
	}
	else if( input.HM  == SMALL )
	{
		if( vhash == small )
		{
			printf("Verification checksum: %lu (Valid)\n", vhash);
			is_invalid = 0;
		}
		else
			printf("Verification checksum: %lu (WARNING - INAVALID CHECKSUM!)\n", vhash);
	}

	return is_invalid;
}

void binary_write( Input in, SimulationData SD )
{
	char * fname = "RS_data.dat";
	printf("Writing all data structures to binary file %s...\n", fname);
	FILE * fp = fopen(fname, "w");

	// Write SimulationData Object. Include pointers, even though we won't be using them.
	fwrite(&SD, sizeof(SimulationData), 1, fp);

	// Write heap arrays in SimulationData Object
        fwrite(SD.n_poles,          sizeof(int),    SD.length_n_poles, fp);
        fwrite(SD.n_windows,        sizeof(int),    SD.length_n_windows, fp);
        fwrite(SD.poles,            sizeof(Pole),   SD.length_poles, fp);
        fwrite(SD.windows,          sizeof(Window), SD.length_windows, fp);
        fwrite(SD.pseudo_K0RS,      sizeof(double), SD.length_pseudo_K0RS, fp);
        fwrite(SD.num_nucs,         sizeof(int),    SD.length_num_nucs, fp);
        fwrite(SD.mats,             sizeof(int),    SD.length_mats, fp);
        fwrite(SD.concs,            sizeof(double), SD.length_concs, fp);

	fclose(fp);
}

SimulationData binary_read( Input in )
{
	SimulationData SD;

	char * fname = "RS_data.dat";
	printf("Reading all data structures from binary file %s...\n", fname);

	FILE * fp = fopen(fname, "r");
	assert(fp != NULL);

	// Read SimulationData Object. Include pointers, even though we won't be using them.
	fread(&SD, sizeof(SimulationData), 1, fp);

	// Allocate space for arrays on heap
        SD.n_poles =          malloc(SD.length_n_poles * sizeof(int));
        SD.n_windows =        malloc(SD.length_n_windows * sizeof(int));
        SD.poles =            malloc(SD.length_poles*sizeof(Pole));
        SD.windows =          malloc(SD.length_windows*sizeof(Window));
        SD.pseudo_K0RS =      malloc(SD.length_pseudo_K0RS*sizeof(double));
        SD.num_nucs =         malloc(SD.length_num_nucs*sizeof(int));
        SD.mats =             malloc(SD.length_mats*sizeof(int));
        SD.concs =            malloc(SD.length_concs*sizeof(double));
        SD.mat_samples =      malloc(SD.length_mat_samples*sizeof(int));

	// Read heap arrays into SimulationData Object
        fread(SD.n_poles,          sizeof(int),    SD.length_n_poles, fp);
        fread(SD.n_windows,        sizeof(int),    SD.length_n_windows, fp);
        fread(SD.poles,            sizeof(Pole),   SD.length_poles, fp);
        fread(SD.windows,          sizeof(Window), SD.length_windows, fp);
        fread(SD.pseudo_K0RS,      sizeof(double), SD.length_pseudo_K0RS, fp);
        fread(SD.num_nucs,         sizeof(int),    SD.length_num_nucs, fp);
        fread(SD.mats,             sizeof(int),    SD.length_mats, fp);
        fread(SD.concs,            sizeof(double), SD.length_concs, fp);

	fclose(fp);

	return SD;
}
