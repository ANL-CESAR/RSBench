#include"multibench.h"

// Prints program logo
void logo(int version)
{
	border_print();
	printf(
"           __  __         _  _    _  ____                      _     \n"
"          |  \\/  | _   _ | || |_ (_)| __ )   ___  _ __    ___ | |__  \n"
"          | |\\/| || | | || || __|| ||  _ \\  / _ \\| '_ \\  / __|| '_ \\ \n"
"          | |  | || |_| || || |_ | || |_) ||  __/| | | || (__ | | | |\n"
"          |_|  |_| \\__,_||_| \\__||_||____/  \\___||_| |_| \\___||_| |_|\n\n"
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
	
	// defaults to max threads on the system	
	input.nthreads = omp_get_num_procs();
	
	// defaults to 355 (corresponding to H-M Large benchmark)
	input.n_nuclides = 355;
	
	// defaults to 1,500,000
	input.lookups = 1500000;
	
	// defaults to H-M Large benchmark
	input.HM = LARGE;
	
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
		// lookups (-l)
		else if( strcmp(arg, "-l") == 0 )
		{
			if( ++i < argc )
				input.lookups = atoi(argv[i]);
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
				else if ( strcmp(argv[i], "XL") == 0 )
					input.HM = XL;
				else if ( strcmp(argv[i], "XXL") == 0 )
					input.HM = XXL;
				else
					print_CLI_error();
			}
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
	
	// Set HM size specific parameters
	// (defaults to large)

	// Return input struct
	return input;
}

void print_CLI_error(void)
{
	printf("Usage: ./XSBench <options>\n");
	printf("Options include:\n");
	printf("  -t <threads>     Number of OpenMP threads to run\n");
	printf("  -s <size>        Size of H-M Benchmark to run (small, large, XL, XXL)\n");
	printf("  -l <lookups>     Number of Cross-section (XS) lookups\n");
	printf("Default is equivalent to: -s large -l 15000000\n");
	printf("See readme for full description of default run values\n");
	exit(4);
}

void print_input_summary(Input input)
{
	// Calculate Estimate of Memory Usage
	size_t mem = get_mem_estimate(input);

	printf("Threads:        %d\n", input.nthreads);
	printf("Nuclides:       %d\n", input.n_nuclides);
	printf("Lookups:        "); fancy_int(input.lookups);
	printf("HM Size:        %d\n", input.HM);
	printf("Ave Resonances: "); fancy_int(input.n_resonances);
	printf("Mem Usage (MB): %.1lf\n", mem / 1024.0 / 1024.0);
}
