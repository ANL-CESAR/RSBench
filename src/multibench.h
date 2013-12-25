#include<stdio.h>
#include<omp.h>
#include<stdlib.h>
#include<time.h>

// io.c

typedef struct{
	int nthreads;
	int n_nuclides;
	int lookups;
	int HM_size;
} Inputs;

void logo(int version);
void center_print(const char *s, int width);
void border_print(void);
void fancy_int( int a );
Inputs read_CLI( int argc, char * argv[] );
void print_CLI_error(void);
