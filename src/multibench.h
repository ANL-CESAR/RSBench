#include<stdio.h>
#include<omp.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>

// io.c

typedef enum __hm{SMALL, LARGE, XL, XXL} HM_size;

typedef struct{
	int nthreads;
	int n_nuclides;
	int lookups;
	HM_size HM;
} Inputs;

void logo(int version);
void center_print(const char *s, int width);
void border_print(void);
void fancy_int( int a );
Inputs read_CLI( int argc, char * argv[] );
void print_CLI_error(void);

// mutlipole.c
int dbl_cmp( const void * a, const void * b );
double * generate_egrid( int n_resonances );
