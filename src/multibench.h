#include<stdio.h>
#include<omp.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>

// typedefs
typedef enum __hm{SMALL, LARGE, XL, XXL} HM_size;

typedef struct{
	int nthreads;
	int n_nuclides;
	int lookups;
	HM_size HM;
	int n_resonances;
	double width;
} Input;

typedef struct{
	int * num_nucs;
	int ** mats;
	double ** concs;
} Materials;

// io.c
void logo(int version);
void center_print(const char *s, int width);
void border_print(void);
void fancy_int( int a );
Input read_CLI( int argc, char * argv[] );
void print_CLI_error(void);
void print_input_summary(Input input);

// mutlipole.c
int dbl_cmp( const void * a, const void * b );
double * generate_nuclide_energies( Input input );
double ** generate_egrid( Input input );

// material.c
int * load_num_nucs(Input input);
int ** load_mats( Input input, int * num_nucs );
double ** load_concs( int * num_nucs );
int pick_mat( unsigned long * seed );
Materials get_materials(Input input);

// utils.c
double rn(unsigned long * seed);
