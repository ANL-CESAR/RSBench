#include<stdio.h>
#include<omp.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<math.h>

#ifdef PAPI
#include "papi.h"
#endif

#define PI 3.14159265359

// typedefs
typedef enum __hm{SMALL, LARGE, XL, XXL} HM_size;

typedef struct{
	int nthreads;
	int n_nuclides;
	int lookups;
	HM_size HM;
	int n_poles;
	int n_windows;
} Input;

typedef struct{
	int * num_nucs;
	int ** mats;
	double ** concs;
} Materials;

typedef struct{
	complex double MP_EA;
	complex double MP_RT;
	complex double MP_RA;
	complex double MP_RF;
	short int l_value;
} Pole;

typedef struct{
	int * n_poles;
	Materials materials;
	double * nuclide_radii;
	Pole ** resonance_params;
	Window ** win_boundaries;
	double ** pseudo_K0RS;
} CalcDataPtrs;

typedef struct{
	double T;
	double A;
	double F;
} FitParams;

typedef struct{
	int start;
	int end;
} Window;

// io.c
void logo(int version);
void center_print(const char *s, int width);
void border_print(void);
void fancy_int( int a );
Input read_CLI( int argc, char * argv[] );
void print_CLI_error(void);
void print_input_summary(Input input);

// init.c
int * generate_n_poles( Input input );
double * generate_nuclide_radii( Input input );
Pole ** generate_resonance_params( Input input, int * n_poles );

// material.c
int * load_num_nucs(Input input);
int ** load_mats( Input input, int * num_nucs );
double ** load_concs( int * num_nucs );
int pick_mat( unsigned long * seed );
Materials get_materials(Input input);

// utils.c
double rn(unsigned long * seed);
size_t get_mem_estimate( Input input );

// xs_kernel.c
void calculate_macro_xs( double * macro_xs, int mat, double E, Input input, CalcDataPtrs data ); 
void calculate_micro_xs( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data);
complex double * calculate_sig_T( double E, Input input, CalcDataPtrs data );

// papi.c
void counter_init( int *eventset, int *num_papi_events );
void counter_stop( int * eventset, int num_papi_events );
