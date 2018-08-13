#include<stdio.h>
#include<omp.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<math.h>
#include<complex.h>

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
	int avg_n_poles;
	int avg_n_windows;
	int numL;
	int doppler;
	int particles;
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
	double T;
	double A;
	double F;
	int start;
	int end;
} Window;

typedef struct{
	int * n_poles;
	int * n_windows;
	Materials materials;
	Pole ** poles;
	Window ** windows;
	double ** pseudo_K0RS;
} CalcDataPtrs;


// io.c
void logo(int version);
void center_print(const char *s, int width);
void border_print(void);
void fancy_int( int a );
Input read_CLI( int argc, char * argv[] );
void print_CLI_error(void);
void print_input_summary(Input input);

// init.c
int * generate_n_poles( Input input,  unsigned long * seed );
int * generate_n_windows( Input input ,  unsigned long * seed);
Pole ** generate_poles( Input input, int * n_poles ,  unsigned long * seed);
Window ** generate_window_params( Input input, int * n_windows, int * n_poles ,  unsigned long * seed);
double ** generate_pseudo_K0RS( Input input ,  unsigned long * seed);

// material.c
int * load_num_nucs(Input input);
int ** load_mats( Input input, int * num_nucs );
double ** load_concs( int * num_nucs, unsigned long * seed );
int pick_mat( unsigned long * seed );
Materials get_materials(Input input, unsigned long * seed);

// utils.c
double rn(unsigned long * seed);
unsigned long rn_i(unsigned long * seed);
size_t get_mem_estimate( Input input );
unsigned int hash(char *str, int nbins);
complex double fast_cexp( double complex z );
double fast_exp(double x);

// xs_kernel.c
double complex fast_nuclear_W( double complex Z );
void calculate_macro_xs( double * macro_xs, int mat, double E, Input input, CalcDataPtrs data, complex double * sigTfactors, long * abrarov, long * alls ); 
void calculate_micro_xs( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, complex double * sigTfactors);
void calculate_micro_xs_doppler( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, complex double * sigTfactors, long * abrarov, long * alls);
void calculate_sig_T( int nuc, double E, Input input, CalcDataPtrs data, complex double * sigTfactors );

// papi.c
void counter_init( int *eventset, int *num_papi_events );
void counter_stop( int * eventset, int num_papi_events );
