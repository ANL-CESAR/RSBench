#include<stdio.h>
#include<omp.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<math.h>
#include<stdint.h>
#include<float.h>

#ifdef PAPI
#include "papi.h"
#endif

#define PI 3.14159265359

// typedefs
typedef enum __hm{SMALL, LARGE, XL, XXL} HM_size;

#define HISTORY_BASED 1
#define EVENT_BASED 2

#define STARTING_SEED 1070
#define INITIALIZATION_SEED 42

typedef struct{
	double r;
	double i;
} RSComplex;

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
	int simulation_method;
} Input;

typedef struct{
	int * num_nucs;
	int * mats;
	double * concs;
	int max_num_nucs;
} Materials;

typedef struct{
	RSComplex MP_EA;
	RSComplex MP_RT;
	RSComplex MP_RA;
	RSComplex MP_RF;
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
	Pole * poles;
	Window ** windows;
	double ** pseudo_K0RS;
	int max_num_poles;
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
int * generate_n_poles( Input input,  uint64_t * seed );
int * generate_n_windows( Input input ,  uint64_t * seed);
Pole * generate_poles( Input input, int * n_poles, uint64_t * seed, int * max_num_poles );
Window ** generate_window_params( Input input, int * n_windows, int * n_poles ,  uint64_t * seed);
double ** generate_pseudo_K0RS( Input input ,  uint64_t * seed);

// material.c
int * load_num_nucs(Input input);
int * load_mats( Input input, int * num_nucs, int * max_num_nucs );
double * load_concs( int * num_nucs, uint64_t * seed, int max_num_nucs );
int pick_mat( uint64_t * seed );
Materials get_materials(Input input, uint64_t * seed);

// utils.c
size_t get_mem_estimate( Input input );
RSComplex fast_cexp( RSComplex z );

// xs_kernel.c
RSComplex fast_nuclear_W( RSComplex Z );
void calculate_macro_xs( double * macro_xs, int mat, double E, Input input, CalcDataPtrs data, RSComplex * sigTfactors, long * abrarov, long * alls ); 
void calculate_micro_xs( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, RSComplex * sigTfactors);
void calculate_micro_xs_doppler( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, RSComplex * sigTfactors, long * abrarov, long * alls);
void calculate_sig_T( int nuc, double E, Input input, CalcDataPtrs data, RSComplex * sigTfactors );

// simulation.c
void run_event_based_simulation(Input input, CalcDataPtrs data, long * abrarov_result, long * alls_result, unsigned long * vhash_result );
void run_history_based_simulation(Input input, CalcDataPtrs data, long * abrarov_result, long * alls_result, unsigned long * vhash_result );
double LCG_random_double(uint64_t * seed);
uint64_t LCG_random_int(uint64_t * seed);
uint64_t fast_forward_LCG(uint64_t seed, uint64_t n);

// rscomplex.c
RSComplex c_add( RSComplex A, RSComplex B);
RSComplex c_sub( RSComplex A, RSComplex B);
RSComplex c_mul( RSComplex A, RSComplex B);
RSComplex c_div( RSComplex A, RSComplex B);
double c_abs( RSComplex A);

// papi.c
void counter_init( int *eventset, int *num_papi_events );
void counter_stop( int * eventset, int num_papi_events );
