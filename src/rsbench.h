#ifndef __RSBENCH_HEADER_H__
#define __RSBENCH_HEADER_H__

#include<stdio.h>
#include<omp.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<math.h>
#include<complex.h>
#include<cuda.h>
#include<cuda_runtime.h>
#include <cuComplex.h>
#include<assert.h>

#ifdef PAPI
#include "papi.h"
#endif

#define PI 3.1415926535897932384626433832795

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
	cuDoubleComplex MP_EA;
	cuDoubleComplex MP_RT;
	cuDoubleComplex MP_RA;
	cuDoubleComplex MP_RF;
	short int l_value;
} Pole_CUDA;

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

static int numL = 0;

#ifdef __cplusplus
  #define RESTRICT __restrict
  #define C_LINKAGE extern "C" 
#else
  #define RESTRICT 
  #define C_LINKAGE 
#endif

#ifdef __cplusplus
  extern "C" {
#endif
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
int * generate_n_windows( Input input );
Pole ** generate_poles( Input input, int * n_poles );
Window ** generate_window_params( Input input, int * n_windows, int * n_poles );
double ** generate_pseudo_K0RS( Input input );

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
void calculate_macro_xs( double * macro_xs, int mat, double E, Input input, CalcDataPtrs data, complex double * sigTfactors, int* counter, int* counter2 ); 
void calculate_micro_xs( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, complex double * sigTfactors, int* counter );
void calculate_sig_T( int nuc, double E, Input input, CalcDataPtrs data, complex double * sigTfactors );
void calculate_sig_T_sim ( double E, int num_iter, const double* data, complex double * sigTfactors );

//CUDA_Driver.c
int dotp_driver(int NTPB);
void calculate_micro_xs_driver( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, complex double * sigTfactors);
void calc_sig_driver ( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, complex double * sigTfactors );

// papi.c
void counter_init( int *eventset, int *num_papi_events );
void counter_stop( int * eventset, int num_papi_events );
#ifdef  __cplusplus
  }
#endif

#endif
