#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<sys/time.h>
#include<math.h>
#include<string.h>
#include<stdint.h>
#include<float.h>
#include<assert.h>
#define CL_TARGET_OPENCL_VERSION 120
#include <CL/cl.h>
#define MAX_SOURCE_SIZE (0x100000)

#include "cl_utils.h"

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
	int kernel_id;
  int platform_id;
  int device_id;
} Input;

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
	unsigned long length_n_poles;
	int * n_windows;
	unsigned long length_n_windows;
	Pole * poles;
	unsigned long length_poles;
	Window * windows;
	unsigned long length_windows;
	double * pseudo_K0RS;
	unsigned long length_pseudo_K0RS;
	int * num_nucs;
	unsigned long length_num_nucs;
	int * mats;
	unsigned long length_mats;
	double * concs;
	unsigned long length_concs;
	int max_num_nucs;
	int max_num_poles;
	int max_num_windows;
	double * p_energy_samples;
	unsigned long length_p_energy_samples;
	int * mat_samples;
	unsigned long length_mat_samples;
} SimulationData;

// io.c
void logo(int version);
void center_print(const char *s, int width);
void border_print(void);
void fancy_int( int a );
Input read_CLI( int argc, char * argv[] );
void print_CLI_error(void);
void print_input_summary(Input input);
int validate_and_print_results(Input input, double runtime, unsigned long vhash, double sim_runtime);

// init.c
SimulationData initialize_simulation( Input input );
int * generate_n_poles( Input input,  uint64_t * seed );
int * generate_n_windows( Input input ,  uint64_t * seed);
Pole * generate_poles( Input input, int * n_poles, uint64_t * seed, int * max_num_poles );
Window * generate_window_params( Input input, int * n_windows, int * n_poles, uint64_t * seed, int * max_num_windows );
double * generate_pseudo_K0RS( Input input, uint64_t * seed );

// material.c
int * load_num_nucs(Input input);
int * load_mats( Input input, int * num_nucs, int * max_num_nucs, unsigned long * length_mats );
double * load_concs( int * num_nucs, uint64_t * seed, int max_num_nucs );
SimulationData get_materials(Input input, uint64_t * seed);

// utils.c
size_t get_mem_estimate( Input input );
RSComplex fast_cexp( RSComplex z );
double get_time(void);

// simulation.c
unsigned long long run_event_based_simulation(Input in, SimulationData SD, double * sim_runtime);
double LCG_random_double(uint64_t * seed);
uint64_t LCG_random_int(uint64_t * seed);
RSComplex c_mul( RSComplex A, RSComplex B);

// CLutils.c
const char *getErrorString(cl_int error);
void check(cl_int error);
void printCompilerError( cl_program program, cl_device_id device );
void print_single_info( cl_platform_id platform, cl_device_id device);
void print_opencl_info(void);
