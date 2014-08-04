#include "rsbench.h"

#define cudaCheckErrors(msg) \
	do { \
		cudaError_t __err = cudaGetLastError(); \
		if (__err != cudaSuccess) { \
			fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
					msg, cudaGetErrorString(__err), \
					__FILE__, __LINE__); \
			fprintf(stderr, "*** FAILED - ABORTING\n"); \
			exit(1); \
		} \
	} while (0)

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true){
	if (code != cudaSuccess)
	{
		fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

//	add val to the vlaue stored at address
__device__ double atomicAdd(double* address, double val) {
	unsigned long long int* address_as_ull =
		(unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed,
				__double_as_longlong(val +
					__longlong_as_double(assumed)));
	} while (assumed != old);
	return __longlong_as_double(old);
}

// Updated distribution built from actual OpenMC look-ups 
__constant__ double dist[12] = {
		0.207834,	// fuel
		0.381401,	// cladding
		0.207763,	// cold, borated water
		0.198185,	// hot, borated water
		0.000036,	// RPV
		0.000032,	// Lower, radial reflector
		0.000039,	// Upper reflector / top plate
		0.000231,	// bottom plate
		0.000406,	// bottom nozzle
		0.000140,	// top nozzle
		0.002414,	// top of fuel assemblies
		0.001519 	// bottom of fuel assemblies
	 /* 0.140,	// fuel
	0.052,	// cladding
	0.275,	// cold, borated water
	0.134,	// hot, borated water
	0.154,	// RPV
	0.064,	// Lower, radial reflector
	0.066,	// Upper reflector / top plate
	0.055,	// bottom plate
	0.008,	// bottom nozzle
	0.015,	// top nozzle
	0.025,	// top of fuel assemblies
	0.013*/ 	// bottom of fuel assemblies
};	

__device__ double devRn(unsigned long * seed) {
	unsigned long m = 2147483647;
	unsigned long n1 = ( 16807ul * (*seed) ) % m;
	(*seed) = n1;
	return (double) n1 / (double) m;
}

// picks a material based on a probabilistic distribution
__device__ int devPick_mat( unsigned long * seed ) {
	double roll = devRn(seed);
	// makes a pick based on the distro
	double running = 0;
	for( int i = 0; i < 12; i++ ) {
		running += dist[i];
		if( roll < running )
			return i;
	}
	return 11;/*
		     for( int i = 0; i < 12; i++ ){
		     double running = 0;
		     for( int j = i; j > 0; j-- )
		     running += dist[j];
		     if( roll < running )
		     return i;
		     }

		     return 0;*/
}

__device__ void devCalc_micro_xs( double * micro_xs, int nuc, double E, Input input, const CalcDataPtrs_d* data, cuDoubleComplex * sigTfactors) {
	// MicroScopic XS's to Calculate
	double sigT;	double sigA;	double sigF;	double sigE;

	// Calculate Window Index
	double spacing = 1.0 / data->n_windows[nuc];
	int window = (int) ( E / spacing );
	if( window == data->n_windows[nuc] )
		window--;

	// Calculate sigTfactors
	double * d_ptr = &(data->pseudo_K0RS_2d [nuc * input.numL]);
	double phi;
	double sqrt_E = sqrt(E);
	for( int i = 0; i < input.numL; i++ ) {
		phi = d_ptr[i] * sqrt_E;

		if( i == 1 )
			phi -= - atan( phi );
		else if( i == 2 )
			phi -= atan( 3.0 * phi / (3.0 - phi*phi));
		else if( i == 3 )
			phi -= atan(phi*(15.0-phi*phi)/(15.0-6.0*phi*phi));

		phi += phi;

		sigTfactors[i].x= cos(phi);
		sigTfactors[i].y= - sin(phi);
	}
	// Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
	Window w = data->windows_2d[nuc * data->pitch_windows + window];
	sigT = E * w.T;	sigA = E * w.A;	sigF = E * w.F;
	// Loop over Poles within window, add contributions
	cuDoubleComplex const1 = make_cuDoubleComplex(0, 1/E), const2 = make_cuDoubleComplex( sqrt_E, 0);
	for( int i = w.start; i < w.end; i++ ){
		Pole pole = data->poles_2d[ nuc * data->pitch_poles + i];
		cuDoubleComplex CDUM = cuCdiv( const1, cuCsub( pole.MP_EA, const2 ) );
		sigT += cuCreal( cuCmul( pole.MP_RT, cuCmul( CDUM, sigTfactors[pole.l_value] ) ) );
		sigA += cuCreal( cuCmul( pole.MP_RA, CDUM) );
		sigF += cuCreal( cuCmul( pole.MP_RF, CDUM) );
	}

	sigE = sigT - sigA;
	micro_xs[0] = sigT; micro_xs[1] = sigA;	micro_xs[2] = sigF; micro_xs[3] = sigE;
}

__device__ void devCalc_macro_xs( double * macro_xs, int mat, double E, Input input, const CalcDataPtrs_d* data, cuDoubleComplex * sigTfactors) {
	// zero out macro vector
	for( int i = 0; i < 4; i++ )
		macro_xs[i] = 0;

	// for nuclide in mat
	for( int i = 0; i < data->materials.num_nucs[mat]; i++ ){
		double micro_xs[4];
		int nuc = data->materials.mats_2d[mat* data->materials.pitch + i];

		devCalc_micro_xs( micro_xs, nuc, E, input, data, sigTfactors);

		for( int j = 0; j < 4; j++ ){
			macro_xs[j] += micro_xs[j] * data->materials.concs_2d[mat* data->materials.pitch+i];
		}
	}
}

//	top level kernel - 4th version
__global__ void calc_kernel (const CalcDataPtrs_d* data, int lookups, int numL/*, int* ints_d*/)   {
	// going to be dynamic
	int tmp = threadIdx.x + blockIdx.x * blockDim.x;
	unsigned long seed = tmp;
	if ( tmp >= lookups)
		return;
	int mat = devPick_mat( &seed );
	double E = devRn( &seed );
	double sqrt_E = sqrt(E);
	double macro_xs[4];
	cuDoubleComplex sigTfactors[4];
	cuDoubleComplex const1 = make_cuDoubleComplex(0, 1/E), const2 = make_cuDoubleComplex( sqrt_E, 0);
	// zero out macro vector
	int i;
#pragma unroll 4
	for( i = 0; i < 4; i++ ){	macro_xs[i] = 0;	}
	// for nuclide in mat;
	for( i = 0; i < data->materials.num_nucs[mat]; i++ ){
		double micro_xs[4];
		int nuc = data->materials.mats_2d[mat* data->materials.pitch + i];

		// MicroScopic XS's to Calculate
		double sigT;	double sigA;	double sigF;	double sigE;

		// Calculate Window Index
		double spacing = 1.0 / data->n_windows[nuc];
		int window = (int) ( E / spacing );
		if( window == data->n_windows[nuc] )
			window--;

		// Calculate sigTfactors
		double * d_ptr = &(data->pseudo_K0RS_2d [nuc * numL]);
		double phi;
		for( int k = 0; k < numL; k++ ) {
			phi = d_ptr[k] * sqrt_E;

			if( k == 1 )
				phi -= - atan( phi );
			else if( k == 2 )
				phi -= atan( 3.0 * phi / (3.0 - phi*phi));
			else if( k == 3 )
				phi -= atan(phi*(15.0-phi*phi)/(15.0-6.0*phi*phi));

			phi += phi;

			sigTfactors[k].x= cos(phi);
			sigTfactors[k].y= - sin(phi);
		}
		// Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
		Window w = data->windows_2d[nuc * data->pitch_windows + window];
		sigT = E * w.T;	sigA = E * w.A;	sigF = E * w.F;
		// Loop over Poles within window, add contributions
		for( int k = w.start; k < w.end; k ++ ){
			Pole pole = data->poles_2d[ nuc * data->pitch_poles + k];
			cuDoubleComplex CDUM = cuCdiv( const1, cuCsub( pole.MP_EA, const2 ) );
			sigT += cuCreal( cuCmul( pole.MP_RT, cuCmul( CDUM, sigTfactors[pole.l_value] ) ) );
			sigA += cuCreal( cuCmul( pole.MP_RA, CDUM) );
			sigF += cuCreal( cuCmul( pole.MP_RF, CDUM) );
		}

		sigE = sigT - sigA;
		micro_xs[0] = sigT; micro_xs[1] = sigA;	micro_xs[2] = sigF; micro_xs[3] = sigE;

#pragma unroll 4
		for( int j = 0; j < 4; j++ ){	macro_xs[j] += micro_xs[j] * data->materials.concs_2d[mat* data->materials.pitch+i];}
	}
}

//	top level driver - 4th version
void top_calc_driver (const CalcDataPtrs_d* data, int ntpb, Input input){
	int num_blocs = input.lookups/ntpb;
	if ( ntpb * num_blocs < input.lookups )
		num_blocs ++;
	printf ("%i %i\n", num_blocs, ntpb);	
	calc_kernel<<<num_blocs, ntpb>>> ( data, input.lookups, input.numL/*, ints_d*/);
	cudaDeviceSynchronize();
}

//	compute SigT (the finest)
__device__  void calc_sig_T ( int i, double phi, cuDoubleComplex* rslt) {
	if ( i == 1 )
		phi -= atan ( phi );
	else if( i == 2 )
		phi -= atan ( 3.0 * phi / (3.0 - phi*phi));
	else if( i == 3 )
		phi -= atan (phi*(15.0-phi*phi)/(15.0-6.0*phi*phi));

	phi += phi;

	rslt -> x = cos(phi); 
	rslt -> y = sin(phi);
}

__global__ void calc_sig_T_sim_kernel ( double E, int num_iter, 
		const double* data, cuDoubleComplex* gTfactors) {
	calc_sig_T (threadIdx.x, data[threadIdx.x] * sqrt(E), &gTfactors[threadIdx.x]);	
}

__global__ void calc_sig_T_sim_dd_kernel ( double E, int num_iter, 
		const CalcDataPtrs_d* data, int drift, cuDoubleComplex* gTfactors) {
	calc_sig_T (threadIdx.x, data->pseudo_K0RS_2d[drift+threadIdx.x] * sqrt(E), &gTfactors[threadIdx.x]);	
}

//	come after the third kernel macro_kernel
__global__ void add_macro_xs (double* macro_xs_d, double* macro_xs_result ) {
	for ( int i = 0;i < 4; i ++)
		atomicAdd( &(macro_xs_d[i]), macro_xs_result[threadIdx.x + i] );
} 

//	come after the third kernel macro_kernel
__global__ void add_macro_xs_single_thread (double* macro_xs_d, double* macro_xs_result, int num ) {
	for ( int j = 0; j < num; j ++)
		for ( int i = 0;i < 4; i ++)
			atomicAdd( &(macro_xs_d[i]), macro_xs_result[j + i] );
} 

//	come after the third kernel macro_kernel
__global__ void add_macro_xs_four_threads (double* macro_xs_d, double* macro_xs_result, int num ) {
	for ( int j = 0; j < num; j ++)
		macro_xs_d[threadIdx.x] += macro_xs_result[j + threadIdx.x];
} 

//	third kernel utilizing CalcDataPtrs_d allocated on device
__global__ void macro_kernel (double * macro_xs, const CalcDataPtrs_d* data, int mat, double E, int numL) {
	cuDoubleComplex sigTfactors[4];
	int nuc = data->materials.mats_2d[mat* data->materials.pitch + threadIdx.x ];
	double * d_ptr = &(data->pseudo_K0RS_2d [nuc * numL]);

	// MicroScopic XS's to Calculate
	double sigT, sigA, sigF, sigE;

	// Calculate Window Index
	double spacing = 1.0 / data->n_windows[nuc];
	int window = (int) ( E / spacing );
	if( window == data->n_windows[nuc] )
		window--;

	// Calculate sigTfactors
	double phi;
	double sqrt_E = sqrt(E);
	for( int i = 0; i < numL; i++ ){
		phi = d_ptr[i] * sqrt_E;

		if( i == 1 )
			phi -= - atan( phi );
		else if( i == 2 )
			phi -= atan( 3.0 * phi / (3.0 - phi*phi));
		else if( i == 3 )
			phi -= atan(phi*(15.0-phi*phi)/(15.0-6.0*phi*phi));

		phi += phi;

		sigTfactors[i].x= cos(phi);
		sigTfactors[i].y= - sin(phi);
	}
	// Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
	Window w = data->windows_2d[nuc * data->pitch_windows + window];
	sigT = E * w.T;	sigA = E * w.A;	sigF = E * w.F;
	// Loop over Poles within window, add contributions
	cuDoubleComplex const1 = make_cuDoubleComplex(0, 1/E), const2 = make_cuDoubleComplex( sqrt_E, 0);
	for( int i = w.start; i < w.end; i++ )	{
		Pole pole = data->poles_2d[nuc * data->pitch_poles + i];
		cuDoubleComplex CDUM = cuCdiv( const1, cuCsub( pole.MP_EA, const2 ) );
		sigT += cuCreal( cuCmul( pole.MP_RT, cuCmul( CDUM, sigTfactors[pole.l_value] ) ) );
		sigA += cuCreal( cuCmul( pole.MP_RA, CDUM) );
		sigF += cuCreal( cuCmul( pole.MP_RF, CDUM) );
	}
	sigE = sigT - sigA;

	int base = 4 * threadIdx.x, idx = mat * data->materials.pitch + threadIdx.x;
	macro_xs[base] += sigT * data->materials.concs_2d[idx ];
	macro_xs[base+1] += sigA * data->materials.concs_2d[idx ];
	macro_xs[base+2] += sigF * data->materials.concs_2d[idx ];
	macro_xs[base+3] += sigE * data->materials.concs_2d[idx ];
}

// 	third driver:	use data_d
void calc_macro_xs_driver ( double * macro_xs, int mat, double E, Input input, CalcDataPtrs* data, CalcDataPtrs_d* data_d, cuDoubleComplex * sigTfactors ) {
	// zero out macro vector
	int num = data->materials.num_nucs[mat];
	double *macro_xs_d, *macro_xs_results;
	//	cuDoubleComplex * sigTfactors_d;
	assert (cudaMalloc((void **) &macro_xs_d, 4*sizeof(double)) == cudaSuccess);
	assert (cudaMalloc((void **) &macro_xs_results, 4*num*sizeof(double)) == cudaSuccess);
	//	assert (cudaMalloc((void **) &sigTfactors_d, input.numL*sizeof(cuDoubleComplex)) == cudaSuccess);
	assert(cudaMemset( macro_xs_d, 0, 4*sizeof(double) ) == cudaSuccess);
	assert(cudaMemset( macro_xs_results, 0, num*4*sizeof(double) ) == cudaSuccess);
	macro_kernel<<<1, num>>> (macro_xs_results, data_d, mat, E, input.numL);
	//	gpuErrchk( cudaGetLastError() );
	//        gpuErrchk( cudaDeviceSynchronize() );
	// for nuclide in mat
	add_macro_xs_four_threads<<<1,4 >>> (macro_xs_d, macro_xs_results, num ); 
	assert(cudaMemcpy( macro_xs, macro_xs_d, 4*sizeof(double),cudaMemcpyDeviceToHost) == cudaSuccess);
	//	assert(cudaMemcpy( sigTfactors, sigTfactors_d, 4*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost) == cudaSuccess);
	cudaFree(macro_xs_d);
	cudaFree(macro_xs_results);
}

__global__ void calc_sig_kernel ( cuDoubleComplex const1, cuDoubleComplex const2, Pole* poles, 
		int base, double* sigT, double* sigA, double* sigF, cuDoubleComplex* sigTfactors) {
	Pole pole = poles[blockIdx.x + base ];
	cuDoubleComplex CDUM = cuCdiv( const1, cuCsub( pole.MP_EA, const2 ) );
	sigT[blockIdx.x + base] = cuCreal( cuCmul( pole.MP_RT, cuCmul( CDUM, sigTfactors[pole.l_value] ) ) );
	sigA[blockIdx.x + base] = cuCreal( cuCmul( pole.MP_RA, CDUM) );
	sigF[blockIdx.x + base] = cuCreal( cuCmul( pole.MP_RF, CDUM) );
}

__global__ void calc_sig_dd_kernel ( cuDoubleComplex const1, cuDoubleComplex const2, CalcDataPtrs_d* data_d, int nuc, 
		int base, double* sigT, double* sigA, double* sigF, cuDoubleComplex* sigTfactors) {
	Pole pole = data_d->poles_2d[ nuc * data_d->pitch_poles +  blockIdx.x + base ];
	cuDoubleComplex CDUM = cuCdiv( const1, cuCsub( pole.MP_EA, const2 ) );
	sigT[blockIdx.x + base] = cuCreal( cuCmul( pole.MP_RT, cuCmul( CDUM, sigTfactors[pole.l_value] ) ) );
	sigA[blockIdx.x + base] = cuCreal( cuCmul( pole.MP_RA, CDUM) );
	sigF[blockIdx.x + base] = cuCreal( cuCmul( pole.MP_RF, CDUM) );
}

__global__ void sum_sigs ( double* sigTs, double* sigAs, double* sigFs, double* sigT, double* sigA, double* sigF, 
		int tpb, int len ) {
	int i, j;
	for ( i = 0; i < len; i += tpb ) {
		if ( ( j = i + threadIdx.x ) < len ) {
			atomicAdd(sigA, sigAs[j]);
			atomicAdd(sigT, sigTs[j]);
			atomicAdd(sigF, sigFs[j]);
		}
	}
}

//	second try of CUDA adaptation of RSBench; data is copied to device every time called
void calc_sig_driver ( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, cuDoubleComplex * sigTfactors ) {
	// MicroScopic XS's to Calculate
	double sigT, sigA, sigF, sigE;

	// Calculate Window Index
	double spacing = 1.0 / data.n_windows[nuc];
	int window = (int) ( E / spacing );
	if( window == data.n_windows[nuc] )
		window--;

	// Calculate sigTfactors
	calculate_sig_T_sim ( E, input.numL, data.pseudo_K0RS[nuc], sigTfactors );

	// Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
	Window w = data.windows[nuc][window];
	sigT = E * w.T;	sigA = E * w.A;	sigF = E * w.F;
	// Loop over Poles within window, add contributions
	cuDoubleComplex const1 = make_cuDoubleComplex(0, 1/E), const2 = make_cuDoubleComplex( sqrt(E), 0);
	int num = w.end - w.start + 1;
	double* sigTs, *sigAs, *sigFs, *sigT_d, *sigA_d, *sigF_d; 
	Pole* poles_d; 
	cuDoubleComplex* sigTfactors_d;
	/* allocate memory on device */
	assert (cudaMalloc((void **) &sigTs, num*sizeof(double)) == cudaSuccess);
	assert (cudaMalloc((void **) &sigAs, num*sizeof(double)) == cudaSuccess);
	assert (cudaMalloc((void **) &sigFs, num*sizeof(double)) == cudaSuccess);
	assert (cudaMalloc((void **) &sigT_d, sizeof(double)) == cudaSuccess);
	assert (cudaMalloc((void **) &sigA_d, sizeof(double)) == cudaSuccess);
	assert (cudaMalloc((void **) &sigF_d, sizeof(double)) == cudaSuccess);
	assert (cudaMalloc((void **) &poles_d, num*sizeof(Pole)) == cudaSuccess);
	assert (cudaMalloc((void **) &sigTfactors_d, num*sizeof(cuDoubleComplex)) == cudaSuccess);
	//	assert(cudaMemcpy( sigTfactors_d, sigTfactors, input.numL*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice) == cudaSuccess);
	assert(cudaMemcpy( poles_d, data.poles[nuc], sizeof(Pole),cudaMemcpyHostToDevice) == cudaSuccess);
	calc_sig_kernel<<<num, 1>>> ( const1, const2, poles_d, 
			w.start, sigTs, sigAs, sigFs, sigTfactors_d );
	sum_sigs<<<1, 512>>> ( sigTs, sigAs, sigFs, sigT_d, sigA_d, sigF_d, 512, num );
	assert(cudaMemcpy( &sigT, sigT_d, sizeof(double),cudaMemcpyDeviceToHost) == cudaSuccess);
	assert(cudaMemcpy( &sigA, sigA_d, sizeof(double),cudaMemcpyDeviceToHost) == cudaSuccess);
	assert(cudaMemcpy( &sigF, sigF_d, sizeof(double),cudaMemcpyDeviceToHost) == cudaSuccess);

	cudaFree( sigTs );  cudaFree( sigAs );  cudaFree( sigFs );  cudaFree( poles_d); cudaFree( sigTfactors_d);
	cudaFree( sigT_d );  cudaFree( sigA_d );  cudaFree( sigF_d );  

	sigE = sigT - sigA;

	micro_xs[0] = sigT;	micro_xs[1] = sigA;	micro_xs[2] = sigF;	micro_xs[3] = sigE;
}

//	second try of CUDA adaptation of RSBench; data is copied to device every time called
void calc_sig_dd_driver ( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, 
		CalcDataPtrs_d* data_d, cuDoubleComplex * sigTfactors ) {
	// MicroScopic XS's to Calculate
	double sigT, sigA, sigF, sigE;

	// Calculate Window Index
	double spacing = 1.0 / data.n_windows[nuc];
	int window = (int) ( E / spacing );
	if( window == data.n_windows[nuc] )
		window--;

	// Calculate sigTfactors
	calculate_sig_T_sim ( E, input.numL, data.pseudo_K0RS[nuc], sigTfactors );

	// Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
	Window w = data.windows[nuc][window];
	sigT = E * w.T;	sigA = E * w.A;	sigF = E * w.F;
	// Loop over Poles within window, add contributions
	cuDoubleComplex const1 = make_cuDoubleComplex(0, 1/E), const2 = make_cuDoubleComplex( sqrt(E), 0);
	int num = w.end - w.start + 1;
	double* sigTs, *sigAs, *sigFs, *sigT_d, *sigA_d, *sigF_d; 
	//	Pole* poles_d; 
	cuDoubleComplex* sigTfactors_d;
	/* allocate memory on device */
	assert (cudaMalloc((void **) &sigTs, num*sizeof(double)) == cudaSuccess);
	assert (cudaMalloc((void **) &sigAs, num*sizeof(double)) == cudaSuccess);
	assert (cudaMalloc((void **) &sigFs, num*sizeof(double)) == cudaSuccess);
	assert (cudaMalloc((void **) &sigT_d, sizeof(double)) == cudaSuccess);
	assert (cudaMalloc((void **) &sigA_d, sizeof(double)) == cudaSuccess);
	assert (cudaMalloc((void **) &sigF_d, sizeof(double)) == cudaSuccess);
	//	assert (cudaMalloc((void **) &poles_d, num*sizeof(Pole)) == cudaSuccess);
	assert (cudaMalloc((void **) &sigTfactors_d, num*sizeof(cuDoubleComplex)) == cudaSuccess);
	//	assert(cudaMemcpy( poles_d, data.poles[nuc], sizeof(Pole),cudaMemcpyHostToDevice) == cudaSuccess);
	calc_sig_dd_kernel<<<num, 1>>> ( const1, const2, data_d, nuc, 
			w.start, sigTs, sigAs, sigFs, sigTfactors_d );
	sum_sigs<<<1, 512>>> ( sigTs, sigAs, sigFs, sigT_d, sigA_d, sigF_d, 512, num );
	assert(cudaMemcpy( &sigT, sigT_d, sizeof(double),cudaMemcpyDeviceToHost) == cudaSuccess);
	assert(cudaMemcpy( &sigA, sigA_d, sizeof(double),cudaMemcpyDeviceToHost) == cudaSuccess);
	assert(cudaMemcpy( &sigF, sigF_d, sizeof(double),cudaMemcpyDeviceToHost) == cudaSuccess);

	cudaFree( sigTs );  cudaFree( sigAs );  cudaFree( sigFs );  //cudaFree( poles_d); 
	cudaFree( sigTfactors_d); cudaFree( sigT_d );  cudaFree( sigA_d );  cudaFree( sigF_d );  

	sigE = sigT - sigA;

	micro_xs[0] = sigT;	micro_xs[1] = sigA;	micro_xs[2] = sigF;	micro_xs[3] = sigE;
}

//	first version of CUDA addption - only the four-iteration loop; each time called, data is copied to the device
void calculate_micro_xs_driver( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, cuDoubleComplex * sigTfactors) {
	// MicroScopic XS's to Calculate
	double sigT, sigA, sigF, sigE;
	double* data_d;
	cuDoubleComplex* cudcomp_d;

	// Calculate Window Index
	double spacing = 1.0 / data.n_windows[nuc];
	int window = (int) ( E / spacing );
	if( window == data.n_windows[nuc] )
		window--;

	/* allocate memory on device */
	assert (cudaMalloc((void **) &data_d, input.numL*sizeof(double)) == cudaSuccess);
	assert (cudaMalloc((void **) &cudcomp_d, input.numL*sizeof(cuDoubleComplex)) == cudaSuccess);
	/* copy host data to device pointers */
	assert(cudaMemcpy(data_d, data.pseudo_K0RS[nuc], input.numL*sizeof(double),cudaMemcpyHostToDevice) == cudaSuccess);
	// Calculate sigTfactors
	calc_sig_T_sim_kernel<<<1, input.numL>>> ( E, input.numL, data_d, cudcomp_d);
	assert(cudaMemcpy( sigTfactors, cudcomp_d, input.numL*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost) == cudaSuccess);
	cudaFree( cudcomp_d );  
	cudaFree( data_d );  
	// Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
	Window w = data.windows[nuc][window];
	sigT = E * w.T;	sigA = E * w.A;	sigF = E * w.F;
	// Loop over Poles within window, add contributions
	cuDoubleComplex const1 = make_cuDoubleComplex(0, 1/E), const2 = make_cuDoubleComplex( sqrt(E), 0);
	for( int i = w.start; i < w.end; i++ ) {
		Pole pole = data.poles[nuc][i];
		cuDoubleComplex CDUM = cuCdiv( const1, cuCsub( pole.MP_EA, const2 ) );
		sigT += cuCreal( cuCmul( pole.MP_RT, cuCmul( CDUM, sigTfactors[pole.l_value] ) ) );
		sigA += cuCreal( cuCmul( pole.MP_RA, CDUM) );
		sigF += cuCreal( cuCmul( pole.MP_RF, CDUM) );
	}

	sigE = sigT - sigA;

	micro_xs[0] = sigT; micro_xs[1] = sigA; micro_xs[2] = sigF; micro_xs[3] = sigE;
}

//	improved first version of CUDA addption - only the four-iteration loop; data is copied to the device in the beginning
void calculate_micro_xs_dd_driver( double * micro_xs, int nuc, double E, Input input, 
		CalcDataPtrs data, CalcDataPtrs_d* data_d, cuDoubleComplex * sigTfactors) {
	// MicroScopic XS's to Calculate
	double sigT, sigA, sigF, sigE;
	cuDoubleComplex* cudcomp_d;

	// Calculate Window Index
	double spacing = 1.0 / data.n_windows[nuc];
	int window = (int) ( E / spacing );
	if( window == data.n_windows[nuc] )
		window--;

	/* allocate memory on device */
	assert (cudaMalloc((void **) &cudcomp_d, input.numL*sizeof(cuDoubleComplex)) == cudaSuccess);
	// Calculate sigTfactors
	calc_sig_T_sim_dd_kernel<<<1, input.numL>>> ( E, input.numL, data_d, nuc * input.numL, cudcomp_d);
	//	assert(cudaMemcpy( sigTfactors, cudcomp_d, input.numL*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost) == cudaSuccess);
	cudaFree( cudcomp_d );  
	// Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
	Window w = data.windows[nuc][window];
	sigT = E * w.T;
	sigA = E * w.A;
	sigF = E * w.F;
	// Loop over Poles within window, add contributions
	cuDoubleComplex const1 = make_cuDoubleComplex(0, 1/E), const2 = make_cuDoubleComplex( sqrt(E), 0);
	for( int i = w.start; i < w.end; i++ ) {
		Pole pole = data.poles[nuc][i];
		cuDoubleComplex CDUM = cuCdiv( const1, cuCsub( pole.MP_EA, const2 ) );
		sigT += cuCreal( cuCmul( pole.MP_RT, cuCmul( CDUM, sigTfactors[pole.l_value] ) ) );
		sigA += cuCreal( cuCmul( pole.MP_RA, CDUM) );
		sigF += cuCreal( cuCmul( pole.MP_RF, CDUM) );
	}

	sigE = sigT - sigA;

	micro_xs[0] = sigT;
	micro_xs[1] = sigA;
	micro_xs[2] = sigF;
	micro_xs[3] = sigE;
}

int dotp_driver(int NTPB){
	int nDevices;
	cudaGetDeviceCount(&nDevices);
	for (int i = 0; i < nDevices; i++) {
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, i);
		printf("Device Number: %d\n", i);
		printf("  Device name: %s\n", prop.name);
		printf("  Memory Clock Rate (KHz): %d\n",
				prop.memoryClockRate);
		printf("  Memory Bus Width (bits): %d\n",
				prop.memoryBusWidth);
		printf("  Maximum number of threads per block: %d\n",
				prop.maxThreadsPerBlock);
		printf("  Maximum size of each dimension of a block: %d %d %d\n",
				prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
		printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
				2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
	}
	return 0;
}
