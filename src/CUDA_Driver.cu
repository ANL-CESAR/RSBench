#include "rsbench.h"
#include "My_Stats.h"

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

__global__ void calc_sig_kernel ( double E, int num_iter, 
		const double* data, cuDoubleComplex* gTfactors) {
/*	Pole pole = data.poles[nuc][i];
	cuDoubleComplex CDUM = -(0.0 - 1.0 * _Complex_I ) / ( pole.MP_EA - sqrt(E) ) / E;
	sigT += creal( pole.MP_RT * CDUM * sigTfactors[pole.l_value] );
	sigA += creal( pole.MP_RA * CDUM);
	sigF += creal( pole.MP_RF * CDUM);*/
}

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
	sigT = E * w.T;
	sigA = E * w.A;
	sigF = E * w.F;
	// Loop over Poles within window, add contributions
	cuDoubleComplex const1 = make_cuDoubleComplex(0, 1/E), const2 = make_cuDoubleComplex( sqrt(E), 0);
	for( int i = w.start; i < w.end; i++ )
	{
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

void calculate_micro_xs_driver( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, cuDoubleComplex * sigTfactors)
{
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
	sigT = E * w.T;
	sigA = E * w.A;
	sigF = E * w.F;
	// Loop over Poles within window, add contributions
	cuDoubleComplex const1 = make_cuDoubleComplex(0, 1/E), const2 = make_cuDoubleComplex( sqrt(E), 0);
	for( int i = w.start; i < w.end; i++ )
	{
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
