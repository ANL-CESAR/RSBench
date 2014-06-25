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

void calculate_micro_xs_driver( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, complex double * sigTfactors)
{
	// MicroScopic XS's to Calculate
	double sigT;
	double sigA;
	double sigF;
	double sigE;
	double* data_d;
	cuDoubleComplex* cudcomp_d, *cudcomp_h;

	// Calculate Window Index
	double spacing = 1.0 / data.n_windows[nuc];
	int window = (int) ( E / spacing );
	if( window == data.n_windows[nuc] )
		window--;

	/* allocate memory on device */
	assert (cudaMalloc((void **) &data_d, input.numL*sizeof(double)) == cudaSuccess);
	assert (cudaMalloc((void **) &cudcomp_d, input.numL*sizeof(cuDoubleComplex)) == cudaSuccess);
	assert (cudaMallocHost((void **) &cudcomp_h, input.numL*sizeof(cuDoubleComplex)) == cudaSuccess);
	/* copy host data to device pointers */
	assert(cudaMemcpy(data_d, data.pseudo_K0RS[nuc], input.numL*sizeof(double),cudaMemcpyHostToDevice) == cudaSuccess);
	// Calculate sigTfactors
	calc_sig_T_sim_kernel<<<1, input.numL>>> ( E, input.numL, data_d, cudcomp_d);
	assert(cudaMemcpy( cudcomp_h, cudcomp_d, input.numL*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost) == cudaSuccess);
	for ( int i = 0; i < input.numL; i ++) {
		sigTfactors[i] = cudcomp_h->x + cudcomp_h->y * _Complex_I;
	}
	cudaFree( cudcomp_d );  
	cudaFreeHost( cudcomp_h );  
	//printf ("Chicago!\n");
	// Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
	Window w = data.windows[nuc][window];
	sigT = E * w.T;
	sigA = E * w.A;
	sigF = E * w.F;
	// Loop over Poles within window, add contributions
	for( int i = w.start; i < w.end; i++ )
	{
		complex double PSIIKI;
		complex double CDUM;
		Pole pole = data.poles[nuc][i];
		PSIIKI = -(0.0 - 1.0 * _Complex_I ) / ( pole.MP_EA - sqrt(E) );
		CDUM = PSIIKI / E;
		sigT += creal( pole.MP_RT * CDUM * sigTfactors[pole.l_value] );
		sigA += creal( pole.MP_RA * CDUM);
		sigF += creal( pole.MP_RF * CDUM);
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
