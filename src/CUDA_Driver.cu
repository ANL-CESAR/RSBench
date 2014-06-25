#include "rsbench.h"
#include "My_Stats.h"

const unsigned int NUM_ITER = 100;
/* vector lengths */
unsigned int sizes[] = {16384, 65536, 262144, 1048576, 4194304, 16777216, 33554432, 50331648, 67108864, 268435456};
//10000, 100000, 1000000, 5000000, 10000000, 50000000, 80000000, 100000000};//, 500000000, 1000000000};

/* performs the dot product of the vectors a and b
   on the gpu and returns the result in global memory
   location indicated by c. n is vector length
 */
__global__ void dotp(int *a, int *b, int *c, int n){
	int i;
	int iglob = threadIdx.x + blockIdx.x*blockDim.x; 
	int iloc  = threadIdx.x                        ;
	/* each block's pairwise products are stored in this temporary array */
	//  __shared__ int block_cache[NTPB]; 
	//note that our block_cache array can be sized statically with a constant
	//or dynamically with the size in bytes passed in as the third argument
	//to the kernel launch specificer
	extern __shared__ int block_cache[]; 

	if (iglob < n)
		block_cache[iloc] = a[iglob]*b[iglob];
	else
		block_cache[iloc] = 0;

	__syncthreads();

	/* on the "master thread" of each block" sum the pairwise products
	   on that block into the block's portion of the global sum */
	if (iloc == 0){
		int sum = 0;
		//		for (i=0;i<NTPB;++i)
		for (i=0;i< blockDim.x;++i)
			sum += block_cache[i];
		atomicAdd(c,sum);  /* now write safely to global memory */
		// *c += sum;         // don't do this -- not safe! 
	}

}

int dotp_driver(int NTPB){
	int *a,   *b,   *c;       /* host pointers */
	int *a_d, *b_d, *c_d;     /* device pointers */
	int i, j, k;			/* vector length */
	/* Number of Threads Per Block */
	assert( NTPB > 0 );
	cudaEvent_t start, stop;  /* timers */
	float time;

  	double trials [NUM_ITER];
	double* stats=NULL;

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

	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	printf("Number of trials: %u\t Num of threads per block: %i\n", NUM_ITER, NTPB);
	printf("Vector Size,DotProd,Min,25%%,Median,75%%,Max,Mean,Variance,Stdev\n");
	for ( j = 0; j < sizeof (sizes)/sizeof(unsigned int); j ++ ) {
		/* allocate host memory */
		assert (cudaMallocHost((void **) &a, sizes[j]*sizeof(int)) == cudaSuccess);
		assert (cudaMallocHost((void **) &b, sizes[j]*sizeof(int)) == cudaSuccess);
		assert (cudaMallocHost((void **) &c, 1*sizeof(int)) == cudaSuccess);

		/* initialize vectors to some simple values and initialize
		   scalar result of dot product to zero */
		for (i=0;i<sizes[j];++i){
			a[i] =  1;
			b[i] = -2;
		}

		/* allocate memory on device */
		assert (cudaMalloc((void **) &a_d, sizes[j]*sizeof(int)) == cudaSuccess);
		assert (cudaMalloc((void **) &b_d, sizes[j]*sizeof(int)) == cudaSuccess);
		assert (cudaMalloc((void **) &c_d, 1*sizeof(int)) == cudaSuccess);

		/* copy host data to device pointers */
		assert(cudaMemcpy(a_d,a,sizes[j]*sizeof(int),cudaMemcpyHostToDevice) == cudaSuccess);
		assert(cudaMemcpy(b_d,b,sizes[j]*sizeof(int),cudaMemcpyHostToDevice) == cudaSuccess);

		for (k = 0; k < NUM_ITER; k++ ){
			*c = 0;
			assert(cudaMemcpy(c_d,c,1*sizeof(int),cudaMemcpyHostToDevice) == cudaSuccess);

			/* one way to set kernel launch values is to fix Number of Threads Per Block (NTPB)
			   as a constant and then calculate # of blocks based on input problem size n. Below
			   is a simple formula that carries out this calculation. THere are other ways to do
			   this that are just as good */

			/* launch and time kernel code */
			cudaEventRecord( start, 0 );  

			dotp<<<(sizes[j])/NTPB,NTPB,NTPB*sizeof(int)>>>(a_d,b_d,c_d,sizes[j]);

			cudaEventRecord( stop, 0 );
			cudaEventSynchronize( stop );
			cudaEventElapsedTime( &time, start, stop );
			trials[k] = time;
			assert(cudaMemcpy(c,c_d,1*sizeof(int),cudaMemcpyDeviceToHost) == cudaSuccess);
		}
		my_stats ( trials, &NUM_ITER, &stats);
		printf ("%u,%i,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", sizes[j], *c, stats[0], stats[1], 
			stats[2], stats[3], stats[4], stats[5], stats[6], stats[7]);
//		printf("size of vectors: %u\tvalue on device:%d\ttime elapsed: %f(ms)\n", sizes[j], *c, time);

		cudaFree(a_d);  cudaFree(b_d);  cudaFree(c_d);
		cudaFreeHost(a);  cudaFreeHost(b);  cudaFreeHost(c);
	}
	cudaEventDestroy( start );
	cudaEventDestroy( stop );
	return 0;
}
