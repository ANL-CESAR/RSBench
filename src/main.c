#include "rsbench.h"

int main(int argc, char * argv[])
{
	// =====================================================================
	// Initialization & Command Line Read-In
	// =====================================================================

	int version = 5;
	unsigned long long vhash = 0;
	struct timeval start, end;
	double wall_time;
	
	// These are the sum of the function values, evaluated in kernel.
	// They are the cumulative result of the random lookups.

	// Vectors for sums of F(x_i).  Dimensions will be V_sums[0:outer_dim].
	// In kernel, Each outer unit j will reduce is results to V_sum[i].
	// In main, we will need to reduce V_sums to get a single V_sum
	double *V_sums;

	// Sum of all F(x_i) from kernel.  Outside of kernel, V_sums will be reduced
	// to get V_sum
	double V_sum[4] = {0, 0, 0, 0};

	// =====================================================================
	// OCCA declarations
	// =====================================================================

	occaKernelInfo lookupInfo;
	occaKernel lookup_kernel, lookup_kernel_doppler;
	occaDevice device;

	occaMemory dev_poles, dev_windows, dev_pseudo_K0RS, dev_n_poles,
	     dev_n_windows, dev_num_nucs, dev_mats, dev_mats_idx, dev_concs,
	     dev_poles_idx, dev_windows_idx;

	occaMemory dev_V_sums;

	// =====================================================================
	// Read command-line input
	// =====================================================================
	Input input = read_CLI( argc, argv );

	logo(version);
	center_print("INPUT SUMMARY", 79);
	border_print();
	print_input_summary(input);

	// =====================================================================
	// Initialize OCCA
	// =====================================================================
	lookupInfo = occaCreateKernelInfo();
	occaKernelInfoAddDefine(lookupInfo, "inner_dim", occaLong(input.inner_dim));
	occaKernelInfoAddDefine(lookupInfo, "outer_dim", occaLong(input.outer_dim));

	#ifdef VERIFICATION
	// occaKernelInfoAddDefine(lookupInfo, "VERIFICATION", occaInt(1));
	#endif

	device = occaGetDevice(input.device_info);

	//lookup_touch = occaBuildKernelFromSource(device,
	//     "lookup_kernel.okl", "lookup_touch", lookupInfo);
	lookup_kernel = occaBuildKernelFromSource(device,
	     input.kernel, "lookup_kernel", lookupInfo);
	lookup_kernel_doppler = occaBuildKernelFromSource(device,
	     input.kernel, "lookup_kernel_doppler", lookupInfo);

	// =====================================================================
	// Initialize RNG
	// =====================================================================
	#ifdef VERIFICATION
	srand(26);
	#else
	srand(time(NULL));
	#endif

	// =====================================================================
	// Prepare Pole Paremeter Grids
	// =====================================================================
	border_print();
	center_print("INITIALIZATION", 79);
	border_print();
	
	gettimeofday(&start, NULL);
	
	// Allocate & fill energy grids
	printf("Generating resonance distributions...\n");
	int * n_poles = generate_n_poles( input );

	// Allocate & fill Window grids
	printf("Generating window distributions...\n");
	int * n_windows = generate_n_windows( input );

	// Get material data
	printf("Loading Hoogenboom-Martin material data...\n");
	Materials materials = get_materials( input ); 

	// Prepare full resonance grid
	printf("Generating resonance parameter grid...\n");
	Pole ** poles = generate_poles( input, n_poles );

	// Prepare full Window grid
	printf("Generating window parameter grid...\n");
	Window ** windows = generate_window_params( input, n_windows, n_poles);

	// Prepare 0K Resonances
	printf("Generating 0K l_value data...\n");
	double ** pseudo_K0RS = generate_pseudo_K0RS( input );

	// Get poles indices
	int * poles_idx = load_poles_idx(input.n_nuclides, n_poles);

	// Get windows indices
	int * windows_idx = load_windows_idx(input.n_nuclides, n_windows);

	CalcDataPtrs data;
	data.n_poles = n_poles;
	data.n_windows = n_windows;
	data.materials = materials;
	data.poles = poles;
	data.windows = windows;
	data.pseudo_K0RS = pseudo_K0RS;
	data.poles_idx = poles_idx;
	data.windows_idx = windows_idx;

	// Prepare verification arrays
	V_sums = (double *) calloc( 4 * input.lookups, sizeof(double) );

	gettimeofday(&end, NULL);
	wall_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1000000.;
	printf("Initialization Complete. (%.3lf seconds)\n", wall_time);
	
	// =====================================================================
	// OCCA mallocs and memcopies
	// =====================================================================

	printf("Allocating and copying to device memory...\n");
	// REMEMBER: memcopy is part of malloc (last arg gets copied to device)

	if (strcasecmp(input.mode, "OpenMP") == 0) {
	dev_poles = occaDeviceWrapMemory(device, data.poles[0], input.n_nuclides*input.avg_n_poles*sizeof(Pole));
	dev_windows = occaDeviceWrapMemory(device, data.windows[0], input.n_nuclides*input.avg_n_windows*sizeof(Window));
	dev_pseudo_K0RS = occaDeviceWrapMemory(device, data.pseudo_K0RS[0], input.n_nuclides*input.numL*sizeof(double));
	dev_n_poles = occaDeviceWrapMemory(device, data.n_poles, input.n_nuclides*sizeof(int));
	dev_n_windows = occaDeviceWrapMemory(device, data.n_windows, input.n_nuclides*sizeof(int));
	dev_poles_idx = occaDeviceWrapMemory(device, data.poles_idx, input.n_nuclides*sizeof(int));
	dev_windows_idx = occaDeviceWrapMemory(device, data.windows_idx, input.n_nuclides*sizeof(int));
	dev_num_nucs = occaDeviceWrapMemory(device, materials.num_nucs, 12*sizeof(int));
	dev_mats = occaDeviceWrapMemory(device, materials.mats, materials.mats_sz*sizeof(int));
	dev_mats_idx = occaDeviceWrapMemory(device, materials.mats_idx, 12*sizeof(int));
	dev_concs = occaDeviceWrapMemory(device, materials.concs, materials.mats_sz*sizeof(double));
	dev_V_sums = occaDeviceWrapMemory(device, V_sums, 4*input.lookups*sizeof(double));
	}
	else{
//	dev_poles = occaDeviceMalloc(device, input.n_nuclides*input.avg_n_poles*sizeof(Pole), NULL);
//	dev_windows = occaDeviceMalloc(device, input.n_nuclides*input.avg_n_windows*sizeof(Window), NULL);
//	dev_pseudo_K0RS = occaDeviceMalloc(device, input.n_nuclides*input.numL*sizeof(double), NULL);
	dev_poles = occaDeviceMalloc(device, input.n_nuclides*input.avg_n_poles*sizeof(Pole), data.poles[0]);
	dev_windows = occaDeviceMalloc(device, input.n_nuclides*input.avg_n_windows*sizeof(Window), data.windows[0]);
	dev_pseudo_K0RS = occaDeviceMalloc(device, input.n_nuclides*input.numL*sizeof(double), data.pseudo_K0RS[0]);
	dev_n_poles = occaDeviceMalloc(device, input.n_nuclides*sizeof(int), data.n_poles);
	dev_n_windows = occaDeviceMalloc(device, input.n_nuclides*sizeof(int), data.n_windows);
	dev_poles_idx = occaDeviceMalloc(device, input.n_nuclides*sizeof(int), data.poles_idx);
	dev_windows_idx = occaDeviceMalloc(device, input.n_nuclides*sizeof(int), data.windows_idx);
	dev_num_nucs = occaDeviceMalloc(device, 12*sizeof(int), materials.num_nucs);
	dev_mats = occaDeviceMalloc(device, materials.mats_sz*sizeof(int), materials.mats);
	dev_mats_idx = occaDeviceMalloc(device, 12*sizeof(int), materials.mats_idx);
	dev_concs = occaDeviceMalloc(device, materials.mats_sz*sizeof(double), materials.concs);
	dev_V_sums = occaDeviceMalloc(device, 4*input.lookups*sizeof(double), V_sums);

	// Call kernel to apply "proper" first-touch on large arrays
//	occaKernelRun(lookup_touch,
//	     dev_poles,
//	     dev_windows,
//	     dev_pseudo_K0RS,
//	     occaInt(input.n_nuclides),
//	     occaInt(input.avg_n_poles),
//	     occaInt(input.avg_n_windows),
//	     occaInt(input.numL));

	// Properly initialize arrays
//	occaCopyPtrToMem(dev_poles, poles[0], occaAutoSize, occaNoOffset);
//	occaCopyPtrToMem(dev_windows, windows[0], occaAutoSize, occaNoOffset);
//	occaCopyPtrToMem(dev_pseudo_K0RS, pseudo_K0RS[0], occaAutoSize, occaNoOffset);
	}

	// =====================================================================
	// Cross Section (XS) Parallel Lookup Simulation Begins
	// =====================================================================
	border_print();
	center_print("SIMULATION", 79);
	border_print();
	
	printf("Beginning Kernel...\n");
	printf("Calculating XS's...\n");
	
	occaDeviceFinish(device);
	gettimeofday(&start, NULL);

	if(input.doppler == 1){
		occaKernelRun(lookup_kernel_doppler,
		     dev_poles,
		     dev_windows,
		     dev_pseudo_K0RS,
		     dev_n_poles,
		     dev_n_windows,
		     dev_poles_idx,
		     dev_windows_idx,
		     dev_num_nucs,
		     dev_mats,
		     dev_mats_idx,
		     dev_concs,
		     occaInt(input.lookups),
		     occaInt(input.n_nuclides),
		     occaInt(input.avg_n_poles),
		     occaInt(input.avg_n_windows),
		     occaInt(input.numL),
		     dev_V_sums);
	}
	else{
		occaKernelRun(lookup_kernel,
		     dev_poles,
		     dev_windows,
		     dev_pseudo_K0RS,
		     dev_n_poles,
		     dev_n_windows,
		     dev_poles_idx,
		     dev_windows_idx,
		     dev_num_nucs,
		     dev_mats,
		     dev_mats_idx,
		     dev_concs,
		     occaInt(input.lookups),
		     occaInt(input.n_nuclides),
		     occaInt(input.avg_n_poles),
		     occaInt(input.avg_n_windows),
		     occaInt(input.numL),
		     dev_V_sums);
	}

	occaDeviceFinish(device);
	gettimeofday(&end, NULL);
	wall_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1000000.;

	printf("Kernel complete.\n" );

	// Device-to-host memcopy
	printf("Copying from device memory...\n");
	occaCopyMemToPtr(V_sums, dev_V_sums, 4*input.lookups*sizeof(double), 0);

	// Reduce sums
	for(int i=0; i<input.lookups; i++){
		V_sum[0] += V_sums[4*i + 0];
		V_sum[1] += V_sums[4*i + 1];
		V_sum[2] += V_sums[4*i + 2];
		V_sum[3] += V_sums[4*i + 3];
	}

	const double V_total_sum = V_sum[0] + V_sum[1] + V_sum[2] + V_sum[3];

	printf("\nSimulation Complete.\n");

	// =====================================================================
	// Print / Save Results and Exit
	// =====================================================================
	border_print();
	center_print("RESULTS", 79);
	border_print();

	printf("Threads:     %d\n", input.nthreads);
	printf("Runtime:     %.3lf seconds\n", wall_time);
	printf("Lookups:     "); fancy_int(input.lookups);
	printf("Lookups/s:   "); fancy_int((double) input.lookups / wall_time);
	printf("Verification: %f\n", V_total_sum / input.lookups);
	#ifdef VERIFICATION
	printf("Verification checksum: %llu\n", vhash);
	#endif

	border_print();

	return 0;
}
