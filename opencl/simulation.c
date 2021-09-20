#include "rsbench.h"

////////////////////////////////////////////////////////////////////////////////////
// BASELINE FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////
// All "baseline" code is at the top of this file. The baseline code is a simple
// implementation of the algorithm, with only minor CPU optimizations in place.
// Following these functions are a number of optimized variants,
// which each deploy a different combination of optimizations strategies. By
// default, RSBench will only run the baseline implementation. Optimized variants
// must be specifically selected using the "-k <optimized variant ID>" command
// line argument.
////////////////////////////////////////////////////////////////////////////////////

unsigned long long run_event_based_simulation(Input in, SimulationData SD, double * sim_runtime)
{

	double start = get_time();

	int * verification_array_host = (int *) malloc( in.lookups * sizeof(int));
	
	////////////////////////////////////////////////////////////////////////////////
	// OpenCL Boilerplate Setup
	////////////////////////////////////////////////////////////////////////////////

	// Let's start setting up our openCL boilerplate
	// Load the kernel source code into the array source_str
	FILE *fp;
	char *source_str;
	size_t source_size;

	fp = fopen("kernel.cl", "r");
	if (!fp) {
		fprintf(stderr, "Failed to load kernel.\n");
		exit(1);
	}
	source_str = (char*) malloc(MAX_SOURCE_SIZE);
	source_size = fread( source_str, 1, MAX_SOURCE_SIZE, fp);
	fclose( fp );

  OpenCLInfo CL = initialize_device(in.platform_id, in.device_id);
  cl_device_id device_id = CL.device_id;
  cl_context context = CL.context;
  cl_command_queue command_queue = CL.command_queue;
	
  printf("Initializing OpenCL data structures and JIT compiling kernel...\n");

  ////////////////////////////////////////////////////////////////////////////////
  // OpenCL Move Memory To Device Buffers
  ////////////////////////////////////////////////////////////////////////////////

  // Create memory buffers on the device for each vector and move data over
  size_t sz = SD.length_num_nucs * sizeof(int);
  cl_int ret;
  cl_mem num_nucs_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  sz, NULL, &ret);
  check(ret);
  ret = clEnqueueWriteBuffer(command_queue, num_nucs_d, CL_TRUE, 0, sz, SD.num_nucs, 0, NULL, NULL);
  check(ret);

  sz = SD.length_concs * sizeof(double);
  cl_mem concs_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  sz, NULL, &ret);
  check(ret);
  ret = clEnqueueWriteBuffer(command_queue, concs_d, CL_TRUE, 0, sz, SD.concs, 0, NULL, NULL);
  check(ret);

  sz = SD.length_mats * sizeof(int);
  cl_mem mats_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  sz, NULL, &ret);
  check(ret);
  ret = clEnqueueWriteBuffer(command_queue, mats_d, CL_TRUE, 0, sz, SD.mats, 0, NULL, NULL);
  check(ret);

  sz = in.lookups * sizeof(int);
  cl_mem verification_array = clCreateBuffer(context, CL_MEM_READ_WRITE,  sz, NULL, &ret);
  check(ret);

  sz = SD.length_n_windows * sizeof(int);
  cl_mem n_windows_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  sz, NULL, &ret);
  check(ret);
  ret = clEnqueueWriteBuffer(command_queue, n_windows_d, CL_TRUE, 0, sz, SD.n_windows, 0, NULL, NULL);
  check(ret);

  sz = SD.length_poles * sizeof(Pole);
  cl_mem poles_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  sz, NULL, &ret);
  check(ret);
  ret = clEnqueueWriteBuffer(command_queue, poles_d, CL_TRUE, 0, sz, SD.poles, 0, NULL, NULL);
  check(ret);

  sz = SD.length_windows * sizeof(Window);
  cl_mem windows_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  sz, NULL, &ret);
  check(ret);
  ret = clEnqueueWriteBuffer(command_queue, windows_d, CL_TRUE, 0, sz, SD.windows, 0, NULL, NULL);
  check(ret);

  sz = SD.length_pseudo_K0RS * sizeof(double);
  cl_mem pseudo_K0RS_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  sz, NULL, &ret);
  check(ret);
  ret = clEnqueueWriteBuffer(command_queue, pseudo_K0RS_d, CL_TRUE, 0, sz, SD.pseudo_K0RS, 0, NULL, NULL);
  check(ret);

  ////////////////////////////////////////////////////////////////////////////////
  // OpenCL Build and Intiailize Kernel Program
  ////////////////////////////////////////////////////////////////////////////////

  // Create a program from the kernel source
  cl_program program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);
  check(ret);

  // Build the program
  ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
  check(ret);

  printCompilerError( program, device_id );

  // Create the OpenCL kernel
  cl_kernel kernel = clCreateKernel(program, "macro_xs_lookup_kernel", &ret);
  check(ret);

  // Set the arguments of the kernel
  ret = clSetKernelArg(kernel, 0, sizeof(Input), (void *)&in);
  check(ret);
  ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&num_nucs_d);
  check(ret);
  ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mats_d);
  check(ret);
  ret = clSetKernelArg(kernel, 3, sizeof(int), (void *)&SD.max_num_nucs);
  check(ret);
  ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&concs_d);
  check(ret);
  ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&n_windows_d);
  check(ret);
  ret = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&pseudo_K0RS_d);
  check(ret);
  ret = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&windows_d);
  check(ret);
  ret = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&poles_d);
  check(ret);
  ret = clSetKernelArg(kernel, 9, sizeof(int), (void *)&SD.max_num_windows);
  check(ret);
  ret = clSetKernelArg(kernel, 10, sizeof(int), (void *)&SD.max_num_poles);
  check(ret);
  ret = clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *)&verification_array);
  check(ret);

  double stop = get_time();
  printf("OpenCL initialization time: %.3lf seconds\n", stop-start);
  start = stop;

  ////////////////////////////////////////////////////////////////////////////////
  // Run Simulation Kernel
  ////////////////////////////////////////////////////////////////////////////////
	
  border_print();
	center_print("SIMULATION", 79);
	border_print();

  printf("Running event based simulation...\n");

	// Execute the OpenCL kernel on the list
	size_t global_item_size = in.lookups; // Process the entire lists
	size_t local_item_size = 256; // Divide work items into groups

  // Add extra work items if global size not evenly divisible by local size
  if( in.lookups % local_item_size != 0 && in.lookups > local_item_size )
    global_item_size = ((in.lookups / local_item_size) + 1) * local_item_size;

  ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL);
  check(ret);

  ////////////////////////////////////////////////////////////////////////////////
  // Retrieve verification data from device and reduce it
  ////////////////////////////////////////////////////////////////////////////////

  // Read the memory buffer C on the device to the local variable C
  ret = clEnqueueReadBuffer(command_queue, verification_array, CL_TRUE, 0, in.lookups * sizeof(int), verification_array_host, 0, NULL, NULL);
  check(ret);

  printf("Reducing verification value...\n");

  unsigned long long verification = 0;

  for( int l = 0; l < in.lookups; l++ )
    verification += verification_array_host[l];

  stop = get_time();
  *sim_runtime = stop-start;
  printf("Simulation + Verification Reduction Runtime: %.3lf seconds\n", *sim_runtime);

  ////////////////////////////////////////////////////////////////////////////////
  // OpenCL cleanup
  ////////////////////////////////////////////////////////////////////////////////

  ret = clFlush(command_queue);
  check(ret);
  ret = clFinish(command_queue);
  check(ret);
  ret = clReleaseKernel(kernel);
  check(ret);
  ret = clReleaseProgram(program);
  check(ret);
  ret = clReleaseMemObject(num_nucs_d);
  check(ret);
  ret = clReleaseMemObject(mats_d);
  check(ret);
  ret = clReleaseMemObject(n_windows_d);
  check(ret);
  ret = clReleaseMemObject(poles_d);
  check(ret);
  ret = clReleaseMemObject(windows_d);
  check(ret);
  ret = clReleaseMemObject(pseudo_K0RS_d);
  check(ret);
  ret = clReleaseMemObject(verification_array);
  check(ret);
  ret = clReleaseCommandQueue(command_queue);
  check(ret);
  ret = clReleaseContext(context);
  check(ret);

  return verification;
}

double LCG_random_double(uint64_t * seed)
{
  const uint64_t m = 9223372036854775808ULL; // 2^63
  const uint64_t a = 2806196910506780709ULL;
  const uint64_t c = 1ULL;
  *seed = (a * (*seed) + c) % m;
  return (double) (*seed) / (double) m;
}	

uint64_t LCG_random_int(uint64_t * seed)
{
  const uint64_t m = 9223372036854775808ULL; // 2^63
  const uint64_t a = 2806196910506780709ULL;
  const uint64_t c = 1ULL;
  *seed = (a * (*seed) + c) % m;
  return *seed;
}	

RSComplex c_mul( RSComplex A, RSComplex B)
{
  double a = A.r;
  double b = A.i;
  double c = B.r;
  double d = B.i;
  RSComplex C;
  C.r = (a*c) - (b*d);
  C.i = (a*d) + (b*c);
  return C;
}
