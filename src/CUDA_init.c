#include "rsbench.h"

// generate CalcDataPtrs_d data on device from CalcDataPtrs* data on host
CalcDataPtrs_d* init_data ( Input input, CalcDataPtrs* data) {//, Input* input_d ) {
	CalcDataPtrs_d* data_d;
	assert(cudaMalloc((void**) &data_d, sizeof(CalcDataPtrs_d))==cudaSuccess);

	int  max;

	// n_poles
	int * n_poles;
	assert(cudaMalloc((void**) &(n_poles),	input.n_nuclides *sizeof(int))==cudaSuccess);
	assert(cudaMemcpy( n_poles, data->n_poles, 
		input.n_nuclides*sizeof(int),cudaMemcpyHostToDevice) == cudaSuccess);
	assert(cudaMemcpy( &(data_d->n_poles), &n_poles, 
		sizeof(int*),cudaMemcpyHostToDevice) == cudaSuccess);
	// n_windows
	int* n_windows;
	assert(cudaMalloc((void**) &(n_windows),input.n_nuclides *sizeof(int))==cudaSuccess);
	assert(cudaMemcpy( n_windows, data->n_windows, 
		input.n_nuclides*sizeof(int),cudaMemcpyHostToDevice) == cudaSuccess);
	assert(cudaMemcpy( &(data_d->n_windows), &n_windows, 
		sizeof(int*),cudaMemcpyHostToDevice) == cudaSuccess);

	max = 0;
	for ( int i = 0; i < input.n_nuclides; i++)
		if ( data->n_poles[i] > max )
			max = data->n_poles[i];
	// 2D Pole array
	Pole* poles, *poles_h;
	assert(cudaMalloc((void**) &(poles), max*sizeof(Pole) * input.n_nuclides) == cudaSuccess);
	poles_h = (Pole*) malloc ( input.n_nuclides * max *sizeof(Pole) );
	for ( int i = 0; i < input.n_nuclides; i++) {
		memcpy( poles_h + i * max, data->poles[i], data->n_poles[i] * sizeof(Pole)) ; 
	}
	printf ("max poles: %i\n", max);
	assert(cudaMemcpy( &(data_d->pitch_poles), &max, sizeof(int),cudaMemcpyHostToDevice) == cudaSuccess);
	assert(cudaMemcpy( poles, poles_h, input.n_nuclides * max *sizeof(Pole), cudaMemcpyHostToDevice ) == cudaSuccess);
	assert(cudaMemcpy( &(data_d->poles_2d), &poles, sizeof(Pole*),cudaMemcpyHostToDevice) == cudaSuccess);

	// 2D Window array
	max = 0;
	for ( int i = 0; i < input.n_nuclides; i++)
		if ( data->n_windows[i] > max )
			max = data->n_windows[i];
	Window* windows_d, *windows_h;
	assert(cudaMalloc((void**) &windows_d, max*sizeof(Window) * input.n_nuclides) == cudaSuccess);
	windows_h = (Window*) malloc ( input.n_nuclides * max *sizeof(Window) );
	for ( int i = 0; i < input.n_nuclides; i++) {
		memcpy( windows_h + i * max, data->windows[i], data->n_windows[i] * sizeof(Window)) ; 
	}
	printf ("max windows: %i\n", max);
	assert(cudaMemcpy( &(data_d->pitch_windows), &max, sizeof(int),cudaMemcpyHostToDevice) == cudaSuccess);
	assert(cudaMemcpy( windows_d, windows_h, input.n_nuclides * max *sizeof(Window), cudaMemcpyHostToDevice ) == cudaSuccess);
	assert(cudaMemcpy( &(data_d->windows_2d), &windows_d, sizeof(Window*),cudaMemcpyHostToDevice) == cudaSuccess);

	// 2D double array - pseudo_K0RS_2d
	double * pseudo_K0RS_h, *pseudo_K0RS_d;
	assert(cudaMalloc((void**) &pseudo_K0RS_d, input.numL * sizeof(double) * input.n_nuclides) == cudaSuccess);
	pseudo_K0RS_h = (double*) malloc ( input.n_nuclides * input.numL *sizeof(double) );
	for ( int i = 0; i < input.n_nuclides; i++) {
		memcpy( pseudo_K0RS_h + i * input.numL, data->pseudo_K0RS[i], input.numL * sizeof(double)) ; 
	}
	assert(cudaMemcpy( pseudo_K0RS_d, pseudo_K0RS_h, input.n_nuclides * input.numL *sizeof(double), cudaMemcpyHostToDevice ) == cudaSuccess);
	assert(cudaMemcpy( &(data_d->pseudo_K0RS_2d), &pseudo_K0RS_d, sizeof(double*),cudaMemcpyHostToDevice) == cudaSuccess);
	
	// materials
	int  *num_nucs_d;
	assert(cudaMalloc((void**) &num_nucs_d, 12 * sizeof(int))==cudaSuccess);
	assert(cudaMemcpy( num_nucs_d, data->materials.num_nucs, 12*sizeof(int), cudaMemcpyHostToDevice) == cudaSuccess);
	assert(cudaMemcpy( &(data_d->materials.num_nucs), &num_nucs_d, sizeof(int*), cudaMemcpyHostToDevice) == cudaSuccess);

	max =0;
	for ( int i = 0; i < 12; i++)
		if ( data->materials.num_nucs[i] > max )
			max = data->materials.num_nucs[i];	
//	printf ( "max_num_nucs: %u\n", max);
	int* mats_d, *mats_h;
	assert(cudaMalloc((void**) &mats_d, 12 * max * sizeof(int))==cudaSuccess);
	mats_h = (int*) malloc ( 12*max*sizeof(int));
	for (int i=0; i < 12; i ++)
		memcpy( mats_h + i * max, data->materials.mats[i], data->materials.num_nucs[i]* sizeof(int) );
	assert(cudaMemcpy( &(data_d->materials.pitch), &max, sizeof(int), cudaMemcpyHostToDevice) == cudaSuccess);
	assert(cudaMemcpy( mats_d, mats_h, 12 * max *sizeof(int), cudaMemcpyHostToDevice ) == cudaSuccess);
	assert(cudaMemcpy( &(data_d->materials.mats_2d), &mats_d, sizeof(int*),cudaMemcpyHostToDevice) == cudaSuccess);

	double* concs_h, *concs_d;
	assert(cudaMalloc((void**) &concs_d, 12 * max * sizeof(double))==cudaSuccess);
	concs_h = (double*) malloc ( 12*max*sizeof(double));
	for (int i=0; i < 12; i ++)
		memcpy( concs_h + i * max, data->materials.concs[i], data->materials.num_nucs[i]* sizeof(double) );
	assert(cudaMemcpy( concs_d, concs_h, 12 * max *sizeof(double), cudaMemcpyHostToDevice ) == cudaSuccess);
	assert(cudaMemcpy( &(data_d->materials.concs_2d), &concs_d, sizeof(double*),cudaMemcpyHostToDevice) == cudaSuccess);

//	assert(cudaMalloc((void**) &input_d, sizeof(Input))==cudaSuccess);
//	assert(cudaMemcpy( input_d, &input, sizeof(Input), cudaMemcpyHostToDevice ) == cudaSuccess);
	return data_d;
}

void free_CalcDataPtrs_d ( CalcDataPtrs_d* data_d ) {
//	int * n_poles, *n_windows;
//	assert(cudaMemcpy( &n_poles, &(data_d->n_poles), 
//		sizeof(int*),cudaMemcpyDeviceToHost) == cudaSuccess);
	cudaFree ( data_d->n_poles );
//	assert(cudaMemcpy( &n_windows, &(data_d->n_windows), 
//		sizeof(int*),cudaMemcpyDeviceToHost) == cudaSuccess);
	cudaFree ( data_d->n_windows );
//	Pole* poles;
//	assert(cudaMemcpy( &poles, &(data_d->poles_2d), 
//		sizeof(Pole*),cudaMemcpyDeviceToHost) == cudaSuccess);
	cudaFree ( data_d->poles_2d );
//	Window* windows;
//	assert(cudaMemcpy( &windows, &(data_d->windows_2d), 
//		sizeof(Window*),cudaMemcpyDeviceToHost) == cudaSuccess);
	cudaFree ( data_d->windows_2d );
//	double* pseudo_K0RS;
//	assert(cudaMemcpy( &pseudo_K0RS, &(data_d->pseudo_K0RS_2d), 
//		sizeof(double*),cudaMemcpyDeviceToHost) == cudaSuccess);
	cudaFree ( data_d->pseudo_K0RS_2d );
//	int * num_nucs;
//	assert(cudaMemcpy( &num_nucs, &(data_d->materials.num_nucs), 
//		sizeof(int*), cudaMemcpyDeviceToHost) == cudaSuccess);
	cudaFree( data_d->materials.num_nucs);
//	int * mats;
//	assert(cudaMemcpy( &mats, &(data_d->materials.mats_2d), 
//		sizeof(int*), cudaMemcpyDeviceToHost) == cudaSuccess);
	cudaFree( data_d->materials.mats_2d);
//	double* concs;
//	assert(cudaMemcpy( &concs, &(data_d->materials.concs_2d), 
//		sizeof(double*), cudaMemcpyDeviceToHost) == cudaSuccess);
	cudaFree( data_d->materials.concs_2d);
	
	cudaFree ( data_d );
}

