#include "rsbench.h"

// Moves all required data structures to the GPU's memory space
SimulationData move_simulation_data_to_device( Input in, SimulationData SD )
{
	printf("Allocating and moving simulation data to GPU memory space...\n");

	size_t sz;
	size_t total_sz = 0;

	// Shallow copy of CPU simulation data to GPU simulation data
	SimulationData GSD = SD;

	auto& rm = umpire::ResourceManager::getInstance();
	umpire::Allocator allocator = rm.getAllocator("DEVICE");

	// Move data to GPU memory space
	sz = GSD.length_num_nucs * sizeof(int);
	GSD.num_nucs = static_cast<int*>(allocator.allocate(sz));
	rm.copy(GSD.num_nucs, SD.num_nucs);
	total_sz += sz;

	sz = GSD.length_concs * sizeof(double);
	GSD.concs = static_cast<double*>(allocator.allocate(sz));
	rm.copy(GSD.concs, SD.concs);
	total_sz += sz;

	sz = GSD.length_mats * sizeof(int);
	GSD.mats = static_cast<int*>(allocator.allocate(sz));
	rm.copy(GSD.mats, SD.mats);
	total_sz += sz;

	sz = GSD.length_n_poles * sizeof(int);
	GSD.n_poles = static_cast<int*>(allocator.allocate(sz));
	rm.copy(GSD.n_poles, SD.n_poles);
	total_sz += sz;

	sz = GSD.length_n_windows * sizeof(int);
	GSD.n_windows = static_cast<int*>(allocator.allocate(sz));
	rm.copy(GSD.n_windows, SD.n_windows);
	total_sz += sz;

	sz = GSD.length_poles * sizeof(Pole);
	GSD.poles = static_cast<Pole*>(allocator.allocate(sz));
	rm.copy(GSD.poles, SD.poles);
	total_sz += sz;

	sz = GSD.length_windows * sizeof(Window);
	GSD.windows = static_cast<Window*>(allocator.allocate(sz));
	rm.copy(GSD.windows, SD.windows);
	total_sz += sz;

	sz = GSD.length_pseudo_K0RS * sizeof(double);
	GSD.pseudo_K0RS = static_cast<double*>(allocator.allocate(sz));
	rm.copy(GSD.pseudo_K0RS, SD.pseudo_K0RS);
	total_sz += sz;
	
	// Allocate verification array on device. This structure is not needed on CPU, so we don't
	// have to copy anything over.
	sz = in.lookups * sizeof(unsigned long);
	GSD.verification = static_cast<unsigned long*>(allocator.allocate(sz));
	total_sz += sz;

	GSD.length_verification = in.lookups;
	
	printf("GPU Intialization complete. Allocated %.0lf MB of data on GPU.\n", total_sz/1024.0/1024.0 );

	return GSD;

}

SimulationData initialize_simulation(Input input)
{
	uint64_t seed = INITIALIZATION_SEED;

	auto& rm = umpire::ResourceManager::getInstance();
	umpire::Allocator allocator = rm.getAllocator("HOST");
	
	// Get material data
	printf("Loading Hoogenboom-Martin material data...\n");
	SimulationData SD = get_materials(input, &seed); 
	
	// Allocate & fill energy grids
	printf("Generating resonance distributions...\n");
	SD.n_poles = generate_n_poles( input, &seed );
	SD.length_n_poles = input.n_nuclides;

	// Allocate & fill Window grids
	printf("Generating window distributions...\n");
	SD.n_windows = generate_n_windows(input, &seed);
	SD.length_n_windows = input.n_nuclides;

	// Prepare full resonance grid
	printf("Generating resonance parameter grid...\n");
	SD.poles = generate_poles( input, SD.n_poles, &seed, &SD.max_num_poles );
	SD.length_poles = input.n_nuclides * SD.max_num_poles;

	// Prepare full Window grid
	printf("Generating window parameter grid...\n");
	SD.windows = generate_window_params(input, SD.n_windows, SD.n_poles, &seed, &SD.max_num_windows);
	SD.length_windows = input.n_nuclides * SD.max_num_windows;

	// Prepare 0K Resonances
	printf("Generating 0K l_value data...\n");
	SD.pseudo_K0RS = generate_pseudo_K0RS(input, &seed);
	SD.length_pseudo_K0RS = input.n_nuclides * input.numL;

	size_t sz = input.lookups * sizeof(unsigned long);
	SD.verification = static_cast<unsigned long*>(allocator.allocate(sz));

	return SD;
}

void release_memory(SimulationData SD) {
	auto& rm = umpire::ResourceManager::getInstance();
	umpire::Allocator allocator = rm.getAllocator("HOST");

	allocator.deallocate(SD.num_nucs);
	allocator.deallocate(SD.concs);
	allocator.deallocate(SD.mats);
	allocator.deallocate(SD.n_poles);
	allocator.deallocate(SD.n_windows);
	allocator.deallocate(SD.poles);
	allocator.deallocate(SD.windows);
	allocator.deallocate(SD.pseudo_K0RS);
}

void release_device_memory(SimulationData GSD) {
	auto& rm = umpire::ResourceManager::getInstance();
	umpire::Allocator allocator = rm.getAllocator("DEVICE");

	allocator.deallocate(GSD.num_nucs);
	allocator.deallocate(GSD.concs);
	allocator.deallocate(GSD.mats);
	allocator.deallocate(GSD.n_poles);
	allocator.deallocate(GSD.n_windows);
	allocator.deallocate(GSD.poles);
	allocator.deallocate(GSD.windows);
	allocator.deallocate(GSD.pseudo_K0RS);
	allocator.deallocate(GSD.verification);
}

int * generate_n_poles( Input input, uint64_t * seed )
{
	int total_resonances = input.avg_n_poles * input.n_nuclides;


	auto& rm = umpire::ResourceManager::getInstance();
	umpire::Allocator allocator = rm.getAllocator("HOST");

	int * R = static_cast<int*>(allocator.allocate(input.n_nuclides * sizeof(int)));
	
	// Ensure all nuclides have at least 1 resonance
	for( int i = 0; i < input.n_nuclides; i++ )
		R[i] = 1;

	// Sample the rest
	for( int i = 0; i < total_resonances - input.n_nuclides; i++ )
		R[LCG_random_int(seed) % input.n_nuclides]++;
	
	/* Debug	
	for( int i = 0; i < input.n_nuclides; i++ )
		printf("R[%d] = %d\n", i, R[i]);
	*/
	
	return R;
}

int * generate_n_windows( Input input, uint64_t * seed )
{
	int total_resonances = input.avg_n_windows * input.n_nuclides;

	auto& rm = umpire::ResourceManager::getInstance();
	umpire::Allocator allocator = rm.getAllocator("HOST");

	int * R = static_cast<int*>(allocator.allocate(input.n_nuclides * sizeof(int)));
	
	// Ensure all nuclides have at least 1 resonance
	for( int i = 0; i < input.n_nuclides; i++ )
		R[i] = 1;

	// Sample the rest
	for( int i = 0; i < total_resonances - input.n_nuclides; i++ )
		R[LCG_random_int(seed) % input.n_nuclides]++;
	
	/* Debug	
	for( int i = 0; i < input.n_nuclides; i++ )
		printf("R[%d] = %d\n", i, R[i]);
	*/
	
	return R;
}

Pole * generate_poles( Input input, int * n_poles, uint64_t * seed, int * max_num_poles )
{
	// Pole Scaling Factor -- Used to bias hitting of the fast Faddeeva
	// region to approximately 99.5% (i.e., only 0.5% of lookups should
	// require the full eval).
	double f = 152.5;
	RSComplex f_c = {f, 0};

	int max_poles = -1;
	
	for( int i = 0; i < input.n_nuclides; i++ )
	{
		if( n_poles[i] > max_poles)
			max_poles = n_poles[i];
	}
	*max_num_poles = max_poles;

	// Allocating 2D matrix as a 1D contiguous vector
	auto& rm = umpire::ResourceManager::getInstance();
	umpire::Allocator allocator = rm.getAllocator("HOST");

	Pole * R = static_cast<Pole*>(allocator.allocate(input.n_nuclides * max_poles * sizeof(Pole)));
	
	// fill with data
	for( int i = 0; i < input.n_nuclides; i++ )
		for( int j = 0; j < n_poles[i]; j++ )
		{
			double r = LCG_random_double(seed);
			double im = LCG_random_double(seed);
			RSComplex t1 = {r, im};
			R[i * max_poles + j].MP_EA = c_mul(f_c,t1);
			r = LCG_random_double(seed);
			im = LCG_random_double(seed);
			RSComplex t2 = {f*r, im};
			R[i * max_poles + j].MP_RT = t2;
			r = LCG_random_double(seed);
			im = LCG_random_double(seed);
			RSComplex t3 = {f*r, im};
			R[i * max_poles + j].MP_RA = t3;
			r = LCG_random_double(seed);
			im = LCG_random_double(seed);
			RSComplex t4 = {f*r, im};
			R[i * max_poles + j].MP_RF = t4;
			R[i * max_poles + j].l_value = LCG_random_int(seed) % input.numL;
		}
	
	return R;
}

Window * generate_window_params( Input input, int * n_windows, int * n_poles, uint64_t * seed, int * max_num_windows )
{
	int max_windows = -1;
	
	for( int i = 0; i < input.n_nuclides; i++ )
	{
		if( n_windows[i] > max_windows)
			max_windows = n_windows[i];
	}
	*max_num_windows = max_windows;

	// Allocating 2D contiguous matrix
	auto& rm = umpire::ResourceManager::getInstance();
	umpire::Allocator allocator = rm.getAllocator("HOST");

	Window * R = static_cast<Window*>(allocator.allocate(input.n_nuclides * max_windows * sizeof(Window)));
	
	// fill with data
	for(int i = 0; i < input.n_nuclides; i++)
	{
		int space = n_poles[i] / n_windows[i];
		int remainder = n_poles[i] - space * n_windows[i];
		int ctr = 0;
		for(int j = 0; j < n_windows[i]; j++)
		{
			R[i * max_windows + j].T = LCG_random_double(seed);
			R[i * max_windows + j].A = LCG_random_double(seed);
			R[i * max_windows + j].F = LCG_random_double(seed);
			R[i * max_windows + j].start = ctr; 
			R[i * max_windows + j].end = ctr + space - 1;

			ctr += space;

			if (j < remainder)
			{
				ctr++;
				R[i * max_windows + j].end++;
			}
		}
	}

	return R;
}

double * generate_pseudo_K0RS( Input input, uint64_t * seed )
{
	auto& rm = umpire::ResourceManager::getInstance();
	umpire::Allocator allocator = rm.getAllocator("HOST");

	double * R = static_cast<double*>(allocator.allocate(input.n_nuclides * input.numL * sizeof(double)));

	for( int i = 0; i < input.n_nuclides; i++)
		for( int j = 0; j < input.numL; j++ )
			R[i * input.numL + j] = LCG_random_double(seed);

	return R;
}
