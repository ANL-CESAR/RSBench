#include "rsbench.h"

// JRT 
int * generate_n_poles( Input input )
{
	int total_resonances = input.avg_n_poles * input.n_nuclides;

	int * R = (int *) calloc( input.n_nuclides, sizeof(int));

	for( int i = 0; i < total_resonances; i++ )
		R[rand() % input.n_nuclides]++;

	// Ensure all nuclides have at least 1 resonance
	for( int i = 0; i < input.n_nuclides; i++ )
		if( R[i] == 0 )
			R[i] = 1;
	
	/* Debug	
	for( int i = 0; i < input.n_nuclides; i++ )
		printf("R[%d] = %d\n", i, R[i]);
	*/
	
	return R;
}

int * generate_n_windows( Input input )
{
	int total_resonances = input.avg_n_windows * input.n_nuclides;

	int * R = (int *) calloc( input.n_nuclides, sizeof(int));

	for( int i = 0; i < total_resonances; i++ )
		R[rand() % input.n_nuclides]++;

	// Ensure all nuclides have at least 1 resonance
	for( int i = 0; i < input.n_nuclides; i++ )
		if( R[i] == 0 )
			R[i] = 1;
	
	/* Debug	
	for( int i = 0; i < input.n_nuclides; i++ )
		printf("R[%d] = %d\n", i, R[i]);
	*/
	
	return R;
}

// 
Pole ** generate_resonance_params( Input input, int * n_poles )
{
	// Allocating 2D contiguous matrix
	Pole ** R = (Pole **) malloc( input.n_nuclides * sizeof( Pole *));
	Pole * contiguous = (Pole *) malloc( input.n_nuclides * input.avg_n_poles * sizeof(Pole));

	int k = 0;
	for( int i = 0; i < input.n_nuclides; i++ )
	{
		R[i] = &contiguous[k];
		k += n_poles[i];
	}
	
	// fill with data
	for( int i = 0; i < input.n_nuclides; i++ )
		for( int j = 0; j < n_poles[i]; j++ )
		{
			R[i][j].EA = (double) rand() / RAND_MAX + (double) rand() / RAND_MAX * _Complex_I;
			R[i][j].RT = (double) rand() / RAND_MAX + (double) rand() / RAND_MAX * _Complex_I;
			R[i][j].RA = (double) rand() / RAND_MAX + (double) rand() / RAND_MAX * _Complex_I;
			R[i][j].RF = (double) rand() / RAND_MAX + (double) rand() / RAND_MAX * _Complex_I;
			R[i][j].l_value = rand() % input.numL;
		}
	
	/* Debug
	for( int i = 0; i < input.n_nuclides; i++ )
		for( int j = 0; j < n_poles[i]; j++ )
			printf("R[%d][%d]: Eo = %lf lambda_o = %lf Tn = %lf Tg = %lf Tf = %lf\n", i, j, R[i][j].Eo, R[i][j].lambda_o, R[i][j].Tn, R[i][j].Tg, R[i][j].Tf);
	*/

	return R;
}

Window ** generate_window_params( Input input, int * n_windows )
{
	// Allocating 2D contiguous matrix
	Pole ** R = (Pole **) malloc( input.n_nuclides * sizeof( Pole *));
	Pole * contiguous = (Pole *) malloc( input.n_nuclides * input.avg_windows * sizeof(Pole));

	int k = 0;
	for( int i = 0; i < input.n_nuclides; i++ )
	{
		R[i] = &contiguous[k];
		k += n_windows[i];
	}
	
	// fill with data
	for( int i = 0; i < input.n_nuclides; i++ )
		for( int j = 0; j < n_windows[i]; j++ )
		{
			R[i][j].EA = (double) rand() / RAND_MAX + (double) rand() / RAND_MAX * _Complex_I;
			R[i][j].RT = (double) rand() / RAND_MAX + (double) rand() / RAND_MAX * _Complex_I;
			R[i][j].RA = (double) rand() / RAND_MAX + (double) rand() / RAND_MAX * _Complex_I;
			R[i][j].RF = (double) rand() / RAND_MAX + (double) rand() / RAND_MAX * _Complex_I;
			R[i][j].l_value = rand() % input.numL;
		}
	return R;
}

double ** generate_pseudo_K0RS( Input input )
{
	double ** R = (double **) malloc( input.n_nuclides * sizeof( double * ));
	double * contiguous = (double *) malloc( input.n_nuclides * input.numL * sizeof(double));

	for( int i = 0; i < input.n_nuclides; i++ )
		R[i] = &contiguous[i*input.numL];

	return R;
}
