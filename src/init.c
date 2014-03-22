#include "rsbench.h"

// Reviewed 
int * generate_n_poles( Input input )
{
	int total_resonances = input.n_poles * input.n_nuclides;

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

// Reviewed
double * generate_nuclide_radii( Input input )
{
	double * A = (double *) malloc( input.n_nuclides * sizeof(double) );

	for( int i = 0; i < input.n_nuclides; i++ )
		A[i] = (double) rand() / RAND_MAX;
	
	/* Debug
	for( int i = 0; i < input.n_nuclides; i++ )
		printf("R[%d] = %.2lf\n", i, A[i]);
	*/
	
	return A;
}

// Reviewed
Resonance ** generate_resonance_params( Input input, int * n_poles )
{
	// Allocating 2D contiguous matrix
	Resonance ** R = (Resonance **) malloc( input.n_nuclides * sizeof( Resonance *));
	Resonance * contiguous = (Resonance *) malloc( input.n_nuclides * input.n_poles * sizeof(Resonance));

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
			R[i][j].Eo = (double) j / n_poles[i] + 0.5 / n_poles[i];
			R[i][j].lambda_o = (double) rand() / RAND_MAX;
			R[i][j].Tn = (double) rand() / RAND_MAX;
			R[i][j].Tg = (double) rand() / RAND_MAX;
			R[i][j].Tf = (double) rand() / RAND_MAX;
		}
	
	/* Debug
	for( int i = 0; i < input.n_nuclides; i++ )
		for( int j = 0; j < n_poles[i]; j++ )
			printf("R[%d][%d]: Eo = %lf lambda_o = %lf Tn = %lf Tg = %lf Tf = %lf\n", i, j, R[i][j].Eo, R[i][j].lambda_o, R[i][j].Tn, R[i][j].Tg, R[i][j].Tf);
	*/

	return R;
}
