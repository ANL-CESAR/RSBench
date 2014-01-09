#include"multibench.h"

void calculate_macro_xs( double * macro_xs, double mat, double E, Input input, Materials materials )
{
	// zero out macro vector
	for( int i = 0; i < 4; i++ )
		macro_xs[i] = 0;

	// for nuclide in mat
	for( int i = 0; i < materials.num_nucs[mat]; i++ )
	{
		double macro_xs_vector[4];

		calculate_micro_xs( double * macro_xs_vector, double 

	
}
