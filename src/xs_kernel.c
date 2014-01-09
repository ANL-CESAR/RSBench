#include"multibench.h"

void calculate_macro_xs( double * macro_xs, int mat, double E, Input input, Materials materials )
{
	// zero out macro vector
	for( int i = 0; i < 4; i++ )
		macro_xs[i] = 0;

	// for nuclide in mat
	for( int i = 0; i < materials.num_nucs[mat]; i++ )
	{
		double micro_xs[4];
		int nuc = materials.mats[mat][i];

		calculate_micro_xs( double * micro_xs, int nuc, double E, Input input, Materials materials);

		for( int j = 0; j < 4; j++ )
		{
			macro_xs[j] += micro_xs[j] * materials.concs[mat][i];
		}
	}
	
}
calculate_micro_xs( double * micro_xs, int nuc, double E, Input input, Materials materials)
{
	// need to determine if it's in resonance, or is interpolable

}
