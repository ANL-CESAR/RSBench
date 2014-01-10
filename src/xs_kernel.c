#include"multibench.h"

void calculate_macro_xs( double * macro_xs, int mat, double E, Input input, CalcDataPtrs data ) 
{
	// zero out macro vector
	for( int i = 0; i < 4; i++ )
		macro_xs[i] = 0;

	// for nuclide in mat
	for( int i = 0; i < (data.materials).num_nucs[mat]; i++ )
	{
		double micro_xs[4];
		int nuc = (data.materials).mats[mat][i];

		calculate_micro_xs( micro_xs, nuc, E, input, data);

		for( int j = 0; j < 4; j++ )
		{
			macro_xs[j] += micro_xs[j] * data.materials.concs[mat][i];
		}
	}
	
}

void calculate_micro_xs( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data)
{
	int idx = (int) (E * data.n_resonances[nuc] );

	Resonance ** R = data.resonance_params;

	Resonance * r = &R[nuc][idx];

	double T = 1.0 / data.n_resonances[nuc];
	double radius = data.nuclide_radii[nuc];

	// Reaction Cross Sections
	double theta_o = 4 * PI * r->lambda_o * r->Tn * r->Eo / T; 
	double term1 = theta_o * T * T * sqrt(r->Eo / E) / ( 4.0 * (E - r->Eo) * (E - r->Eo) + T * T);

	double XS_g = ( r->Tg / T ) * term1;
	double XS_f = ( r->Tf / T ) * term1;

	// Scattering Cross Section
	double XS_s = theta_o * term1 * ( r->Tn / T + ( 4.0*(E - r->Eo) * radius ) / ( T * theta_o * sqrt( r->Eo / E ) ) ) + 4.0 * PI * radius * radius;

	// Total Cross Section
	double XS_t = XS_g + XS_f + XS_s;

	micro_xs[0] = XS_g;
	micro_xs[1] = XS_f;
	micro_xs[2] = XS_s;
	micro_xs[3] = XS_t;
}