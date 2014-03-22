#include "rsbench.h"

// Reviewed
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

	/* Debug
	printf("E = %.2lf, mat = %d, macro_xs[0] = %.2lf, macro_xs[1] = %.2lf, macro_xs[2] = %.2lf, macro_xs[3] = %.2lf\n",
	E, mat, macro_xs[0], macro_xs[1], macro_xs[2], macro_xs[3] );
	*/
	
}

void calculate_micro_xs( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data)
{
	// MicroScopic XS's to Calculate
	double sigT;
	double sigA;
	double sigF;
	double sigE;

	// Calculate Window Index
	double spacing = 1.0 / input.n_windows;
	int window = (int) ( E / spacing );
	if( window == input.n_windows )
		window--;

	// Calculate sigTfactors
	complex double * sigTfactors = calculate_sig_T( E, input, data );

	// Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
	sigT = E * CalcDataPtrs.fitting[window].T;
	sigA = E * CalcDataPtrs.fitting[window].A;
	sigF = E * CalcDataPtrs.fitting[window].F;

	// Loop over Poles within window, add contributions
	int start = data.win_boundaries[window].start;
	int end = data.win_boundaries[window].end;
	for( int i = start; i < end; i++ )
	{
		complex double PSIIKI;
		complex double CDUM;
		complex double pole = data.mpdata[i];
		PSIIKI = -(0.0 - 1.0 * complex_I ) / ( pole.MP_EA - sqrt(E) );
		CDUM = PSIIKI / E;
		sigT += creal( pole.MP_RT * CDUM * sigTfactors[i] );
		sigA += creal( pole.MP_RA * CDUM);
		sigF += creal( pole.MP_RF * CDUM);
	}

	sigE = sigT - sigA;

	micro_xs[0] = sigT;
	micro_xs[1] = sigA;
	micro_xs[2] = sigF;
	micro_xs[3] = sigE;
	
}

complex double * calculate_sig_T( double E, Input input, CalcDataPtrs data )
{
	double phi;
	int num_L = input.num_L;
	complex double * sigTfactors = (complex double *) malloc( num_L * sizeof(complex double) );

	for( int i = 0; i < num_L; i++ )
	{
		phi = data.pseudo_K0RS[i] * sqrt(E);

		if( i == 1 )
			phi -= - atan( phi );
		else if( i == 2 )
			phi -= atan( 3.0 * phi / (3.0 - phi*phi));
		else if( i == 3 )
			phi -= atan(phi*(15.0-phi*phi)/(15.0-6.0*phi*phi));

		phi *= 2.0;

		sigTfactors[i] = cos(phi) - sin(phi) * complex_I;
	}

	return sigTfactors;
}
/*
// Reviewed
void old_calculate_micro_xs( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data)
{
	int idx = (int) (E * data.n_poles[nuc] );
	Resonance ** R = data.resonance_params;
	Resonance * r = &R[nuc][idx];
	double T = 1.0 / data.n_poles[nuc];
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

	// Store in output array
	micro_xs[0] = XS_g;
	micro_xs[1] = XS_f;
	micro_xs[2] = XS_s;
	micro_xs[3] = XS_t;

}
*/
