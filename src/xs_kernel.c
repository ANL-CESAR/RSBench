#include "rsbench.h"

// Reviewed
void calculate_macro_xs( double * macro_xs, int mat, double E, Input input, CalcDataPtrs data, complex double * sigTfactors ) 
{
	// zero out macro vector
	for( int i = 0; i < 4; i++ )
		macro_xs[i] = 0;

	// for nuclide in mat
	for( int i = 0; i < (data.materials).num_nucs[mat]; i++ )
	{
		double micro_xs[4];
		int nuc = (data.materials).mats[mat][i];

		calculate_micro_xs( micro_xs, nuc, E, input, data, sigTfactors);

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

void calculate_micro_xs( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, complex double * sigTfactors)
{
	// MicroScopic XS's to Calculate
	double sigT;
	double sigA;
	double sigF;
	double sigE;

	// Calculate Window Index
	double spacing = 1.0 / data.n_windows[nuc];
	int window = (int) ( E / spacing );
	if( window == data.n_windows[nuc] )
		window--;

	// Calculate sigTfactors
	calculate_sig_T(nuc, E, input, data, sigTfactors );

	// Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
	Window w = data.windows[nuc][window];
	sigT = E * w.T;
	sigA = E * w.A;
	sigF = E * w.F;
//	printf ("w.size: %i\n", w.end-w.start+1);
	if ( numL < w.end-w.start+1)
		numL = w.end-w.start+1;
	// Loop over Poles within window, add contributions
	for( int i = w.start; i < w.end; i++ )
	{
		complex double PSIIKI;
		complex double CDUM;
		Pole pole = data.poles[nuc][i];
		PSIIKI = -(0.0 - 1.0 * _Complex_I ) / ( pole.MP_EA - sqrt(E) );
		CDUM = PSIIKI / E;
		sigT += creal( pole.MP_RT * CDUM * sigTfactors[pole.l_value] );
		sigA += creal( pole.MP_RA * CDUM);
		sigF += creal( pole.MP_RF * CDUM);
	}

	sigE = sigT - sigA;

	micro_xs[0] = sigT;
	micro_xs[1] = sigA;
	micro_xs[2] = sigF;
	micro_xs[3] = sigE;
}

void calculate_sig_T( int nuc, double E, Input input, CalcDataPtrs data, complex double * sigTfactors )
{
	double phi;
	for( int i = 0; i < input.numL; i++ )
	{
		phi = data.pseudo_K0RS[nuc][i] * sqrt(E);

		if( i == 1 )
			phi -= - atan( phi );
		else if( i == 2 )
			phi -= atan( 3.0 * phi / (3.0 - phi*phi));
		else if( i == 3 )
			phi -= atan(phi*(15.0-phi*phi)/(15.0-6.0*phi*phi));

		phi *= 2.0;

		sigTfactors[i] = cos(phi) - sin(phi) * _Complex_I;
	}
}

void calculate_sig_T_sim ( double E, int num_iter, const double* data, complex double * sigTfactors )
{
	double phi;
	for( int i = 0; i < num_iter; i++ )
	{
		phi = data[i] * sqrt(E);

		if( i == 1 )
			phi -= - atan( phi );
		else if( i == 2 )
			phi -= atan( 3.0 * phi / (3.0 - phi*phi));
		else if( i == 3 )
			phi -= atan(phi*(15.0-phi*phi)/(15.0-6.0*phi*phi));

		phi *= 2.0;

		sigTfactors[i] = cos(phi) - sin(phi) * _Complex_I;
	}
}
