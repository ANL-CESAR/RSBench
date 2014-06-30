#include "rsbench.h"

// Reviewed
void calculate_macro_xs( double * macro_xs, int mat, double E, Input input, CalcDataPtrs data, cuDoubleComplex * sigTfactors, int* counter, int* counter2 ) 
{
	// zero out macro vector
	for( int i = 0; i < 4; i++ )
		macro_xs[i] = 0;

//	if ( *counter2 < (data.materials).num_nucs[mat] )
//		*counter2 = (data.materials).num_nucs[mat];
	// for nuclide in mat
	for( int i = 0; i < (data.materials).num_nucs[mat]; i++ )
	{
		double micro_xs[4];
		int nuc = (data.materials).mats[mat][i];

//		calculate_micro_xs( micro_xs, nuc, E, input, data, sigTfactors, counter);
//		calculate_micro_xs_driver( micro_xs, nuc, E, input, data, sigTfactors);
		calc_sig_driver( micro_xs, nuc, E, input, data, sigTfactors);

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

void calculate_micro_xs( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, cuDoubleComplex * sigTfactors, int* counter)
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
//	calculate_sig_T(nuc, E, input, data, sigTfactors );
	calculate_sig_T_sim ( E, input.numL, data.pseudo_K0RS[nuc], sigTfactors );
	// Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
	Window w = data.windows[nuc][window];
	sigT = E * w.T;
	sigA = E * w.A;
	sigF = E * w.F;
	// Loop over Poles within window, add contributions
//	if ( *counter < ( w.end - w.start + 1) )
//		*counter = w.end - w.start + 1;
	cuDoubleComplex const1 = make_cuDoubleComplex(0, 1/E), const2 = make_cuDoubleComplex( sqrt(E), 0);
	for( int i = w.start; i < w.end; i++ )
	{
		Pole pole = data.poles[nuc][i];
		cuDoubleComplex CDUM = cuCdiv( const1, cuCsub( pole.MP_EA, const2 ) );
		sigT += cuCreal( cuCmul( pole.MP_RT, cuCmul( CDUM, sigTfactors[pole.l_value] ) ) );
		sigA += cuCreal( cuCmul( pole.MP_RA, CDUM) );
		sigF += cuCreal( cuCmul( pole.MP_RF, CDUM) );
	}

	sigE = sigT - sigA;

	micro_xs[0] = sigT;
	micro_xs[1] = sigA;
	micro_xs[2] = sigF;
	micro_xs[3] = sigE;
}

void calculate_sig_T( int nuc, double E, Input input, CalcDataPtrs data, cuDoubleComplex * sigTfactors )
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

		sigTfactors[i].x= cos(phi);
		sigTfactors[i].y= - sin(phi) ;
	}
}

void calculate_sig_T_sim ( double E, int num_iter, const double* data, cuDoubleComplex * sigTfactors )
{
	double phi;
	double sqrt_E = sqrt(E);
	for( int i = 0; i < num_iter; i++ )
	{
		phi = data[i] * sqrt_E;

		if( i == 1 )
			phi -= - atan( phi );
		else if( i == 2 )
			phi -= atan( 3.0 * phi / (3.0 - phi*phi));
		else if( i == 3 )
			phi -= atan(phi*(15.0-phi*phi)/(15.0-6.0*phi*phi));

		phi += phi;

		sigTfactors[i].x= cos(phi);
		sigTfactors[i].y= - sin(phi);
	}
}
