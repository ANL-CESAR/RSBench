#include "rsbench.h"

#define N 10
#define Tm 12.0
#define Tm2 144.0
// This one works!
double complex Abrarov_W( double complex Z )
{
	// Precomputed parts for speeding things up
	// (N = 10, Tm = 12.0)
	double complex prefactor = 8.124330e+01 * I;
	const double an[N] = {
		2.758402e-01,
		2.245740e-01,
		1.594149e-01,
		9.866577e-02,
		5.324414e-02,
		2.505215e-02,
		1.027747e-02,
		3.676164e-03,
		1.146494e-03,
		3.117570e-04
	};
	const double neg_1n[N] = {
		-1.0,
		1.0,
		-1.0,
		1.0,
		-1.0,
		1.0,
		-1.0,
		1.0,
		-1.0,
		1.0
	};

	const double denominator_left[N] = {
		9.869604e+00,
		3.947842e+01,
		8.882644e+01,
		1.579137e+02,
		2.467401e+02,
		3.553058e+02,
		4.836106e+02,
		6.316547e+02,
		7.994380e+02,
		9.869604e+02
	};

	double complex W = I * ( 1 - cexp(I*Tm*Z) ) / (Tm * Z );
	double complex sum = 0;
	#pragma simd
	for( int n = 1; n <= N; n++ )
	{
		int idx = n-1;
		complex double top = neg_1n[idx] * cexp(I*Tm*Z) - 1.0;
		complex double bot = denominator_left[idx] - Tm2*Z*Z;
		sum += an[idx] * (top/bot);
	}
	W += prefactor * Z  * sum;
	return W;
}
// This one works!
double complex tramm_faddeeva2( double complex Z )
{
	double complex W = I * ( 1 - cexp(I*Tm*Z) ) / (Tm * Z );
	double complex sum = 0;
	for( int n = 1; n <= N; n++ )
	{
		double an = 2.0 * sqrt(M_PI) / Tm * exp( - n*n * M_PI * M_PI / (Tm * Tm));
		complex double top = pow(-1,n) * cexp(I*Tm*Z) - 1.0;
		complex double bot = n*n * M_PI*M_PI - Tm*Tm*Z*Z;
		sum += an * (top/bot);
	}
	W += I * Tm*Tm * Z / sqrt(M_PI) * sum;
	return W;
}


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

		if( input.doppler == 1 )
			calculate_micro_xs_doppler( micro_xs, nuc, E, input, data, sigTfactors);
		else
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

// No Temperature dependence (i.e., 0K evaluation)
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

// Temperature Dependent Variation of Kernel
// (This involves using the Complex Faddeeva function to
// Doppler broaden the poles within the window)
void calculate_micro_xs_doppler( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, complex double * sigTfactors)
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

	double dopp = 0.5;

	// Loop over Poles within window, add contributions
	for( int i = w.start; i < w.end; i++ )
	{
		Pole pole = data.poles[nuc][i];

		// Prep Z
		double complex Z = (E - pole.MP_EA) * dopp;

		// Evaluate Fadeeva Function
		//double faddeeva = Faddeeva_w(Z, 0.1);
		//double faddeeva = tramm_faddeeva2( Z );
		double faddeeva = Abrarov_W( Z );

		// Update W
		sigT += creal( pole.MP_RT * faddeeva * sigTfactors[pole.l_value] );
		sigA += creal( pole.MP_RA * faddeeva);
		sigF += creal( pole.MP_RF * faddeeva);
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
