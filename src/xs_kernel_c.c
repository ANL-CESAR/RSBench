#include "rsbench.h"

// Reviewed
void calculate_macro_xs( double * macro_xs, int mat, double E, Input input, CalcDataPtrs data, Complex * sigTfactors ) 
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

	/*// Debug
	printf("E = %.2lf, mat = %d, macro_xs[0] = %.2lf, macro_xs[1] = %.2lf, macro_xs[2] = %.2lf, macro_xs[3] = %.2lf\n",
	E, mat, macro_xs[0], macro_xs[1], macro_xs[2], macro_xs[3] );
	*/
	
}

Complex cadd(Complex a, Complex b)
{
	Complex c;
	c.real = a.real + b.real;
	c.imag = a.imag + b.imag;
	return c;
}

Complex csubtract(Complex a, Complex b)
{
	Complex c;
	c.real = a.real - b.real;
	c.imag = a.imag - b.imag;
	return c;
}

Complex cmultiply(Complex a, Complex b)
{
	Complex c;
	c.real = a.real*b.real - a.imag*b.imag;
	c.imag = a.imag*b.real + a.real*b.imag;
	return c;
}

Complex cdivide(Complex a, Complex b)
{
	Complex c;
	c.real = (a.real*b.real + a.imag*b.imag)/(b.real*b.real + b.imag*b.imag);
	c.imag = (a.imag*b.real - a.real*b.imag)/(b.real*b.real + b.imag*b.imag);
	return c;
}

void calculate_micro_xs( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, Complex * sigTfactors)
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
		Complex PSIIKI;
		Complex CDUM;
		Pole pole = data.poles[nuc][i];
		//PSIIKI = -(0.0 - 1.0 * _Complex_I ) / ( pole.MP_EA - sqrt(E) );
		//CDUM = PSIIKI / E;
		//sigT += creal( pole.MP_RT * CDUM * sigTfactors[pole.l_value] );
		//sigA += creal( pole.MP_RA * CDUM);
		//sigF += creal( pole.MP_RF * CDUM);
/*		double d = (pole.MP_EA.real - sqrt(E)) * (pole.MP_EA.real - sqrt(E)) \
		     + (pole.MP_EA.imag * pole.MP_EA.imag);
		PSIIKI.real = pole.MP_EA.imag / d;
		PSIIKI.imag = sqrt(E) - pole.MP_EA.real / d;
		CDUM.real = PSIIKI.real / E;
		CDUM.imag = PSIIKI.imag / E;
		sigT += (pole.MP_RT.real * CDUM.real - pole.MP_RT.imag * CDUM.imag) * sigTfactors[pole.l_value].real \
		     - (pole.MP_RT.imag * CDUM.real + pole.MP_RT.real * CDUM.imag) * sigTfactors[pole.l_value].imag;
		sigA += pole.MP_RA.real * CDUM.real - pole.MP_RA.imag * CDUM.imag;
		sigF += pole.MP_RF.real * CDUM.real - pole.MP_RF.imag * CDUM.imag;
*/
		Complex tmp, t, a, f;
		PSIIKI.real = 0.0;
		PSIIKI.imag = 1.0;
		tmp.real = pole.MP_EA.real - sqrt(E);
		tmp.imag = pole.MP_EA.imag;
		PSIIKI = cdivide(PSIIKI, tmp);
		tmp.real = E;
		tmp.imag = 0;
		CDUM = cdivide(PSIIKI, tmp);
		tmp = cmultiply(pole.MP_RT, CDUM);
		t = cmultiply(tmp, sigTfactors[pole.l_value]);
		sigT += t.real;
		a = cmultiply(pole.MP_RA, CDUM);
		sigA += a.real;
		f = cmultiply(pole.MP_RF, CDUM);
		sigF += f.real;
	}

	sigE = sigT - sigA;

	micro_xs[0] = sigT;
	micro_xs[1] = sigA;
	micro_xs[2] = sigF;
	micro_xs[3] = sigE;
}

// Temperature Dependent Variation of Kernel
// (This involves using the Faddeeva function to
// Doppler broaden the poles within the window)
void calculate_micro_xs_doppler( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, Complex * sigTfactors)
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
//		double Z = (E - creal(pole.MP_EA)) * dopp;

		// Evaluate Fadeeva Function
//		double faddeeva = exp(-1.0 * creal(Z * Z)) * erfc(-1.0 * creal(Z * I));

		// Update W
//		sigT += creal( pole.MP_RT * faddeeva * sigTfactors[pole.l_value] );
//		sigA += creal( pole.MP_RA * faddeeva);
//		sigF += creal( pole.MP_RF * faddeeva);
	}

	sigE = sigT - sigA;

	micro_xs[0] = sigT;
	micro_xs[1] = sigA;
	micro_xs[2] = sigF;
	micro_xs[3] = sigE;
}

void calculate_sig_T( int nuc, double E, Input input, CalcDataPtrs data, Complex * sigTfactors )
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

		//sigTfactors[i] = cos(phi) - sin(phi) * _Complex_I;
		sigTfactors[i].real = cos(phi);
		sigTfactors[i].imag =  -sin(phi);
	}
}
