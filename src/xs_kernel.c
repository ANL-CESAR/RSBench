#include "rsbench.h"

// This function uses a combination of the Abrarov Approximation
// and the QUICK_W three term asymptotic expansion.
// Only expected to use Abrarov ~0.5% of the time.
RSComplex fast_nuclear_W( RSComplex Z )
{
	// Abrarov 
	if( c_abs(Z) < 6.0 )
	{
		// Precomputed parts for speeding things up
		// (N = 10, Tm = 12.0)
		//RSComplex prefactor = 8.124330e+01 * I;
		RSComplex prefactor = {0, 8.124330e+01};
		double an[10] = {
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
		double neg_1n[10] = {
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

		double denominator_left[10] = {
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

		RSComplex t1 = {0, 12};
		RSComplex t2 = {12, 0};
		RSComplex i = {0,1};
		RSComplex one = {1, 0};
		//RSComplex W = i * ( 1 - fast_cexp(t1*Z) ) / (12. * Z );
		//RSComplex W = i * ( one - fast_cexp(c_mul(t1, Z)) ) / c_mul(t2, Z);
		//RSComplex W = i * ( c_sub(one, fast_cexp(c_mul(t1, Z))) ) / c_mul(t2, Z);
		//RSComplex W = c_mul(i, ( c_sub(one, fast_cexp(c_mul(t1, Z))) )) / c_mul(t2, Z);
		RSComplex W = c_div(c_mul(i, ( c_sub(one, fast_cexp(c_mul(t1, Z))) )) , c_mul(t2, Z));
		RSComplex sum = {0,0};
		for( int n = 0; n < 10; n++ )
		{
			//RSComplex top = neg_1n[n] * fast_cexp(t1*Z) - 1.;
			//RSComplex top = c_sub(neg_1n[n] * fast_cexp(t1*Z), 1.);
			//RSComplex top = c_sub(neg_1n[n] * fast_cexp(c_mul(t1, Z)), 1.);
			RSComplex t3 = {neg_1n[n], 0};
			RSComplex top = c_sub(c_mul(t3, fast_cexp(c_mul(t1, Z))), one);
			RSComplex t4 = {denominator_left[n], 0};
			//RSComplex bot = c_sub(t4, 144.*Z*Z);
			RSComplex t5 = {144, 0};
			RSComplex bot = c_sub(t4, c_mul(t5,c_mul(Z,Z)));
			RSComplex t6 = {an[n], 0};
			sum = c_add(sum, c_mul(t6, c_div(top,bot)));
		}
		//W += prefactor * Z  * sum;
		W = c_add(W, c_mul(prefactor, c_mul(Z, sum)));
		return W;
	}
	else
	{
		// QUICK_2 3 Term Asymptotic Expansion (Accurate to O(1e-6)).
		// Pre-computed parameters
		RSComplex a = {0.512424224754768462984202823134979415014943561548661637413182,0};
		RSComplex b = {0.275255128608410950901357962647054304017026259671664935783653, 0};
		RSComplex c = {0.051765358792987823963876628425793170829107067780337219430904, 0};
		RSComplex d = {2.724744871391589049098642037352945695982973740328335064216346, 0};

		RSComplex i = {0,1};
		RSComplex Z2 = c_mul(Z, Z);
		// Three Term Asymptotic Expansion
		//RSComplex W = I * Z * (a/(Z*Z - b) + c/(Z*Z - d));
		RSComplex W = c_mul(c_mul(Z,i), (c_add(c_div(a,(c_sub(Z2, b))) , c_div(c,(c_sub(Z2, d))))));

		return W;
	}
}

void calculate_macro_xs( double * macro_xs, int mat, double E, Input input, CalcDataPtrs data, RSComplex * sigTfactors, long * abrarov, long * alls ) 
{
	// zero out macro vector
	for( int i = 0; i < 4; i++ )
		macro_xs[i] = 0;

	// Pre-processing state to determine poles and indexing
	// Determine total number of poles for this lookup
	int num_nucs = data.materials.num_nucs[mat];
	int npoles = 0;
	//int npoles_per_nuc[500]; // We are assuming maximum 500 nuclides being tracked
	int nuc_of_pole[5000]; // We are assuming 10 windows per nuclide, so 5k poles max
	int offset_of_nuc[500]; // Offsets where nuclides begin
	int w_starts[500]; // window start offset locations
	int p_idx = 0;

	// The final macro values we want
	double sigT = 0;
	double sigA = 0;
	double sigF = 0;
	double sigE = 0;

	for( int nuclide = 0; nuclide < num_nucs; nuclide++ )
	{
		int n = data.materials.mats[mat][nuclide];
		offset_of_nuc[nuclide] = p_idx;
		// Calculate Window Index
		double spacing = 1.0 / data.n_windows[n];
		//printf("spacing = %lf\n", spacing);
		int window = (int) ( E / spacing );
		if( window == data.n_windows[n] )
			window--;

		//printf("window = %d\n", window);
		Window w = data.windows[n][window];

		// Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
		double conc = data.materials.concs[mat][nuclide];
		sigT += E * w.T * conc;
		sigA += E * w.A * conc;
		sigF += E * w.F * conc;

		// Compute & store number of poles in the window
		int np = w.end - w.start;
		w_starts[nuclide] = w.start;
		//npoles_per_nuc[n] = np;
		npoles += np;

		// Compute & store the nuclide for each pole
		for( int p = 0; p < np; p++ )
		{
			nuc_of_pole[p_idx] = nuclide;
			p_idx++;
		}
		// Debug
		/*
		printf("E = %.2lf, mat = %d, macro_xs[0] = %.2lf, macro_xs[1] = %.2lf, macro_xs[2] = %.2lf, macro_xs[3] = %.2lf\n",
		E, mat, sigT, sigA, sigF, sigE );
		*/
	}
	//printf("npoles = %d\n", npoles);

	// By this point, we have added in background contributions, and know all our indices

	// Loop over poles
	for( int p = 0; p < npoles; p++ )
	{
		// Determine nuclide index
		int nuclide = nuc_of_pole[p];
		int n = data.materials.mats[mat][nuclide];
		
		double conc = data.materials.concs[mat][nuclide];

		// Determine pole index
		int p_idx = p - offset_of_nuc[nuclide];
		p_idx += w_starts[nuclide];

		double dopp = 0.5;

		Pole pole = data.poles[n][p_idx];
		//printf("pole: n = %d  p_idx = %d\n", n, p_idx);

		// Prep Z
		RSComplex E_c = {E, 0};
		RSComplex dopp_c = {dopp, 0};
		RSComplex Z = c_mul(c_sub(E_c, pole.MP_EA), dopp_c);

		// Evaluate Fadeeva Function
		RSComplex faddeeva = fast_nuclear_W( Z );
		//printf("real fad = %lf\n", faddeeva.r);

		int l_value = pole.l_value;

		// Compute SigT Factor
		double phi = data.pseudo_K0RS[n][l_value] * sqrt(E);

		if( l_value == 1 )
			phi -= - atan( phi );
		else if( l_value == 2 )
			phi -= atan( 3.0 * phi / (3.0 - phi*phi));
		else if( l_value == 3 )
			phi -= atan(phi*(15.0-phi*phi)/(15.0-6.0*phi*phi));

		phi *= 2.0;

		RSComplex sigTfactor;
		sigTfactor.r = cos(phi);
		sigTfactor.i = -sin(phi);

		// Update W
		sigT += ((c_mul( pole.MP_RT, c_mul(faddeeva, sigTfactor) )).r ) * conc;
		sigA += ((c_mul( pole.MP_RA , faddeeva)).r) * conc;
		sigF += ((c_mul( pole.MP_RF , faddeeva)).r) * conc;
	}
	
	// Compute SigE
	sigE = sigT - sigA;

	macro_xs[0] = sigT;
	macro_xs[1] = sigA;
	macro_xs[2] = sigF;
	macro_xs[3] = sigE;

	// Debug
	//printf("E = %.2lf, mat = %d, macro_xs[0] = %.2lf, macro_xs[1] = %.2lf, macro_xs[2] = %.2lf, macro_xs[3] = %.2lf\n",
	//E, mat, macro_xs[0], macro_xs[1], macro_xs[2], macro_xs[3] );
	
}

// No Temperature dependence (i.e., 0K evaluation)
void calculate_micro_xs( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, RSComplex * sigTfactors)
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
		RSComplex PSIIKI;
		RSComplex CDUM;
		Pole pole = data.poles[nuc][i];
		RSComplex t1 = {0, 1};
		RSComplex t2 = {sqrt(E), 0 };
		PSIIKI = c_div( t1 , c_sub(pole.MP_EA,t2) );
		RSComplex E_c = {E, 0};
		CDUM = c_div(PSIIKI, E_c);
		sigT += (c_mul(pole.MP_RT, c_mul(CDUM, sigTfactors[pole.l_value])) ).r;
		sigA += (c_mul( pole.MP_RA, CDUM)).r;
		sigF += (c_mul(pole.MP_RF, CDUM)).r;
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
void calculate_micro_xs_doppler( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data, RSComplex * sigTfactors, long * abrarov, long * alls)
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
		RSComplex E_c = {E, 0};
		RSComplex dopp_c = {dopp, 0};
		RSComplex Z = c_mul(c_sub(E_c, pole.MP_EA), dopp_c);
		if( c_abs(Z) < 6.0 )
			(*abrarov)++;
		(*alls)++;

		// Evaluate Fadeeva Function
		RSComplex faddeeva = fast_nuclear_W( Z );

		// Update W
		sigT += (c_mul( pole.MP_RT, c_mul(faddeeva, sigTfactors[pole.l_value]) )).r;
		sigA += (c_mul( pole.MP_RA , faddeeva)).r;
		sigF += (c_mul( pole.MP_RF , faddeeva)).r;
	}

	sigE = sigT - sigA;

	micro_xs[0] = sigT;
	micro_xs[1] = sigA;
	micro_xs[2] = sigF;
	micro_xs[3] = sigE;
}

void calculate_sig_T( int nuc, double E, Input input, CalcDataPtrs data, RSComplex * sigTfactors )
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
		sigTfactors[i].r = cos(phi);
		sigTfactors[i].i = -sin(phi);

	}
}
