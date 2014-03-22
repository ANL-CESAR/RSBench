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
	int idx = (int) (E * data.n_resonances[nuc] );

	// Calculate sigTfactors
	Complex * sigTfactors = calculate_sig_T( E, input, data );


	
Complex * calculate_sig_T( int nuc, double E, Input input, CalcDataPtrs data )
{
	int max_L = input.num_L;
	Complex * sigTfactors = (Complex *) malloc( num_L * sizeof(Complex) );
	double phi;
	
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

		sigTfactors[i].r = cos(phi);
		sigTfactors[i].i = -sin(phi);
	}
	return sigTfactors;
}

 do iL = 1, multipoles % numL
      twophi(iL) = multipoles%pseudo_k0RS(iL)*sqrtE
      if(iL == 2) then
        twophi(iL) = twophi(iL) - atan(twophi(iL))
      else if (iL == 3) then
        arg = 3.0_8*twophi(iL)/(3.0_8-twophi(iL)**2)
        twophi(iL) = twophi(iL) - atan(arg)
      else if (iL == 4) then
        arg = twophi(iL)*(15.0_8-twophi(iL)**2)/(15.0_8-6.0_8*twophi(iL)**2)
        twophi(iL) = twophi(iL) - atan(arg)
      end if

      twophi(iL) = 2.0_8 * twophi(iL)
      sigT_factor(iL) = cmplx(cos(twophi(iL)),-sin(twophi(iL)), KIND=8)
    end do

// Reviewed
void old_calculate_micro_xs( double * micro_xs, int nuc, double E, Input input, CalcDataPtrs data)
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

	// Store in output array
	micro_xs[0] = XS_g;
	micro_xs[1] = XS_f;
	micro_xs[2] = XS_s;
	micro_xs[3] = XS_t;

	/* Debug
	printf("nuc = %d, E = %.2lf :: idx = %d, T = %.2lf, radius = %.2lf, theta_o = %.2lf, term1 = %.2lf, XS_g = %.2lf, XS_f = %.2lf, XS_s = %.2lf, XS_t = %.2lf\n", nuc, E, idx, T, radius, theta_o, term1, XS_g, XS_f, XS_s, XS_t);
	*/
}
