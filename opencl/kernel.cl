#define PI 3.14159265359

typedef enum __hm{SMALL, LARGE, XL, XXL} HM_size;

#define HISTORY_BASED 1
#define EVENT_BASED 2

#define STARTING_SEED 1070
#define INITIALIZATION_SEED 42

typedef struct{
	double r;
	double i;
} RSComplex;

typedef struct{
	int nthreads;
	int n_nuclides;
	int lookups;
	HM_size HM;
	int avg_n_poles;
	int avg_n_windows;
	int numL;
	int doppler;
	int particles;
	int simulation_method;
	int kernel_id;
  int platform_id;
  int device_id;
} Input;

typedef struct{
	RSComplex MP_EA;
	RSComplex MP_RT;
	RSComplex MP_RA;
	RSComplex MP_RF;
	short int l_value;
} Pole;

typedef struct{
	double T;
	double A;
	double F;
	int start;
	int end;
} Window;
void calculate_macro_xs( double * macro_xs, int mat, double E, Input input, __global const int * num_nucs, __global const int * mats, int max_num_nucs, __global const double * concs, __global const int * n_windows, __global const double * pseudo_K0Rs, __global const Window * windows, __global const Pole * poles, int max_num_windows, int max_num_poles ) ;
void calculate_micro_xs( double * micro_xs, int nuc, double E, Input input, __global const int * n_windows, __global const double * pseudo_K0RS, __global const Window * windows, __global const Pole * poles, int max_num_windows, int max_num_poles);
void calculate_micro_xs_doppler( double * micro_xs, int nuc, double E, Input input, __global const int * n_windows, __global const double * pseudo_K0RS, __global const Window * windows, __global const Pole * poles, int max_num_windows, int max_num_poles );
int pick_mat( unsigned long * seed );
void calculate_sig_T( int nuc, double E, Input input, __global const double * pseudo_K0RS, RSComplex * sigTfactors );
RSComplex fast_nuclear_W( RSComplex Z );
double LCG_random_double(unsigned long * seed);
unsigned long LCG_random_int(unsigned long * seed);
unsigned long fast_forward_LCG(unsigned long seed, unsigned long n);
RSComplex c_add( RSComplex A, RSComplex B);
RSComplex c_sub( RSComplex A, RSComplex B);
RSComplex c_mul( RSComplex A, RSComplex B);
RSComplex c_div( RSComplex A, RSComplex B);
double c_abs( RSComplex A);
double fast_exp(double x);
RSComplex fast_cexp( RSComplex z );

__kernel void macro_xs_lookup_kernel(		Input in,
											__global const int *num_nucs,
											__global const int * mats,
											int max_num_nucs,
											__global const double *concs,
											__global const int *n_windows,
											__global const double *pseudo_K0RS,
											__global const Window * windows,
											__global const Pole * poles,
											int max_num_windows,
											int max_num_poles,
											__global int * verification_array )
{
	// Get the index of the current element to be processed
	int i = get_global_id(0);
  
  if( i >= in.lookups )
    return;
		
	// Set the initial seed value
	unsigned long seed = STARTING_SEED;	

	// Forward seed to lookup index (we need 2 samples per lookup)
	seed = fast_forward_LCG(seed, 2*i);
		
	// Randomly pick an energy and material for the particle
	double p_energy = LCG_random_double(&seed);
	int mat         = pick_mat(&seed); 
	
	double macro_xs_vector[4] = {0};

	// Perform macroscopic Cross Section Lookup
	calculate_macro_xs( macro_xs_vector, mat, p_energy, in, num_nucs, mats, max_num_nucs, concs, n_windows, pseudo_K0RS, windows, poles, max_num_windows, max_num_poles );
		
	// For verification, and to prevent the compiler from optimizing
	// all work out, we interrogate the returned macro_xs_vector array
	// to find its maximum value index, then increment the verification
	// value by that index. In this implementation, we prevent thread
	// contention by using an OMP reduction on the verification value.
	// For accelerators, a different approach might be required
	// (e.g., atomics, reduction of thread-specific values in large
	// array via CUDA thrust, etc).
	double max = -DBL_MAX;
	int max_idx = 0;
	for(int j = 0; j < 4; j++ )
	{
		if( macro_xs_vector[j] > max )
		{
			max = macro_xs_vector[j];
			max_idx = j;
		}
	}
	verification_array[i] = max_idx+1;

}

void calculate_macro_xs( double * macro_xs, int mat, double E, Input input, __global const int * num_nucs, __global const int * mats, int max_num_nucs, __global const double * concs, __global const int * n_windows, __global const double * pseudo_K0Rs, __global const Window * windows, __global const Pole * poles, int max_num_windows, int max_num_poles ) 
{
	// zero out macro vector
	for( int i = 0; i < 4; i++ )
		macro_xs[i] = 0;

	// for nuclide in mat
	for( int i = 0; i < num_nucs[mat]; i++ )
	{
		double micro_xs[4];
		int nuc = mats[mat * max_num_nucs + i];

		if( input.doppler == 1 )
			calculate_micro_xs_doppler( micro_xs, nuc, E, input, n_windows, pseudo_K0Rs, windows, poles, max_num_windows, max_num_poles);
		else
			calculate_micro_xs( micro_xs, nuc, E, input, n_windows, pseudo_K0Rs, windows, poles, max_num_windows, max_num_poles);

		for( int j = 0; j < 4; j++ )
		{
			macro_xs[j] += micro_xs[j] * concs[mat * max_num_nucs + i];
		}
		// Debug
		/*
		printf("E = %.2lf, mat = %d, macro_xs[0] = %.2lf, macro_xs[1] = %.2lf, macro_xs[2] = %.2lf, macro_xs[3] = %.2lf\n",
		E, mat, macro_xs[0], macro_xs[1], macro_xs[2], macro_xs[3] );
		*/
	}

	// Debug
	/*
	printf("E = %.2lf, mat = %d, macro_xs[0] = %.2lf, macro_xs[1] = %.2lf, macro_xs[2] = %.2lf, macro_xs[3] = %.2lf\n",
	E, mat, macro_xs[0], macro_xs[1], macro_xs[2], macro_xs[3] );
	*/
}

// No Temperature dependence (i.e., 0K evaluation)
void calculate_micro_xs( double * micro_xs, int nuc, double E, Input input, __global const int * n_windows, __global const double * pseudo_K0RS, __global const Window * windows, __global const Pole * poles, int max_num_windows, int max_num_poles)
{
	// MicroScopic XS's to Calculate
	double sigT;
	double sigA;
	double sigF;
	double sigE;

	// Calculate Window Index
	double spacing = 1.0 / n_windows[nuc];
	int window = (int) ( E / spacing );
	if( window == n_windows[nuc] )
		window--;

	// Calculate sigTfactors
	RSComplex sigTfactors[4]; // Of length input.numL, which is always 4
	calculate_sig_T(nuc, E, input, pseudo_K0RS, sigTfactors );

	// Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
	Window w = windows[nuc * max_num_windows + window];
	sigT = E * w.T;
	sigA = E * w.A;
	sigF = E * w.F;

	// Loop over Poles within window, add contributions
	for( int i = w.start; i < w.end; i++ )
	{
		RSComplex PSIIKI;
		RSComplex CDUM;
		Pole pole = poles[nuc * max_num_poles + i];
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
void calculate_micro_xs_doppler( double * micro_xs, int nuc, double E, Input input, __global const int * n_windows, __global const double * pseudo_K0RS, __global const Window * windows, __global const Pole * poles, int max_num_windows, int max_num_poles )
{
	// MicroScopic XS's to Calculate
	double sigT;
	double sigA;
	double sigF;
	double sigE;

	// Calculate Window Index
	double spacing = 1.0 / n_windows[nuc];
	int window = (int) ( E / spacing );
	if( window == n_windows[nuc] )
		window--;

	// Calculate sigTfactors
	RSComplex sigTfactors[4]; // Of length input.numL, which is always 4
	calculate_sig_T(nuc, E, input, pseudo_K0RS, sigTfactors );

	// Calculate contributions from window "background" (i.e., poles outside window (pre-calculated)
	Window w = windows[nuc * max_num_windows + window];
	sigT = E * w.T;
	sigA = E * w.A;
	sigF = E * w.F;

	double dopp = 0.5;

	// Loop over Poles within window, add contributions
	for( int i = w.start; i < w.end; i++ )
	{
		Pole pole = poles[nuc * max_num_poles + i];

		// Prep Z
		RSComplex E_c = {E, 0};
		RSComplex dopp_c = {dopp, 0};
		RSComplex Z = c_mul(c_sub(E_c, pole.MP_EA), dopp_c);

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

// picks a material based on a probabilistic distribution
int pick_mat( unsigned long * seed )
{
	// I have a nice spreadsheet supporting these numbers. They are
	// the fractions (by volume) of material in the core. Not a 
	// *perfect* approximation of where XS lookups are going to occur,
	// but this will do a good job of biasing the system nonetheless.

	double dist[12];
	dist[0]  = 0.140;	// fuel
	dist[1]  = 0.052;	// cladding
	dist[2]  = 0.275;	// cold, borated water
	dist[3]  = 0.134;	// hot, borated water
	dist[4]  = 0.154;	// RPV
	dist[5]  = 0.064;	// Lower, radial reflector
	dist[6]  = 0.066;	// Upper reflector / top plate
	dist[7]  = 0.055;	// bottom plate
	dist[8]  = 0.008;	// bottom nozzle
	dist[9]  = 0.015;	// top nozzle
	dist[10] = 0.025;	// top of fuel assemblies
	dist[11] = 0.013;	// bottom of fuel assemblies

	double roll = LCG_random_double(seed);

	// makes a pick based on the distro
	for( int i = 0; i < 12; i++ )
	{
		double running = 0;
		for( int j = i; j > 0; j-- )
			running += dist[j];
		if( roll < running )
			return i;
	}

	return 0;
}

void calculate_sig_T( int nuc, double E, Input input, __global const double * pseudo_K0RS, RSComplex * sigTfactors )
{
	double phi;

	for( int i = 0; i < 4; i++ )
	{
		phi = pseudo_K0RS[nuc * input.numL + i] * sqrt(E);

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
}

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
		RSComplex W = c_div(c_mul(i, ( c_sub(one, fast_cexp(c_mul(t1, Z))) )) , c_mul(t2, Z));
		RSComplex sum = {0,0};
		for( int n = 0; n < 10; n++ )
		{
			RSComplex t3 = {neg_1n[n], 0};
			RSComplex top = c_sub(c_mul(t3, fast_cexp(c_mul(t1, Z))), one);
			RSComplex t4 = {denominator_left[n], 0};
			RSComplex t5 = {144, 0};
			RSComplex bot = c_sub(t4, c_mul(t5,c_mul(Z,Z)));
			RSComplex t6 = {an[n], 0};
			sum = c_add(sum, c_mul(t6, c_div(top,bot)));
		}
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
		RSComplex W = c_mul(c_mul(Z,i), (c_add(c_div(a,(c_sub(Z2, b))) , c_div(c,(c_sub(Z2, d))))));

		return W;
	}
}

double LCG_random_double(unsigned long * seed)
{
	const unsigned long m = 9223372036854775808UL; // 2^63
	const unsigned long a = 2806196910506780709UL;
	const unsigned long c = 1UL;
	*seed = (a * (*seed) + c) % m;
	return (double) (*seed) / (double) m;
}	

unsigned long LCG_random_int(unsigned long * seed)
{
	const unsigned long m = 9223372036854775808UL; // 2^63
	const unsigned long a = 2806196910506780709UL;
	const unsigned long c = 1UL;
	*seed = (a * (*seed) + c) % m;
	return *seed;
}	

unsigned long fast_forward_LCG(unsigned long seed, unsigned long n)
{
	const unsigned long m = 9223372036854775808UL; // 2^63
	unsigned long a = 2806196910506780709UL;
	unsigned long c = 1UL;

	n = n % m;

	unsigned long a_new = 1;
	unsigned long c_new = 0;

	while(n > 0) 
	{
		if(n & 1)
		{
			a_new *= a;
			c_new = c_new * a + c;
		}
		c *= (a + 1);
		a *= a;

		n >>= 1;
	}

	return (a_new * seed + c_new) % m;
}

// Complex arithmetic functions

RSComplex c_add( RSComplex A, RSComplex B)
{
	RSComplex C;
	C.r = A.r + B.r;
	C.i = A.i + B.i;
	return C;
}

RSComplex c_sub( RSComplex A, RSComplex B)
{
	RSComplex C;
	C.r = A.r - B.r;
	C.i = A.i - B.i;
	return C;
}

RSComplex c_mul( RSComplex A, RSComplex B)
{
	double a = A.r;
	double b = A.i;
	double c = B.r;
	double d = B.i;
	RSComplex C;
	C.r = (a*c) - (b*d);
	C.i = (a*d) + (b*c);
	return C;
}

RSComplex c_div( RSComplex A, RSComplex B)
{
	double a = A.r;
	double b = A.i;
	double c = B.r;
	double d = B.i;
	RSComplex C;
	double denom = c*c + d*d;
	C.r = ( (a*c) + (b*d) ) / denom;
	C.i = ( (b*c) - (a*d) ) / denom;
	return C;
}

double c_abs( RSComplex A)
{
	return sqrt(A.r*A.r + A.i * A.i);
}


// Fast (but inaccurate) exponential function
// Written By "ACMer":
// https://codingforspeed.com/using-faster-exponential-approximation/
// We use our own to avoid small differences in compiler specific
// exp() intrinsic implementations that make it difficult to verify
// if the code is working correctly or not.
double fast_exp(double x)
{
  x = 1.0 + x * 0.000244140625;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  return x;
}

// Implementation based on:
// z = x + iy
// cexp(z) = e^x * (cos(y) + i * sin(y))
RSComplex fast_cexp( RSComplex z )
{
	double x = z.r;
	double y = z.i;

	// For consistency across architectures, we
	// will use our own exponetial implementation
	//double t1 = exp(x);
	double t1 = fast_exp(x);
	double t2 = cos(y);
	double t3 = sin(y);
	RSComplex t4 = {t2, t3};
	RSComplex t5 = {t1, 0};
	RSComplex result = c_mul(t5, (t4));
	return result;
}	
	
