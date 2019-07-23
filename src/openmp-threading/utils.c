#include "rsbench.h"

// Implementation based on:
// z = x + iy
// cexp(z) = e^x * (cos(y) + i * sin(y))
RSComplex fast_cexp( RSComplex z )
{
	double x = z.r;
	double y = z.i;

	double t1 = exp(x);
	double t2 = cos(y);
	double t3 = sin(y);
	RSComplex t4 = {t2, t3};
	RSComplex t5 = {t1, 0};
	RSComplex result = c_mul(t5, (t4));
	return result;
}	

size_t get_mem_estimate( Input input )
{
	size_t poles = input.n_nuclides * input.avg_n_poles * sizeof(Pole) + input.n_nuclides * sizeof(Pole *);
	size_t windows = input.n_nuclides * input.avg_n_windows * sizeof(Window) + input.n_nuclides * sizeof(Window *);
	size_t pseudo_K0RS = input.n_nuclides * input.numL * sizeof( double ) + input.n_nuclides * sizeof(double);
	size_t other = input.n_nuclides * 2 * sizeof(int);

	size_t total = poles + windows + pseudo_K0RS + other;
	
	return total;
}
