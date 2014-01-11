#include "resobench.h"

// Park & Miller Multiplicative Conguential Algorithm
// From "Numerical Recipes" Second Edition
double rn(unsigned long * seed)
{
	double ret;
	unsigned long n1;
	unsigned long a = 16807;
	unsigned long m = 2147483647;
	n1 = ( a * (*seed) ) % m;
	*seed = n1;
	ret = (double) n1 / m;
	return ret;
}

size_t get_mem_estimate( Input input )
{
	size_t radii = input.n_nuclides * sizeof(double);
	size_t resonances = input.n_nuclides * input.n_resonances * sizeof(Resonance) + input.n_nuclides * sizeof(Resonance *);
	size_t total = radii + resonances;
	
	return total;
}
