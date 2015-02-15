#include "rsbench.h"

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

// RNG Used for Verification Option.
// This one has a static seed (must be set manually in source).
// Park & Miller Multiplicative Conguential Algorithm
// From "Numerical Recipes" Second Edition
double rn_v(void)
{
	static unsigned long seed = 1337;
	double ret;
	unsigned long n1;
	unsigned long a = 16807;
	unsigned long m = 2147483647;
	n1 = ( a * (seed) ) % m;
	seed = n1;
	ret = (double) n1 / m;
	return ret;
}

unsigned int hash(unsigned char *str, int nbins)
{
	unsigned int hash = 5381;

	while(*str != '\0')
		hash = ((hash << 5) + hash) + ((int) *(str++));

	return hash % nbins;
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
