#ifndef MY_STATS_H
#define	MY_STATS_H

#include<stdlib.h>
#include<math.h>

/*	double comparator		*/
int compare_double (const void * a, const void * b) {
	double tmp= ( *(double*)a - *(double*)b );
	if ( tmp < 0 )
		return -1;
	if ( tmp > 0 )
		return 1;
	return 0;
}

/*	calculate mean and standard deviation for the vecotr of doubles buf	*/
void mean_var_stdev (const double* buf, const unsigned int* sz, double* mean, double* variance, double* stdev ) {
	if ( *sz < 2 )
		return;
	double sum=0, sum_sq=0;
	int i;
	for ( i = 0;i < *sz; i ++ ) {
		sum += buf[i];
		sum_sq += buf[i] * buf[i];
	}
	*mean = sum/(*sz);
	/*	unbiased estimate	*/
	*variance = (sum_sq- (*mean)*(*mean)*(*sz))/(*sz-1);
	*stdev= sqrt( *variance );
}

/**/
void my_stats ( double* vec, const unsigned int* sz, double** stats) {
	if ( vec == NULL || *sz < 2 )
		return;

	if ( *stats == NULL )
		*stats = (double*)malloc( 8 * sizeof(double) );
	else
		*stats = (double*)realloc ( *stats, 8 * sizeof(double) );
	qsort ( vec, *sz, sizeof(double), compare_double);
	(*stats)[0] = vec[0];
	(*stats)[1] = vec[*sz/4];
	(*stats)[2] = vec[*sz/2];
	(*stats)[3] = vec[*sz*3/4];
	(*stats)[4] = vec[*sz-1];
	mean_var_stdev( vec, sz, &(*stats)[5], &(*stats)[6], &(*stats)[7] );
}

#endif

