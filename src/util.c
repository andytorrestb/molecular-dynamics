#include "util.h"
#include <math.h>

double util_sumsq(size_t i, double *arr[3]) {
	size_t j;
	double sum = 0.0;
	for(j = 0; j < 3; j += 1) {
		sum += arr[j][i]*arr[j][i];
	}
	return sum;
}

void util_stat(size_t n, double *data, double *mean, double *stdev) {
	size_t i;
	double sum = 0.0, sumsq = 0.0;
	for(i = 0; i < n; i += 1) {
		sum += data[i];
		sumsq += data[i]*data[i];
	}
	*mean = sum/n;
	*stdev = sqrt(sumsq/n - sum*sum/n/n);
}
