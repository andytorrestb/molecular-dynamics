#include "util.h"

double util_sumsq(size_t i, double *arr[3]) {
	size_t j;
	double sum = 0.0;
	for(j = 0; j < 3; j += 1) {
		sum += arr[j][i]*arr[j][i];
	}
	return sum;
}
