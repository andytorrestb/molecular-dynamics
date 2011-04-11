#ifndef _UTIL_H
#define _UTIL_H

#include <stdlib.h>

double util_sumsq(size_t i, double *arr[3]);
void util_stat(size_t n, double *data, double *mean, double *stdev);

#endif
