#ifndef _RAND_H
#define _RAND_H

#include <stdlib.h>

double rand_uniform();
double rand_gaussian();
void rand_maxboltz(size_t n, double *arr[3]);

#endif
