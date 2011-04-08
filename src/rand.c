#include "rand.h"
#include <stdlib.h>
#include <math.h>

// Returns a random number uniformly distributed in the range [1/(RAND_MAX+1), 1] ≈ (0,1].
double rand_uniform() {
	return (rand()+1)/((double)RAND_MAX+1);
}

// Returns a Gaussian-distributed random number using the Box-Muller method:
// http://en.wikipedia.org/wiki/Box_muller#Basic_form (note: second value discarded)
double rand_gaussian() {
	double u1 = rand_uniform(), u2 = rand_uniform();
	return sqrt(-2*log(u1))*cos(2*M_PI*u2);
}

// Given a 3-array of pointers to n-arrays representing the x, y, and z components of vectors,
// place the vectors in a Maxwell-Boltzmann distribution.
void rand_maxboltz(size_t n, double *arr[3]) {
	size_t i, j;
	for(i = 0; i < 3; i += 1) {
		for(j = 0; j < n; j += 1) {
			// In a Maxwell-Boltzmann distribution, the components of the velocity vectors are Gaussian distributed:
			// http://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution#Distribution_for_the_velocity_vector
			arr[i][j] = rand_gaussian();
		}
	}
}
