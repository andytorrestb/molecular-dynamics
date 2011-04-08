#include "sys.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "rand.h"
#include "util.h"

// Set up a new system with n particles.
sys_t *sys_alloc(size_t n) {
	size_t i;
	sys_t *s = malloc(sizeof(sys_t));
	s->n = n;
	for(i = 0; i < 3; i += 1) {
		s->x[i] = calloc(n, sizeof(double));
		s->v[i] = calloc(n, sizeof(double));
		s->f[i] = calloc(n, sizeof(double));
	}
	return s;
}

// Free the memory used by a system.
void sys_free(sys_t *s) {
	size_t i;
	for(i = 0; i < 3; i += 1) {
		free(s->x[i]);
		free(s->v[i]);
		free(s->f[i]);
	}
	free(s);
}

// Put the particle velocities in a Maxwell-Boltzmann distribution with temperature temp.
void sys_maxboltz(sys_t *s, double temp) {
	size_t i, j;
	double scale;
	rand_maxboltz(s->n, s->v);
	scale = sqrt(temp/sys_temp(s));
	for(i = 0; i < 3; i += 1) {
		for(j = 0; j < s->n; j += 1) {
			s->v[i][j] *= scale;
		}
	}
}

// Find the positive integer k such that n = 4*k^3. Returns -1 if there is not such integer.
int sys_fcc_n(size_t n) {
	int k, k3;
	if(n < 0 || n % 4 != 0) {
		return -1;
	}
	k3 = n/4;
	for(k = 0; k < k3; k += 1) {
		if(k*k*k == k3) {
			return k;
		}
	}
	return -1;
}

// Put the particle positions in an face-centered cubic configuration. Returns false if the number
// of particles in the system is an invalid number (not of the form 4*k^3).
int sys_fcc(sys_t *s) {
	int dim = sys_fcc_n(s->n);
	if(dim == -1) {
		return 0;
	}
	// TODO
}

// Get the temperature of the system.
double sys_temp(sys_t *s) {
	size_t i;
	double v2 = 0.0;
	for(i = 0; i < s->n; i += 1) {
		v2 += util_sumsq(i, s->v);
	}
	return v2/3/(s->n-1);
}
