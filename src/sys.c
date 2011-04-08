#include "sys.h"
#include <stdlib.h>
#include "rand.h"

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
	rand_maxboltz(s->n, s->v);
}
