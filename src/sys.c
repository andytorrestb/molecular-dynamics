#include "sys.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "rand.h"
#include "util.h"

int sys_fcc_n(size_t n);

// Set up a new system with n particles.
sys_t *sys_alloc(size_t n, double width) {
	size_t i;
	sys_t *s = malloc(sizeof(sys_t));
	s->n = n;
	s->width = width;
	for(i = 0; i < 3; i += 1) {
		s->x[i] = calloc(n, sizeof(double));
		s->disp[i] = calloc(n, sizeof(double));
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
		free(s->disp[i]);
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
	sys_zero_cmvel(s);
}

// Set the center of mass velocity to 0
void sys_zero_cmvel(sys_t *s) {
	size_t i, axis;
	double vel[3];
	sys_cmvel(s, &vel[0], &vel[1], &vel[2]);
	for(i = 0; i < s->n; i += 1) {
		for(axis = 0; axis < 3; axis += 1) {
			s->v[axis][i] -= vel[axis];
		}
	}
}

// Get the center of mass velocity
void sys_cmvel(sys_t *s, double *xout, double *yout, double *zout) {
	size_t i, axis;
	double vel[3] = {0.0, 0.0, 0.0};
	for(i = 0; i < s->n; i += 1) {
		for(axis = 0; axis < 3; axis += 1) {
			vel[axis] += s->v[axis][i];
		}
	}
	*xout = vel[0] / s->n;
	*yout = vel[1] / s->n;
	*zout = vel[2] / s->n;
}

// Find the positive integer k such that n = 4*k^3. Returns -1 if there is not such integer.
int sys_fcc_n(size_t n) {
	int k, k3;
	if(n < 0 || n % 4 != 0) {
		return -1;
	}
	k3 = n/4;
	for(k = 0; k <= k3; k += 1) {
		if(k*k*k == k3) {
			return k;
		}
	}
	return -1;
}

// Put the particle positions in an face-centered cubic configuration. Returns false if the number
// of particles in the system is an invalid number (not of the form 4*k^3).
int sys_fcc(sys_t *s) {
	int i, j, k, atom, idx = 0, dim;
	double pos[3][4] = {
		{0.0, 0.0, 0.5, 0.5},
		{0.0, 0.5, 0.0, 0.5},
		{0.0, 0.5, 0.5, 0.0}
	};
	double scale;
	
	// Check if the system has a valid number of particles
	dim = sys_fcc_n(s->n);
	if(dim == -1) {
		return 0;
	}
	scale = s->width/dim;
	
	// Set the position of each atom in each unit cell
	for(i = 0; i < dim; i += 1) {
		for(j = 0; j < dim; j += 1) {
			for(k = 0; k < dim; k += 1) {
				for(atom = 0; atom < 4; atom += 1) {
					s->x[0][idx] = (i+pos[0][atom])*scale;
					s->x[1][idx] = (j+pos[1][atom])*scale;
					s->x[2][idx] = (k+pos[2][atom])*scale;
					idx += 1;
				}
			}
		}
	}
	return 1;
}

// Get the temperature of the system.
double sys_temp(sys_t *s) {
	return 2.0*sys_kinetic(s)/3.0/(s->n-1);
}

// Get the kinetic energy of the system.
double sys_kinetic(sys_t *s) {
	size_t i;
	double v2 = 0.0;
	for(i = 0; i < s->n; i += 1) {
		v2 += util_sumsq(i, s->v);
	}
	return v2/2.0;
}

// Return the shortest vector from particle i to particle j, taking out to the
// hypertoroidal topology of the system (periodic boundaries).
void sys_dist(sys_t *s, int i, int j, double *xout, double *yout, double *zout) {
	size_t axis;
	double r, *out[3] = {xout, yout, zout};
	for(axis = 0; axis < 3; axis += 1) {
		r = s->x[axis][j] - s->x[axis][i];
		*out[axis] = r - (int)(r/s->width + (r<0?-0.5:0.5))*s->width;
	}
}

// Step the system by dt time units using the provided force accumulation function.
// The force function should modify the `f` variable in the provided sys_t structure.
void sys_step(sys_t *s, void (*force)(sys_t*), double dt) {
	size_t i, axis;
	double tmp;
	for(i = 0; i < s->n; i += 1) {
		for(axis = 0; axis < 3; axis += 1) {
			tmp = s->v[axis][i]*dt + s->f[axis][i]*dt*dt/2.0;
			s->x[axis][i] += tmp;
			s->disp[axis][i] += tmp;
			s->v[axis][i] += s->f[axis][i]*dt/2.0;
		}
	}
	force(s);
	for(i = 0; i < s->n; i += 1) {
		for(axis = 0; axis < 3; axis += 1) {
			s->v[axis][i] += s->f[axis][i]*dt/2.0;
			s->x[axis][i] = fmod(s->width + fmod(s->x[axis][i], s->width), s->width);
		}
	}
}

// Find the mean square displacement of the particles in the system.
double sys_msd(sys_t *s) {
	size_t i, axis;
	double sum_dsq = 0.0;
	for(i = 0; i < s->n; i += 1) {
		for(axis = 0; axis < 3; axis += 1) {
			sum_dsq += s->disp[axis][i] * s->disp[axis][i];
		}
	}
	return sum_dsq / s->n;
}
