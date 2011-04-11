#include "lj.h"
#include "sys.h"
#include <stdlib.h>
#include "util.h"
#include <math.h>

void lj_force(sys_t *s) {
	size_t i, j, axis;
	double r2;
	for(i = 0; i < s->n-1; i += 1) {
		for(j = i+1; j < s->n; j += 1) {
			sys_dist(s, i, j, &s->f[0][i], &s->f[1][i], &s->f[2][i]);
			r2 = util_sumsq(i, s->f);
			for(axis = 0; axis < 3; axis += 1) {
				s->f[axis][i] *= -24.0*(2.0/pow(r2,7) - 1.0/pow(r2,4));
				s->f[axis][j] = -s->f[axis][i];
			}
		}
	}
}

double lj_potential(sys_t *s) {
	size_t i, j;
	double x, y, z, r2, poten = 0.0;
	for(i = 0; i < s->n-1; i += 1) {
		for(j = i+1; j < s->n; j += 1) {
			sys_dist(s, i, j, &x, &y, &z);
			r2 = x*x + y*y + z*z;
			poten += 4.0*(1.0/pow(r2,6) - 1.0/pow(r2,3));
		}
	}
	return poten;
}
