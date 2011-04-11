#ifndef _SYS_H
#define _SYS_H

#include <stdlib.h>

typedef struct sys {
	size_t n;
	double width;
	double *x[3], *v[3], *f[3];
} sys_t;

sys_t *sys_alloc(size_t n, double width);
void sys_free(sys_t *s);
void sys_maxboltz(sys_t *s, double temp);
int sys_fcc(sys_t *s);
double sys_temp(sys_t *s);
double sys_kinetic(sys_t *s);
void sys_dist(sys_t *s, int i, int j, double *xout, double *yout, double *zout);
void sys_step(sys_t *s, void (*force)(sys_t*), double dt);

#endif
