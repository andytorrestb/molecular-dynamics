#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "rand.h"
#include "sys.h"
#include "util.h"
#include "lj.h"

int print_vel_dist(int n, double temp);
int print_fcc_pos(int dim);
int print_lj_pos(int iterations);
int print_temp_energy(double temp);

int main(int argc, char **argv) {
	srand(time(NULL));
	
	if(argc > 1) {
		if(argc > 3 && strcmp(argv[1], "--vel-dist") == 0) {
			return print_vel_dist(atoi(argv[2]), atof(argv[3]));
		}
		if(argc > 2 && strcmp(argv[1], "--fcc-pos") == 0) {
			return print_fcc_pos(atoi(argv[2]));
		}
		if(argc > 2 && strcmp(argv[1], "--lj-pos") == 0) {
			return print_lj_pos(atoi(argv[2]));
		}
	}
	
	fprintf(stderr, "Usage: %s [--vel-dist n T | --fcc-pos dim | --lj-pos iter]\n\n", argv[0]);
	fprintf(stderr,
		"--vel-dist n T: Generate a set of n particles with Maxwell-Boltzmann distributed "
		"velocities at temperature T and print out their magnitudes.\n\n");
	fprintf(stderr,
		"--fcc-pos dim: Set up dim FCC unit cells in each direction and print out "
		"particle positions.\n\n");
	fprintf(stderr,
		"--lj-pos iter: Set up a system with two atoms at (-1,0,0) and (1,0,0) in a LJ "
		"potential and step the system ahead by 0.4 time units (step size of 0.004), "
		"printing the x coordinate of the each atom at every step. Runs for the given "
		"number of iterations (`iter`).\n");
	return 0;
}

int print_vel_dist(int n, double temp) {
	size_t i;
	sys_t *s = sys_alloc(n, 1.0);
	sys_maxboltz(s, temp);
	for(i = 0; i < s->n; i += 1) {
		printf("%12f\n", sqrt(util_sumsq(i, s->v)));
	}
	sys_free(s);
	return 0;
}

int print_fcc_pos(int dim) {
	size_t i;
	sys_t *s = sys_alloc(4*dim*dim*dim, 1.0);
	if(!sys_fcc(s)) {
		printf("Error generating FCC lattice.\n");
		return -1;
	}
	for(i = 0; i < s->n; i += 1) {
		printf("%12f%12f%12f\n", s->x[0][i], s->x[1][i], s->x[2][i]);
	}
	sys_free(s);
	return 0;
}

int print_lj_pos(int iterations) {
	int i;
	double energy, dt = 0.001;
	sys_t *s = sys_alloc(2, 10.0);
	s->x[0][0] = 4;
	s->x[0][1] = 6;
	printf("%12s %12s %12s %12s\n", "Time", "Atom 1", "Atom 2", "Energy");
	for(i = 0; i < iterations; i += 1) {
		sys_step(s, lj_force, dt);
		energy = sys_kinetic(s) + lj_potential(s);
		printf("%12f %12f %12f %12f\n", i*dt, s->x[0][0] - 5.0, s->x[0][1] - 5.0, energy);
	}
	sys_free(s);
	return 0;
}
