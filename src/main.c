#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "rand.h"
#include "sys.h"

void print_vel_dist(int n);

int main(int argc, char **argv) {
	srand(time(NULL));
	
	if(argc > 1) {
		if(argc > 2 && strcmp(argv[1], "--vel-dist") == 0) {
			print_vel_dist(atoi(argv[2]));
			return 0;
		}
	}
	
	fprintf(stderr, "Usage: %s [--vel-dist n]\n", argv[0]);
	fprintf(stderr,
		"--vel-dist n: Generate a set of n particles with Maxwell-Boltzmann distributed "
		"velocities and print out their magnitudes.\n");
	return 0;
}

void print_vel_dist(int n) {
	size_t i;
	double r2;
	sys_t *s = sys_alloc(n);
	sys_maxboltz(s, 1.0);
	for(i = 0; i < s->n; i += 1) {
		r2 = s->v[0][i]*s->v[0][i] + s->v[1][i]*s->v[1][i] + s->v[2][i]*s->v[2][i];
		printf("%12f\n", sqrt(r2));
	}
	sys_free(s);
}
