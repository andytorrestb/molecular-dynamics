#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "rand.h"
#include "sys.h"
#include "util.h"

void print_vel_dist(int n, double temp);

int main(int argc, char **argv) {
	srand(time(NULL));
	
	if(argc > 1) {
		if(argc > 3 && strcmp(argv[1], "--vel-dist") == 0) {
			print_vel_dist(atoi(argv[2]), atof(argv[3]));
			return 0;
		}
	}
	
	fprintf(stderr, "Usage: %s [--vel-dist n T]\n", argv[0]);
	fprintf(stderr,
		"--vel-dist n: Generate a set of n particles with Maxwell-Boltzmann distributed "
		"velocities at temperature T and print out their magnitudes.\n");
	return 0;
}

void print_vel_dist(int n, double temp) {
	size_t i;
	sys_t *s = sys_alloc(n);
	sys_maxboltz(s, temp);
	for(i = 0; i < s->n; i += 1) {
		printf("%12f\n", sqrt(util_sumsq(i, s->v)));
	}
	sys_free(s);
}
