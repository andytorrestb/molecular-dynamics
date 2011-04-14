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
int print_temp_energy(double width, double temp);
double avg_energy(double density, double temp);
int print_temp_response(double density, double temp_min, double temp_max, double temp_step);
int print_avg_msd(double density, double temp_min, double temp_max, double temp_step);
int print_cmvel(double density, double temp);

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
		if(argc > 3 && strcmp(argv[1], "--temp-energy") == 0) {
			return print_temp_energy(atof(argv[2]), atof(argv[3]));
		}
		if(argc > 5 && strcmp(argv[1], "--temp-response") == 0) {
			return print_temp_response(atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]));
		}
		if(argc > 5 && strcmp(argv[1], "--avg-msd") == 0) {
			return print_avg_msd(atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]));
		}
		if(argc > 3 && strcmp(argv[1], "--cm-vel") == 0) {
			return print_cmvel(atof(argv[2]), atof(argv[3]));
		}
	}
	
	fprintf(stderr, "Usage: %s [\e[4mMODE\e[0m]\n\n", argv[0]);
	fprintf(stderr, "Where \e[4mMODE\e[0m is one of:\n\n");
	fprintf(stderr,
		"\e[1m--vel-dist \e[4mn\e[0;1m \e[4mT\e[0m\n"
		"Generate a set of \e[4mn\e[0m particles with Maxwell-Boltzmann distributed "
		"velocities at temperature \e[4mT\e[0m and print out their magnitudes.\n\n");
	fprintf(stderr,
		"\e[1m--fcc-pos \e[4mdim\e[0m\n"
		"Set up \e[4mdim\e[0m FCC unit cells in each direction and print out "
		"particle positions.\n\n");
	fprintf(stderr,
		"\e[1m--lj-pos \e[4miter\e[0m\n"
		"Set up a system with two atoms at (-1,0,0) and (1,0,0) in a LJ "
		"potential and step the system ahead by 0.4 time units (step size of 0.004), "
		"printing the x coordinate of the each atom at every step. Runs for \e[4miter\e[0m "
		"iterations.\n\n");
	fprintf(stderr,
		"\e[1m--temp-energy \e[4mL\e[0;1m \e[4mT\e[0m\n"
		"Build an FCC lattice with 32 particles, place the particles in a Maxwell-Boltzmann "
		"distribution with temperature \e[4mT\e[0m, then evolve the system for 20 time units.\n\n");
	fprintf(stderr,
		"\e[1m--temp-response \e[4mdensity\e[0;1m \e[4mTmin\e[0;1m \e[4mTmax\e[0;1m \e[4mTstep\e[0m\n"
		"Starting with a FCC lattice at a given density, find the response of the total energy to the "
		"starting temperature.\n\n");
	fprintf(stderr,
		"\e[1m--avg-msd \e[4mdensity\e[0;1m \e[4mTmin\e[0;1m \e[4mTmax\e[0;1m \e[4mTstep\e[0m\n"
		"For each temperature in the range given by [\e[4mTmin\e[0m,\e[4mTmax\e[0m,\e[4mTstep\e[0m], "
		"set up a FCC lattice at \e[4mdensity\e[0m and find the average over all atoms of the square "
		"of the atom's displacement from its original position.\n\n");
	fprintf(stderr,
		"\e[1m--cm-vel \e[4mdensity\e[0;1m \e[4mT\e[0m\n"
		"Center of mass velocity (in each direction) as a function of time, starting in a FCC "
		"lattice with \e[4mdensity\e[0m and temperature \e[4mT\e[0m.\n\n");
	return 0;
}

int print_vel_dist(int n, double temp) {
	size_t i;
	sys_t *s = sys_alloc(n, 1.0);
	sys_maxboltz(s, temp);
	printf("%12s\n", "Velocity");
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
		fprintf(stderr, "Error generating FCC lattice.\n");
		return -1;
	}
	printf("%12s %12s %12s\n", "X", "Y", "Z");
	for(i = 0; i < s->n; i += 1) {
		printf("%12f %12f %12f\n", s->x[0][i], s->x[1][i], s->x[2][i]);
	}
	sys_free(s);
	return 0;
}

int print_lj_pos(int iterations) {
	int i;
	double energy, dt = 0.004;
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

int print_temp_energy(double width, double temp) {
	int i, iterations = 5000, dim = 2;
	double dt = 0.004, energy = 0.0, ke, pe;
	sys_t *s = sys_alloc(4*dim*dim*dim, width);
	if(!sys_fcc(s)) {
		fprintf(stderr, "Error generating FCC lattice.\n");
		return -1;
	}
	sys_maxboltz(s, temp);
	printf("%12s %12s %12s %12s\n", "Time", "Kinetic", "Potential", "Total");
	for(i = 0; i < iterations; i += 1) {
		sys_step(s, lj_force, dt);
		ke = sys_kinetic(s);
		pe = lj_potential(s);
		energy += log(ke+pe);
		printf("%12f %12e %12e %12e\n", i*dt, ke, pe, ke+pe);
	}
	fprintf(stderr, "Average energy = %12e\n", exp(energy/iterations));
	sys_free(s);
	return 0;
}

double avg_energy(double density, double temp) {
	int i, iterations = 2500, dim = 2;
	double dt = 0.004, energy = 0.0;
	sys_t *s = sys_alloc(4*dim*dim*dim, 4*dim*dim*dim/density);
	sys_fcc(s);
	sys_maxboltz(s, temp);
	for(i = 0; i < iterations; i += 1) {
		sys_step(s, lj_force, dt);
		energy += log(sys_kinetic(s) + lj_potential(s));
	}
	sys_free(s);
	return energy/iterations;
}

int print_temp_response(double density, double temp_min, double temp_max, double temp_step) {
	int i, ensemble = 1;
	double temp, energy;
	printf("%12s %12s\n", "Temperature", "Energy");
	for(temp = temp_min; temp <= temp_max; temp += temp_step) {
		energy = 0.0;
		for(i = 0; i < ensemble; i += 1) {
			energy += avg_energy(density, temp);
		}
		energy /= ensemble;
		printf("%12f %12e\n", temp, energy);
	}
	return 0;
}

int print_avg_msd(double density, double temp_min, double temp_max, double temp_step) {
	int i, j, iterations = 250, dim = 2, ensemble = 5;
	double temp, dt = 0.004, msd;
	sys_t **s = malloc(ensemble * sizeof(sys_t*));
	printf("%12s %12s %18s\n", "Temperature", "Time", "Avg. displacement");
	for(temp = temp_min; temp <= temp_max; temp += temp_step) {
		for(j = 0; j < ensemble; j += 1) {
			s[j] = sys_alloc(4*dim*dim*dim, 4*dim*dim*dim/density);
			sys_fcc(s[j]);
			sys_maxboltz(s[j], temp);
		}
		for(i = 0; i < iterations; i += 1) {
			for(j = 0; j < ensemble; j += 1) {
				sys_step(s[j], lj_force, dt);
			}
			if(i % 10 == 0) {
				msd = 0.0;
				for(j = 0; j < ensemble; j += 1) {
					msd += sys_msd(s[j]);
				}
				msd /= ensemble;
				printf("%12f %12f %12f\n", temp, i*dt, msd);
			}
		}
		for(j = 0; j < ensemble; j += 1) {
			sys_free(s[j]);
		}
	}
	return 0;
}

int print_cmvel(double density, double temp) {
	int i, iterations = 2500, dim = 2;
	double dt = 0.004, vx, vy, vz;
	sys_t *s = sys_alloc(4*dim*dim*dim, 4*dim*dim*dim/density);
	sys_fcc(s);
	sys_maxboltz(s, temp);
	printf("%12s %15s %12s %12s\n", "Time", "CM velocity (X)", "(Y)", "(Z)");
	for(i = 0; i < iterations; i += 1) {
		sys_step(s, lj_force, dt);
		sys_cmvel(s, &vx, &vy, &vz);
		printf("%12f %15g %12g %12g\n", i*dt, vx, vy, vz);
	}
	sys_free(s);
	return 0;
}
