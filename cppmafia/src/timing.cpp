/** @file timing.cpp implementation of simple timing */

#include "timing.h"

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** the phase times */
static double phase_times_g[PhaseEnd + 1];

/** the names of the phases */
static const char *phase_names_g[] = {
	"reading points",
	"building histogram",
	"building windows",
	"building bitmaps",
	"finding and deduplicating CDUs",
	"finding dense CDUs",
	"finding unjoined DUs",
	"building graph",
	"building clusters",
	"writing data"
};

/** the current timing phase */
static int cur_phase_g = PhaseNone;

/** the other time moment*/
static double t1_g = 0.0;

void start_phase(int phase) {
	double t2 = omp_get_wtime();
	if(cur_phase_g != PhaseNone) {
		// normal operation
		phase_times_g[cur_phase_g] += t2 - t1_g;
	}
	t1_g = t2;
	cur_phase_g = phase;
}  // start_phase

double phase_time(int phase) {
	return phase_times_g[phase];
}

const char *phase_name(int phase) {
	return phase_names_g[phase];
}

void print_timing_info() {
	// print individual timing info
	for(int phase = 0; phase < PhaseEnd; phase++)
		printf("%s: %lf s\n", phase_name(phase), phase_time(phase));
	// total timing and breaking to I/O-non-I/O parts
	double total_time = 0;
	for(int phase = 0; phase < PhaseEnd; phase++)
		total_time += phase_time(phase);
	double io_time = phase_time(PhaseReadData) + phase_time(PhaseWriteData);
	double compute_time = total_time - io_time;
	printf("total time: %lf s\n", total_time);
	printf("I/O time: %lf s\n", io_time);
	printf("compute time: %lf s\n", compute_time);
}  // print_timing_info
