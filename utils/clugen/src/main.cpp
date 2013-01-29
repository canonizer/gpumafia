/** @file main.c implementation of the cluster generation utility 
		usage:
		-n<number-of-points> -d<dataset-dimensionality> -k<cluster-dimensionality>
		<output-file> 
 */

#include "cluster.h"
#include "gen-context.h"

#include <algorithm>
#include <assert.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

const char *out_path = "./my-cluster.dat";

GenContext *gc = 0;

/** parses an integer >= specific integer value (inclusive) */
int parse_int(const char *name, int min_val);

/** parses a double in the specific interval (inclusive) */
double parse_double(const char *name, double min_val, double max_val);

/** prints the usage of the program and exits */
void print_usage(int exit_code);

/** parse command line options and fill in the global variables */
void parse_cmdline(int argc, char **argv);

/** fill in cluster dimension numbers */
void fill_cluster_dims(void);

/** randomizes the random generator with the timer */
void time_randomize(void);

/** generate the dataset */
void generate_points(void);

/** write the dataset to the output file */
void write_file(const double *ps, int n, int d);

int main(int argc, char **argv) {
	//fprintf(stderr, "the program started\n");
	gc = new GenContext();
	parse_cmdline(argc, argv);
	//fprintf(stderr, "parsed cmdline\n");
	time_randomize();
	gc->generate_template();
	double *ps = gc->generate_points();
	//fprintf(stderr, "randomized timer\n");
	write_file(ps, gc->n, gc->d);
	gc->print_info();
	//fprintf(stderr, "written to file\n");
	delete[] ps;
}  // main

void print_usage(int exit_code) {
	fprintf(stderr, "usage: clugen [options] [outfile]\n");
	fprintf(stderr, "options:\n");
	fprintf(stderr, "  %-20s %s\n", "-n <npoints>", 
					"number of points to generate (default 20000)");
	fprintf(stderr, "  %-20s %s\n", "-d <ndims>", 
					"point dimensionality (default 2)");
	fprintf(stderr, "  %-20s %s\n", "-k <max_clu_dims>", 
					"maximum cluster dimensionality (default 2)");
	fprintf(stderr, "  %-20s %s\n", "-m <nclusters>",
					"number of clusters (default 1)");
	fprintf(stderr, "  %-20s %s\n", "-f <fraction>", 
					"fraction of points which belong to the clusters (default 0.9)");
	fprintf(stderr, "  %-20s %s\n", "-s", 
					"generate ball clusters, not boxes");
	fprintf(stderr, "  %-20s %s\n", "-l", 
					"also generate clusters of lesser dimensionality than "
					"the maximum");
	fprintf(stderr, "  %-20s %s\n", "-N", 
					"try to generate non-intersecting clusters; an error is produced if "
					"cluster generation cannot be finished in fixed time");
	fprintf(stderr, "  %-20s %s\n", "-h", "print this message and exit\n");
	exit(exit_code);
}  // print_usage

int parse_int(const char *name, int min_val) {
	int val;
	if(!sscanf(optarg, "%d", &val)) {
		// can't parse integer
		fprintf(stderr, "expected integer for argument %s\n", name);
		print_usage(-1);
	} else if(val < min_val) {
		fprintf(stderr, "value for %s must be >= %d\n", name, min_val);
		print_usage(-1);
	}
	return val;
}  // parse_int

double parse_double(const char *name, double min_val, double max_val) {
	double val;
	if(!sscanf(optarg, "%lf", &val)) {
		// can't parse double
		fprintf(stderr, "expected double for argument %s\n", name);
		print_usage(-1);
	} else if(val < min_val || val > max_val) {
		fprintf(stderr, "value for %s must be between %lf and %lf, "
						"both ends inclusive", name, min_val, max_val);
		print_usage(-1);
	}
	return val;
}  // parse_double

void parse_cmdline(int argc, char **argv) {
	const char *optstr = "shlNn:d:k:m:f:";
	int cur_opt;
	optind = 1;
	bool kmin_is_kmax = true;
	while((cur_opt = getopt(argc, argv, optstr)) > 0) {
		switch(cur_opt) {
		case 's':
			gc->ct = SPHERE;
			break;
		case 'l':
			kmin_is_kmax = false;
			break;
		case 'N':
			gc->clusters_not_intersect = true;
			break;
		case 'n':
			gc->n = parse_int("n", 1);
			break;
		case 'd':
			gc->d = parse_int("d", 1);
			break;
		case 'k':
			gc->kmax = parse_int("k", 2);
			break;
		case 'm':
			gc->m = parse_int("m", 1);
			break;
		case 'f':
			gc->cf = parse_double("f", 0, 1);
			break;
		case 'h':
			print_usage(0);
			break;
		case '?':
			fprintf(stderr, "unknown option %c\n", optopt);
			print_usage(-1);
			break;
		case ':':
			fprintf(stderr, "option %c requires an argument\n", optopt);
			print_usage(-1);
			break;
		default:
			assert(0);
			break;
		}  // switch(cur_opt)
	}  // while(cur_opt > 0)
	

	// clamp kmax to d
	gc->kmax = min(gc->kmax, gc->d);

	// set kmin
	if(kmin_is_kmax)
		gc->kmin = gc->kmax;

	// set the output file name
	if(optind < argc)
		out_path = argv[optind];
}  // parse_cmdline

void time_randomize(void) {
	srandom((unsigned int)time(0));
}  // time_randomize

void write_file(const double *ps, int n, int d) {
	FILE *file = fopen(out_path, "w");
	if(!file) {
		fprintf(stderr, "cannot open file %s for writing", out_path);
		exit(-1);
	}
	for(int i = 0; i < n; i++) {
		for(int idim = 0; idim < d; idim++) {
			fprintf(file, "%lf", ps[i * d + idim]);
			if(idim < d - 1) 
				fprintf(file, " ");
		}
		fprintf(file, "\n");
	}
	fclose(file);
}  // write_file
