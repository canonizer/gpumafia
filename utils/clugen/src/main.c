/** @file main.c implementation of the cluster generation utility 
		usage:
		-n<number-of-points> -d<dataset-dimensionality> -k<cluster-dimensionality>
		<output-file> 
 */

#include <assert.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/** the total number of points */
int n = 50000;
/** the dimensionality of the dataset */
int d = 4;
/** the dimensionality of the cluster */
int k = 2;
/** the fraction of points belonging to the cluster(s) */
double cf = 0.9;
/** the extents of the region to which the points belong */
double pmin = 0, pmax = 100;
/** the extent of the cluster */
double cmin = 42, cmax = 58;
/** the name of the output file */
char *out_path = "./my-cluster.dat";

/** the points*/
double *ps = 0;

/** generate a random double value within the specific range */
double drandom(double a, double b);

/** parses an integer >= specific integer value (inclusive) */
int parse_int(const char *name, int min_val);

/** prints the usage of the program and exits */
void print_usage(int exit_code);

/** parse command line options and fill in the global variables */
void parse_cmdline(int argc, char **argv);

/** randomizes the random generator with the timer */
void time_randomize();

/** generate the dataset */
void generate_points();

/** write the dataset to the output file */
void write_file();

int main(int argc, char **argv) {
	//fprintf(stderr, "the program started\n");
	parse_cmdline(argc, argv);
	//fprintf(stderr, "parsed cmdline\n");
	time_randomize();
	//fprintf(stderr, "randomized timer\n");
	generate_points();
	//fprintf(stderr, "generated points\n");
	write_file();
	//fprintf(stderr, "written to file\n");
}  // main

double drandom(double a, double b) {
	return a + (b - a) * ((double)random() / RAND_MAX);
}  // drandom

void print_usage(int exit_code) {
	fprintf(stderr, "usage:\n");
	fprintf(stderr, "\tclugen [-n<npoints>] [-d<ndims>] [-k<cluster_ndims>] " 
					"[<outfile>]\n");
	fprintf(stderr, "\tclugen -h\n");
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

void parse_cmdline(int argc, char **argv) {
	const char *optstr = "n:d:k:h:";
	int cur_opt;
	optind = 1;
	while((cur_opt = getopt(argc, argv, optstr)) > 0) {
		switch(cur_opt) {
		case 'n':
			n = parse_int("n", 1);
			break;
		case 'd':
			d = parse_int("d", 1);
			break;
		case 'k':
			k = parse_int("k", 1);
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

	// set the output file name
	if(optind < argc)
		out_path = argv[optind];
	
}  // parse_cmdline

void time_randomize() {
	srandom((unsigned int)time(0));
}  // time_randomize

void generate_points() {
	ps = (double *)malloc(sizeof(*ps) * n * d);
	for(int i = 0; i < n; i++) {
		int in_cluster = drandom(0, 1) < cf;
		for(int idim = 0; idim < d; idim++) {
			// TODO: select cluster dimensions randomly
			if(in_cluster && idim < k) {
				// generate point inside the cluster
				ps[i * d + idim] = drandom(cmin, cmax);
			} else {
				// generate a random point inside the region
				ps[i * d + idim] = drandom(pmin, pmax);
			}
		}
	}  // for(each point)
}  // generate_points

void write_file() {
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
