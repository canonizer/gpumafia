/** @file options.cpp implementations of parsing options */

#include "options.h"
#include "utils.h"

#include <assert.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// used to parse options
static struct option long_opts[] = {
	{"alpha", required_argument, 0, 'a'},
	{"beta", required_argument, 0, 'b'},
	{"bins", required_argument, 0, 'n'},
	{"min-wins", required_argument, 0, 'u'},
	{"max-wins", required_argument, 0, 'M'},
	{"no-set-dedup", no_argument, 0, 'D'},
	{"no-bitmap", no_argument, 0, 'P'}, 
	{"no-bitmaps", no_argument, 0, 'P'}, 
	{"verbose", no_argument, 0, 'V'},
	{"timing", no_argument, 0, 't'},
	{"device", no_argument, 0, 'D'},
	{"output-points", no_argument, 0, 'p'},
	{"help", no_argument, 0, 'h'},
	{0, 0, 0, 0}
};

// global options storage
Options *Options::opts = 0;

Options::Options(int argc, char **argv) 
	: in_path(0), out_path(0), min_nbins(200), min_nwindows(5), max_nwindows(20),
		alpha(1.5), beta(0.25), flags(OptionSetDedup | OptionUseBitmaps) {
	int cur_opt = 0;
	optind = 1;
	int opt_ind = -1;
	while((cur_opt = getopt_long
				 (argc, argv, ":a:b:n:u:M:Vhdp", long_opts, &opt_ind)) > 0) {
		switch(cur_opt) {
		case 'a':
			// alpha argument
			parse_double("alpha", &alpha);
			break;
		case 'b':
			// beta 
			parse_double("beta", &beta);
			break;
		case 'n':
			// number of bins
			parse_int("nbins", &min_nbins, 1);
			break;
		case 'u': 
			// minimum number of windows, = #windows when uniform
			parse_int("min-wins", &min_nwindows, 1);
			break;
		case 'M':
			// maximum number of windows
			parse_int("max-wins", &max_nwindows, 1);
			break;
		case 'V':
			// verbosity
			flags |= OptionVerbose;
			break;
		case 'D':
			// disable set deduplication
			flags &= ~OptionSetDedup;
			break;
		case 'p':
			// whether to output points
			flags |= OptionOutputPoints;
			break;
		case 'P':
			// disable bitmaps during point counting
			flags &= ~OptionUseBitmaps;
			break;
		case 'd':
			// use a GPU
			flags |= OptionUseDevice;
#ifndef MAFIA_USE_DEVICE
			// this option is an error if the program was not compiled for device
			// support 
			fprintf(stderr, "-d: option not supported, cppmafia was compiled "
							"without device (GPU) support\n");
			print_usage(-1);
#endif
			break;
		case 't':
			// timing information
			flags |= OptionTiming;
			break;
		case 'h':
			// print help message and exit
			print_usage(0);
			break;
		case '?':
			// unknown option
			fprintf(stderr, "%s: unknown option\n", long_opts[opt_ind].name);
			print_usage(-1);
			break;
		case ':':
			// missing argument
			fprintf(stderr, "argument missing for option %s\n", 
							long_opts[opt_ind].name);
			print_usage(-1);
			break;
		default:
			assert(0);
			break;
		}  // switch()
	}  // while(cur_opt)
	
	// optind contains the last argument, the name of the file
	if(optind >= argc) {
		// file path argument missing
		fprintf(stderr, "file with data for clustering expected\n");
		print_usage(-1);
	}
	in_path = (char *)malloc((strlen(argv[optind]) + 1) * sizeof(char));
	strcpy(in_path, argv[optind]);
	fprintf(stderr, "%s\n", in_path);
	
	// form the base output path; do this by removing the extension from the input
	// file name
	char *pdot_char = strrchr(in_path, '.');
	size_t base_len = pdot_char ? pdot_char - in_path : strlen(in_path);
	out_path = (char *)malloc((base_len + 1) * sizeof(char));
	memset(out_path, 0, base_len + 1);
	strncpy(out_path, in_path, base_len);
	
} // Options

void Options::parse_cmdline(int argc, char **argv) {
	opts = new Options(argc, argv);
}  // parse_cmdline

void Options::parse_double(const char *opt, double *pval) {
	if(!sscanf(optarg, "%lf", pval)) {
		fprintf(stderr, "%s: invalid %s argument, expecting double\n", 
						optarg, opt);
		print_usage(-1);
	}
} // parse_double
 
void Options::parse_int(const char *opt, int *pval, int min_val) {
	if(!sscanf(optarg, "%d", pval)) {
		fprintf(stderr, "%s: invalid %s argument, expecting int\n", 
						optarg, opt);
		print_usage(-1);
	} else if(*pval < min_val) {
		fprintf(stderr, "%d: value is less that minimum allowed %d\n", 
						*pval, min_val);
		print_usage(-1);
	}
} // parse_double

void Options::print_usage(int exit_code) {
	fprintf(stderr, "\n");
	fprintf(stderr, "usage: cppmafia [options] input-file\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "options (GNU option syntax):\n");
	fprintf(stderr, "\t-a,--alpha alpha - dense threshold (alpha)\n");
	fprintf(stderr, "\t-b,--beta beta - window merge threshold (alpha)\n");
	fprintf(stderr, "\t-n,--bins nbins - minimum number of bins\n");
	fprintf(stderr, "\t-u,--min-wins nwins - number of windows for uniform " 
					"dimensions\n");
	fprintf(stderr, "\t-M,--max-wins nwins - maximum number of windows for " 
					"a dimension\n");
	fprintf(stderr, "\t-h,--help - print this message and exit\n");
	fprintf(stderr, "\t-V,--verbose - print intermediate data as the algorithm " 
					"runs\n");
#ifdef MAFIA_USE_DEVICE
	fprintf(stderr, "\t-d,--device - offload some work to device, i.e. GPU\n");
#endif
	fprintf(stderr, "\t-p,--output-points - output cluster points in "
					"addition to indices \n");
	fprintf(stderr, "\t--no-set-dedup - turns off set deduplication; can degrade "
					"speed\n");
	fprintf(stderr, "\t--no-bitmap,--no-bitmaps - disables point counting with "
					"bitmaps;\n\t\tcan degrade speed extremely, and does not work on " 
					"devices\n");
	fprintf(stderr, "\t--timing - prints out additional timing info\n");
	exit(exit_code);
}  // print_usage

Options::~Options() {
	free(in_path);
	free(out_path);
} // ~Options
