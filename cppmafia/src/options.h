#ifndef MAFIA_OPTIONS_H_
#define MAFIA_OPTIONS_H_

/** @file options.h options for MAFIA algorithm and program */

/** option flags for MAFIA */
enum {
	/** no option flags */
	OptionNone = 0x0,
	/** whether to produce verbose output */
	OptionVerbose = 0x1,
	/** whether to use set cluster deduplication */
	OptionSetDedup = 0x2,
	/** whether to use bitmaps to count point membership */
	OptionUseBitmaps = 0x4,
	/** whether to output timing information as well */
	OptionTiming = 0x8
};

/** options for MAFIA algorithm and program */
struct Options {
	/** parses the options from the command line */
	Options(int argc, char **argv);
	/** destructor for the options class */
	~Options();	
private:
	/** prints the program's usage and exists with a specific code */
	static void print_usage(int exit_code);
	/** parses the double option, exits on failure */
	static void parse_double(const char *opt, double *pval);
	/** parses the unsigned integer option, exits on failure to parse the integer
			or if the integer is <= min_val */
	static void parse_int(const char *opt, int *pval, int min_val);
public:
	/** the input file for the program
			@note this class is responsible for releasing memory
	*/
	char* in_path;
	/** the output path for the program 
			@note this class is responsible for releasing memory
	 */
	char* out_path;
	/** the minimum number of bins per dimension */
	int min_nbins;
	/** the number of windows for a uniform dimension, sometimes also called the
	"minimum number of windows" */
	int min_nwindows;
	/** the maximum number of windows per dimension */
	int max_nwindows;
	/** alpha, the dense unit threshold */
	double alpha;
	/** beta, the threshold for maximums for merging neighbouring windows */
	double beta;
	/** the option flags*/
	int flags;
};  // struct Options

#endif
