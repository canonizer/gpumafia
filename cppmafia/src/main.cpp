/** @file main.cpp the main program */

#include <stdio.h>
#include <vector>

#include "mafia-io.h"
#include "mafia-solver.h"
#include "options.h"
#include "timing.h"

using namespace std;

int main(int argc, char **argv) {

  int npoints, ndims;
  double *points;

	// parse optinos
	Options opts(argc, argv);
	
  // load the points
	start_phase(PhaseReadData);
	read_points(opts.in_path, &points, &npoints, &ndims);

  // run the solver
  vector<vector<int> > cluster_idxs = mafia_solve(points, npoints, ndims, opts);
  
  // write the clusters
	start_phase(PhaseWriteData);
	write_clusters(opts.out_path, points, npoints, ndims, cluster_idxs);
	start_phase(PhaseEnd);

	// print timing info
	if(opts.flags & OptionTiming)
		print_timing_info();

  return 0;
}  // main
