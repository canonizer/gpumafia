/** @file mafia-solver.h interface for C++ MAFIA solver */

#ifndef MAFIA_SOLVER_H_
#define MAFIA_SOLVER_H_

#include <vector>

#include "cdu.h"
#include "options.h"
#include "ref.h"
#include "utils.h"
#include "window.h"

using namespace std;

/** The main class for MAFIA solving */
template<class T>
class MafiaSolver {
 public:
  /** creates a MAFIA solver; the parameters are the same as with
      mafia_solve()  */
  MafiaSolver
    (const T *points, int n, int d, const Options &opts);
  /** finds clusters using MAFIA algorithm; the return value is the
      same as with  mafia_solve() */
  vector<vector<int> > find_clusters();
  /** the destructor of the object */
  virtual ~MafiaSolver();

 private:
  // internal functions used to build clusters
  /** access the point coordinate */
  //inline const T & coord(int i, int idim) { return ps[i * d + idim]; }
  /** build the histograms */
  void build_histos();
	/** computes the histogram on host for a specific dimension */
	void compute_histo_host(int idim);
	/** compute point limits (pmin/pmax) on host */
	void compute_limits_host();
  /** build the windows */
  void build_windows();
  /** build the uniform windows 
      @param idim the dimension along which to build the windows
      @param nwindows the number of windows to build
      @param[in, out] the vector to which the resulting windows are appended
   */
  void build_uniform_windows(int idim, int nwindows, vector<Window> &ws);
	/** builds point membership bitmaps for dense windows */
	void build_bitmaps();
	/** computes a bitmap for a single window on host */
	void compute_bitmap_host(int idim, int iwin);
  /** finds the candidate dense units */
  void find_cdus();
	/** deduplicate the CDUs found */
	void dedup_cdus();
	/** naive N**2 CDU deduplication */
	void naive_dedup_cdus();
	/** find the dense units from current CDUs */
	void find_dense_cdus();
	/** clears away dense units which are new or unassimilated */
	void find_unjoined_dus();
	/** builds the graph of DUs */
	void build_du_graph();
	/** find connected components of the DU graph; 
			these components are what will become the current clusters
	 */
	void build_du_clusters();
	/** build the final clusters, in terms of point lists */
	void build_clusters();

	// inline getter functions for various properties
	/** check whether the output needs to be verbose */
	inline bool is_verbose() const { return flags & OptionVerbose; }
	/** check whether to use set deduplication */
	inline bool use_set_dedup() const { return flags & OptionSetDedup; }
	/** checks whether to use bitmaps */
	inline bool use_bitmaps() const { return flags & OptionUseBitmaps; }
	/** checks whether to use the device (i.e., GPU) */
	inline bool use_device() const { return flags & OptionUseDevice; }

	// functions which use the device (GPU)
#ifdef MAFIA_USE_DEVICE
	/** makes a simple on-device operation to ensure that the device is
			initialized before any further time measurements are performed */
	void touch_dev();
	/** allocates point data and copies them to device; this must be performed
			before any on-device operations with points can be performed */
	void copy_ps_to_device();
	/** compute point limits (pmin/pmax) on device */
	void compute_limits_dev();
	/** computes the histogram on the device for a specific dimension; note that
			it does only histogram computation proper, and not the other related
			tasks
			*/
	void compute_histo_dev(int idim);
	/** computes a bitmap for a single window on device */
	void compute_bitmap_dev(int idim, int iwin);
	/** frees the resources on the device */
	void free_dev_resources();
#endif

  // debug functions for printing content after specific stages
  /** print the histograms */
  void print_histos();
  /** print the windows */
  void print_windows();
	/** prints the bitmaps */
	void print_bitmaps();
	/** prints the terminal dense units found */
	void print_terminal_dus();
	/** print just a list of CDUs */
	void print_dus(vector<ref<Cdu> > &dus);

 private:
	// initial data of the algorithm
  /** data dimensionality */
  int d; 
  /** number of points in the data */
  int n;
  /** the data points, stored in dim-first order */
  const T *ps;
	/** the data points stored on device */
	T *d_ps;
  // parameters
  /** the minimum number of bins */
  int min_nbins;
  /** the minimum number of windows */
  int min_nwindows;
  /** the maximum number of windows */
  int max_nwindows;
  /** dense window threshold */
  T alpha;
  /** separate window threshold */
  T beta;
  /** by how much to extend the domain to avoid special handling of
      the right border */
  T eps;
	/** the flags */
	int flags;
  
  // temporary working data used for finding clusters;
	// some data also have device versions
  /** data minimums, one per dimension */
  T *pmins;
  /** data maximums, one per dimension */
  T *pmaxs;
  /** bin histograms along each dimension */
  //vector<vector<int> > histos;
	/** an array containing the data for all histograms */
	int *histo_data;
	/** the starts of individual histograms in the histo_data array */
	int **histos;
	/** the number of bins in the histogram for each dimension */
	int *nbinss;
  /** windows along each dimension */
  vector<vector<Window> > windows;
  /** terminal dense units which will be connected into a graph;
   each array corresponds to a single dimensionality */
  vector<vector<ref<Cdu> > > terminal_dus;
  /** current list of dense units; of current dimension if currently
      identifying them, or of previous dimension if not */
  vector<ref<Cdu> > cur_dus;
	/** the list of new dense units; will be copied back to cur_dus at the end of
  the iteration */
	vector<ref<Cdu> > new_dus;
  /** current candidate dense units */
  vector<ref<Cdu> > cdus;
	/** the clusters consisting of DUs */
	vector<vector<ref<Cdu> > > du_clusters;
	/** the clusters consisting of point indices */
	vector<vector<int> > clusters;	
  /** current dimension for DU search */
  int cur_dim;

};  // MafiaSolver

/** a wrapper function for the MAFIA solver */
template <class T>
vector<vector<int> > mafia_solve
(const T *points, int npoints, int ndims, const Options &opts);

#endif
