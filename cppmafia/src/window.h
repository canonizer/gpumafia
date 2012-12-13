#ifndef MAFIA_WINDOW_H_
#define MAFIA_WINDOW_H_

/** @file window.h definition of a single adaptive window  */

#include "bitmap.h"

struct Window {
  /** creates a new window; threshold is initialized later */
  Window(int left, int width, int max);
  /** checks whether this window can merge with the other one 
      @note only successive windows will be merged 
   */
  bool can_merge_with(const Window &w2, double beta) const;
  /** merge the window with another one 
      @note only successive windows will be merged
   */
  Window merge_with(const Window &w2) const;
  /** computes the threshold for the window, based on its width
      @param alpha the alpha parameter of the MAFIA algorithm
      @param npoints the total number of points fed to the algorithm
      @param nbins the number of bins along window's dimension
   */
  void compute_threshold(double alpha, int npoints, int nbins);
	/** checks whether the window is dense */
	bool is_dense();
  /** non-inclusive right border of the window */
  inline int right() const { return left + width - 1; }
  
  /** the window width, i.e. the number of bins in the window */
  int width;
  /** the window left bound, i.e. the start of the window */
  int left;
  /** the window threshold */
  int threshold;
  /** the window maximum, the maximum number of elements in a bin
      belonging to a window */ 
  int max;
  /** the left boundary of the window in point coordinates */
  double pleft;
  /** the right boundary of the window in point coordinates */
  double pright;
	/** the window's bitmap, i.e. the set of points belonging to the window; it is
  set for dense windows only */
	bitmap pset;
};

#endif
