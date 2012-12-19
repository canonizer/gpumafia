/** @file window.cpp implementation of window functions */
#include <algorithm>

#include "window.h"

using namespace std;

Window::Window(int idim, int left, int width, int max) {
	this->idim = idim;
  this->left = left;
  this->width = width;
  this->max = max;
}

bool Window::can_merge_with(const Window &w2, double beta) const {
  if(!max || !w2.max) {
    // one is zero, merge only if both are zero
    return !max && !w2.max;
  } else {
    // both non-zero, thresholding
    double max_ratio = double(max) / w2.max;
    return 1 - min(max_ratio, 1 / max_ratio) <= beta;
	}
}

Window Window::merge_with(const Window &w2) const {
  Window w = Window(idim, min(left, w2.left), width + w2.width, ::max(max, w2.max));
	w.idim = idim;
  w.pleft = ::min(pleft, w2.pleft);
  w.pright = ::max(pright, w2.pright);
	return w;
}

void Window::compute_threshold(double alpha, int npoints, int nbins) {
  threshold = int(alpha * width * npoints / nbins);
}

bool Window::is_dense() {
	return max * width >= threshold;
}
