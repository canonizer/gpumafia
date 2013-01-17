/** @file cdu.cpp implementation of (candidate) dense unit */


#include <stdio.h>
#include <stdlib.h>

#include "cdu.h"
#include "utils.h"
#include "window.h"

bool Cdu::operator==(const Cdu& cdu2) const {
  if(coords.size() != cdu2.coords.size())
    return false;
  for(int i = 0; i < coords.size(); i++) 
    if(coords[i] != cdu2.coords[i])
      return false;
  return true;
}

bool Cdu::operator<(const Cdu& cdu2) const {
	int n1 = coords.size(), n2 = cdu2.coords.size();
	if(n1 < n2)
		return true;
	else if(n1 > n2)
		return false;
	else {
		// n1 == n2
		// compare dimension numbers
		int n = n1;
		for(int i = 0; i < n; i++) {
			if(coords[i].dim < cdu2.coords[i].dim)
				return true;
			else if(coords[i].dim > cdu2.coords[i].dim)
				return false;
		}
		// compare window numbers
		for(int i = 0; i < n; i++) {
			if(coords[i].win < cdu2.coords[i].win)
				return true;
			else if(coords[i].win > cdu2.coords[i].win)
				return false;
		}
		// both CDUs are equal
		return false;
	}
}  // operator<

bool Cdu::can_merge_with(const Cdu& cdu2) {
  int n = coords.size();
  if(cdu2.coords.size() != n)
    return false;
  int ndifferent = 0;
  // simulate merging of two coordinate vectors for CDUs
  for(int i1 = 0, i2 = 0; ndifferent <= 2 && (i1 < n || i2 < n); ) {
    if(i1 < n && i2 < n) {
      dimpair_t dp1 = coords[i1], dp2 = cdu2.coords[i2];
      if(dp1.dim == dp2.dim) {
				if(dp1.win == dp2.win) {
					// same dimpair, move forward
					i1++, i2++;
				} else {
					// different coordinate for the same dimension, can't merge
					return false;
				}
      } else {
				// different, move one with smaller dimension
				ndifferent++;
				if(dp1.dim < dp2.dim)
					i1++;
				else
					i2++;
      }
    } else {
      // one is at the end, just more the other one
      ndifferent++;
      if(i2 == n)
				i1++;
      else
				i2++;
    }
  }  // for (i1, i2)
	return ndifferent == 2;
}  // can_merge_with

Cdu *Cdu::merge_with(const Cdu &cdu2) {
  Cdu *cdu = new Cdu();
  int n1 = coords.size(), n2 = cdu2.coords.size();
  for(int i1 = 0, i2 = 0; i1 < n1 || i2 < n2; ) {
    if(i1 < n1 && i2 < n2) {
      dimpair_t dp1 = coords[i1], dp2 = cdu2.coords[i2];
      if(dp1.dim == dp2.dim) {
				cdu->coords.push_back(dp1);
				i1++;
				i2++;
      } else if(dp1.dim < dp2.dim) {
				cdu->coords.push_back(dp1);
				i1++;
      } else {
				cdu->coords.push_back(dp2);
				i2++;
      }
    } else {
      // just add the coordinates
      if(i1 < n1) {
				cdu->coords.push_back(coords[i1++]);
      } else {
				cdu->coords.push_back(cdu2.coords[i2++]);
			}
    }
  }  // for(i1, i2)
  return cdu;
}  // merge_with

bool Cdu::is_assimilated_into(const Cdu& cdu2) const {
	int n1 = coords.size(), n2 = cdu2.coords.size();
	// check for subset with the same algorithm as doing merges
	// 
	for(int i1 = 0, i2 = 0; i1 < n1 || i1 < n2;) {
		// check if finished with first one only
		if(i1 == n1)
			return true;
		// check if finished with second one only
		if(i2 == n2)
			return false;
		// else not finished, check dimensions
		dimpair_t dp1 = coords[i1], dp2 = cdu2.coords[i2];
		if(dp1.dim == dp2.dim) {
			// same dimension, check windows; if different, not a subset
			if(dp1.win == dp2.win) {
				i1++; 
				i2++;
			} else 
				return false;
		} else if(dp1.dim < dp2.dim) {
			// this has an element cdu2 hasn't, not a subset
			return false;
		} else {
			// cdu2 has an element this hasn't, OK
			i2++;
		}
	}  // for(i1, i2)
	return true;
}  // is_assimilated_into

template<class T> bool Cdu::contains_point
(const T *ps, int n, int d, int i, const vector<Window> &ws) const {
  for(int icoord = 0; icoord < coords.size(); icoord++) {
    dimpair_t dp = coords[icoord];
    int dim = dp.dim;
    const Window &w = ws[dp.win];
    if(PS(i, dim) < w.pleft || PS(i, dim) >= w.pright)
      return false;
  }
  return true;
}  // contains_point()

template<class T> int Cdu::count_points_direct
(const T* ps, int n, int d, const vector<Window> &ws) {
	if(coords.size() == 1) {
		dimpair_t dp = coords[0];
		const Window &w = ws[dp.win];
		// 1-d, npoints = window.width * window.max
		npoints = w.width * w.max;
	} else {
		// multi-d, do real point counting
		// new variable for points, to satisfy OpenMP
		int local_npoints = 0;
		#pragma omp parallel for reduction(+:local_npoints)
		for(int i = 0; i < n; i++) {
			if(this->contains_point(ps, n, d, i, ws))
				local_npoints++;
		}
		npoints = local_npoints;
	}
  return npoints;
}  // count_points_direct

// moved __builtin_popcount into a separate function, to satisfy GCC
// and its OpenMP implementation
static inline int popcnt(unsigned word) {
	return __builtin_popcount(word);
}

int Cdu::count_points_bitmaps
(int nwords, unsigned *bmps, const vector<Window> & ws) {
	if(coords.size() == 1) {
		dimpair_t dp = coords[0];
		const Window &w = ws[dp.win];
		// 1-d, npoints = window.width * window.max
		npoints = w.width * w.max;
	} else {
		// multi-d, do real point counting
		// new variable for points, to satisfy OpenMP
		int local_npoints = 0;
		#pragma omp parallel for reduction(+:local_npoints)
		for(int iword = 0; iword < nwords; iword++) {
			unsigned word = ~0u;
			for(int icoord = 0; icoord < coords.size(); icoord++)
				word &= BMPS(coords[icoord].win, iword);
			local_npoints += popcnt(word);
		}  // for(iword)
		npoints = local_npoints;
	}
  return npoints;
}  // count_points

bool Cdu::is_dense(const vector<Window> &ws) const {
	for(int icoord = 0; icoord < coords.size(); icoord++) {
		dimpair_t dp = coords[icoord];
		const Window &w = ws[dp.win];
		if(npoints < w.threshold)
			return false;
	}
	return true;
}  // is_dense

bool Cdu::has_common_face_with(const Cdu &cdu2, const vector<Window> &ws) 
	const {
	// for a common face, dimensions must be the same, 
	// and the difference in window numbers exactly 1
	int n = coords.size();
	if(cdu2.coords.size() != n)
		return false;
	int sum = 0;
	for(int icoord = 0; icoord < n && sum <= 1; icoord++) {
		dimpair_t dp1 = coords[icoord], dp2 = cdu2.coords[icoord];
		if(dp1.dim != dp2.dim)
			return false;
		sum += abs(ws[dp1.win].iwin - ws[dp2.win].iwin);
	}
	return sum == 1;
}  // has_common_face_with()

// explicit instantiations
template bool Cdu::contains_point
(const float *ps, int n, int d, int i, const vector<Window> &ws) const;
template int Cdu::count_points_direct
(const float* p, int n, int d, const vector<Window> &ws);
template bool Cdu::contains_point
(const double *p, int n, int d, int i, const vector<Window> &ws) const;
template int Cdu::count_points_direct
(const double* p, int n, int d, const vector<Window> &ws);
