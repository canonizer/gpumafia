#ifndef MAFIA_CDU_H_
#define MAFIA_CDU_H_

/** @file cdu.h a single (candidate) dense unit in MAFIA algorithm */

#include <functional>
#include <vector>

#include "ref.h"
#include "window.h"

using namespace std;

/** a coordinate of a single window for a single dimension */
struct dimpair_t {
  dimpair_t(int dim, int win) {
    this->dim = dim;
    this->win = win;
  }
  bool operator==(const dimpair_t& dp2) const {
    return dim == dp2.dim && win == dp2.win;
  }
  bool operator!=(const dimpair_t& dp2) const {
    return !(*this == dp2);
  }
  /** the dimension number */
  int dim;
  /** the number of window inside a dimension */
  int win;
}; 

/** a single (candidate) dense unit */
class Cdu : public ref_object {
  /** creates an empty CDU */
  Cdu() : npoints(0), flag(false) {}
 public:
  /** creates a 1D CDU */
  Cdu(int dim, int win) : npoints(0), flag(false) {
    coords.push_back(dimpair_t(dim, win));
  }
	/** the dimensionality (length) of the CDU */
	inline int len() const { return coords.size(); }
  /** checks whether the CDU is the same as another one */
  bool operator==(const Cdu &cdu2) const;
  /** checks whether the CDU is distinct from another one */
  bool operator!=(const Cdu &cdu2) const { return !(*this == cdu2); }
	/** checks whether this CDU is 'less' than the other one. The less comparison
			is somewhat artificial, but ensures strict order on the set of all
			possible CDUs.
			CDU a is less than CDU b iff one of the following holds:
			- dimensionality of a is less than dimensionality of b
			- a and b are of the same dimensionality, but the vector of dimension
			numbers for a is lexicographically smaller than the vector of dimension
			numbers for b
			- a and b are of the same dimensionality, and their dimension number
			vectors are equal, but the corresponding window number vector for a is smaller
	*/
	bool operator<(const Cdu &cdu2) const;
	/** gets the subsequence of the CDU by omitting one of the components 
			@param ic the component to omit
			@remarks the caller is responsible for deallocation of the returned value
	 */
	Cdu* subsequence(int ic) const;
  /** checks whether the CDU can merge with another one */
  bool can_merge_with(const Cdu &cdu2);
  /** merge two units and create a third one; the new CDU is
      allocated in dynamic memory and is owned by the caller */
  Cdu *merge_with(const Cdu &cdu2);
	/** checks whether this unit is assimilated into another one, i.e. its set of
			coordinates is a subset of cdu2's set of coordinates */
	bool is_assimilated_into(const Cdu &cdu2) const;
  /** checks whether the CDU contains the point */
  template<class T> bool contains_point
	(const T *ps, int n, int d, int i, const vector<Window> &ws) const;
  /** counts the points inside the CDU (the naive way), also sets the point count */
  template<class T> int count_points_direct
	(const T* ps, int n, int d, const vector<Window> &ws);
	/** conts the points inside the CDU using bitmaps from windows, and also sets
			the point count */
	int count_points_bitmaps
	(int nwords, unsigned *bmps, const vector<Window> &ws);
	/** checks whether CDU is dense with respect to window thresholds */
	bool is_dense(const vector<Window> &ws) const;
	/** checks whether this DU shares a face with another one; 
			@notes it returns false if two DUs are the same
	 */
	bool has_common_face_with(const Cdu &cdu2, const vector<Window> &ws) const;
  /** the coordinates of the dense unit, sorted by ascending dimension
      number */
  vector<dimpair_t> coords;
  /** number of points in the CDU */
  int npoints;
	/** the flag which is used for various purposes connected with CDUs; in
			particular, it can mark it as joined or already visited */
  bool flag;
};

/** the comparator for ref<Cdu>, so that we can put it into a set */
struct CduCmp : binary_function<ref<Cdu>, ref<Cdu>, bool> {
	inline bool operator() (const ref<Cdu> &p, const ref<Cdu> &q) const {
		// TODO: handle the 0-pointer in some way
		return *p < *q;
	}  // operator()
};  // CduCmp

#endif
