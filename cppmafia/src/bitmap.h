#ifndef MAFIA_BITMAP_H_
#define MAFIA_BITMAP_H_

/** @file bitmap.h header file for a simple bitmap set */

#include "ref.h"

#include <stdlib.h>
#include <vector>

using namespace std;

// forward references
class bitmap;

/** reference to a single bit; this is a transient class */
struct BitRef {
public:
	/** creates a new bit reference */	
	inline BitRef(unsigned *data, int ibit) {
		this->data = data;
		this->ibit = ibit;
	}
	/** gets the bitref's value */
	operator bool() const { return (*data >> ibit) & 1u; } 
	/** boolean assignment; surely non-atomic */
	inline bool operator=(bool b) {
		unsigned mask = 1u << ibit;
		return *data = *data & ~mask | (b << ibit);
	}
	/** the data referenced */
	unsigned *data;
	/** the bit referenced */
	int ibit;
}; // struct BitRef

/** internal implementation of the bitmap */
struct BitmapImpl : ref_object {
	
	/** creates a new bitmap implementation with the specific number of bits */
	BitmapImpl(int nbits, bool set_zero);
	/** destroys the bitmap implementation */
	virtual ~BitmapImpl();

	/** the number of elements (uints) to represent the bitmap */
	size_t n;
	/** the size of the bitmap, in bits; bits outside of this size have no meaning */
	size_t nbits;
	/** the bitmap data */
	unsigned *data;
	
};  // class BitmapImpl

/** the public handle to the bitmap implementation; provides built-in reference counting */
class bitmap : public ref<BitmapImpl> {

private:
	/** gets index of an element inside the array from the bit index */
	inline int elt_index(int ibit) const {return ibit / (sizeof(unsigned) * 8);}
	/** gets index of a bit inside the element from bit inside array index */
	inline int elt_bit_index(int ibit) const {return ibit % (sizeof(unsigned) * 8); }
	
public:
	/** creates a NULL bitset; note that it's not the same as bitmap(0) */
	bitmap();
	/** creates a new bitmap with the specified number of bits, and optionally
	zeroes out the elements */
	bitmap(size_t nbits, bool set_zero = true);
	/** makes a clone of this bitmap */
	bitmap clone();
	/** fills the vector with the indices of 1-bits in the array */
	void expand_into(vector<int>& inds);
	// inline functions
	/** gets the size of the bitmap (in bits) */
	inline size_t size() const { return ptr->nbits; }
	/** indexing of a specific bit */
	inline BitRef operator[](int ibit) {
		return BitRef(ptr->data + elt_index(ibit), elt_bit_index(ibit));
	}
	inline bool operator[](int ibit) const {
		return (ptr->data[elt_index(ibit)] >> elt_bit_index(ibit)) & 1;
	}

	/** counts the number of bits in the sets */
	int count() const;	
	// friend operator overloads
	/** the intersection of two sets; null & anything = anything */
	friend bitmap operator &(const bitmap &a, const bitmap &b);
	/** the union of two sets; null | anything = anything 
			@notes doesn't work when the number of bits are different
	 */
	friend bitmap operator |(const bitmap &a, const bitmap &b);
	/** the mutual difference of two sets; null ^ anything = anything */
	//friend bitmap operator ^(const bitmap &a, const bitmap &b);
	/** checks two sets for equality; two null references are equal, and two sets
	of different length can't be equal */
	//friend bool operator ==(const bitmap &a, const bitmap &b);
	/** checks two sets for inequality; null is not equal to any other set, and
	two sets of different length are never equal */
	//friend bool operator !=(const bitmap &a, const bitmap &b);
	/** checks whether the first set is a subset of the second one; always returns
	false for sets of different lengths */
	//friend bool operator <(const bitmap &a, const bitmap &b);
	/** checks whether the second set is the subset of the first one; always
	returns false if the sets are of different lengths; null is a subset of any set */
	//friend bool operator >(const bitmap &a, const bitmap &b);
	/** checks whether the first set is an (improper) subset of the second one;
	always returns false for sets of different lengths; null is a subset of any
	set */
	//friend bool operator <=(const bitmap &a, const bitmap &b);
	/** checks whether the second subset is an (improper) subset of the first one;
	always returns false for sets of different lengths; null is a subset of any
	set */
	//friend bool operator >=(const bitmap &a, const bitmap &b);

	/** assignment with intersection */
	const bitmap& operator&=(const bitmap &a);
	/** assignment with union */
	const bitmap& operator|=(const bitmap &a);
	/** assignment with mutual difference (xor) */
	//const bitmap& operator^=(const bitmap &a);
	
}; // class bitmap

#endif
