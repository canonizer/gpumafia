/** @file bitmap.cpp bitmap implementation */

#include "bitmap.h"
#include "utils.h"

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

// BitmapImpl methods
BitmapImpl::BitmapImpl(int nbits, bool set_zero) : data(0) {
	this->nbits = nbits;
	n = divup(nbits, sizeof(unsigned) * 8);
	data = (unsigned *)malloc(sizeof(*data) * n);
	// last data element always zeroed out, to avoid leftovers
	if(n && data)
		data[n - 1] = 0;
	if(set_zero)
		memset(data, 0, sizeof(*data) * n);
} // BitmapImpl

BitmapImpl::~BitmapImpl() {
	free(data);
} // ~BitmapImpl

// bitmap methods
bitmap::bitmap() {}

bitmap::bitmap(size_t nbits, bool set_zero) : 
ref<BitmapImpl>(new BitmapImpl(nbits, set_zero)) {}

bitmap bitmap::clone() {
	if(!this->ptr)
		return *this;
	bitmap new_bmp = bitmap((*this)->nbits, false);
	memcpy(new_bmp->data, (*this)->data, (*this)->n * sizeof(unsigned));
	return new_bmp;
}  // clone

void bitmap::expand_into(vector<int>& inds) {
	if(!this->ptr)
		return;
	unsigned *data = (*this)->data;
	int n = (*this)->n;
	for(int i = 0; i < n; i++) {
		unsigned d = data[i];
		for(int shift = 0; shift < sizeof(d) * 8; shift++)
			if((d >> shift) & 1)
				inds.push_back(i * sizeof(d) * 8 + shift);
	}
}  // expand_into

int bitmap::count() const {
	if(!ptr)
		return 0;
	int cnt = 0;
	for(size_t i = 0; i < ptr->n; i++)
		cnt += __builtin_popcount(ptr->data[i]);
	return cnt;
}  // count

void bitmap::print() const {
	if(!this->ptr)
		return;
	printf("[");
	for(int i = 0; i < (*this)->n; i++) {
		printf("%0x\n", (*this)->data[i]);
		if(i < (*this)->n - 1)
			printf(" ");
	}
	printf("]");
}

const bitmap &bitmap::operator &=(const bitmap &a) {
	if(this->ptr && a.ptr && (*this)->nbits == a->nbits) {
		// fast-track: compute in-place
		int n = (*this)->n;
		unsigned *data1 = (*this)->data, *data2 = a->data;
		for(int i = 0; i < n; i++)
			data1[i] &= data2[i];
	} else {
		// do the usual way
		*this = *this & a;
	}
	return *this;
}  // operator &= 

const bitmap &bitmap::operator |=(const bitmap &a) {
	return *this = *this | a;
}  // operator |=

bitmap operator &(const bitmap &a, const bitmap &b) {
	if(!a.ptr || !b.ptr)
		if(a.ptr)
			return a;
		else
			return b;
	// allocate new bitmap, set extra bits to zero
	size_t nbits = max(a->nbits, b->nbits);
	bitmap c(nbits, false);
	size_t n_min = min(a->n, b->n), n_max = max(a->n, b->n);
	memset(c->data + n_min, (n_max - n_min) * sizeof(unsigned), 0);
	// compute the intersection
	for(int i = 0; i < n_min; i++)
		c->data[i] = a->data[i] & b->data[i];
	return c;
}  // operator&

bitmap operator|(const bitmap &a, const bitmap &b) {
	if(!a.ptr || !b.ptr)
		if(a.ptr)
			return a;
		else 
			return b;
	size_t nbits = max(a->nbits, b->nbits);
	bitmap c(nbits, false);
	size_t n_min = min(a->n, b->n), n_max = max(a->n, b->n);
	const bitmap &a_max = n_max == a-> n ? a : b;
	// TODO: handle the case when the number of bits are not equal 
	memcpy(c->data + n_min, a_max->data + n_min, (n_max - n_min) * sizeof(unsigned));
	// compute the union
	for(int i = 0; i < n_min; i++)
		c->data[i] = a->data[i] | b->data[i];
	return c;
}  // operator|
