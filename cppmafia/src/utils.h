#ifndef MAFIA_UTILS_H_
#define MAFIA_UTILS_H_

/** @file utils.h some utility functions and definitions */

/** compiled with device support */ 
#define MAFIA_USE_DEVICE

#include <stdlib.h>

/** a macro to check whether CUDA calls are successful */
#define CHECK(call) \
	{\
	cudaError_t res = (call);\
	if(res != cudaSuccess) {\
	const char* err_str = cudaGetErrorString(res);\
	fprintf(stderr, "%s (%d): %s in %s", __FILE__, __LINE__, err_str, #call);	\
	exit(-1);\
	}\
	}

/** division with rounding upwards */
inline int divup(int a, int b) { return a / b + (a % b ? 1 : 0); }

/** bulk allocation of memory of specified size */
void *bulk_alloc(size_t nbytes);

/** freeing of bulk-allocated memory */
void bulk_free(void *ptr);

/** accessing points; variables ps, n and d need be defined */
//#define PS(i, idim) ps[i * d + idim]
#define PS(i, idim) ps[(idim) * n + (i)]

/** accessing bitmap words; variables bmps, nwords and nwidnows need be defined
		*/ 
#define BMPS(iw, iword) bmps[(iw) * nwords + iword]

#endif
