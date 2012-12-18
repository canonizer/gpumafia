/** @file utils.cpp implementation of some utility functions */

#include "options.h"
#include "utils.h"

#ifdef MAFIA_USE_DEVICE
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#include <stdio.h>
#include <stdlib.h>

void *bulk_alloc(size_t n) {
#ifdef MAFIA_USE_DEVICE
	//#ifdef XXX_YYY
	if(Options::options().use_device()) {
		void *ptr;
		CHECK(cudaMallocHost(&ptr, n));
		return ptr;
	} else
#endif
			return malloc(n);
}

void bulk_free(void *ptr) {
#ifdef MAFIA_USE_DEVICE
	//#ifdef XXX_YYY
	if(Options::options().use_device())
		cudaFreeHost(ptr);
	else
#endif
		free(ptr);
}
