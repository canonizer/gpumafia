/** @file mafia-solver-device.cpp device-related parts of MAFIA solver
		this file is not compiled if compiling without GPU support
 */

#include "mafia-solver.h"

#ifdef MAFIA_USE_DEVICE

#include <cuda.h>
#include <limits>
#include <stdio.h>
#include <thrust/reduce.h>

using namespace thrust;

template<class T> void MafiaSolver<T>::touch_dev() {
	char h_arr[1], *d_arr;
	CHECK(cudaMalloc((void**)&d_arr, 1));
	CHECK(cudaMemcpy(d_arr, h_arr, 1, cudaMemcpyHostToDevice));
	CHECK(cudaMemcpy(h_arr, d_arr, 1, cudaMemcpyDeviceToHost));
	cudaFree(d_arr);
}  // touch_dev()

template<class T> void MafiaSolver<T>::copy_ps_to_device() {
	CHECK(cudaMalloc((void**)&d_ps, sizeof(*d_ps) * n * d));
	CHECK(cudaMemcpy(d_ps, ps, sizeof(*d_ps) * n * d, cudaMemcpyHostToDevice));
}

template<class T> void MafiaSolver<T>::compute_limits_dev() {
	//T *d_pmins, *d_pmaxs;
	//size_t limits_sz = sizeof(*d_pmins) * d;
	//CHECK(cudaMalloc((void**)&d_pmins, limits_sz));
	//CHECK(cudaMalloc((void**)&d_pmaxs, limits_sz));
	device_ptr<T> d_psptr(d_ps);
	for(int idim = 0; idim < d; idim++) {
		pmins[idim] = reduce(d_psptr + idim * n, d_psptr + (idim + 1) * n, 
												 std::numeric_limits<T>::infinity(), minimum<T>());
		pmaxs[idim] = reduce(d_psptr + idim * n, d_psptr + (idim + 1) * n, 
												 -std::numeric_limits<T>::infinity(), maximum<T>());
	}	
	//CHECK(cudaMemcpy(pmins, d_pmins, limits_sz, cudaMemcpyDeviceToHost));
	//CHECK(cudaMemcpy(pmaxs, d_pmaxs, limits_sz, cudaMemcpyDeviceToHost));
	//cudaFree(pmins);
	//cudaFree(pmaxs);
}  // compute_limits_dev

/** histogram computation kernel 
		@param psd a single dimension of the points on device
		@param n the number of points on the device
		@param pmin the minimum coordinate of the point along the current dimension
		@param piwidth the inverted width of the bin along the current dimension
		@param histo the histogram along the current dimension
		@param nbins the number of bins along the current dimension
*/
template<class T>
__global__ void histo_kernel
(T* psd, int n, double pmin, double piwidth, int *histo, int nbins) {
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if(i >= n)
		return;
	int ibin = (int)floor((psd[i] - pmin) * piwidth);
	ibin = min(max(ibin, 0), nbins - 1);
	atomicAdd(histo + ibin, 1);
}  // histo_kernel
// instantiations
template __global__ void histo_kernel
(float* psd, int n, double pmin, double piwidth, int *histo, int nbins);
template __global__ void histo_kernel
(double* psd, int n, double pmin, double piwidth, int *histo, int nbins);

/** bins in shared memory
		@param nels number of input elements to be processed by each input thread
 */
extern __shared__ int lhisto[];
template<class T>
__global__ void local_histo_kernel
(T* psd, int n, double pmin, double piwidth, int *histo, int nbins, int nels) {
	// zero shared memory bins
	int ii = threadIdx.x;
	int bs = blockDim.x;
	for(int ibin = ii; ibin < nbins; ibin += bs) 
		lhisto[ibin] = 0;
	__syncthreads();
	// compute the histogram in shared memory
	int istart = ii + blockIdx.x * bs * nels;
	int iend = min(ii + (blockIdx.x + 1) * bs * nels, n);
	for(int i = istart; i < iend; i += bs) {
		int ibin = (int)floor((psd[i] - pmin) * piwidth);
		ibin = min(max(ibin, 0), nbins - 1);
		atomicAdd(lhisto + ibin, 1);		
	}
	__syncthreads();
	// accumulate the histogram into the global memory
	for(int ibin = ii; ibin < nbins; ibin += bs)
		atomicAdd(histo + ibin, *(lhisto + ibin));
}  // local_histo_kernel

template<class T> void MafiaSolver<T>::compute_histo_dev(int idim) {
	int nbins = nbinss[idim];
	// on-device data for a histogram
	int *d_histo;
	size_t histo_sz = sizeof(*d_histo) * nbins;
	CHECK(cudaMalloc((void**)&d_histo, histo_sz));
	CHECK(cudaMemset(d_histo, 0, histo_sz));
	// kernel call
	// TODO: handle more than 256M points
	size_t bs = 256;
	//histo_kernel<<<divup(n, bs), bs>>>
	//	(d_ps + idim * n, n, pmins[idim], nbins / (pmaxs[idim] - pmins[idim]), 
	//	 d_histo, nbins);	
	int nels = 16;
	local_histo_kernel<<<divup(n, bs * nels), bs, nbins * sizeof(int)>>>
		(d_ps + idim * n, n, pmins[idim], nbins / (pmaxs[idim] - pmins[idim]), 
		 d_histo, nbins, nels);	
	CHECK(cudaDeviceSynchronize());
	
	// copy the data back
	CHECK(cudaMemcpy(histos[idim], d_histo, histo_sz, cudaMemcpyDeviceToHost));
	cudaFree(d_histo);
}  // compute_histo_dev

/** kernel to compute the bitmaps on device 
		@param bmp bitmap data on device
		@param nwords number of 32-bit words in the bitmap
		@param ps device point data for thebitmap's dimension
		@param n the number of points
		@param pleft the left boundary of the window's range (inclusive)
		@param pright the right boundary of the window's range (non-inclusive)
 */
template<class T> __global__ void bitmap_kernel
(unsigned *bmp, int nwords, T *ps, int n, T pleft, T pright) {
	//int iword = threadIdx.x + blockIdx.x * blockDim.x;
	// TODO: use local memory
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if(i >= n)
		return;
	int bit = pleft <= ps[i] && ps[i] < pright ? 1 : 0;
	int shift = i % (sizeof(int) * 8);
	atomicOr(bmp + i / (sizeof(int) * 8), bit << shift);
	/* int iword = threadIdx.x + blockIdx.x * blockDim.x;
	if(iword >= nwords)
		return;
	// accumulate the data word over the points
	int word = 0; 
	int istart = iword * sizeof(int) * 8; 
	int iend = min(n, (iword + 1) * (int)sizeof(int) * 8);
	for(int i = istart; i < iend; i++) {
		int bit = pleft <= ps[i] && ps[i] < pright ? 1 : 0;
		word |= bit << (i - istart);
	}  // for(i)
	bmp[iword] = word; */
}  // bitmap_kernel
// instantiations
template __global__ void bitmap_kernel
(unsigned *bmp, int nwords, float *ps, int n, float pleft, float pright);
template __global__ void bitmap_kernel
(unsigned *bmp, int nwords, double *ps, int n, double pleft, double pright);

template<class T> 
void MafiaSolver<T>::compute_bitmap_dev(int idim, int iwin) {
	Window &w = windows[idim][iwin];
	// allocate memory on device
	unsigned *d_bmp;
	int nwords = w.pset->n, bmp_sz = nwords * sizeof(*d_bmp);
	CHECK(cudaMalloc((void**)&d_bmp, bmp_sz));
	CHECK(cudaMemset(d_bmp, 0, bmp_sz));
	// call the kernel
	size_t bs = 256;
	//bitmap_kernel<<<divup(nwords, bs), bs>>>
	//	(d_bmp, nwords, d_ps + n * idim, n, (T)w.pleft, (T)w.pright);
	bitmap_kernel<<<divup(n, bs), bs>>>
		(d_bmp, nwords, d_ps + n * idim, n, (T)w.pleft, (T)w.pright);
	cudaDeviceSynchronize();

	// copy data back
	CHECK(cudaMemcpy(w.pset->data, d_bmp, bmp_sz, cudaMemcpyDeviceToHost));
	cudaFree(d_bmp);
}  // compute_bitmap_dev

template<class T>
void MafiaSolver<T>::free_dev_resources() {
	cudaFree(d_ps);
}  // free_dev_resources

// explicit instantiations
template class MafiaSolver<float>;
template class MafiaSolver<double>;
  
#endif
