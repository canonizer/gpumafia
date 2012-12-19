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
	// compute limits using thrust
	device_ptr<T> d_psptr(d_ps);
	for(int idim = 0; idim < d; idim++) {
		pmins[idim] = reduce(d_psptr + idim * n, d_psptr + (idim + 1) * n, 
												 std::numeric_limits<T>::infinity(), minimum<T>());
		pmaxs[idim] = reduce(d_psptr + idim * n, d_psptr + (idim + 1) * n, 
												 -std::numeric_limits<T>::infinity(), maximum<T>());
	}	
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

template<class T> 
void MafiaSolver<T>::alloc_bitmaps_dev() {
	int nwindows = dense_ws.size();
	CHECK(cudaMalloc((void**)&d_bmps, sizeof(*d_bmps) * nwindows * nwords));
	CHECK(cudaMemset(d_bmps, 0, sizeof(*d_bmps) * nwindows * nwords));
}  // alloc_bitmaps_dev

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
void MafiaSolver<T>::compute_bitmap_dev(int iwin) {
	Window &w = dense_ws[iwin];
	int idim = w.idim;
	// allocate memory on device
	unsigned *h_bmp = bmps + iwin * nwords, *d_bmp = d_bmps + iwin * nwords;
	// call the kernel
	size_t bs = 256;
	//bitmap_kernel<<<divup(nwords, bs), bs>>>
	//	(d_bmp, nwords, d_ps + n * idim, n, (T)w.pleft, (T)w.pright);
	bitmap_kernel<<<divup(n, bs), bs>>>
		(d_bmp, nwords, d_ps + n * idim, n, (T)w.pleft, (T)w.pright);
	CHECK(cudaDeviceSynchronize());

	// copy data back
	int bmp_sz = nwords * sizeof(*d_bmp);
	CHECK(cudaMemcpy(h_bmp, d_bmp, bmp_sz, cudaMemcpyDeviceToHost));
}  // compute_bitmap_dev

__global__ void point_count_kernel
(int *pcounts, int *iwins, int ncdus, int ncoords, unsigned *bmps, int nwords,
 int nwords_pthr) {
	int icdu = threadIdx.y + blockIdx.y * blockDim.y;
	if(icdu >= ncdus)
		return;
	int bs = blockDim.x;
	int istart = threadIdx.x + blockIdx.x * bs * nwords_pthr;
	int iend = min(istart + nwords_pthr * bs, nwords);
	int pcount = 0;
	for(int iword = istart; iword < iend; iword += bs) {
		unsigned word = ~0u;
		for(int icoord = 0; icoord < ncoords; icoord++) {
			int iwin = iwins[icdu * ncoords + icoord];
			word &= bmps[iwin * nwords + iword];
		}
		pcount += __popc(word);
	}
	atomicAdd(pcounts + icdu, pcount);
}  // point_count_kernel

template<class T>
void MafiaSolver<T>::count_points_dev() {
	// copy CDU data to device, first aggregate on host; just window indices will do
	int ncdus = cdus.size();
	int ncoords = cur_dim + 1;
	// window numbers, in CDU-major order
	size_t win_sz = sizeof(int) * ncdus * ncoords;
	int *h_iwins = (int*)bulk_alloc(win_sz);
	for(int icdu = 0; icdu < ncdus; icdu++) {
		Cdu &cdu = *cdus[icdu];
		for(int icoord = 0; icoord < ncoords; icoord++)
			h_iwins[icdu * ncoords + icoord] = cdu.coords[icoord].win;
	}
	int *d_iwins; // window numbers, on device
	CHECK(cudaMalloc((void**)&d_iwins, win_sz));
	CHECK(cudaMemcpy(d_iwins, h_iwins, win_sz, cudaMemcpyHostToDevice));
	bulk_free(h_iwins);
	
	// run the kernel on device
	int *h_pcounts, *d_pcounts;  // point counts on host and device
	size_t pcount_sz = sizeof(*h_pcounts) * ncdus;
	CHECK(cudaMalloc((void**)&d_pcounts, pcount_sz));
	CHECK(cudaMemset(d_pcounts, 0, pcount_sz));
	int nwords_pthr = min(max(nwords / 64, 2), 64); // number of words per thread
	dim3 bs(64, 8);  // block size
	// iterate over CDU parts
	int ncdus_ppart = 32768 * bs.y;
	int ncdu_parts = divup(ncdus, ncdus_ppart);
	for(int icdu_part = 0; icdu_part < ncdu_parts; icdu_part++) {
		int cur_ncdus = min(ncdus_ppart, ncdus - icdu_part * ncdus_ppart);
		dim3 grid(divup(nwords, nwords_pthr * bs.x), divup(cur_ncdus, bs.y));
		// TODO: support more than 2**18-4 CDUs
		point_count_kernel<<<grid, bs>>>
			(d_pcounts + icdu_part * ncdus_ppart, d_iwins + icdu_part * ncdus_ppart, cur_ncdus, ncoords, 
			 d_bmps, nwords, nwords_pthr);
	}
	CHECK(cudaThreadSynchronize());

	// copy data back
	h_pcounts = (int*)bulk_alloc(pcount_sz);
	CHECK(cudaMemcpy(h_pcounts, d_pcounts, pcount_sz, cudaMemcpyDeviceToHost));
	for(int icdu = 0; icdu < ncdus; icdu++)
		cdus[icdu]->npoints = h_pcounts[icdu];

	// free everything
	cudaFree(d_iwins);
	cudaFree(d_pcounts);
	bulk_free(h_pcounts);
}  // count_points_dev

template<class T>
void MafiaSolver<T>::free_dev_resources() {
	cudaFree(d_ps);
	cudaFree(d_bmps);
}  // free_dev_resources

// explicit instantiations
template class MafiaSolver<float>;
template class MafiaSolver<double>;
  
#endif
