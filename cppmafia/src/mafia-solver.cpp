#include "cdu.h"
#include "mafia-solver.h"
#include "timing.h"
#include "utils.h"

#include <algorithm>
#include <limits>
#include <vector>
#include <queue>
#include <set>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

template<class T>
MafiaSolver<T>::MafiaSolver
(const T *points, int n, int d, const Options &opts) : 
	d_ps(0), pmins(0), pmaxs(0), nbinss(0), histo_data(0), histos(0), bmps(0), 
	d_bmps(0) {
  // initialize the class
  this->ps = points;
  this->n = n;
  this->d = d;
  this->min_nbins = opts.min_nbins;
  this->min_nwindows = opts.min_nwindows;
  this->max_nwindows = opts.max_nwindows;
  this->alpha = opts.alpha;
  this->beta = opts.beta;
  this->eps = 1e-8;
	this->flags = opts.flags;
  
#ifdef MAFIA_USE_DEVICE
	// initialize the device if it is used
	if(use_device())
		touch_dev();
#endif
}  // MafiaSolver

template<class T>
vector<vector<int> > MafiaSolver<T>::find_clusters() {
  vector<vector<int> > cluster_idxs;
	start_phase(PhaseBuildHisto);
  build_histos();
	if(is_verbose())
	 	print_histos();
	start_phase(PhaseBuildWindows);
  build_windows();
	if(is_verbose())
	 	print_windows();
	start_phase(PhaseBuildBitmaps);
	if(use_bitmaps()) {
		build_bitmaps();
		// if(verbose())
		//   print_bitmaps();
	}
  cur_dim = 0;
  do {
		if(is_verbose())
			printf("dimension %d\n", cur_dim);
		start_phase(PhaseFindCdus);
		find_cdus();
		if(is_verbose()) {
		 	//printf("CDUs: ");
		 	//print_dus(cdus);
		}
		if(cur_dim > 0) {
			dedup_cdus();
			if(is_verbose()) {
			 	//printf("dedupped CDUs: ");
			 	//print_dus(cdus);
			}
		}
		start_phase(PhaseFindDense);
		find_dense_cdus();
		if(is_verbose()) {
			//printf("found DUs: ");
		 	//print_dus(new_dus);
		}
		if(cur_dim > 0) {
			start_phase(PhaseFindUnjoined);
			find_unjoined_dus();
			if(is_verbose()) {
				printf("new terminal DUs: ");
				print_dus(terminal_dus[cur_dim - 1]);
			}
		}
		cur_dus = new_dus;
    cur_dim++;
  } while(cdus.size() > 0);
	// find connected components in the graph
	start_phase(PhaseBuildGraph);
	build_du_graph();
	build_du_clusters();
	start_phase(PhaseBuildClusters);
	// build clusters as index lists
	build_clusters();
	return clusters;
}  // find_clusters()

template <class T>
MafiaSolver<T>::~MafiaSolver() {
#ifdef MAFIA_USE_DEVICE
	if(use_device())
		free_dev_resources();
#endif
  // free the memory
  free(pmins);
  free(pmaxs);
	free(nbinss);
	free(histos);
	free(histo_data);
	bulk_free(bmps);
}  // ~MafiaSolver

template<class T>
void MafiaSolver<T>::build_histos() {
#ifdef MAFIA_USE_DEVICE
	if(use_device())
		copy_ps_to_device();
#endif

  // compute minima and maxima for per-dimension coordinates
  pmins = (T*)malloc(sizeof(*pmins) * d);
  pmaxs = (T*)malloc(sizeof(*pmaxs) * d);
	nbinss = (int *)malloc(sizeof(*nbinss) * d);
#ifdef MAFIA_USE_DEVICE
	if(use_device())
		compute_limits_dev();
	else
#endif
		compute_limits_host();
	
	// count the number of bins in each histogram
	for(int idim = 0; idim < d; idim++) {
    // increase max by eps, to avoid special handling of the right corner, 
    // and make all bins and windows right-exclusive, i.e. [lb, rb)
    pmaxs[idim] += eps;
		nbinss[idim] = max(min_nbins, (int)ceilf(pmaxs[idim] - pmins[idim]));
	}  // for(idim)

	// allocate the data for the entire histogram
	int total_nbins = 0;
	for(int idim = 0; idim < d; idim++) 
		total_nbins += nbinss[idim];
	histo_data = (int*)malloc(sizeof(*histo_data) * total_nbins);
	memset(histo_data, 0, sizeof(*histo_data) * total_nbins);

	// set the pointer for each histogram
	histos = (int**)malloc(sizeof(*histos) * d);
	int sum_nbins = 0;
	for(int idim = 0; idim < d; idim++) {
		histos[idim] = histo_data + sum_nbins;
		sum_nbins += nbinss[idim];
	}
  
	// compute the point histograms
  for(int idim = 0; idim < d; idim++) {
#ifdef MAFIA_USE_DEVICE
		if(use_device())
			compute_histo_dev(idim);
		else
#endif
			compute_histo_host(idim);
  } 
}  // build_histos

template<class T>
void MafiaSolver<T>::compute_limits_host() {
 for(int idim = 0; idim < d; idim++) {
    pmins[idim] = numeric_limits<T>::infinity(); 
    pmaxs[idim] = -numeric_limits<T>::infinity();
    for(int i = 0; i < n; i++) {
      pmins[idim] = min(pmins[idim], PS(i, idim));
      pmaxs[idim] = max(pmaxs[idim], PS(i, idim));
    }
  }  // for(idim)
}  // compute_limits_host

template<class T>
void MafiaSolver<T>::compute_histo_host(int idim) {
	// TODO: support unnormalized data; for now, assume that maximum
	// width of a bin is 1
	int nbins = nbinss[idim];
	T bin_iwidth = nbins / (pmaxs[idim] - pmins[idim]);
	int *histo = histos[idim];
	for(int i = 0; i < n; i++) {
		int ibin = (int)floor((PS(i, idim) - pmins[idim]) * bin_iwidth);
		ibin = min(max(ibin, 0), nbins - 1);
		histo[ibin]++;
	}
}  // compute_histo_host

template<class T>
void MafiaSolver<T>::build_uniform_windows
(int idim, int nwindows, vector<Window > &ws) {
  int *histo = histos[idim];
  int nbins = nbinss[idim];;
  int wwidth = divup(nbins, nwindows);
  for(int left = 0; left < nbins; left += wwidth) {
    int right = min(left + wwidth, nbins);
    int width = right - left;
    // find the maximum number of elements in a window bin
    int wmax = 0;
    for(int ibin = left; ibin < right; ibin++) 
      wmax = max(wmax, histo[ibin]);
    Window w = Window(idim, left, width, wmax);
    w.pleft = pmins[idim] + left * (pmaxs[idim] - pmins[idim]) / nbins;
    w.pright = pmins[idim] + right * (pmaxs[idim] - pmins[idim]) / nbins;
    ws.push_back(w);
  }  // for each uniform window width interval
}  // build_uniform_windows

template<class T>
void MafiaSolver<T>::build_bitmaps() {
	int nwindows = dense_ws.size();
	nwords = divup(n, sizeof(*bmps) * 8);
	size_t bmp_sz = sizeof(*bmps) * nwindows * nwords;
	bmps = (unsigned *)bulk_alloc(bmp_sz);
	memset(bmps, 0, bmp_sz);
#ifdef MAFIA_USE_DEVICE
	if(use_device())
		alloc_bitmaps_dev();
#endif
	for(int iwin = 0; iwin < nwindows; iwin++) {
#ifdef MAFIA_USE_DEVICE
		if(use_device())
			compute_bitmap_dev(iwin);
		else
#endif
			compute_bitmap_host(iwin);
	}
}  // build_bitmaps

template<class T>
void MafiaSolver<T>::compute_bitmap_host(int iwin) {
	// check all points for membership (only the coordinate in this 
	// dimension)
	Window &w = dense_ws[iwin];
	int idim = w.idim;
	for(int i = 0; i < n; i++) {
		int iword = i / (sizeof(*bmps) * 8);
		int ibit = i % (sizeof(*bmps) * 8);
		unsigned bit = w.pleft <= PS(i, idim) && PS(i, idim) < w.pright ? 1u : 0u;
		BMPS(iwin, iword) |= bit << ibit;
	}
}  // compute_bitmap_host()

template<class T>
void MafiaSolver<T>::build_windows() {

	// build windows for each dimension
  for(int idim = 0; idim < d; idim++) {
    if(pmaxs[idim] - pmins[idim] > eps) {

      // start with non-adaptive windows
      vector<Window> initial_windows;
      build_uniform_windows(idim, max_nwindows, initial_windows);

      // merge into adaptive windows
      windows.push_back(vector<Window>());
      vector<Window> &dim_windows = windows[idim];
      dim_windows.push_back(initial_windows[0]);
      for(int iwin = 1; iwin < initial_windows.size(); iwin++) {
				Window &w2 = initial_windows[iwin];
				Window w = dim_windows.back();
				if(w.can_merge_with(w2, beta)) {
					dim_windows.pop_back();
					dim_windows.push_back(w.merge_with(w2));
				} else 
					dim_windows.push_back(w2);
      }  // for each initial window

      // handle uniform case (i.e., 1 window)
      if(dim_windows.size() == 1) {
				dim_windows.clear();
				build_uniform_windows(idim, min_nwindows, dim_windows);
      }

      // compute thresholds, also set window numbers
      int nbins = nbinss[idim];
      for(int iwin = 0; iwin < dim_windows.size(); iwin++) {
				dim_windows[iwin].compute_threshold(alpha, n, nbins);
				dim_windows[iwin].iwin = iwin;
			}

    } else {
      // no spread in dimension, just ignore
      windows.push_back(vector<Window>());
    }
  } // for each dimension

	// now separate windows which are dense
	for(int idim = 0; idim < d; idim++) {
		vector<Window> &dim_windows = windows[idim];
		for(int iwin = 0; iwin < dim_windows.size(); iwin++) {
			Window &w = dim_windows[iwin];
			if(w.is_dense())
				dense_ws.push_back(w);
		}
	}  // for(idim)
}  // build_windows

template<class T>
void MafiaSolver<T>::find_cdus() {
	if(cur_dim > 0) {
		// do search & merging of CDUs
		// for the case of set deduplication, deduplication is done during generation
		cdus.clear();
		set<ref<Cdu>, CduCmp> new_cdus;
		for(int i1 = 0; i1 < cur_dus.size(); i1++) {
			for(int i2 = i1 + 1; i2 < cur_dus.size(); i2++) {
				Cdu &du1 = *cur_dus[i1], &du2 = *cur_dus[i2];
				if(du1.can_merge_with(du2)) {
					// merge two dus
					ref<Cdu> new_cdu = du1.merge_with(du2);
					du1.flag = du2.flag = true;
					if(use_set_dedup())
						new_cdus.insert(new_cdu);
					else
						cdus.push_back(new_cdu);
				}
			}
		}
		// if set deduplication was used, everything must be copied back into the
		// vector 
		for(set<ref<Cdu>, CduCmp>::iterator pcdu = new_cdus.begin(); 
				pcdu !=	new_cdus.end(); pcdu++) {
			cdus.push_back(*pcdu);
		}
	} else {
		// just add all windows as CDUs
		for(int iwin = 0; iwin < dense_ws.size(); iwin++)
			cdus.push_back(new Cdu(dense_ws[iwin].idim, iwin));
	}  // if(cur_dim > 0)
}  // find_cdus

template<class T>
void MafiaSolver<T>::dedup_cdus() {
	if(!use_set_dedup())
		naive_dedup_cdus();
}  // dedup_cdus

template<class T>
void MafiaSolver<T>::naive_dedup_cdus() {
	// TODO: here and in some other places, avoid unnecessary vector copying
	// njoined field is used to mark duplicate CDUs
	vector<ref<Cdu> > new_cdus;
	for(int i1 = 0; i1 < cdus.size(); i1++) {
		Cdu &du1 = *cdus[i1];
		// check if marked as duplicate already
		if(du1.flag)
			continue;
		// no previous duplicates, push back this node
		new_cdus.push_back(cdus[i1]);
		for(int i2 = i1 + 1; i2 < cdus.size(); i2++) {
			Cdu &du2 = *cdus[i2];
			if(du2.flag)
				continue;
			du2.flag = du1 == du2;
		}  // for(cdu2)
	}  // for(cdu1)
	cdus = new_cdus;
}  // dedup_cdus

template<class T>
void MafiaSolver<T>::find_dense_cdus() {
	new_dus.clear();
	// for dim 0, point counting is always done on host
#ifdef MAFIA_USE_DEVICE
	if(use_device() && cur_dim > 0)
		count_points_dev();
	else
#endif
		count_points_host();
	for(int icdu = 0; icdu < cdus.size(); icdu++) {
		Cdu &cdu = *cdus[icdu];
		if(cdu.is_dense(dense_ws))
			new_dus.push_back(&cdu);
	}  // for(cdu)
}  // find_dense_cdus

template<class T>
void MafiaSolver<T>::count_points_host() {
	for(int icdu = 0; icdu < cdus.size(); icdu++) {
		Cdu &cdu = *cdus[icdu];
		if(use_bitmaps())
			cdu.count_points_bitmaps(nwords, bmps, dense_ws);
		else 
			cdu.count_points_direct(ps, n, d, dense_ws);
	}  // for(cdu)
}  // count_points_host

template<class T>
void MafiaSolver<T>::find_unjoined_dus() {
	terminal_dus.push_back(vector<ref<Cdu> > ());
	vector<ref<Cdu> > &dim_terminal_dus = terminal_dus[cur_dim - 1];
	for(int idu = 0; idu < cur_dus.size(); idu++) {
		// check whether du is unjoined or unassimilated
		Cdu &du = *cur_dus[idu];
		bool joined = true;
		if(du.flag) {
			joined = false;
			// check if assimilated into a new DU
			for(int inew_du = 0; inew_du < new_dus.size(); inew_du++) {
				if(du.is_assimilated_into(*new_dus[inew_du])) {
					joined = true;
					break;
				}
			}
		} else {
			// DU hasn't been even joined
			joined = false;
		}
		// add to terminal dense units
		if(!joined) {
			dim_terminal_dus.push_back(cur_dus[idu]);
		}
	}  // for(current du)
}  // find_unjoined_dus

template<class T>
void MafiaSolver<T>::build_du_graph() {
	for(int idim = 0; idim < terminal_dus.size(); idim++) {
		vector<ref<Cdu> > &dus = terminal_dus[idim];
		for(int idu1 = 0; idu1 < dus.size(); idu1++) {
			Cdu* pdu1 = dus[idu1];
			for(int idu2 = idu1 + 1; idu2 < dus.size(); idu2++) {
				Cdu* pdu2 = dus[idu2];
				if(pdu1->has_common_face_with(*pdu2, dense_ws)) {
					// add an edge in both directions
					pdu1->neighbours.push_back(pdu2);
					pdu2->neighbours.push_back(pdu1);
				}
			}  // for(idu2)
		}  // for(idu1)
	}  // for(idim)
}  // build_du_graph

template<class T>
void MafiaSolver<T>::build_du_clusters() {
	for(int idim = 0; idim < terminal_dus.size(); idim++) {
		vector<ref<Cdu> > &dus = terminal_dus[idim];
		// reset flags on DUs; it will be used to indicate "visited"
		for(int idu = 0; idu < dus.size(); idu++)
			dus[idu]->flag = false;
		// until there are unvisited DUs
		for(int idu = 0; idu < dus.size(); idu++) {
			if(dus[idu]->flag)
				continue;
			du_clusters.push_back(vector<ref<Cdu> > ());
			vector<ref<Cdu> > &cluster = du_clusters.back();
			// find connected component with breadth-first search
			queue<Cdu *> q;
			q.push(dus[idu]);
			while(!q.empty()) {
				Cdu *pdu = q.front();
				q.pop();
				pdu->flag = true;
				cluster.push_back(pdu);
				// push unattended neighbours to the queue
				for(int ineigh = 0; ineigh < pdu->neighbours.size(); ineigh++) {
					Cdu *pneigh = pdu->neighbours[ineigh];
					if(!pneigh->flag) 
						q.push(pneigh);
				} // for()
			}  // while(queue not empty)
		}  // for(idu unvisited)
	}  // for(idim)
}  // build_du_clusters

template<class T>
void MafiaSolver<T>::build_clusters() {
	for(int iclu = 0; iclu < du_clusters.size(); iclu++) {
		vector<ref<Cdu> > &du_cluster = du_clusters[iclu];
		clusters.push_back(vector<int>());
		vector<int> &cluster = clusters.back();
		if(use_bitmaps()) {
			for(int iword = 0; iword < nwords; iword++) {
				unsigned cword = 0u;
				// find point set, then transform it into a set of indices
				for(int idu = 0; idu < du_cluster.size(); idu++) {
					Cdu &du = *du_cluster[idu];
					unsigned duword = ~0u;
					for(int icoord = 0; icoord < du.coords.size(); icoord++) {
						dimpair_t dp = du.coords[icoord];
						duword &= BMPS(dp.win, iword);
					}
					cword |= duword;
				}  // for(idu)
				// now iterate through 1-bits in the word
				for(int ibit = 0; ibit < sizeof(cword) * 8; ibit++)
					if((cword >> ibit) & 1u)
						cluster.push_back(iword * (int)sizeof(cword) * 8 + ibit);
			}  // for(iword)
		} else {
			for(int i = 0; i < n; i++) {
				for(int idu = 0; idu < du_cluster.size(); idu++) {
					if(du_cluster[idu]->contains_point(ps, n, d, i, dense_ws)) {
						cluster.push_back(i);
						break;
					}					
				}  // for(idu)
			}  // for(i(point)) 
		}  // if(use_bitmaps())
	}  // for(iclu)
}  // build_clusters

template<class T>
void MafiaSolver<T>::print_histos() {
  for(int idim = 0; idim < d; idim++) {
    int *histo = histos[idim];
    printf("dimension %d: [", idim);
    for(int ibin = 0; ibin < nbinss[idim]; ibin++) {
      printf("%d", histo[ibin]);
      if(ibin != nbinss[idim] - 1)
				printf(" ");
    }
    printf("]\n");
  }
}  // print_histos

template<class T>
void MafiaSolver<T>::print_windows() {
  // print the windows
  for(int idim = 0; idim < d; idim++) {
    vector<Window> &dim_windows = windows[idim];
    printf("dimension %d: [", idim);
    for(int iwin = 0; iwin < dim_windows.size(); iwin++) {
      Window &w = dim_windows[iwin];      
      printf("(%d..%d max=%d t=%d)", w.left, w.right(), w.max, w.threshold);
      if(iwin != dim_windows.size() - 1)
				printf(" ");
    }
    printf("]\n");
  }
}  // print_windows

template<class T>
void MafiaSolver<T>::print_bitmaps() {
	printf("[\n");
	for(int iwin = 0; iwin < dense_ws.size(); iwin++) {
		Window &w = dense_ws[iwin];
		printf("[");
		for(int iword = 0; iword < nwords; iword++) {
			printf("%0x\n", BMPS(iwin, iword));
			if(iword < nwords - 1)
				printf(" ");
		}
		printf("]");
		if(iwin < dense_ws.size() - 1)
			printf("\n");
	}
	printf("]\n");
}  // print_bitmaps

template<class T>
void MafiaSolver<T>::print_terminal_dus() {
	for(int idim = 0; idim < terminal_dus.size(); idim++) {
		vector<ref<Cdu> > &dim_terminal_dus = terminal_dus[idim];
		printf("dimension %d: ", idim);
		print_dus(dim_terminal_dus);
	}
}  // print_terminal_dus()

template<class T>
void MafiaSolver<T>::print_dus(vector<ref<Cdu> > &dus) {
	printf("[ ");	
	for(int idu = 0; idu < dus.size(); idu++) {
		Cdu &du = *dus[idu];
		printf("[");
		for(int idp = 0; idp < du.coords.size(); idp++) {
			printf("%d:%d", du.coords[idp].dim, dense_ws[du.coords[idp].win].iwin);
			if(idp < du.coords.size() - 1)
				printf(" ");
		}
		printf("] ");
	}
	printf("]\n");	
}  // print_dus()

template<class T>
vector<vector<int> > mafia_solve
(const T *points, int npoints, int ndims, const Options &opts) {
  MafiaSolver<T> solver(points, npoints, ndims, opts);
  return solver.find_clusters();
}

// explicit instantiations
template class MafiaSolver<float>;
template class MafiaSolver<double>;

template
vector<vector<int> > mafia_solve<float>
(const float *points, int npoints, int ndims, const Options &opts);

template
vector<vector<int> > mafia_solve<double>
(const double *points, int npoints, int ndims, const Options &opts);
