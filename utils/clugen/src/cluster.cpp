/** @file cluster.cpp implementation of cluster generation */

#include "cluster.h"
#include "gen-context.h"

#include <algorithm>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

const double NO_INTERS_DISTANCE = 3;
const int NO_INTERS_TRIES = 10001;

Cluster::Cluster(GenContext *ctx, int k, double cf, int ct) : 
	ctx(ctx), dims(0), cmin(0), cmax(0), cmid(0), crad(0), cinvr2(0),
  k(k), cf(cf), ct(ct) 
{}  // Cluster

Cluster::~Cluster() {
	delete[] dims;
	delete[] cmin;
	delete[] cmax;
	delete[] cmid;
	delete[] crad;
	delete[] cinvr2;
}  // ~Cluster

void Cluster::generate_template(void) {
	// dimension numbers
	dims = new int[k];
	if(ctx->m == 1) {
		// one cluster, no randomness
		for(int iidim = 0; iidim < k; iidim++)
			dims[iidim] = iidim;
	} else {
		// random dimensions
		for(int iidim = 0; iidim < k; iidim++) {
			bool is_different;
			do {
				dims[iidim] = random() % ctx->d;
				is_different = true;
				for(int jjdim = 0; jjdim < iidim; jjdim++) 
					if(dims[iidim] == dims[jjdim]) {
						is_different = false;
						break;
					}
			} while(!is_different);
		}
		// sort the dimensions
		sort(dims, dims + k);
	}  // if(single cluster)-else

	// cluster position 
	cmin = new double[k];
	cmax = new double[k];
	for(int iidim = 0; iidim < k; iidim++) {
		cmin[iidim] = 0;
		cmax[iidim] = -1;
	}
	cmid = new double[k];
	crad = new double[k];
	cinvr2 = new double[k];
	for(int iidim = 0; iidim < k; iidim++) {
		int ntries = 0;
		do {
			cmin[iidim] = drandom(ctx->pmin, ctx->pmax - ctx->csize);
			cmax[iidim] = cmin[iidim] + ctx->csize;
		} while(ctx->clusters_not_intersect && ctx->intersects_with(this) && 
						++ntries < NO_INTERS_TRIES);
		if(ntries >= NO_INTERS_TRIES) {
			fprintf(stderr, "cannot generate a non-intersecting cluster\n");
			exit(-1);
		}
		cmid[iidim] = (cmax[iidim] + cmin[iidim]) / 2;
		crad[iidim] = (cmax[iidim] - cmin[iidim]) / 2;
		cinvr2[iidim] = 1 / (crad[iidim] * crad[iidim]);
	}  // for(iidim)
	
}  // generate_template

void Cluster::generate_point(double *p) const {
	ctx->generate_other_point_dims(p, dims, k);
	switch(ct) {
	case BOX:
		generate_box_point(p);
		break;
	case SPHERE:
		generate_sphere_point(p);
		break;
	default:
		assert(0);
		break;
	}  // switch(ct)
}  // generate_point

bool Cluster::intersects_with(const Cluster &clu2) const {
	int iidim1, iidim2, k1 = k, k2 = clu2.k;
	for(iidim1 = 0, iidim2 = 0; iidim1 < k1 && iidim2 < k2; ) {
		int idim1 = dims[iidim1], idim2 = clu2.dims[iidim2];
		if(idim1 == idim2) {
			double cmin1 = cmin[iidim1], cmax1 = cmax[iidim1], 
				cmin2 = clu2.cmin[iidim2], cmax2 = clu2.cmax[iidim2];
			if(cmin1 < cmax1 && cmin2 < cmax2) {
				// valid ranges
				double min_cmax = min(cmax1, cmax2), max_cmin = max(cmin1, cmin2);
				if(max_cmin - min_cmax < NO_INTERS_DISTANCE)
					return true;
			}
			iidim1++, iidim2++;
		} else if(idim1 < idim2)
			iidim1++;
		else
			iidim2++;
	}  // for(iidim)
	return false;
}  // intersects_with

void Cluster::print_info(void) const {
	// print probability, dimensions and extents
	printf("%s (%d, f=%.3lf) [", ct == BOX ? "box" : "sphere", k, cf);
	for(int iidim = 0; iidim < k; iidim++)
		printf("%d:%.2lf..%.2lf%s", dims[iidim], cmin[iidim], 
					 cmax[iidim], iidim < k - 1 ? " " : "]");
}  // print_info

void Cluster::generate_box_point(double *p) const {
	for(int iidim = 0; iidim < k; iidim++) {
		int idim = dims[iidim];
		p[idim] = drandom(cmin[iidim], cmax[iidim]);
	}
}  // generate_box_point

void Cluster::generate_sphere_point(double *p) const {
	// just generate point inside a box until it falls inside the cluster
	double dist2;
	do {
		generate_box_point(p);
		dist2 = 0;
		for(int iidim = 0; iidim < k; iidim++) {
			int idim = dims[iidim];
			double dcoord = p[idim] - cmid[idim];
			dist2 += dcoord * dcoord * cinvr2[idim];
		}
	} while(dist2 > 1);
}  // generate_sphere_point
