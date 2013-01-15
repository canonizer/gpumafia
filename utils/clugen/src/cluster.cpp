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
	cmid = new double[k];
	crad = new double[k];
	cinvr2 = new double[k];
	for(int iidim = 0; iidim < k; iidim++) {
		cmin[iidim] = drandom(ctx->pmin, ctx->pmax - ctx->csize);
		cmax[iidim] = cmin[iidim] + ctx->csize;
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

void Cluster::print_info(void) const {
	// print probability, dimensions and extents
	printf("%s (%.3lf) [", ct == BOX ? "box" : "sphere", cf);
	for(int iidim = 0; iidim < k; iidim++)
		printf("%d:%.2lf..%.2lf%s", dims[iidim], cmin[iidim], 
					 cmax[iidim], iidim < k - 1 ? " " : "]");
}  // print_info

void Cluster::generate_box_point(double *p) const {
	for(int iidim = 0; iidim < k; iidim++) {
		int idim = dims[iidim];
		p[idim] = drandom(cmin[idim], cmax[idim]);
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
