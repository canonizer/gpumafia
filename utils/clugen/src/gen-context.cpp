/** @file gen-context.cpp implementation of GenContext class */

#include "cluster.h"
#include "gen-context.h"

#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

GenContext::GenContext() :
	n(10000), d(2), kmin(2), kmax(2), m(1), ct(BOX), cf(0.9), pmin(0), 
	pmax(100), csize(10), clusters(0), clusters_not_intersect(false)
{}  // GenContext

GenContext::~GenContext() {
	if(clusters) {
		for(int i = 0; i < m; i++)
			delete clusters[i];
	}
	delete[] clusters;
}  // ~GenContext

void GenContext::generate_template(void) {
	clusters = new Cluster*[m];
	memset(clusters, 0, sizeof(*clusters) * m);
	for(int i = 0; i < m; i++) {
		int k1 = irandom(kmin, kmax);
		clusters[i] = new Cluster(this, k1, cf / m, ct);
		clusters[i]->generate_template();
	}
}  // generate_template()

void GenContext::generate_point(double *p) const {
	// "cluster or background noise"
	if(drandom(0, 1) < cf) {
		// cluster; select one
		int ic = min(max((int)floor(drandom(0, 1) * m), 0), m - 1);
		clusters[ic]->generate_point(p);
	} else {
		// noise
		generate_other_point_dims(p, 0, 0);
	}
}  // generate_point

void GenContext::generate_other_point_dims
(double *p, const int *cluster_dims, int k)	const {
	int jdim = 0;
	for(int idim = 0; idim < d; idim++) {
		if(cluster_dims && jdim < k && cluster_dims[jdim] == idim)
			jdim++;
		else
			p[idim] = drandom(pmin, pmax);
	}
}  // generate_other_point_dims

double *GenContext::generate_points() const {
	double *ps = new double[n * d];
	for(int i = 0; i < n; i++)
		generate_point(ps + i * d);
	return ps;
}  // generate_points

void GenContext::print_info(void) const {
	// print the number of clusters, and then the clusters
	printf("%d clusters, %d%% points noise\n", m, (int)round((1 - cf) * 100));
	for(int i = 0; i < m; i++) {
		printf("cluster %d: ", i);
		clusters[i]->print_info();
		printf("\n");
	}
}  // print_info

bool GenContext::intersects_with(const Cluster* clu2) const {
	for(int iclu = 0; iclu < m; iclu++) {
		Cluster *clu = clusters[iclu];
		if(clu && clu != clu2 && clu->intersects_with(*clu2))
			return true;
	}
	return false;
}  // intersects_with

double drandom(double a, double b) {
	return a + (b - a) * ((double)random() / RAND_MAX);
}  // drandom

int irandom(int a, int b) {
	return a + (random() % (b - a + 1));
}  // irandom
