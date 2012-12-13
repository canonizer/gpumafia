#ifndef MAFIA_IO_H_
#define MAFIA_IO_H_

#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;

template<class T> void read_points
(const char* path, T** points, int* npoints, int* ndims);

template<class T> void write_clusters
(const char* path, T* points, int npoints, int ndims, 
 const vector<vector<int> > &cluster_idxs);

#endif
