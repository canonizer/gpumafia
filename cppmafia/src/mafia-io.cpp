/** @file mafia-io.cpp implementation of MAFIA helper I/O routines */

#include "mafia-io.h"
#include "utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

template<class T> void read_points
(const char* path, T** points, int* npoints, int* ndims) {
  // open file
  FILE *file = fopen(path, "r");
  if(!file) {
    fprintf(stderr, "cannot open file %s\n", path);
    exit(-1);
  }
  
  // read the file, line-by-line
  int line_num = 0;
  size_t buf_size = 100;
  char *line_buf = (char*)malloc(sizeof(char) * buf_size);
  vector<double> vpoints;
  *npoints = 0;
  *ndims = 0;
  while(getline(&line_buf, &buf_size, file) > 0) {
    line_num++;
    double coord;
    int cur_ndims = 0;
    char *buf_ptr = line_buf, *next_buf_ptr = 0;
    while(1) {
      coord = strtod(buf_ptr, &next_buf_ptr);
      if(buf_ptr == next_buf_ptr)
				break;
      buf_ptr = next_buf_ptr;
      cur_ndims++;
      vpoints.push_back(coord);
    }
    if(cur_ndims) {
      ++*npoints;
      // non-empty line
      if(*ndims && *ndims != cur_ndims) {
				// point of wrong dimensionality
				fprintf
					(stderr, "point at line %d has wrong dimensionality: %d instead of %d\n",
					 line_num, cur_ndims, *ndims);
				exit(-1);
      } else if(!*ndims) {
				// first non-empty line, set the number of dimensions
				*ndims = cur_ndims;
      }
    }
  }  // while(lines)

  // check the vector length
  if(vpoints.size() != *ndims * *npoints) {
    fprintf(stderr, "wrong number of floating-point values in file\n");
    exit(-1);
  }

  // just fill in the points
	int n = *npoints, d = *ndims;
  *points = (T*)malloc(sizeof(**points) * *ndims * *npoints);
	T *ps = *points;
  //for(size_t icoord = 0; icoord < *ndims * *npoints; icoord++)
  //  (*points)[icoord] = vpoints[icoord];
	for(int i = 0; i < n; i++)
		for(int idim = 0; idim < d; idim++)
			PS(i, idim) = vpoints[i * d + idim];

  // close the file
  fclose(file);
}  // read_points

template<class T> void write_clusters
(const char* path, T* points, int npoints, int ndims, 
 const vector<vector<int> > &cluster_idxs) {
  // number length
  int num_len = sizeof(cluster_idxs.size()) * 3;
  // buffer for file names
  int name_len = strlen(path) + num_len + strlen("-.dat");
  char dat_name[name_len + 1], idx_name[name_len + 1];
  for(int icluster = 0; icluster < cluster_idxs.size(); icluster++) {
    // open files
    //sprintf(dat_name, "%s-%d.dat", path, icluster);
    sprintf(idx_name, "%s-%d.idx", path, icluster);
    //FILE *dat_file = fopen(dat_name, "w");
    FILE *idx_file = fopen(idx_name, "w");
    // if(!dat_file || !idx_file) {
		if(!idx_file) {
      fprintf(stderr, "cannot write cluster data to file\n");
      exit(-1);
    }
    
    // write the cluster data
    const vector<int> & cluster = cluster_idxs[icluster];
    for(int ipidx = 0; ipidx < cluster.size(); ipidx++) {
      int pidx = cluster[ipidx];
      //T *point = points + pidx * ndims;
      // for(int idim = 0; idim < ndims; idim++) {
			// 	fprintf(dat_file, "%.9lf", (double)point[idim]);
			// 	if(idim < ndims - 1)
			// 		fprintf(dat_file, "\t");
      // }  // for each point coordinate
      fprintf(idx_file, "%d\n", pidx);
      // fprintf(dat_file, "\n");
    }  // for each point in cluster

    // close the files
    // fflush(dat_file);
    // fclose(dat_file);
    fflush(idx_file);
    fclose(idx_file);
  }  // for each cluster

}  // write_clusters

// explicit instantiations, required to compile
template void read_points<float>
(const char* path, float** points, int* npoints, int* ndims);

template void read_points<double>
(const char* path, double** points, int* npoints, int* ndims);

template void write_clusters<float>
(const char* path, float* points, int npoints, int ndims, 
 const vector<vector<int> > &cluster_idxs);

template void write_clusters<double>
(const char* path, double* points, int npoints, int ndims, 
 const vector<vector<int> > &cluster_idxs);
