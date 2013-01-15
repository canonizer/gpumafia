#ifndef MAFIA_CLUSTER_H_
#define MAFIA_CLUSTER_H_

/** @file cluster.h a class which generates a single cluster */

// forward references
struct GenContext;

/** cluster type, box or sphere */
enum {BOX, SPHERE};

/** generates a single cluster */
struct Cluster {

	/** the generator context to which the cluster belongs */	
	GenContext *ctx;
	/** cluster dimensionality */
	int k;
	/** fraction of points falling inside the cluster */
	double cf;
	/** cluster type, box or sphere */
	int ct;
	/** dimension numbers of the cluster, [k] */
	int *dims;
	/** the minimum extents of the cluster, [k] */
	double *cmin;
	/** the maximum extents of the cluster, [k] */
	double *cmax;
	/** the midpoint of the cluster, [k] */
	double *cmid;
	/** the half-axes of the cluster in each dimension, [k] */
	double *crad;
	/** the inverted square radii in each dimension, [k] */
	double *cinvr2;
	
	/** the constructor for the cluster */
	Cluster(GenContext *ctx, int k, double cf, int ct);

	/** generates the template for the cluster */
	void generate_template(void);

	/** generates a single point */
	void generate_point(double *p) const;

	/** prints information about the cluster */
	void print_info(void) const;

	/** the destructor */
	~Cluster();

private:

	/** generates a point for cluster dimensions in a box */
	void generate_box_point(double *p) const;
	
	/** generates a point for cluster dimensions in a sphere */
	void generate_sphere_point(double *p) const;

};  // struct Cluster

#endif
