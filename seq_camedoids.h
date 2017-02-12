#ifndef SEQ_CAMEDOIDS
#define SEQ_CAMEDOIDS

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

void init_clusters(int nclusters, int nelements, int centroids[]);
double euclid(int n, double data[][n], int index1, double *index2);
double euclid_dynamic(double **data, int n, int index1, double *index2);
int find_nearest_cluster(int numClusters, int ncoords, int object,
                         double objects[][ncoords], double cluster[][ncoords]);
int find_nearest_clusterdd(int numClusters, int ncoords, int object,
                           double **objects, double **clstr);
double *mean_distance(int numObjs, int numCoords, double objects[][numCoords],
                      int *membership, double means[][numCoords]);
double **seq_camedoids2(int numObjs, int numCoords, double objects[][numCoords],
                        int numClusters, int *membership, int *clusterSize,
                        double clusters[][numCoords], int *delta, int *flag);
double seq_camedoids(int numObjs, int numCoords, double **objects,
                     int numClusters, int nCands, int *membership,
                     double clusters[][numCoords], double *d, int delta);

#endif
