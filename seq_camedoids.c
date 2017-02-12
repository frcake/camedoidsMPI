#include "seq_camedoids.h"

/* Euclid distance between two multi-dimensional points  */
/*min_dist = euclid(ncoords, objects, object, cluster[0]);*/
double euclid(int n,            /* no. coords / features */
              double data[][n], /*[numObjs][numCoords] objects */
              int index1,       /*position of the first element*/
              double *index2 /*attributes's vector of the second element*/) {
  double result = 0.;
  int i;

  for (i = 0; i < n; i++) {
    double term = data[index1][i] - index2[i];
    result += term * term;
  }

  result = sqrt(result);

  return result;
}

double euclid_dynamic(double **data, int n,

                      int index1, double *index2) {
  double result = 0.;
  int i;

  for (i = 0; i < n; i++) {
    double term = data[index1][i] - index2[i];
    result += term * term;
  }

  result = sqrt(result);

  return result;
}

/* Initialize medoids with random elements of the dataset */

/* find the cluster id that has min distance to object */
/*index =find_nearest_cluster(numClusters, numCoords, i, objects, clusters);*/
int find_nearest_cluster(
    int numClusters, /* no. clusters */
    int ncoords,     /* no. coordinates */
    int object,      /* position of the specific object */
    double objects[][ncoords] /*[numObjs][numCoords] */,
    double cluster[][ncoords] /*[numClusters][numCoords] medoids*/) {
  int index, i;
  double dist, min_dist;
  double **clstr;
  index = 0;
  /*ncoords = numCoords*/
  clstr = malloc(sizeof(double *) * numClusters);
  for (int qq = 0; qq < numClusters; qq++) {
    clstr[qq] = malloc(sizeof(double) * ncoords);
  }

  for (int qq = 0; qq < numClusters; qq++) {
    for (int ww = 0; ww < ncoords; ww++) {
      clstr[qq][ww] = cluster[qq][ww];
      // printf("%f\n", clstr[qq][ww]);
    }
  }

  min_dist = euclid(ncoords, objects, object, clstr[0]);
  // printf("%f \t", min_dist);

  for (i = 1; i < numClusters; i++) {

    dist = euclid(ncoords, objects, object, clstr[i]);

    if (dist < min_dist) {
      min_dist = dist;
      index = i;
    }
  }
  for (int qq = 0; qq < numClusters; qq++) {
    free(clstr[qq]);
  }
  free(clstr);
  return (index);
}

int find_nearest_clusterdd(
    int numClusters, /* no. clusters */
    int ncoords,     /* no. coordinates */
    int object,      /* position of the specific object */
    double **objects /*[numObjs][numCoords] */,
    double **clstr /*[numClusters][numCoords] medoids*/) {
  int index, i;
  double dist, min_dist;

  index = 0;
  min_dist = euclid_dynamic(objects, ncoords, object, clstr[0]);

  for (i = 1; i < numClusters; i++) {

    dist = euclid_dynamic(objects, ncoords, object, clstr[i]);

    if (dist < min_dist) {
      min_dist = dist;
      index = i;
    }
  }

  return (index);
}

double *mean_distance(
    int numObjs /*eg 100k no. objects */,
    int numCoords /* no. features / dimensions */,
    double objects[][numCoords] /* in: [numObjs][numCoords]*/,
    int *membership /* out: [numObjs] */, double means[][numCoords]
    /* out: [numClusters][numCoords] medoids (position of medoids into the dataset) */) {
  int i, j;
  double term, d;
  double *distance;
  distance = calloc(numObjs, sizeof(double));

  for (i = 0; i < numObjs; i++) {
    /*calculate the distance of the current element of its cluster mean*/
    d = 0;
    for (j = 0; j < numCoords; j++) {
      double term = objects[i][j] - means[membership[i]][j];

      d += term * term;
    }
    /*euclidian distance from each object to the calculated cluster mean*/

    distance[i] = sqrt(d);
    // printf("\n %f ", distance[i]);
  }
  return distance;
}

double **
seq_camedoids2(int numObjs /*eg 100k no. objects */,
               int numCoords /* no. features / dimensions */,
               double objects[][numCoords] /* in: [numObjs][numCoords]*/,
               int numClusters /* no. clusters = 10*/,
               // int nCands /* no. candidates */,
               int *membership /* out: [numObjs] */, int *clusterSize,
               double clusters[][numCoords],
               /* out: [numClusters][numCoords] medoids (position of
                  medoids into the dataset) */
               int *delta, int *flag) {
  int i, j, k, rep, index, loop = 0, pos;
  // int *clusterSize;
  /* [numClusters]: no. objects assigned in each new cluster */

  // double **clusters_means;
  double(*cluster_distance_sum)[numCoords] =
      malloc(sizeof(double[(numClusters)][numCoords]));
  /* [numClusters][numCoords]  cluster's means */

  double **OBJ, **CLUSTOR;

  /* [numClusters][nCands][numCoords] the final medoids'
            candidates for every cluster*/
  OBJ = malloc(sizeof(double *) * numObjs);
  for (int qq = 0; qq < numObjs; qq++) {
    OBJ[qq] = malloc(sizeof(double) * numCoords);
  }

  for (int qq = 0; qq < numObjs; qq++) {
    for (int ww = 0; ww < numCoords; ww++) {
      OBJ[qq][ww] = objects[qq][ww];
      // printf("%f\n", OBJ[qq][ww]);
    }
  }

  CLUSTOR = malloc(sizeof(double *) * numClusters);
  for (int qq = 0; qq < numClusters; qq++) {
    CLUSTOR[qq] = malloc(sizeof(double) * numCoords);
  }

  for (int qq = 0; qq < numClusters; qq++) {
    for (int ww = 0; ww < numCoords; ww++) {
      CLUSTOR[qq][ww] = clusters[qq][ww];
      // printf("%f\n", OBJ[qq][ww]);
    }
  }
  (*delta) = 0;

  if (*flag == 0) {
    for (i = 0; i < numObjs; i++) {
      membership[i] = -1;

      if (i == 0) {
        *flag = 1;
      }
    }
  }

  for (i = 0; i < numObjs; i++) {
    /* find the array index of nearest cluster center */
    /* printf("\n numObjs : %d numClusters : %d , numCoords : %d , i : %d\n",
            numObjs, numClusters, numCoords, i);*/

    index = find_nearest_cluster(numClusters, numCoords, i, objects, clusters);
    // printf("%d,", index);
    /* if membership changes, increase delta by 1 */
    if (membership[i] != index) {
      *delta += 1.0;
    }

    /* assign the membership to object i */
    membership[i] = index;

    /*update new cluster center : sum of objects located within */
    clusterSize[index]++;

    /*calculate each feature sum for each medoid featu
    e*/
    for (j = 0; j < numCoords; j++) {
      cluster_distance_sum[index][j] += objects[i][j];
    }
  }
  for (int qq = 0; qq < numObjs; qq++) {
    free(OBJ[qq]);
  }
  free(OBJ);
  for (int qq = 0; qq < numClusters; qq++) {
    free(CLUSTOR[qq]);
  }
  free(CLUSTOR);
  // free(objects);
  return cluster_distance_sum;
}

/***********************************************************************************************************************/

double seq_camedoids(

    int numObjs, int numCoords, double **objects, int numClusters, int nCands,
    int *membership, double clusters[][numCoords], double *d, int delta) {

  int i, j, k, rep, index, loop = 0, pos;

  int *clusterSize;
  /* [numClusters]: no. objects assigned in each new cluster */

  double **clusters_means;
  /* [numClusters][numCoords]  cluster's means */

  double error = 0, *errors;

  int **medoid_cand;
  /* [numClusters][nCands]  medoid candidates for each
              cluster */
  double **cand_sum_dis;
  /* [numClusters][nCands]  accumulative distance of each
          medoid candidate to all other elements of its cluster*/
  double **cand_dis;
  /* [numClusters][nCands]  distance of each medoid candidate to
the mean of its cluster*/
  double *max_cand;
  /* [numClusters] candidate with the max distance of every
         other candidate in the same cluster*/
  int *cand;
  /* [numClusters] counter of candidates that have already been
 found for every cluster*/
  int *max_cand_pos;
  /* [numClusters] posistion of the candidate with the max
            distance of every other candidate in the same cluster*/
  double ***final_cand;
  /* [numClusters][nCands][numCoords] the final medoids'
            candidates for every cluster*/

  /*ALLOCATIONS*/
  cand = (int *)calloc(numClusters, sizeof(int));

  max_cand_pos = (int *)calloc(numClusters, sizeof(int));

  medoid_cand = (int **)malloc(numClusters * sizeof(int *));
  for (j = 0; j < numClusters; j++)
    medoid_cand[j] = (int *)calloc(nCands, sizeof(int));

  final_cand = (double ***)malloc(numClusters * sizeof(double **));
  for (j = 0; j < numClusters; j++) {
    final_cand[j] = (double **)malloc(nCands * sizeof(double *));
    for (i = 0; i < nCands; i++)
      final_cand[j][i] = (double *)malloc(numCoords * sizeof(double));
  }

  cand_dis = (double **)malloc(numClusters * sizeof(double *));
  for (j = 0; j < numClusters; j++)
    cand_dis[j] = (double *)calloc(nCands, sizeof(double));

  max_cand = (double *)calloc(numClusters, sizeof(double));

  cand_sum_dis = (double **)malloc(numClusters * sizeof(double *));
  for (i = 0; i < numClusters; i++) {
    cand_sum_dis[i] = (double *)calloc(nCands, sizeof(double));
  }

  errors = (double *)malloc(numClusters * sizeof(double));

  clusters_means = (double **)malloc(numClusters * sizeof(double *));
  for (i = 0; i < numClusters; i++)
    clusters_means[i] = (double *)calloc(numCoords, sizeof(double));

  clusterSize = (int *)calloc(numClusters, sizeof(int));
  error = 0;
  /*ALLOCATIONS END*/

  /*for (i = 0; i < numObjs; i++)
    membership[i] = -1;*/

  /*do {

    delta = 0;
    error = 0;*/

  /*PARALLEL PART ENDS HERE*/
  // printf("\nDELTA IN SEQ %d", delta);
  /* find nCands candidates for every cluster  */

  for (i = 0; i < numObjs; i++) {

    // if candidates of the cluster < nCands

    if (cand[membership[i]] < nCands) {
      // add the object to the cluster's candidate
      medoid_cand[membership[i]][cand[membership[i]]] = i;
      cand_dis[membership[i]][cand[membership[i]]] = d[i];
      if (d[i] > max_cand[membership[i]]) {

        max_cand[membership[i]] = d[i];
        max_cand_pos[membership[i]] = cand[membership[i]];
      }
      cand[membership[i]]++;
    } else if (d[i] < max_cand[membership[i]])
    // or distance of the current element smaller than the biggest //distance
    // of
    //   cluster's candidates so far
    {

      // add the object to the cluster's candidate
      medoid_cand[membership[i]][max_cand_pos[membership[i]]] = i;
      cand_dis[membership[i]][max_cand_pos[membership[i]]] = d[i];

      max_cand[membership[i]] = cand_dis[membership[i]][0];
      max_cand_pos[membership[i]] = 0;
      for (j = 1; j < nCands; j++)

        if (cand_dis[membership[i]][j] > max_cand[membership[i]]) {

          max_cand[membership[i]] = cand_dis[membership[i]][j];
          max_cand_pos[membership[i]] = j;
        }
    }
  }

  /* 3D */
  for (i = 0; i < numClusters; i++)
    for (j = 0; j < nCands; j++)
      for (k = 0; k < numCoords; k++)
        final_cand[i][j][k] = objects[medoid_cand[i][j]][k];

  /* calculate the accumulative distance of every candidate to the other
   * objects of its cluster*/
  for (i = 0; i < numObjs; i++) {
    for (j = 0; j < nCands; j++) {

      cand_sum_dis[membership[i]][j] +=
          euclid_dynamic(objects, numCoords, i, final_cand[membership[i]][j]);
    }
  }

  /* set the candidate with the smallest accumulative distance as the medoid
   * of its cluster */

  for (i = 0; i < numClusters; i++) {

    errors[i] = cand_sum_dis[i][0];
    pos = 0;

    for (j = 1; j < nCands; j++) {
      if (cand_sum_dis[i][j] < errors[i]) {
        errors[i] = cand_sum_dis[i][j];
        pos = j;
      }
    }

    for (j = 0; j < numCoords; j++)
      clusters[i][j] = final_cand[i][pos][j];
  }
  /* printf("\n");
   for (i = 0; i < numClusters; i++) {
     for (j = 0; j < numCoords; j++) {
       printf("%f,", clusters[i][j]);
     }
     printf("\n");
   }*/

  /* initialize for the next loop*/
  /*for (i = 0; i < numClusters; i++) {
    for (j = 0; j < numCoords; j++)
      clusters_means[i][j] = 0;
  }*/

  for (j = 0; j < numClusters; j++) {
    cand[j] = 0;
    max_cand[j] = 0;
    for (k = 0; k < nCands; k++)
      cand_dis[j][k] = 0;
  }

  for (i = 0; i < numClusters; i++) {
    clusterSize[i] = 0;
    for (j = 0; j < nCands; j++)
      cand_sum_dis[i][j] = 0;
  }

  // to ypologizei kathe fora ara if delta == 0 den ypologizei tipota.

  for (i = 0; i < numClusters; i++) {
    error += errors[i];
    printf("\n ERROR : %f\n", error);
  }

  for (j = 0; j < numClusters; j++) {
    for (i = 0; i < nCands; i++)
      free(final_cand[j][i]);
    free(final_cand[j]);
  }
  free(final_cand);
  for (i = 0; i < numClusters; i++)
    free(clusters_means[i]);
  free(clusters_means);
  free(clusterSize);
  for (i = 0; i < numClusters; i++)
    free(medoid_cand[i]);
  free(medoid_cand);
  for (i = 0; i < numClusters; i++)
    free(cand_sum_dis[i]);
  free(cand_sum_dis);
  free(errors);
  free(max_cand_pos);
  for (i = 0; i < numClusters; i++)
    free(cand_dis[i]);
  free(cand_dis);
  free(max_cand);
  free(cand);

  return error;
}
