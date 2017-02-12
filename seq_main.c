#include "file_io.h"
#include "seq_camedoids.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define FLAG 0

static void usage(char *argv0) {
  char *help =
      "Usage: %s [switches] -i filename -n num_clusters -k num_cands\n"
      "       -f filename      : file containing data to be clustered\n"
      "       -c num_clusters  : number of clusters (K must > 1)\n"
      "       -i filename      : file containing the initial centroids \n"
      "       -l lines_to_skip : lines to be ignored from the beggining of the "
      "file\n"
      "       -a attr_to_skip  : attributes to be ignored (options 1:first, "
      "2:last, 3: first & last)\n"
      "       -k num_cands     : number of candidates\n"
      "       -h               : print this help information\n";
  fprintf(stderr, help, argv0);
  exit(-1);
}

/******************************************************************/
void p() { printf("\n"); }
void print_array_1d_double(double *A, int height, int width, int rank) {
  int count = 0;
  for (int i = 0; i < height * width; i++) {
    printf("R:%d - %f\t", rank, A[i]);
    count++;
    if (count == width) {
      printf("\n");
      count = 0;
    }
  }
  p();
}

void print_array_1d_int(int *A, int height, int width, int rank) {
  int count = 0;
  for (int i = 0; i < height * width; i++) {
    printf("R:%d - %d\t", rank, A[i]);
    count++;
    if (count == width) {
      printf("\n");
      count = 0;
    }
  }
  p();
}

void print_array(double *A, size_t height, size_t width, int process) {
  int count = 0;
  for (size_t i = 0; i < height; ++i) {
    for (size_t j = 0; j < width; ++j) {
      // printf("A[%zu][%zu] = %d\t", i, j, A[i * width + j]);
      printf("%f,", A[i * width + j]);
    }
  }
  p();
  p();
}

/*************************************************************************/

int main(int argc, char **argv) {

  // MY VAR initial
  int extraObjs, size, my_rank, sliceSize, sliceSizeExt, target, startSlice,
      endSlice;
  int *a;

  int i, j, opt;
  extern char *optarg;
  int numCoords, numObjs, numclusters = 2;
  int numCands = 100;
  double **objects; /* [numObjs][numCoords] data objects */

  int *clusters; /* [numClusters] cluster medoids */
  // double **medoids;
  int *clusterid, *clusteridFext; /* [numObjs] membership */
  double time2 = 0, time3;
  char *input_file, *clusters_file;
  unsigned long start, end;
  struct timeval tv1, tv2;
  int lines_to_skip = 0, attr_to_skip = 0, clustersfromfile = 0;
  // Command line arguments
  while ((opt = getopt(argc, argv, "f:c:i:l:a:t:k:h")) != EOF) {
    switch (opt) {
    case 'f':
      input_file = optarg;
      break;
    case 'c':
      numclusters = atoi(optarg);
      break;
    case 'i':
      clustersfromfile = 1;
      clusters_file = optarg;
      break;
    case 'l':
      lines_to_skip = atoi(optarg);
      break;
    case 'a':
      attr_to_skip = atoi(optarg);
      break;
    case 'k':
      numCands = atoi(optarg);
      break;
    case 'h':
    default:
      usage(argv[0]);
      break;
    }
  }

  if (numclusters <= 1) {
    printf("Too few clusters\n");
    return -1;
  }

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Status status;

  // dataset slicing!
  a = malloc(2 * sizeof(int));
  int *clusterSize, *clusteridAll, *REDUCEBUFSIZE, delta, *dlt, *displs,
      *recvcount, sumdelta = -2, k;
  int *flag = FLAG, loop;
  double *distance, *distanceFext, *distanceAll, error;

  /* read data points from file ------------------------------------------*/
  if (my_rank == 0) {

    objects = file_read(input_file, &numObjs, &numCoords, lines_to_skip,
                        attr_to_skip);

    // printf("\n Rank %d , Calculating Slice Size: %d", my_rank, sliceSize);
    if (objects == NULL)
      return -1;

    printf("::Objects loaded::\n");

    a[0] = numObjs;
    a[1] = numCoords;
  }

  MPI_Bcast(a, 2, MPI_INT, 0, MPI_COMM_WORLD);
  numObjs = a[0];
  numCoords = a[1];

  extraObjs = numObjs % size;
  sliceSize = numObjs / size;
  sliceSizeExt = (numObjs / size) + extraObjs;
  /*printf("\n-------Slicing statistics---------\n\n\n\nExtraObjects : %d \t "
         "sliceSize: %d \t sliceSizeExt: %d \n---------End of statistics "
         "-------\n\n\n\n",
         extraObjs, sliceSize, sliceSizeExt);*/
  double(*subArray)[numCoords] = malloc(sizeof(double[(sliceSize)][numCoords]));
  double(*medoids)[numCoords] =
      malloc(sizeof(double[(numclusters)][numCoords]));
  clusterid = (int *)malloc(sliceSize * sizeof(int));
  clusteridFext = (int *)malloc(sliceSizeExt * sizeof(int));
  clusteridAll = (int *)malloc(numObjs * sizeof(int));
  clusterSize = (int *)calloc(numclusters, sizeof(int));
  distance = calloc(sliceSize, sizeof(double));
  distanceFext = calloc(sliceSizeExt, sizeof(double));
  distanceAll = calloc(numObjs, sizeof(double));
  displs = calloc(size, sizeof(int));
  recvcount = calloc(size, sizeof(int));
  REDUCEBUFSIZE = (int *)calloc(numclusters, sizeof(int));
  /* medoids = (double **)malloc(numclusters * sizeof(double *));
   for (i = 0; i < numclusters; i++)
     medoids[i] = (double *)malloc(numCoords * sizeof(double));*/
  double(*cluster_distance_sum)[numCoords] =
      malloc(sizeof(double[(numclusters)][numCoords]));
  double(*clusters_means)[numCoords] =
      malloc(sizeof(double[(numclusters)][numCoords]));
  double(*REDUCEBUF)[numCoords] =
      malloc(sizeof(double[(numclusters)][numCoords]));

  if (my_rank == 0) {
    clusters = clusters_read(clusters_file, numclusters);
    if (clusters == NULL)
      return -1;
    for (i = 0; i < numclusters; i++) {
      for (j = 0; j < numCoords; j++) {
        medoids[i][j] = objects[clusters[i]][j];
      }
    }
  }
  double(*subArrayFext)[numCoords] =
      malloc(sizeof(double[sliceSizeExt][numCoords]));

  if (my_rank == 0) {

    // RANDOM MEDOIDS - FIRST TIME
    // MASTER PROCESS CODE - SLicing!
    for (target = 0; target < size; target++) {
      if (target == 0) {
        startSlice = target * sliceSize;
        endSlice = startSlice + sliceSizeExt - 1;
      } else if (target < size) {
        startSlice = target * sliceSize + extraObjs;
        endSlice = startSlice + sliceSize - 1;
      }

      for (i = startSlice; i <= endSlice; i++) {
        for (j = 0; j < numCoords; j++) {
          if (target == 0) {

            subArrayFext[i - startSlice][j] = objects[i][j];
          } else if (target > 0 && target < size) {
            subArray[i - startSlice][j] = objects[i][j];
          }
        }
      }

      if (target > 0 && target < size) {
        MPI_Send(subArray, sliceSize * numCoords, MPI_DOUBLE, target, 10,
                 MPI_COMM_WORLD);
      }
    }

  } else if (my_rank > 0 && my_rank < size) {

    MPI_Recv(subArray, sliceSize * numCoords, MPI_DOUBLE, target, 10,
             MPI_COMM_WORLD, &status);

    /* [numClusters][numCoords]  cluster's distance sum */
  }

  // COLLECTIVE COMMUNICATION!
  // for (int duo; duo < 2; duo++) {

  MPI_Barrier(MPI_COMM_WORLD);
  time_t t0, t1;
  if (my_rank == 0) {
    gettimeofday(&t0, 0);
  }
  for (int REC = 0; REC < 500; REC++) {
    MPI_Bcast(medoids, numclusters * numCoords, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {

      // RANDOM MEDOIDS - FIRST TIME
      // MASTER PROCESS CODE - SLicing!
      sumdelta = 0;

      cluster_distance_sum =
          seq_camedoids2(sliceSizeExt, numCoords, subArrayFext, numclusters,
                         clusteridFext, clusterSize, medoids, &dlt, &flag);
      delta = dlt;
      /*printf("\n********************* RANK %d DELTA: %d**************\n",
             my_rank, delta);*/

      /*MPI_Allreduce(cluster_distance_sum, REDUCEBUF, numclusters * numCoords,
                    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);*/

    } else if (my_rank > 0 && my_rank < size) {

      /* [numClusters][numCoords]  cluster's distance sum */

      cluster_distance_sum =
          seq_camedoids2(sliceSize, numCoords, subArray, numclusters, clusterid,
                         clusterSize, medoids, &dlt, &flag);
      // print_array(cluster_distance_sum, numclusters, numCoords, my_rank);
      delta = dlt;
      /*printf("\n********************* RANK %d DELTA: %d**************\n",
             my_rank, delta);*/
    }
    MPI_Allreduce(cluster_distance_sum, REDUCEBUF, numclusters * numCoords,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(clusterSize, REDUCEBUFSIZE, numclusters, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    /*if (my_rank == 0) {
      print_array(REDUCEBUF, numclusters, numCoords, my_rank);
    }*/
    /*if (my_rank == 0) {
      print_array_1d_int(REDUCEBUFSIZE, numclusters, 1, my_rank);
    }*/
    // print_array_1d_int(REDUCEBUFSIZE, numclusters, 1, my_rank);
    if (my_rank == 0) {
      //  print_array(cluster_distance_sum, numclusters, numCoords, my_rank);
      for (int qq = 0; qq < numclusters; qq++) {
        if (REDUCEBUFSIZE[qq] > 1) {
          for (k = 0; k < numCoords; k++)
            REDUCEBUF[qq][k] /= REDUCEBUFSIZE[qq];
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&delta, &sumdelta, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Bcast(&sumdelta, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(REDUCEBUF, numclusters * numCoords, MPI_DOUBLE, 0,
              MPI_COMM_WORLD);

    if (my_rank == 0) {
      distanceFext = mean_distance(sliceSizeExt, numCoords, subArrayFext,
                                   clusteridFext, REDUCEBUF);
      /*for (i = 0; i < numclusters; i++) {
        for (j = 0; j < numCoords; j++)
          REDUCEBUF[i][j] = 0;
      }*/
      for (int dcounter = 1; dcounter < size; dcounter++) {
        recvcount[dcounter] = sliceSize;
        displs[dcounter] = dcounter * sliceSize + extraObjs;
      }
      recvcount[0] = sliceSizeExt;
      displs[0] = 0;
    } else if (my_rank > 0 && my_rank < size) {
      /* for (i = 0; i < numclusters; i++) {
         for (j = 0; j < numCoords; j++)
           REDUCEBUF[i][j] = 0;
       }*/
      distance =
          mean_distance(sliceSize, numCoords, subArray, clusterid, REDUCEBUF);
    }

    MPI_Bcast(recvcount, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(displs, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == 0) {
      MPI_Gatherv(distanceFext, recvcount[my_rank], MPI_DOUBLE, distanceAll,
                  recvcount, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      MPI_Gatherv(clusteridFext, recvcount[my_rank], MPI_INT, clusteridAll,
                  recvcount, displs, MPI_INT, 0, MPI_COMM_WORLD);
    } else if (my_rank > 0 && my_rank < size) {
      MPI_Gatherv(distance, recvcount[my_rank], MPI_DOUBLE, distanceAll,
                  recvcount, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv(clusterid, recvcount[my_rank], MPI_INT, clusteridAll,
                  recvcount, displs, MPI_INT, 0, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // print_array_1d(REDUCEBUFSIZE, numclusters, 1, my_rank);

    // print_array(REDUCEBUF, numclusters, numCoords, my_rank);

    if (my_rank == 0) {
      // print_array_1d_int(clusteridAll, numObjs, 1, my_rank);
      seq_camedoids(numObjs, numCoords, objects, numclusters, numCands,
                    clusteridAll, medoids, distanceAll, sumdelta);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (i = 0; i < numclusters; i++) {
      clusterSize[i] = 0;
      REDUCEBUFSIZE[i] = 0;
      for (j = 0; j < numCoords; j++) {
        REDUCEBUF[i][j] = 0;
        cluster_distance_sum[i][j] = 0;
      }
    }

    // MPI_Bcast(&sumdelta, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /*MPI_Bcast(&sumdelta, 1, MPI_INT, 0, MPI_COMM_WORLD);*/
    if (sumdelta == 0) {
      break;
    }
  }
  /*if (my_rank == 0) {
    for (i = 0; i < numObjs; i++) {
      printf("%d \n", clusteridAll[i]);
    }
  }*/
  if (my_rank == 0) {

    gettimeofday(&t1, 0);
    long elapsed = (t1 - t0) * 1000000 + t1 - t0;
    printf("\n\n\n\n**********************************************************"
           "********* \n=REPORTING\n=time elapsed : %ld or approx. "
           "%ldseconds\n=histograms are saved in the code DIR. histogram 0 / "
           "output 0\n"
           "**********************************************************"
           "*********\n",
           elapsed, elapsed / 1000000);
  }
  if (my_rank == 0) {
    print_results(input_file, numCoords, objects, clusteridAll, medoids,
                  numObjs, numclusters);
  }
  /* output performance and clean
     ------------------------------------------*/
  if (my_rank == 0) {
    printf("::Clustering done::\n\n\n");
    printf("--- Clustering info ---\n");
    printf("File: %s\n", input_file);
    printf("Initial clusters: ");
    if (clustersfromfile == 1)
      printf("inserted from file\n");
    else
      printf("random\n");
    printf("Skipped lines: %d\n", lines_to_skip);
    printf("Ignored attribute: ");
    switch (attr_to_skip) {
    case NONE:
      printf("none\n");
      break;
    case FIRST:
      printf("first\n");
      break;
    case LAST:
      printf("last\n");
      break;
    case BOTH:
      printf("first and last\n");
      break;
    default:
      printf("\n");
    }
    printf("Objects: %d\n", numObjs);
    printf("Attributes: %d\n", numCoords);
    printf("Clusters: %d\n", numclusters);
    printf("Number of candidates: %d\n\n", numCands);

    printf("--- Results: ---\n");
    printf("Error: %f\n", error);
    printf("Time for clusters' initialization: %lf\n", time2);
    printf("Time for clustering: %lf\n", time3);
    printf("Total time: %lf\n", time2 + time3);

    free(objects);
    free(clusters);
    free(clusterid);
    free(clusteridFext);
    free(subArray);
    free(subArrayFext);
    free(clusteridAll);
    free(REDUCEBUF);
    free(REDUCEBUFSIZE);
    free(cluster_distance_sum);
    free(clusters_means);
    free(distance);
    free(distanceFext);
    free(distanceAll);
  }

  MPI_Finalize();

  return (0);
}
