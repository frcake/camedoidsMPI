#ifndef FILE_IO
#define FILE_IO

#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#define MAX_CHAR_PER_LINE 128

#define NONE 0
#define FIRST 1
#define LAST 2
#define BOTH 3

int *clusters_read(char *filename, int numclusters);
double **file_read(char *filename, int *numObjs, int *numCoords,
                   int lines_to_skip, int attr_to_skip);
void print_results(char *filename, int numCoords, double objects[][numCoords],
                   int *clusterid, double clusters[][numCoords], int numObjs,
                   int numclusters);

#endif
