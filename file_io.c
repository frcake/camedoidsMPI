#include "file_io.h"

int *clusters_read(char *filename, int numclusters) {
  // printf("\n FILENAME: %s", filename);
  FILE *infile;
  char *line;
  int lineLen, i;
  int *clusters = (int *)malloc(numclusters * sizeof(int));
  lineLen = MAX_CHAR_PER_LINE;

  if ((infile = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "Error: no such file (%s)\n", filename);
    return NULL;
  }

  for (i = 0; i < numclusters; i++) {
    /*1 2 3 4 5 6 7 8 9 10*/
    fscanf(infile, "%d", &clusters[i]);
  }

  fclose(infile);
  return clusters;
}

double **file_read(char *filename, int *numObjs, int *numCoords,
                   int lines_to_skip, int attr_to_skip) {
  double **objects;
  int i, j, len;
  ssize_t numBytesRead;
  int done = 0;
  FILE *infile;
  char *line, *ret;
  int lineLen;

  if ((infile = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "Error: no such file (%s)\n", filename);
    return NULL;
  }

  /* first find the number of objects */
  lineLen = MAX_CHAR_PER_LINE;
  line = (char *)malloc(lineLen);
  assert(line != NULL);

  (*numObjs) = 0;

  while (fgets(line, lineLen, infile) != NULL) {
    /* check each line to find the max line length */

    while (strlen(line) == lineLen - 1) {
      /* this line read is not complete */
      len = strlen(line);
      fseek(infile, -len, SEEK_CUR);

      /* increase lineLen */
      lineLen += MAX_CHAR_PER_LINE;
      line = (char *)realloc(line, lineLen);
      assert(line != NULL);

      ret = fgets(line, lineLen, infile);
      assert(ret != NULL);
    }
    // COUNT OBJECTS! BASED ON LINES
    if (strtok(line, " \t\n") != 0)
      (*numObjs)++;
  }

  (*numObjs) -= lines_to_skip;

  if ((*numObjs) <= 0) {
    fprintf(stderr, "Error: No objects found\n");
    return NULL;
  }

  rewind(infile);

  /*find the number of attributes*/
  (*numCoords) = 0;

  fgets(line, lineLen, infile);

  char *pch;
  pch = strtok(line, ",;");

  while (pch != NULL) {

    pch = strtok(NULL, ",;");
    (*numCoords)++;
  }

  if (attr_to_skip != NONE) {
    (*numCoords)--;
    if (attr_to_skip == BOTH)
      (*numCoords)--;
  }

  rewind(infile);

  /* allocate space for objects[][] and read all objects */
  len = (*numObjs) * (*numCoords);
  objects = (double **)malloc((*numObjs) * sizeof(double *));
  assert(objects != NULL);
  objects[0] = (double *)malloc(len * sizeof(double));
  assert(objects[0] != NULL);
  for (i = 1; i < (*numObjs); i++)
    objects[i] = objects[i - 1] + (*numCoords);

  /* read all objects */
  /*de mpenei*/
  for (i = 0; i < lines_to_skip; i++) {
    fgets(line, lineLen, infile);
  }

  i = 0;
  j = 0;

  while (fgets(line, lineLen, infile) != NULL) {
    pch = strtok(line, ",;");
    while (pch != NULL && j < (*numCoords)) {
      if (attr_to_skip % 2 == 1 && j == 0 && done == 0) {
        done = 1;
        pch = strtok(NULL, ",;");
        continue;
      }
      objects[i][j] = atof(pch);
      pch = strtok(NULL, ",;");
      j++;
    }
    i++;
    j = 0;
    done = 0;
  }

  assert(i == *numObjs);
  // printf("\n ASSERT %d %d", i, (*numObjs));

  fclose(infile);
  free(line);

  return objects;
}

void print_results(char *filename, int numCoords, double objects[][numCoords],
                   int *clusterid, double clusters[][numCoords], int numObjs,
                   int numclusters) {
  FILE *fptr;
  int i, j;
  char outFileName[1024];

  /* output: the coordinates of the clusters' medoids ----------------------*/
  sprintf(outFileName, "%s.medoids", filename);

  fptr = fopen(outFileName, "w");
  for (i = 0; i < numclusters; i++) {
    for (j = 0; j < numCoords; j++)
      fprintf(fptr, "%f ", clusters[i][j]);
    fprintf(fptr, "\n");
  }

  fclose(fptr);

  /* output: the cluster of each object ----------------------*/
  sprintf(outFileName, "%s.membership", filename);

  fptr = fopen(outFileName, "w");

  for (i = 0; i < numObjs; i++)
    fprintf(fptr, "%d\n", clusterid[i]);

  fclose(fptr);
}
