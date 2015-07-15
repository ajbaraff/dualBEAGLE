#include <stdio.h>
#include <stdlib.h>
#define MAX 1048576

int main(int argc, char* argv[])
{
  FILE *infile, *outfile;
  char *filename, line[MAX], ***position, ***allele;
  int i, j, k, n, nsamps, seed, *nsites;

  /* check for too few arguments */
  if(argc < 2) {
    printf("Error: wrong number of arguments - \"ms2beagle <input file> <output file 1> ... <output file n>\".\n");
    return 1;
  }

  /* open and read the input file */
  filename = argv[1];
  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return 1;
  }

  /* command line and seed */
  fgets(line, MAX, infile);
  sscanf(line, "%*s %d %d", &n, &nsamps);
  fgets(line, MAX, infile);
  fgets(line, MAX, infile);

  /* check for too few arguments */
  if(argc != 2+nsamps) {
    printf("Error: wrong number of arguments - \"ms2beagle <input file> <output file 1> ... <output file n>\".\n");
    return 1;
  }

  nsites = malloc(nsamps*sizeof(int));
  position = malloc(nsamps*sizeof(char*));
  allele = malloc(nsamps*sizeof(int*));
  for(i = 0; i < nsamps; i++) {
    fgets(line, MAX, infile);
    fgets(line, MAX, infile);
    sscanf(line, "%*s %d", nsites+i);

    fgets(line, MAX, infile);
    sscanf(line, "%*s %[^\n]", line);
    position[i] = malloc(nsites[i]*sizeof(char*));
    for(k = 0; k < nsites[i]; k++) {
      position[i][k] = malloc(16*sizeof(char));
      sscanf(line, "%s %[^\n]", position[i][k], line);
    }

    allele[i] = malloc(n*sizeof(int*));
    for(j = 0; j < n; j++) {
      allele[i][j] = malloc(nsites[i]*sizeof(int));
      fgets(line, MAX, infile);
      for(k = 0; k < nsites[i]; k++) allele[i][j][k] = line[k];
    }
    fgets(line, MAX, infile);
  }
  
  fclose(infile);

  for(i = 0; i < nsamps; i++) {
    /* open and read the output files */
    filename = argv[2+i];
    outfile = fopen(filename, "w");
    if(outfile == NULL) {
      printf("Error: could not open file \"%s\".\n", filename);
      return 1;
    }

    fprintf(outfile, "I id");
    for(j = 0; j < n; j++) 
      fprintf(outfile, " s%d", j+1);
    fputc('\n', outfile);
    for(k = 0; k < nsites[i]; k++) {
      fprintf(outfile, "M m%s", position[i][k]);
      for(j = 0; j < n; j++) fprintf(outfile, " %c", allele[i][j][k]+1);
      fputc('\n', outfile);
    }
  
    fclose(outfile);
  }

  return 0;
}


