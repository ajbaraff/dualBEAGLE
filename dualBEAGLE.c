#include "dualBEAGLE.h"

struct tree *read_input(char *filename)
{
  FILE *infile;
  char line[MAX_LINE], *line_ptr;
  int i, j, k, k_max, n, curr_level, curr_index, curr_allele, max_allele, num_samp, num_level = START_SIZE, *num_node, *num_allele, *node_list;
  double curr_p_sum, new_p_sum, p_temp, p_max;
  struct tree_node *curr_node, ***node;
  struct tree *tree;

  /* open the input file */
  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return NULL;
  }

  /* reads opening line and calculates sample size */
  fgets(line, MAX_LINE, infile);
  num_samp = 0;
  line_ptr = line;
  sscanf(line_ptr, "%*s %*s %n", &n);
  line_ptr += n;
  while(sscanf(line_ptr, "%*s %n", &n) == 0) {
    num_samp++;
    line_ptr += n;
  }

  /* read input to find number and size of levels */
  num_allele = malloc(num_level*sizeof(int));
  curr_level = 0;
  max_allele = 2; // temporary fix - what if test has allele training doesn't
  while(fgets(line, MAX_LINE, infile) != NULL) {
    if(line[0] != '\n') {
      line_ptr = line;
      sscanf(line_ptr, "%*s %*s %n", &n);
      line_ptr += n;
      
      while(sscanf(line_ptr, "%d %n", &curr_allele, &n) == 1) {
	if(curr_allele > max_allele) max_allele = curr_allele;
 	line_ptr += n;
      }

      num_allele[curr_level] = max_allele;
      max_allele = 2;

      curr_level++;
      if(curr_level == num_level) {
        num_level *= 2;
        num_allele = realloc(num_allele, num_level*sizeof(int));
      }
    }
  }
  num_level = curr_level+1;
  num_allele = realloc(num_allele, num_level*sizeof(int));
  num_allele[num_level-1] = 0;

  rewind(infile);
  fgets(line, MAX_LINE, infile);

  tree = malloc(sizeof(struct tree));
  tree->num_level = num_level;
  tree->num_allele = num_allele;

  node = malloc(num_level*sizeof(struct tree_node**));
  node[0] = malloc(sizeof(struct tree_node*));
  node[0][0] = new_node(0, 0, num_allele[0]);
  num_node = calloc(num_level, sizeof(int));
  num_node[0] = 1;
  node_list = calloc(num_samp, sizeof(int));
  curr_level = 0;
  while(fgets(line, MAX_LINE, infile) != NULL) {
    if(line[0] != '\n') {
      line_ptr = line;
      sscanf(line_ptr, "%*s %*s %n", &n);
      line_ptr += n;
 
      if(curr_level+1 < num_level-1) {
	node[curr_level+1] = malloc(num_allele[curr_level]*num_node[curr_level]*sizeof(struct tree_node*));
	curr_index = 0;
	while(sscanf(line_ptr, "%d %n", &curr_allele, &n) == 1) {
	  curr_allele--;  // fix to allele key later
	  curr_node = node[curr_level][node_list[curr_index]];
	  if(curr_node->allele_miss[curr_allele]) {
	    node[curr_level+1][num_node[curr_level+1]] = new_node(curr_level+1, num_node[curr_level+1], num_allele[curr_level+1]);
	    add_edge(tree, curr_node, node[curr_level+1][num_node[curr_level+1]], curr_allele, 1);
	    node_list[curr_index] = num_node[curr_level+1];
	    num_node[curr_level+1]++;
	  }
	  else {
	    curr_node->count++;
	    curr_node->count_child[curr_allele]++;
	    node_list[curr_index] = curr_node->child[curr_allele]->id;
	  }
	  curr_index++;
	  line_ptr += n;
	}
      }

      else {
	node[curr_level+1] = malloc(sizeof(struct tree_node*));
	node[curr_level+1][0] = new_node(curr_level+1, 0, num_allele[curr_level+1]);
	num_node[curr_level+1] = 1;
        curr_index = 0;
        while(sscanf(line_ptr, "%d %n", &curr_allele, &n) == 1) {
	  curr_allele--;  // fix to allele key later
          curr_node = node[curr_level][node_list[curr_index]];
          if(curr_node->allele_miss[curr_allele]) 
	    add_edge(tree, curr_node, node[curr_level+1][0], curr_allele, 1);
	  else {
            curr_node->count++;
            curr_node->count_child[curr_allele]++;
          }
          curr_index++;
          line_ptr += n;
        }
      }

      node[curr_level+1] = realloc(node[curr_level+1], num_node[curr_level+1]*sizeof(struct tree_node*));

      curr_level++;
      if(curr_level == num_level) {
        num_level *= 2;
        num_allele = realloc(num_allele, num_level*sizeof(int));
      }
    }
  }

  tree->num_node = num_node;
  tree->node = node;

  /* done with input */
  fclose(infile);

  /* calculate edge and node weights */
  curr_p_sum = 1;
  for(i = 0; i < num_level; i++) {
    new_p_sum = 0;
    for(j = 0; j < num_node[i]; j++) {
      curr_node = node[i][j];
      curr_node->p_node /= curr_p_sum;
      for(k = 0; k < curr_node->num_allele; k++) {
        curr_node->p_child[k] = curr_node->p_child_test[k] = ((double)curr_node->count_child[k] + OFFSET) / ((double)curr_node->count + OFFSET*(double)curr_node->num_allele);
        p_temp = curr_node->p_node * curr_node->p_child[k];
        if(!curr_node->allele_miss[k]) {
          curr_node->child[k]->p_node += p_temp;
          new_p_sum += p_temp;
        }
      }
    }
    curr_p_sum = new_p_sum;
  }

  /* add missing parallel edges */
  for(i = 0; i < tree->num_level; i++)
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      p_max = 0;
      k_max = -1;
      for(k = 0; k < curr_node->num_allele; k++)
	if(curr_node->p_child[k] > p_max) {
	  p_max = curr_node->p_child[k];
	  k_max = k;
	}
      if(k_max >= 0)
	for(k = 0; k < curr_node->num_allele; k++)
	  if(curr_node->allele_miss[k])
	    add_edge(tree, curr_node, curr_node->child[k_max], k, 0);
    }

  return(tree);
}

struct tree *read_bgl(char *filename)
{
  FILE *infile;
  char line[MAX_LINE];
  int i, j, k, k_max, prev_level, curr_level, curr_index, next_index, max_index, curr_allele, max_allele, num_level = START_SIZE, curr_count, *num_node, *num_allele;
  double curr_p_sum, new_p_sum, p_temp, p_max;
  struct tree_node *curr_node, ***node;
  struct tree *tree;

  // open the input file 
  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return NULL;
  }

  // skip first two lines 
  fgets(line, MAX_LINE, infile);
  fgets(line, MAX_LINE, infile);

  // read input to find number and size of levels 
  num_node = malloc(num_level*sizeof(int));
  num_allele = malloc(num_level*sizeof(int));
  prev_level = 0;
  max_index = 0;
  max_allele = 2; // temporary fix - what if test has allele training doesn't
  while(fgets(line, MAX_LINE, infile) != NULL) {
    if(line[0] != '\n') {
      sscanf(line, "%d %*s %d %*d %d", &curr_level, &curr_index, &curr_allele);
      
      if(curr_level == prev_level) {
        if(curr_index > max_index) max_index = curr_index;
	if(curr_allele > max_allele) max_allele = curr_allele;
      }

      else {
	if(curr_level == num_level) {
	  num_level *= 2;
          num_node = realloc(num_node, num_level*sizeof(int));
	  num_allele = realloc(num_allele, num_level*sizeof(int));
	}

	num_node[prev_level] = max_index + 1;
	num_allele[prev_level] = max_allele;
	max_index = 0;
	max_allele = 2;

	prev_level = curr_level;
      }
    }
  }
  num_level = curr_level+2;
  num_node = realloc(num_node, num_level*sizeof(int));
  num_allele = realloc(num_allele, num_level*sizeof(int));
  num_node[prev_level] = max_index + 1;
  num_allele[prev_level] = max_allele;
  num_node[num_level-1] = 1;
  num_allele[num_level-1] = 0;

  rewind(infile);
  fgets(line, MAX_LINE, infile);
  fgets(line, MAX_LINE, infile);

  tree = malloc(sizeof(struct tree));
  tree->num_level = num_level;
  tree->num_node = num_node;
  tree->num_allele = num_allele;

  // make nodes 
  node = malloc(num_level*sizeof(struct tree_node**));
  for(i = 0; i < num_level; i++) {
    node[i] = malloc(num_node[i]*sizeof(struct tree_node*));
    for(j = 0; j < num_node[i]; j++)
      node[i][j] = new_node(i, j, num_allele[i]);
  }

  tree->node = node;

  // read input to get edges 
  while(fgets(line, MAX_LINE, infile) != NULL) {
    if(line[0] != '\n') {
      sscanf(line, "%d %*s %d %d %d %d", &curr_level, &curr_index, &next_index, &curr_allele, &curr_count);
      curr_allele--;  // fix to allele key later
      add_edge(tree, node[curr_level][curr_index], node[curr_level+1][next_index], curr_allele, curr_count);
    }
  }

  // done with input 
  fclose(infile);

  // calculate edge and node weights
  curr_p_sum = 1;
  for(i = 0; i < num_level; i++) {
    new_p_sum = 0;
    for(j = 0; j < num_node[i]; j++) {
      curr_node = node[i][j];
      curr_node->p_node /= curr_p_sum;
      for(k = 0; k < curr_node->num_allele; k++) {
        curr_node->p_child[k] = curr_node->p_child_test[k] = ((double)curr_node->count_child[k] + OFFSET) / ((double)curr_node->count + OFFSET*(double)curr_node->num_allele);
        p_temp = curr_node->p_node * curr_node->p_child[k];
        if(!curr_node->allele_miss[k]) {
          curr_node->child[k]->p_node += p_temp;
          new_p_sum += p_temp;
        }
      }
    }
    curr_p_sum = new_p_sum;
  }

  /* add missing parallel edges */
  /*for(i = 0; i < tree->num_level; i++)
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      p_max = 0;
      k_max = -1;
      for(k = 0; k < curr_node->num_allele; k++)
        if(curr_node->p_child[k] > p_max) {
          p_max = curr_node->p_child[k];
          k_max = k;
        }
      if(k_max >= 0)
        for(k = 0; k < curr_node->num_allele; k++)
          if(curr_node->allele_miss[k])
            add_edge(tree, curr_node, curr_node->child[k_max], k, 0);
	    }*/

  return(tree);
}

char *reverse_input(char *filename) 
{
  FILE *infile, *outfile;
  char line[MAX_LINE], *temp;
  int i, j, num_line = 0;

  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return NULL;
  }

  /* open the output file */
  temp = malloc((strlen(filename)+5)*sizeof(char));
  strcpy(temp, filename);
  strcat(temp, ".rev");
  outfile = fopen(temp, "w");
  if(outfile == NULL) {
    printf("Error: could not open file \"%s\".\n", temp);
    return NULL;
  }

  /* copy opening line */
  fgets(line, MAX_LINE, infile);
  fprintf(outfile, "%s", line);

  while(fgets(line, MAX_LINE, infile) != NULL) num_line++;
  rewind(infile);

  for(i = 0; i < num_line; i++) {
    fgets(line, MAX_LINE, infile);
    for(j = 0; j < num_line-i; j++)
      fgets(line, MAX_LINE, infile);
    fprintf(outfile, "%s", line);
    rewind(infile);
  }

  fclose(infile);
  fclose(outfile);

  return(temp);
}

void split_input(char *filename, int num_split, int seed)
{
  FILE *infile, *outfile1, *outfile2;
  char line[MAX_LINE], *line_ptr, *temp, init, entry[MAX_ENTRY];
  int i, n, num_samp, curr_allele, *selected;
  double rand_u;

  /* open the input file */
  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return;
  }

  /* reads opening line and calculates sample size */
  fgets(line, MAX_LINE, infile);
  num_samp = 0;
  line_ptr = line;
  sscanf(line_ptr, "%*s %*s %n", &n);
  line_ptr += n;
  while(sscanf(line_ptr, "%*s %n", &n) == 0) {
    num_samp++;
    line_ptr += n;
  }

  if(num_split > num_samp) {
    printf("Error: Split size larger than sample size.\n");
    return;
  }
  
  i = n = 0;
  selected = calloc(num_samp, sizeof(int));
  srand(seed);
  while (n < num_split) {
    rand_u = (double)rand() / ((double)RAND_MAX + 1);
    if((num_samp-i)*rand_u >= num_split-n)
      i++;
    else {
      selected[i] = 1;
      i++; 
      n++;
    }
  }

  temp = malloc((strlen(filename)+7)*sizeof(char));
  strcpy(temp, filename);
  strcat(temp, ".train");
  outfile1 = fopen(temp, "w");
  if(outfile1 == NULL) {
    printf("Error: could not open file \"%s\".\n", temp);
    return;
  }
  free(temp);

  temp = malloc((strlen(filename)+6)*sizeof(char));
  strcpy(temp, filename);
  strcat(temp, ".test");
  outfile2 = fopen(temp, "w");
  if(outfile2 == NULL) {
    printf("Error: could not open file \"%s\".\n", temp);
    return;
  }
  free(temp);

  line_ptr = line;
  sscanf(line_ptr, "%c %s %n", &init, entry, &n);
  fprintf(outfile1, "%c %s ", init, entry);
  fprintf(outfile2, "%c %s ", init, entry);
  line_ptr += n;
  i = 0;
  while(sscanf(line_ptr, "%s %n", entry, &n) == 1) {
    if(selected[i]) fprintf(outfile1, "%s ", entry);
    else fprintf(outfile2, "%s ", entry);
    line_ptr += n;
    i++;
  }
  fprintf(outfile1, "\n");
  fprintf(outfile2, "\n");

  while(fgets(line, MAX_LINE, infile) != NULL) {
    if(line[0] != '\n') {
      line_ptr = line;
      sscanf(line_ptr, " %c %s %n", &init, entry, &n);
      fprintf(outfile1, " %c %s ", init, entry);
      fprintf(outfile2, " %c %s ", init, entry);
      line_ptr += n;
      i = 0;
      while(sscanf(line_ptr, "%d %n", &curr_allele, &n) == 1) {
        if(selected[i]) fprintf(outfile1, "%d ", curr_allele);
	else fprintf(outfile2, "%d ", curr_allele);
        line_ptr += n;
	i++;
      }
      fprintf(outfile1, "\n");
      fprintf(outfile2, "\n");
    }
  }

  fclose(infile);
  fclose(outfile1);
  fclose(outfile2);

  free(selected);
}

void boot_input(char *filename, int seed)
{
  FILE *infile, *outfile;
  char line[MAX_LINE], *line_ptr, *temp, init, entry[MAX_ENTRY];
  int i, j, n, num_samp, curr_allele, *selected;
  double rand_u;

  /* open the input file */
  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return;
  }

  /* reads opening line and calculates sample size */
  fgets(line, MAX_LINE, infile);
  num_samp = 0;
  line_ptr = line;
  sscanf(line_ptr, "%*s %*s %n", &n);
  line_ptr += n;
  while(sscanf(line_ptr, "%*s %n", &n) == 0) {
    num_samp++;
    line_ptr += n;
  }

  selected = calloc(num_samp, sizeof(int));
  srand(seed);
  for(i = 0; i < num_samp; i++) {
    rand_u = (double)rand() / ((double)RAND_MAX + 1);
    selected[(int)(num_samp * rand_u)]++;
  }

  temp = malloc((strlen(filename)+6)*sizeof(char));
  strcpy(temp, filename);
  strcat(temp, ".boot");
  outfile = fopen(temp, "w");
  if(outfile == NULL) {
    printf("Error: could not open file \"%s\".\n", temp);
    return;
  }
  free(temp);

  line_ptr = line;
  sscanf(line_ptr, "%c %s %n", &init, entry, &n);
  fprintf(outfile, "%c %s ", init, entry);
  line_ptr += n;
  i = 0;
  while(sscanf(line_ptr, "%s %n", entry, &n) == 1) {
    for(j = 0; j < selected[i]; j++)
      fprintf(outfile, "%s ", entry);
    line_ptr += n;
    i++;
  }
  fprintf(outfile, "\n");

  while(fgets(line, MAX_LINE, infile) != NULL) {
    if(line[0] != '\n') {
      line_ptr = line;
      sscanf(line_ptr, " %c %s %n", &init, entry, &n);
      fprintf(outfile, " %c %s ", init, entry);
      line_ptr += n;
      i = 0;
      while(sscanf(line_ptr, "%d %n", &curr_allele, &n) == 1) {
	for(j = 0; j < selected[i]; j++)
	  fprintf(outfile, "%d ", curr_allele);
        line_ptr += n;
        i++;
      }
      fprintf(outfile, "\n");
    }
  }

  fclose(infile);
  fclose(outfile);

  free(selected);
}

void prep_boot(char *filename, struct tree *tree, int num_level, int start, int length, int reverse) {
  int i, j, k, n;
  double p_temp, curr_p_sum, new_p_sum;
  struct tree_node *curr_node;

  read_test(filename, tree, num_level, start, length, reverse);

  tree->node[0][0]->p_node = 1;
  for(i = 1; i < tree->num_level; i++) 
    for(j = 0; j < tree->num_node[i]; j++) 
      tree->node[i][j]->p_node = 0;
  
  /* calculate edge and node weights */
  curr_p_sum = 1;
  for(i = 0; i < tree->num_level; i++) {
    new_p_sum = 0;
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      curr_node->p_node /= curr_p_sum;
      for(k = 0; k < curr_node->num_allele; k++) {
	if(curr_node->count_test > 0)
	  curr_node->p_child[k] = curr_node->p_child_test[k] = ((double)curr_node->count_test_child[k] + OFFSET) / ((double)curr_node->count_test + OFFSET*(double)curr_node->num_allele);
	else
	  curr_node->p_child[k] = 0;
        p_temp = curr_node->p_node * curr_node->p_child[k];
        if(!curr_node->allele_miss[k]) {
          curr_node->child[k]->p_node += p_temp;
          new_p_sum += p_temp;
        }
      }
    }
    curr_p_sum = new_p_sum;
  }
}

struct tree_node *new_node(int level, int id, int num_allele)
{
  int i;
  struct tree_node *new = malloc(sizeof(struct tree_node));

  new->level = level;
  new->id = id;
  new->count = new->num_back = new->num_haplo = new->marked = 0;
  new->count_test = 0;
  new->num_allele = num_allele;
  new->count_child = calloc(num_allele, sizeof(int));
  new->count_test_child = calloc(num_allele, sizeof(double));
  new->allele_miss = malloc(num_allele*sizeof(int));
  for(i = 0; i < num_allele; i++) new->allele_miss[i] = 1;
  new->p_node = (level==0);
  new->back_weight = NULL;
  new->p_child = malloc(num_allele*sizeof(double));
  new->p_child_test = malloc(num_allele*sizeof(double));
  new->haplo = NULL;
  new->child = malloc(num_allele*sizeof(struct tree_node*));
  for(i = 0; i < num_allele; i++) new->child[i] = NULL;

  return new;
}

void add_edge(struct tree *tree, struct tree_node *parent, struct tree_node *child, int allele, int count)
{
  parent->count += count;
  parent->count_child[allele] += count;

  if(count > 0)
    parent->allele_miss[allele] = 0;
  parent->child[allele] = child;
}

void free_tree(struct tree *tree)
{
  int i, j;

  if(tree == NULL) return;
  
  for(i = 0; i < tree->num_level; i++) {
    for(j = 0; j < tree->num_node[i]; j++) free_node(tree->node[i][j]);
    free(tree->node[i]);
  }

  free(tree->node);
  free(tree->num_node);
  free(tree->num_allele);
    
  free(tree);
}

void free_node(struct tree_node *node)
{
  int i;

  if(node == NULL) return;

  free(node->count_child);
  free(node->count_test_child);
  free(node->allele_miss);
  
  free(node->p_child);
  free(node->p_child_test);
  if(node->num_back > 0) free(node->back_weight);

  /* for(i = 0; i < node->num_haplo; i++) {
    free(node->haplo[i]->allele);
    free(node->haplo[i]);
  }
  free(node->haplo); */
  
  free(node->child);

  free(node);
}

struct tree *copy_tree(struct tree *tree)
{
  int i, j, k;
  struct tree_node *curr_node;
  struct tree *new;

  new = malloc(sizeof(struct tree));
  new->num_level = tree->num_level;
  new->num_node = malloc(tree->num_level*sizeof(int));
  memcpy(new->num_node, tree->num_node, tree->num_level*sizeof(int));
  new->num_allele = malloc(tree->num_level*sizeof(int));
  memcpy(new->num_allele, tree->num_allele, tree->num_level*sizeof(int));
  new->loss = tree->loss;
  new->loss_test = tree->loss_test;
  new->start = tree->start;
  new->length = tree->length;

  new->node = malloc(tree->num_level*sizeof(struct tree_node**));
  for(i = 0; i < tree->num_level; i++) {
    new->node[i] = malloc(tree->num_node[i]*sizeof(struct tree_node*));
    for(j = 0; j < tree->num_node[i]; j++)
      new->node[i][j] = copy_node(tree->node[i][j]);
  }

  for(i = 0; i < tree->num_level; i++)
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = new->node[i][j];
      for(k = 0; k < curr_node->num_allele; k++)
        if(!curr_node->allele_miss[k])
	  curr_node->child[k] = new->node[i+1][tree->node[i][j]->child[k]->id];
    }

  return(new);
}

struct tree_node *copy_node(struct tree_node *node)
{
  int i;
  struct tree_node *new = malloc(sizeof(struct tree_node));

  new->level = node->level;
  new->id = node->id;
  new->count = node->count;
  new->count_test = node->count_test;
  new->num_allele = node->num_allele;
  new->num_back = node->num_back;
  new->num_haplo = node->num_haplo;
  new->marked = node->marked;

  new->count_child = malloc(node->num_allele*sizeof(int));
  memcpy(new->count_child, node->count_child, node->num_allele*sizeof(int));
  new->count_test_child = malloc(node->num_allele*sizeof(double));
  memcpy(new->count_test_child, node->count_test_child, node->num_allele*sizeof(double));
  new->allele_miss = malloc(node->num_allele*sizeof(int));
  memcpy(new->allele_miss, node->allele_miss, node->num_allele*sizeof(int));

  new->p_node = node->p_node;

  new->p_child = malloc(node->num_allele*sizeof(double));
  memcpy(new->p_child, node->p_child, node->num_allele*sizeof(double));
  new->p_child_test = malloc(node->num_allele*sizeof(double));
  memcpy(new->p_child_test, node->p_child_test, node->num_allele*sizeof(double));
  new->back_weight = malloc(node->num_back*sizeof(double)); 
  memcpy(new->back_weight, node->back_weight, node->num_back*sizeof(double));

  new->haplo = NULL;
  /* new->haplo = malloc(node->num_haplo*sizeof(struct haplotype*));
     for(i = 0; i < node->num_haplo; i++) {
     new->haplo[i] = malloc(sizeof(struct haplotype));
     new->haplo[i]->length = node->haplo[i]->length;
     new->haplo[i]->allele = malloc(node->haplo[i]->length*sizeof(int));
     memcpy(new->haplo[i]->allele, node->haplo[i]->allele, node->haplo[i]->length*sizeof(int));
     new->haplo[i]->freq = node->haplo[i]->freq;
     } */

  new->child = malloc(node->num_allele*sizeof(struct tree_node*));
  for(i = 0; i < node->num_allele; i++) 
    new->child[i] = NULL;
  
  return new;
}

struct tree *cut_tree(struct tree *tree, int start, int length, int reverse)
{
  int i, j, k;
  struct tree_node *curr_node;
  struct tree *new;

  if(reverse) start = tree->num_level - start - length + 1;

  new = malloc(sizeof(struct tree));
  new->num_level = length + 1;
  new->num_node = malloc(new->num_level*sizeof(int));
  memcpy(new->num_node, tree->num_node + start - 1, new->num_level*sizeof(int));
  new->num_allele = malloc(new->num_level*sizeof(int));
  memcpy(new->num_allele, tree->num_allele + start - 1, new->num_level*sizeof(int));
  new->num_allele[length] = 0;

  new->node = malloc(new->num_level*sizeof(struct tree_node**));
  for(i = 0; i < new->num_level; i++) {
    new->node[i] = malloc(new->num_node[i]*sizeof(struct tree_node*));
    for(j = 0; j < new->num_node[i]; j++) {
      new->node[i][j] = copy_node(tree->node[i+start-1][j]);
      new->node[i][j]->level -= start - 1;
    }
  }

  for(i = 0; i < new->num_node[length]; i++) {
    curr_node = new->node[length][i];
    for(j = 0; j < curr_node->num_allele; j++)
      curr_node->child[j] = NULL;
    curr_node->num_allele = 0;
  }

  for(i = 0; i < new->num_level; i++) 
    for(j = 0; j < new->num_node[i]; j++) {
      curr_node = new->node[i][j];
      for(k = 0; k < curr_node->num_allele; k++) 
	curr_node->child[k] = new->node[i+1][tree->node[i+start-1][j]->child[k]->id];
    }

  while(new->num_node[0] > 1) 
    merge_node(new, new->node[0][0], new->node[0][1], 0, 1);
  while(new->num_node[length] > 1) 
    merge_node(new, new->node[length][0], new->node[length][1], 0, 1);

  return(new);
}

void calc_test(struct tree *train, struct tree *test)
{
  int i, j, k, *node_list, *next_list;
  struct tree_node *curr_node;

  if(train->node[0][0]->count_test > 0)
    for(i = 0; i < train->num_level; i++)
      for(j = 0; j < train->num_node[i]; j++) {
	curr_node = train->node[i][j];
	curr_node->count_test = 0;
	for(k = 0; k < curr_node->num_allele; k++)
	  curr_node->count_test_child[k] = 0;
      }
  
  node_list = malloc(sizeof(int));
  node_list[0] = 0;

  for(i = 0; i < test->num_level-1; i++) {
    next_list = malloc(test->num_node[i+1]*sizeof(int));
    for(j = 0; j < test->num_node[i]; j++) {
      curr_node = test->node[i][j];
      train->node[i][node_list[j]]->count_test += curr_node->count;
      for(k = 0; k < curr_node->num_allele; k++)
	if(!curr_node->allele_miss[k]) {
	  train->node[i][node_list[j]]->count_test_child[k] += curr_node->count_child[k];
	  next_list[curr_node->child[k]->id] = train->node[i][node_list[j]]->child[k]->id;
	}
    }
    free(node_list);
    node_list = next_list;
  }
  train->node[test->num_level-1][0]->count_test += test->node[test->num_level-1][0]->count;

  free(node_list);
}

void read_test(char *filename, struct tree *tree, int num_level, int start, int length, int reverse)
{
  FILE *infile;
  char line[MAX_LINE], *line_ptr;
  int i, j, k, n, curr_level, curr_allele, curr_samp, num_samp;
  double *temp, **node_list;
  struct tree_node *curr_node;

  if(reverse) start = num_level - start - length + 1;

  if(tree->node[0][0]->count_test > 0)
    for(i = 0; i < tree->num_level; i++)
      for(j = 0; j < tree->num_node[i]; j++) {
        curr_node = tree->node[i][j];
        curr_node->count_test = 0;
        for(k = 0; k < curr_node->num_allele; k++)
          curr_node->count_test_child[k] = 0;
      }

  // open the input file 
  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return;
  }

  // reads opening line and calculates sample size 
  fgets(line, MAX_LINE, infile);
  num_samp = 0;
  line_ptr = line;
  sscanf(line_ptr, "%*s %*s %n", &n);
  line_ptr += n;
  while(sscanf(line_ptr, "%*s %n", &n) == 0) {
    num_samp++;
    line_ptr += n;
  }
  tree->node[0][0]->count_test = num_samp;

  // read input to calculate BIC 
  curr_level = 1 - start;
  node_list = malloc(num_samp*sizeof(double*));
  for(curr_samp = 0; curr_samp < num_samp; curr_samp++) {
    node_list[curr_samp] = malloc(sizeof(double));
    node_list[curr_samp][0] = 1;
  }
  while(fgets(line, MAX_LINE, infile) != NULL && curr_level < length) {
    if(line[0] != '\n') {
      if(curr_level >= 0) {
	line_ptr = line;
	sscanf(line_ptr, "%*s %*s %n", &n);
	line_ptr += n;
	
	curr_samp = 0;
	while(sscanf(line_ptr, "%d %n", &curr_allele, &n) == 1) {
	  curr_allele--;  // fix to allele key later
	  temp = calloc(tree->num_node[curr_level+1], sizeof(double));
	  for(i = 0; i < tree->num_node[curr_level]; i++) {
	    curr_node = tree->node[curr_level][i];
	    curr_node->count_test_child[curr_allele] += node_list[curr_samp][i];
	    /*if(curr_node->allele_miss[curr_allele])
	      for(j = 0; j < tree->num_node[curr_level+1]; j++)
		temp[j] += node_list[curr_samp][i] * tree->node[curr_level+1][j]->p_node;
	    else
	    temp[curr_node->child[curr_allele]->id] += node_list[curr_samp][i]; */
	    temp[curr_node->child[curr_allele]->id] += node_list[curr_samp][i];
	  }
	  for(i = 0; i < tree->num_node[curr_level+1]; i++)
	    tree->node[curr_level+1][i]->count_test += temp[i];
	  free(node_list[curr_samp]);
	  node_list[curr_samp] = temp;
	  curr_samp++;
	  line_ptr += n;
	}
      }
      curr_level++;
    }
  }

  fclose(infile);
}

/*
void read_test(char *filename, struct tree *tree, int num_level, int start, int length, int reverse)
{
  FILE *infile;
  char line[MAX_LINE], *line_ptr;
  int i, j, k, n, curr_level, curr_allele, curr_samp, num_samp;
  double *temp, **node_list;
  struct tree_node *curr_node;

  if(reverse) start = num_level - start - length + 1;

  if(tree->node[0][0]->count_test > 0)
    for(i = 0; i < tree->num_level; i++)
      for(j = 0; j < tree->num_node[i]; j++) {
        curr_node = tree->node[i][j];
        curr_node->count_test = 0;
        for(k = 0; k < curr_node->num_allele; k++)
          curr_node->count_test_child[k] = 0;
      }

  // open the input file 
  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return;
  }

  // reads opening line and calculates sample size 
  fgets(line, MAX_LINE, infile);
  num_samp = 0;
  line_ptr = line;
  sscanf(line_ptr, "%*s %*s %n", &n);
  line_ptr += n;
  while(sscanf(line_ptr, "%*s %n", &n) == 0) {
    num_samp++;
    line_ptr += n;
  }
  tree->node[0][0]->count_test = num_samp;

  // read input to calculate BIC 
  curr_level = 1 - start;
  node_list = malloc(num_samp*sizeof(double*));
  for(curr_samp = 0; curr_samp < num_samp; curr_samp++) {
    node_list[curr_samp] = malloc(sizeof(double));
    node_list[curr_samp][0] = 1;
  }
  while(fgets(line, MAX_LINE, infile) != NULL && curr_level < length) {
    if(line[0] != '\n') {
      if(curr_level >= 0) {
	printf("curr_level = %d\n", curr_level);
	line_ptr = line;
	sscanf(line_ptr, "%*s %*s %n", &n);
	line_ptr += n;
	
	curr_samp = 0;
	while(sscanf(line_ptr, "%d %n", &curr_allele, &n) == 1) {
	  curr_allele--;  // fix to allele key later
	  temp = calloc(tree->num_node[curr_level+1], sizeof(double));
	  for(i = 0; i < tree->num_node[curr_level]; i++) {
	    curr_node = tree->node[curr_level][i];
	    curr_node->count_test_child[curr_allele] += node_list[curr_samp][i];
	    /* if(curr_node->allele_miss[curr_allele])
	      for(j = 0; j < tree->num_node[curr_level+1]; j++)
		temp[j] += node_list[curr_samp][i] * tree->node[curr_level+1][j]->p_node;
	    else
	    temp[curr_node->child[curr_allele]->id] += node_list[curr_samp][i]; 
	    temp[curr_node->child[curr_allele]->id] += node_list[curr_samp][i];
	  }
	  for(i = 0; i < tree->num_node[curr_level+1]; i++)
	    tree->node[curr_level+1][i]->count_test += temp[i];
	  free(node_list[curr_samp]);
	  node_list[curr_samp] = temp;
	  curr_samp++;
	  line_ptr += n;
	}
      }	
      curr_level++;
    }
  }

  fclose(infile);
} */

void print_tree(char *filename, struct tree *tree, int reverse)
{
  int i, j, k;
  struct tree_node *curr_node;
  FILE *outfile;

  outfile = fopen(filename, "w");
  if(outfile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return;
  }

  for(i = 0; i < tree->num_level; i++) {
    //for(i = tree->start - 1; i < tree->start + tree->length; i++) {
    if(reverse) fprintf(outfile, "Level %d: %d\n", tree->num_level-i-1, tree->num_node[i]);
    else fprintf(outfile, "Level %d: %d\n", i, tree->num_node[i]);
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      fprintf(outfile, "\tNode ID = %d, Node Weight = %8.8lf, Node Count = %d, Test Count = %8.8lf\n", curr_node->id+1, curr_node->p_node, curr_node->count, curr_node->count_test);
      for(k = 0; k < curr_node->num_allele; k++) {
	if(curr_node->allele_miss[k])
	  fprintf(outfile, "\t\tChild Node = %d, Child Allele = %d, Child Weight = %8.8lf, Child Count = %d. Test Count = %8.8lf\n", 0, k+1, curr_node->p_child[k], curr_node->count_child[k], curr_node->count_test_child[k]);
	else
	  fprintf(outfile, "\t\tChild Node = %d, Child Allele = %d, Child Weight = %8.8lf, Child Count = %d. Test Count = %8.8lf\n", curr_node->child[k]->id+1, k+1, curr_node->p_child[k], curr_node->count_child[k], curr_node->count_test_child[k]);
      }
    }
  }

  fclose(outfile);
}

double loglik(struct tree *tree)
{
  int i, j, k;
  double loglik = 0;
  struct tree_node *curr_node;
  
  for(i = 0; i < tree->num_level; i++) {
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      for(k = 0; k < curr_node->num_allele; k++) {
	if(curr_node->p_child[k] > 0)
	  loglik += curr_node->count_child[k] * log(curr_node->p_child[k]);
	else
	  loglik += curr_node->count_child[k] * log(1.0/(curr_node->count+1.0));
      }
    }
  }
  
  return loglik;
}


double loss(struct tree *tree, double lambda)
{
  int i, j, k;
  double loss = 0;
  struct tree_node *curr_node;

  /* for(i = 0; i < tree->num_node[tree->start-1]; i++)
     loss -= tree->node[tree->start-1][i]->count * log(tree->node[tree->start-1][i]->p_node); */

  for(i = tree->start - 1; i < tree->start + tree->length - 1; i++) {
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      for(k = 0; k < curr_node->num_allele; k++) {
        if(curr_node->p_child[k] > 0)
          loss -= curr_node->count_child[k] * log(curr_node->p_child[k]);
        else
          loss -= curr_node->count_child[k] * log(1.0/(curr_node->count+1.0));
      }
    }
    loss += lambda * tree->num_node[i] * (tree->num_allele[i]-1);
    //loss += lambda * (tree->num_node[i] * (tree->num_node[i]-1) / 2) * (tree->num_allele[i]-1);
  }
  
  tree->loss = loss;
  return loss;
}

double aic(struct tree *tree)
{
  int i, j, k;
  double loss = 0;
  struct tree_node *curr_node;

  /* for(i = 0; i < tree->num_node[tree->start-1]; i++)
     loss -= tree->node[tree->start-1][i]->count * log(tree->node[tree->start-1][i]->p_node); */

  for(i = tree->start - 1; i < tree->start + tree->length - 1; i++) {
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      for(k = 0; k < curr_node->num_allele; k++) {
        if(curr_node->p_child[k] > 0)
          loss -= curr_node->count_child[k] * log(curr_node->p_child[k]);
        else
          loss -= curr_node->count_child[k] * log(1.0/(curr_node->count+1.0));
      }
    }
    loss += tree->num_node[i] * (tree->num_allele[i]-1);
  }

  return 2*loss;
}

double bic(struct tree *tree)
{
  int i, j, k;
  double loss = 0;
  struct tree_node *curr_node;

  /* for(i = 0; i < tree->num_node[tree->start-1]; i++)
     loss -= tree->node[tree->start-1][i]->count * log(tree->node[tree->start-1][i]->p_node); */

  for(i = tree->start - 1; i < tree->start + tree->length - 1; i++) {
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      for(k = 0; k < curr_node->num_allele; k++) {
        if(curr_node->p_child[k] > 0)
          loss -= curr_node->count_child[k] * log(curr_node->p_child[k]);
        else
          loss -= curr_node->count_child[k] * log(1.0/(curr_node->count+1.0));
      }
    }
    loss += log(tree->node[0][0]->count)/2 * tree->num_node[i] * (tree->num_allele[i]-1);
  }

  return 2*loss;
}

double test_loss(struct tree *tree)
{
  int i, j, k;
  double loss = 0, pen = 0;
  struct tree_node *curr_node;
  FILE *outfile;

  /* for(i = 0; i < tree->num_node[tree->start-1]; i++)
     loss -= tree->node[tree->start-1][i]->count_test * log(tree->node[tree->start-1][i]->p_node); */

  for(i = tree->start - 1; i < tree->start + tree->length - 1; i++) {
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      for(k = 0; k < curr_node->num_allele; k++) {
        if(curr_node->p_child[k] > 0)
          loss -= curr_node->count_test_child[k] * log(curr_node->p_child[k]);
	else {
	  loss -= curr_node->count_test_child[k] * log(1.0/(curr_node->count+1.0)); 
	  //pen -= curr_node->count_test_child[k] * log(1.0/(curr_node->count+1000000.0));
	}
      }
    }
  }
  
  /*
  outfile = fopen("bgl2.pen", "a");
  fprintf(outfile, "%8.8lf %8.8lf\n", loss-pen, pen);
  fclose(outfile);
  */

  return loss;
}

double merge_test(struct tree_node *node1, struct tree_node *node2, double lambda, double diff, double min_diff)
{
  int i, n;
  double temp_diff;
  // double diff = 0;

  if(node1 == node2) return min_diff;

  if(min_diff < 0) return min_diff;

  for(i = 0; i < node1->num_allele; i++) {
    node1->p_child_test[i] = node2->p_child_test[i] = (node1->p_child[i]*node1->p_node + node2->p_child[i]*node2->p_node) / (node1->p_node + node2->p_node);
    if(node1->p_child_test[i] > 0)
      diff += (node1->count_child[i] + node2->count_child[i]) * log(node1->p_child_test[i]);
    else
      diff += (node1->count_child[i] + node2->count_child[i]) * log(1.0/(node1->count+node2->count+1.0));
    if(node1->p_child[i] > 0)
      diff -= node1->count_child[i] * log(node1->p_child[i]);
    else
      diff -= node1->count_child[i] * log(1.0/(node1->count+1.0));
    if(node2->p_child[i] > 0)
      diff -= node2->count_child[i] * log(node2->p_child[i]);
    else
      diff -= node2->count_child[i] * log(1.0/(node2->count+1.0));
  }
  diff += lambda * (node1->num_allele-1);
  //diff += lambda * (num_node[node1->level]-1) * (node1->num_allele-1);

  if(diff < min_diff) 
    min_diff = diff;

  for(i = 0; i < node1->num_allele; i++)
    if(!node1->allele_miss[i] && !node2->allele_miss[i]) {
      temp_diff = merge_test(node1->child[i], node2->child[i], lambda, diff, min_diff);
      if(temp_diff < min_diff)
	min_diff = temp_diff;
    }

  //num_node[node1->level]--;

  return min_diff;
} 

/* double merge_test(struct tree *tree, struct tree_node *node1, struct tree_node *node2, double lambda, int first)
{
  int i;

  if(first) tree->loss_test = tree->loss;

  if(node1 == node2) return tree->loss_test;

  for(i = 0; i < node1->num_allele; i++) {
    node1->p_child_test[i] = node2->p_child_test[i] = (node1->p_child[i]*node1->p_node + node2->p_child[i]*node2->p_node) / (node1->p_node + node2->p_node);
    tree->loss_test -= (node1->count_child[i] + node2->count_child[i]) * log(node1->p_child_test[i]) - node1->count_child[i] * log(node1->p_child[i]) - node2->count_child[i] * log(node2->p_child[i]);
    if(!node1->allele_miss[i] && !node2->allele_miss[i])
      merge_test(tree, node1->child[i], node2->child[i], lambda, 0);
  }
  tree->loss_test -= lambda * (node1->num_allele-1);
  
  return tree->loss_test;
  } */

struct tree_node *merge_node(struct tree *tree, struct tree_node *node1, struct tree_node *node2, double lambda, int first)
{
  int i, j, k, i_max = -1;
  double diff = 0, p_max = 0;
  struct tree_node *curr_node, *new = malloc(sizeof(struct tree_node));

  if(node1->id == node2->id) return node1;

  if(node1->id > node2->id) {
    curr_node = node1;
    node1 = node2;
    node2 = curr_node;
  }

  new->level = node1->level;
  new->id = node1->id;
  new->count = node1->count + node2->count;
  new->count_test = node1->count_test + node2->count_test;
  new->num_allele = node1->num_allele;
  new->num_back = node1->num_back;
  new->num_haplo = node1->num_haplo + node2->num_haplo;
  new->marked = 0;
  new->p_node = node1->p_node + node2->p_node;

  new->count_child = malloc(new->num_allele*sizeof(int));
  new->count_child = malloc(new->num_allele*sizeof(int));
  new->count_child = malloc(new->num_allele*sizeof(int));
  new->count_child = malloc(new->num_allele*sizeof(int));
  new->count_test_child = malloc(new->num_allele*sizeof(double));
  new->allele_miss = malloc(new->num_allele*sizeof(int));
  new->p_child = malloc(new->num_allele*sizeof(double));
  new->p_child_test = malloc(new->num_allele*sizeof(double));
  new->child = malloc(new->num_allele*sizeof(struct tree_node*));
  for(i = 0; i < new->num_allele; i++) {
    new->count_child[i] = node1->count_child[i] + node2->count_child[i];
    new->count_test_child[i] = node1->count_test_child[i] + node2->count_test_child[i];
    new->allele_miss[i] = node1->allele_miss[i] && node2->allele_miss[i];
    new->p_child[i] = new->p_child_test[i] = (node1->p_child[i]*node1->p_node + node2->p_child[i]*node2->p_node) / new->p_node;
    if(new->p_child[i] > 0)
      diff += new->count_child[i] * log(new->p_child[i]);
    else
      diff += new->count_child[i] * log(1.0/(new->count+1.0));
    if(node1->p_child[i] > 0)
      diff -= node1->count_child[i] * log(node1->p_child[i]);
    else
      diff -= node1->count_child[i] * log(1.0/(node1->count+1.0));
    if(node2->p_child[i] > 0)
      diff -= node2->count_child[i] * log(node2->p_child[i]);
    else
      diff -= node2->count_child[i] * log(1.0/(node2->count+1.0));
  }
  diff += lambda * (new->num_allele-1);
  //diff += lambda * (tree->num_node[node1->level]-1) * (new->num_allele-1);

  if(node1->level >= tree->start - 1 && node1->level < tree->start + tree->length)
    tree->loss -= diff;

  new->back_weight = malloc(new->num_back*sizeof(double));
  for(i = 0; i < new->num_back; i++)
    new->back_weight[i] = (node1->back_weight[i]*node1->p_node + node2->back_weight[i]*node2->p_node) / new->p_node; 

  for(i = 0; i < new->num_allele; i++) {
    if(new->allele_miss[i]) 
      new->child[i] = NULL;
    else {
      if(!node1->allele_miss[i] && node2->allele_miss[i])
        new->child[i] = node1->child[i];
      else if(node1->allele_miss[i] && !node2->allele_miss[i])
        new->child[i] = node2->child[i];
      else
        new->child[i] = merge_node(tree, node1->child[i], node2->child[i], lambda, 0);
    }
  }

  for(i = 0; i < new->num_allele; i++)
    if(new->p_child[i] > p_max) {
      p_max = new->p_child[i];
      i_max = i;
    }
  if(i_max >= 0)
    for(i = 0; i < new->num_allele; i++)
      if(new->allele_miss[i])
	new->child[i] = new->child[i_max];

  // haplo merging?

  if(new->level > 0) {
    for(i = 0; i < tree->num_node[new->level-1]; i++) {
      curr_node = tree->node[new->level-1][i];
      for(j = 0; j < curr_node->num_allele; j++)
	if(!curr_node->allele_miss[j])
	  if(curr_node->child[j]->id == node1->id || curr_node->child[j]->id == node2->id)
	    curr_node->child[j] = new;
    }
  }

  tree->node[node1->level][node1->id] = new;
  for(i = node2->id; i < tree->num_node[node1->level]-1; i++) {
    tree->node[node1->level][i] = tree->node[node1->level][i+1];
    tree->node[node1->level][i]->id--;
  }
  tree->num_node[node1->level]--;
  tree->node[node1->level] = realloc(tree->node[node1->level], tree->num_node[node1->level]*sizeof(struct tree_node*));

  free_node(node1);
  free_node(node2);

  return new;
}

/* struct tree_node *merge_node(struct tree *tree, struct tree_node *node1, struct tree_node *node2, double lambda, int first)
{
  int i, j, k;
  double diff = 0;
  struct tree_node *curr_node, *new = malloc(sizeof(struct tree_node));

  printf("Merging (inside) nodes (%d, %d) and (%d, %d)\n", node1->level, node1->id+1, node2->level, node2->id+1);

  if(node1 == node2) return node1;

  if(node1->id > node2->id) {
    curr_node = node1;
    node1 = node2;
    node2 = curr_node;
  }

  tree->node[node1->level][node1->id] = new;
  for(i = node2->id; i < tree->num_node[node1->level]-1; i++) {
    tree->node[node1->level][i] = tree->node[node1->level][i+1];
    tree->node[node1->level][i]->id--;
  }
  tree->num_node[node1->level]--;
  tree->node[node1->level] = realloc(tree->node[node1->level], tree->num_node[node1->level]*sizeof(struct tree_node*));

  new->level = node1->level;
  new->id = node1->id;
  new->count = node1->count + node2->count;
  new->count_test = node1->count_test + node2->count_test;
  new->num_allele = node1->num_allele;
  new->num_back = node1->num_back;
  new->num_haplo = node1->num_haplo + node2->num_haplo;
  new->marked = 0;
  new->p_node = node1->p_node + node2->p_node;

  new->count_child = malloc(new->num_allele*sizeof(int));
  new->count_test_child = malloc(new->num_allele*sizeof(double));
  new->allele_miss = malloc(new->num_allele*sizeof(int));
  new->p_child = malloc(new->num_allele*sizeof(double));
  new->p_child_test = malloc(new->num_allele*sizeof(double));
  new->child = malloc(new->num_allele*sizeof(struct tree_node*));
  for(i = 0; i < new->num_allele; i++) {
    new->count_child[i] = node1->count_child[i] + node2->count_child[i];
    new->count_test_child[i] = node1->count_test_child[i] + node2->count_test_child[i];
    new->allele_miss[i] = node1->allele_miss[i] && node2->allele_miss[i];
    new->p_child[i] = new->p_child_test[i] = (node1->p_child[i]*node1->p_node + node2->p_child[i]*node2->p_node) / new->p_node;
    diff += new->count_child[i] * log(new->p_child[i]) - node1->count_child[i] * log(node1->p_child[i]) - node2->count_child[i] * log(node2->p_child[i]);
  }
  diff += lambda * (new->num_allele-1);

  if(node1->level >= tree->start - 1 && node1->level < tree->start + tree->length)
    tree->loss -= diff;

  new->back_weight = malloc(new->num_back*sizeof(double));
  for(i = 0; i < new->num_back; i++)
    new->back_weight[i] = (node1->back_weight[i]*node1->p_node + node2->back_weight[i]*node2->p_node) / new->p_node;

  for(i = 0; i < new->num_allele; i++) {
    if(new->allele_miss[i])
      new->child[i] = NULL;
    else {
      if(!node1->allele_miss[i] && node2->allele_miss[i])
        new->child[i] = node1->child[i];
      else if(node1->allele_miss[i] && !node2->allele_miss[i])
        new->child[i] = node2->child[i];
      else {
        printf("Merging (outside) nodes (%d, %d) and (%d, %d)\n", node1->child[i]->level, node1->child[i]->id+1, node2->child[i]->level, node2->child[i]->id+1);
        new->child[i] = merge_node(tree, node1->child[i], node2->child[i], lambda, 0);
      }
    }
  }

  // haplo merging?

  if(first) {
    for(i = 0; i < tree->num_node[new->level-1]; i++) {
      curr_node = tree->node[new->level-1][i];
      for(j = 0; j < curr_node->num_allele; j++)
        if(curr_node->child[j] == node1 || curr_node->child[j] == node2)
          curr_node->child[j] = new;
    }
  }

  printf("Freeing nodes (%d, %d) and (%d, %d)\n", node1->level, node1->id+1, node2->level, node2->id+1);
  free_node(node1);
  free_node(node2);
  printf("Freed.\n");

  return new;
} */
/*
struct tree *merge_tree(struct tree *tree, double lambda, int start, int length, int reverse)
{
  int i, j, k, l, curr_label, num_node, *label;
  double min_diff;
  struct tree_node *curr_node1, *curr_node2, **node_list, ***node_ptr;

  if(reverse) start = tree->num_level - start - length + 1;

  tree->start = start;
  tree->length = length;

  loss(tree, lambda);

  for(i = 0; i < tree->num_level; i++) {
    printf("\t\t%d\n", i);
    if(tree->num_node[i] > MAX_NODE) return NULL;

    curr_label = 1;
    label = calloc(tree->num_node[i], sizeof(int));
    node_list = malloc(tree->num_node[i]*sizeof(struct tree_node*));
    node_ptr = malloc(tree->num_node[i]*sizeof(struct tree_node**));

    for(j = 0; j < tree->num_node[i]; j++) {
      if(!label[j]) 
	label[j] = curr_label++;
      curr_node1 = tree->node[i][j];
      node_list[j] = tree->node[i][j];
      node_ptr[j] = &(node_list[j]);
      for(k = j+1; k < tree->num_node[i]; k++) {
        curr_node2 = tree->node[i][k];
	min_diff = merge_test(tree, curr_node1, curr_node2, lambda, 0, lambda); // should be times max of num_allele-1
	if(min_diff >= 0) {
	  if(!label[k]) 
	    label[k] = label[j];
	  else
	    for(l = 0; l < tree->num_node[i]; l++)
	      if(label[l] == label[k]) 
		label[l] = label[j];
	}
      }
    }

    //printf("Level: %d\n\tLabels:", i);
    //for(j = 0; j < tree->num_node[i]; j++)
    //  printf(" %d", label[j]);
    //printf("\n");

    num_node = tree->num_node[i];
    for(j = 0; j < num_node; j++) {
      for(k = j+1; k < num_node; k++)
        if(label[j] == label[k]) {
          *node_ptr[j] = merge_node(tree, *node_ptr[j], *node_ptr[k], lambda, 1);
          node_ptr[k] = node_ptr[j];
        }
    }

    free(node_list);
    free(node_ptr);
  }

  return(tree);
}
*/

struct tree *merge_tree(struct tree *tree, double lambda, int start, int length, int reverse)
{
  int i, j, k, max_j, max_k, merge;
  double min_diff, max_diff;
  struct tree_node *curr_node1, *curr_node2;
  FILE *outfile;

  if(reverse) start = tree->num_level - start - length + 1;

  tree->start = start;
  tree->length = length;

  loss(tree, lambda);
  
  for(i = 0; i < tree->num_level; i++) {
    //printf("\t\t%d\n", i);
    if(tree->num_node[i] > MAX_NODE) return NULL;
    do {
      max_diff = 0;
      merge = 0;
      //printf("curr_loss = %8.8lf\n", curr_loss);
      for(j = 0; j < tree->num_node[i]-1; j++) {
        //printf("Checking (%d, %d)...\n", i, j+1);
        for(k = j+1; k < tree->num_node[i]; k++) {
          curr_node1 = tree->node[i][j];
          curr_node2 = tree->node[i][k];
	  min_diff = merge_test(curr_node1, curr_node2, lambda, 0, lambda); // should be times max of num_allele-1 
          //printf("\twith (%d, %d): min_diff = %8.8lf\n", i, k+1, min_diff);
	  if(min_diff >= 0 && min_diff >= max_diff) {
            max_diff = min_diff;
	    max_j = j;
	    max_k = k;
	    merge = 1;
	  }
          //printf("\n");
        }
      }
      if(merge) {
	//        printf("Merging (%d, %d) and (%d, %d).\n", i, max_j+1, i, max_k+1);
        curr_node1 = tree->node[i][max_j];
        curr_node2 = tree->node[i][max_k];
        merge_node(tree, curr_node1, curr_node2, lambda, 1);
      }
    } while(merge);
  }

  // print size?
  //outfile = fopen("size", "a");
  //if(outfile == NULL) {
  //  printf("Error: could not open file \"size\".\n");
  //  return;
  //}
  //fprintf(outfile, "%8.2lf", lambda);
  //for(i = 0; i < tree->num_level; i++)
  //  fprintf(outfile, " %d", tree->num_node[i]);
  //fprintf(outfile, "\n");
  //fclose(outfile); 

  return tree;
} 

/* struct tree *merge_tree(struct tree *tree, double lambda)
{
  int i, j, k, l, new_test, max_j, max_k;
  double curr_loss, new_loss, forw_loss, max_diff;
  struct tree_node *curr_node1, *curr_node2;
  FILE *outfile;

  loss(tree, lambda);
  
  for(i = 0; i < tree->num_level; i++) {
    if(tree->num_node[i] > MAX_NODE) return NULL;
    do {
      max_diff = 0;
      curr_loss = tree->loss;
      //printf("curr_loss = %8.8lf\n", curr_loss);
      for(j = 0; j < tree->num_node[i]; j++) {
	//printf("Checking (%d, %d)...\n", i, j+1);
	for(k = j+1; k < tree->num_node[i]; k++) {
	  curr_node1 = tree->node[i][j];
	  curr_node2 = tree->node[i][k];
	  new_loss = merge_test(tree, curr_node1, curr_node2, lambda, 1);
	  //printf("\twith (%d, %d): Curr = %8.8lf", i, k+1, curr_loss - new_loss);
	  //printf("\tnew_loss = %8.8lf\n", new_loss);
	  if(curr_loss > new_loss) {
	    new_test = 1;
	    for(l = 0; l < tree->num_allele[i]; l++)
	      if(!curr_node1->allele_miss[l] && !curr_node2->allele_miss[l]) {
		merge_test(tree, curr_node1->child[l], curr_node2->child[l], lambda, new_test);
		new_test = 0;
	      }
	    forw_loss = tree->loss_test;
	    //printf("  Forw = %8.8lf", forw_loss - new_loss);
	    //printf("\t\tforw_loss = %8.8lf\n", forw_loss);
	    if(forw_loss - new_loss > max_diff) {
	      max_diff = forw_loss - new_loss;
	      max_j = j;
	      max_k = k;
	    }
	  }
	  //printf("\n");
	}
      }
      if(max_diff > 0) {
	//printf("Merging (%d, %d) and (%d, %d).\n", i, max_j+1, i, max_k+1);
	curr_node1 = tree->node[i][max_j];
	curr_node2 = tree->node[i][max_k];
	merge_node(tree, curr_node1, curr_node2, lambda, 1);
      }
    } while(max_diff > 0);
  }

  outfile = fopen("size", "a");
  if(outfile == NULL) {
    printf("Error: could not open file \"size\".\n");
    return;
  }
  fprintf(outfile, "%8.2lf", lambda);
  for(i = 0; i < tree->num_level; i++)
    fprintf(outfile, " %d", tree->num_node[i]);
  fprintf(outfile, "\n");
  fclose(outfile);

  return tree;
  } */

void print_view(char *filename, struct tree *tree)
{
  int i, j, k, sum, **node_list;
  char *temp;
  struct tree_node *curr_node;
  FILE *outfile;

  temp = malloc((strlen(filename)+7)*sizeof(char));
  strcpy(temp, filename);
  strcat(temp, ".nodes");
  outfile = fopen(temp, "w");
  if(outfile == NULL) {
    printf("Error: could not open file \"%s\".\n", temp);
    return;
  }
  free(temp);

  fprintf(outfile, "ID LEVEL NUMBER COUNT\n");

  node_list = malloc(tree->num_level*sizeof(int*));
  sum = 0;
  for(i = 0; i < tree->num_level; i++) {
    node_list[i] = malloc(tree->num_node[i]*sizeof(int));
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      fprintf(outfile, "%d.%d %d %d %3.3lf\n", curr_node->level, curr_node->id+1, curr_node->level, curr_node->id+1, curr_node->p_node);
      node_list[i][j] = sum + curr_node->id;
    }
    sum += tree->num_node[i];
  }

  fclose(outfile);

  temp = malloc((strlen(filename)+7)*sizeof(char));
  strcpy(temp, filename);
  strcat(temp, ".edges");
  outfile = fopen(temp, "w");
  if(outfile == NULL) {
    printf("Error: could not open file \"%s\".\n", temp);
    return;
  }
  free(temp);

  fprintf(outfile, "NODE1 NODE2 ALLELE COUNT\n");

  for(i = 0; i < tree->num_level-1; i++)
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      for(k = 0; k < curr_node->num_allele; k++)
	if(!curr_node->allele_miss[k])
	  fprintf(outfile, "%d %d %d %3.3lf\n", node_list[i][j], node_list[curr_node->child[k]->level][curr_node->child[k]->id], k+1, curr_node->p_child[k]);
    }

  fclose(outfile);
  
  for(i = 0; i < tree->num_level; i++)
    free(node_list[i]);
  free(node_list);
}

void calc_tree(struct tree *tree1, struct tree *tree2, int print)
{
  int i, j, k, curr_level, curr_id, back_level, back_id, check;
  double child_sum, back_mult;
  struct tree_node *curr_node, *back_node;
  // struct queue *queue;

  // queue = new_queue();
  curr_node = tree1->node[0][0];
  if(curr_node->num_back != 1) {
    curr_node->num_back = 1;
    curr_node->back_weight = malloc(sizeof(double));
    curr_node->back_weight[0] = 1;
  }
  // enqueue(queue, curr_node);

  for(curr_level = 0; curr_level < tree1->num_level; curr_level++) {
    back_level = tree2->num_level - curr_level - 2;
    for(curr_id = 0; curr_id < tree1->num_node[curr_level]; curr_id++) {
      if(print && curr_level >= 0) printf("Level: %d, ID: %d\n", curr_level, curr_id+1);
      curr_node = tree1->node[curr_level][curr_id];
      curr_node->marked = 0;
      child_sum = 0;
      for(i = 0; i < curr_node->num_allele; i++) {
	if(!curr_node->allele_miss[i] && !curr_node->child[i]->marked) {
	  curr_node->child[i]->marked = 1;
	  curr_node->child[i]->p_node = 0;
	  if(curr_node->child[i]->num_back != tree2->num_node[back_level]) {
	    curr_node->child[i]->num_back = tree2->num_node[back_level];
	    free(curr_node->child[i]->back_weight);
	    curr_node->child[i]->back_weight = calloc(curr_node->child[i]->num_back, sizeof(double));
	  }
	  else
	    for(j = 0; j < tree2->num_node[back_level]; j++)
	      curr_node->child[i]->back_weight[j] = 0;
	}
	curr_node->p_child[i] = 0;
	if(print && curr_level >= 0) printf("\tAllele %d = ", i+1);
	for(j = 0; j < tree2->num_node[back_level]; j++) {
	  back_node = tree2->node[back_level][j];
	  if(back_node->allele_miss[i])
	    back_mult = 0;
	  else 
	    back_mult = curr_node->back_weight[back_node->child[i]->id];
	  if(curr_node->p_node > 0) {
	    curr_node->p_child[i] += back_mult * back_node->p_child[i] * back_node->p_node / curr_node->p_node;
	    if(print && curr_level >= 0) printf("%4.4lf * %4.4lf * %4.4lf / %4.4lf", back_mult, back_node->p_child[i], back_node->p_node, curr_node->p_node); 
	  }
	  else if(print && curr_level >= 0) printf("0");
	  if(print && curr_level >= 0) printf("\n\t         + ");
	  if(print && curr_level >= 0 && j == tree2->num_node[back_level]-1) printf(" = %4.4lf\n", curr_node->p_child[i]);
	  if(curr_node->allele_miss[i] && curr_node->p_child[i] > 0) {
	    if(curr_level+1 == tree1->num_level-1)
	      curr_node->child[i] = tree1->node[tree1->num_level-1][0];
	    else {
	      tree1->num_node[curr_level+1]++;
	      tree1->node[curr_level+1] = realloc(tree1->node[curr_level+1], tree1->num_node[curr_level+1]*sizeof(struct tree_node*));
	      tree1->node[curr_level+1][tree1->num_node[curr_level+1]-1] = new_node(curr_level+1, tree1->num_node[curr_level+1]-1, tree1->num_allele[curr_level+1]);
	      curr_node->child[i] = tree1->node[curr_level+1][tree1->num_node[curr_level+1]-1];
	      curr_node->child[i]->marked = 1;
	      curr_node->child[i]->p_node = 0;
	      curr_node->child[i]->num_back = tree2->num_node[back_level];
	      curr_node->child[i]->back_weight = calloc(curr_node->child[i]->num_back, sizeof(double));
	    }
	    curr_node->allele_miss[i] = 0;
	  }
	  if(!curr_node->allele_miss[i])
	    curr_node->child[i]->back_weight[back_node->id] += back_mult * back_node->p_child[i];
	}
	child_sum += curr_node->p_child[i];
      }
      for(i = 0; i < curr_node->num_allele; i++) {
	if(child_sum > 0)
	  curr_node->p_child[i] /= child_sum;
	if(!curr_node->allele_miss[i]) {
	  //printf("\tAllele %d = , Node ID = %d, p_node += %8.8lf\n", i+1, curr_node->child[i]->id+1, curr_node->p_node * curr_node->p_child[i]);
	  curr_node->child[i]->p_node += curr_node->p_node * curr_node->p_child[i];
	}
      }
    }
  }      

  /*
  while(!empty(queue)) {
    curr_node = dequeue(queue);
    curr_node->marked = 0;
    back_level = tree2->num_level - curr_node->level - 2;
    p_sum = 0;
    for(i = 0; i < curr_node->num_allele; i++) {
      if(!curr_node->allele_miss[i] && !curr_node->child[i]->marked) {
        enqueue(queue, curr_node->child[i]);
        curr_node->child[i]->marked = 1;
        curr_node->child[i]->p_node = 0;
	if(curr_node->child[i]->num_back != tree2->num_node[back_level]) {
	  curr_node->child[i]->num_back = tree2->num_node[back_level];
	  free(curr_node->child[i]->back_weight);
	  curr_node->child[i]->back_weight = calloc(curr_node->child[i]->num_back, sizeof(double));
	}
	else
	  for(j = 0; j < tree2->num_node[back_level]; j++)
	    curr_node->child[i]->back_weight[j] = 0;
      }
      curr_node->p_child[i] = 0;
      for(j = 0; j < tree2->num_node[back_level]; j++) {
	back_node = tree2->node[back_level][j];
	if(back_node->allele_miss[i]) {
	  back_mult = 0;
	  for(k = 0; k < tree2->num_node[back_level+1]; k++) 
	    back_mult += curr_node->back_weight[k] * tree2->node[back_level+1][k]->p_node;
	}
	else back_mult = curr_node->back_weight[back_node->child[i]->id];
	curr_node->p_child[i] += back_mult * back_node->p_child[i] * back_node->p_node / curr_node->p_node;
	if(!curr_node->allele_miss[i])
	  curr_node->child[i]->back_weight[back_node->id] += back_mult * back_node->p_child[i];
      }
      p_sum += curr_node->p_child[i];
    }
    printf("%lf8.8\n", p_sum);
    for(i = 0; i < curr_node->num_allele; i++) {
      curr_node->p_child[i] /= p_sum;
      curr_node->child[i]->p_node += curr_node->p_node * curr_node->p_child[i];
    }
  }
  */
}

int calc_freq(struct haplotype ***results, struct tree* tree, int start, int length, int reverse)
{
  int i, j, k, l, m, index, num_haplo = 0, *temp;
  struct tree_node *curr_node;

  if(reverse) start = tree->num_level - start - length + 1;

  for(i = 0; i < tree->num_node[start-1]; i++) {
    curr_node = tree->node[start-1][i];
    curr_node->num_haplo = 1;
    curr_node->haplo = malloc(sizeof(struct haplotype*));
    curr_node->haplo[0] = malloc(sizeof(struct haplotype));
    curr_node->haplo[0]->length = 0;
    curr_node->haplo[0]->allele = NULL;
    curr_node->haplo[0]->freq = log10(curr_node->p_node);
  }

  for(i = start-1; i < start+length; i++)
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      for(k = 0; k < curr_node->num_allele; k++) {
	if(!curr_node->allele_miss[k]) { // remove?
	  curr_node->child[k]->haplo = realloc(curr_node->child[k]->haplo, (curr_node->child[k]->num_haplo+curr_node->num_haplo)*sizeof(struct haplotype*));
	  for(l = 0; l < curr_node->num_haplo; l++) {
	    index = curr_node->child[k]->num_haplo + l;
	    curr_node->child[k]->haplo[index] = malloc(sizeof(struct haplotype));
	    curr_node->child[k]->haplo[index]->length = curr_node->haplo[l]->length + 1;
	    curr_node->child[k]->haplo[index]->allele = malloc(curr_node->child[k]->haplo[index]->length*sizeof(int));
	    memcpy(curr_node->child[k]->haplo[index]->allele, curr_node->haplo[l]->allele, curr_node->haplo[l]->length*sizeof(int));
	    curr_node->child[k]->haplo[index]->allele[curr_node->haplo[l]->length] = k+1;  /* allele fix */
	    curr_node->child[k]->haplo[index]->freq = curr_node->haplo[l]->freq + log10(curr_node->p_child[k]);
	  }
	  curr_node->child[k]->num_haplo += curr_node->num_haplo;
	}
      }
    }
  /* free old haplos too */

  /* free old results? */
  *results = NULL;
  for(i = 0; i < tree->num_node[start+length-1]; i++) {
    curr_node = tree->node[start+length-1][i];
    *results = realloc(*results, (num_haplo+curr_node->num_haplo)*sizeof(struct haplotype*));
    for(j = 0; j < curr_node->num_haplo; j++)
      (*results)[num_haplo+j] = curr_node->haplo[j];
    num_haplo += curr_node->num_haplo;
    curr_node->num_haplo = 0;
  }

  qsort(*results, num_haplo, sizeof(struct haplotype*), compare);

  j = 0;
  for(i = 1; i < num_haplo; i++) {
    if(compare((*results)+j, (*results)+i) == 0) {
      if((*results)[i]->freq == -1.0/0.0)
	continue;
      if((*results)[j]->freq == -1.0/0.0)
	(*results)[j]->freq = (*results)[i]->freq;
      else
	(*results)[j]->freq += log10(1 + pow(10, (*results)[i]->freq - (*results)[j]->freq));
    }
    else (*results)[++j] = (*results)[i];
  }
  num_haplo = j + 1;
  *results = realloc(*results, num_haplo*sizeof(struct haplotype*));

  if(reverse) {
    for(i = 0; i < num_haplo; i++) {
      temp = malloc((*results)[i]->length*sizeof(int));
      memcpy(temp, (*results)[i]->allele, (*results)[i]->length*sizeof(int));
      for(j = 0; j < (*results)[i]->length; j++)
        (*results)[i]->allele[j] = temp[(*results)[i]->length-j-1];
      free(temp);
    }
    qsort(*results, num_haplo, sizeof(struct haplotype*), compare);
  }

  return num_haplo;
}

void print_haplo(char* filename, struct haplotype **results, int num_haplo)
{
  FILE *outfile;
  int i, j;

  outfile = fopen(filename, "w");
  if(outfile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return;
  }

  for(i = 0; i < num_haplo; i++) {
    for (j = 0; j < results[i]->length; j++)
      fprintf(outfile, "%d", results[i]->allele[j]);
    fprintf(outfile, " %8.8lf\n", results[i]->freq);
  }

  fclose(outfile);
}

int compare(const void* a, const void* b)
{
  int i, n = (*(struct haplotype**)a)->length;

  for(i = 0; i < (*(struct haplotype**)a)->length; i++)
    if((*(struct haplotype**)a)->allele[i] != (*(struct haplotype**)b)->allele[i])
      return (*(struct haplotype**)a)->allele[i] - (*(struct haplotype**)b)->allele[i];

  return 0;
}

