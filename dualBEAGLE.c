#include "dualBEAGLE.h"

struct tree *read_input(char *filename)
{
  FILE *infile;
  char line[MAX_LINE], *line_ptr;
  int i, j, k, n, curr_level, curr_index, curr_allele, max_allele, num_samp, num_level = START_SIZE, *num_node, *num_allele, *node_list;
  double curr_p_sum, new_p_sum, p_temp;
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

  return(tree);
}

struct tree *read_bgl(char *filename)
{
  FILE *infile;
  char line[MAX_LINE];
  int i, j, k, prev_level, curr_level, curr_index, next_index, max_index, curr_allele, max_allele, num_level = START_SIZE, curr_count, *num_node, *num_allele;
  double curr_p_sum, new_p_sum, p_temp;
  struct tree_node *curr_node, ***node;
  struct tree *tree;

  /* open the input file */
  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return NULL;
  }

  /* skip first two lines */
  fgets(line, MAX_LINE, infile);
  fgets(line, MAX_LINE, infile);

  /* read input to find number and size of levels */
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

  /* make nodes */
  node = malloc(num_level*sizeof(struct tree_node**));
  for(i = 0; i < num_level; i++) {
    node[i] = malloc(num_node[i]*sizeof(struct tree_node*));
    for(j = 0; j < num_node[i]; j++)
      node[i][j] = new_node(i, j, num_allele[i]);
  }

  tree->node = node;

  /* read input to get edges */
  while(fgets(line, MAX_LINE, infile) != NULL) {
    if(line[0] != '\n') {
      sscanf(line, "%d %*s %d %d %d %d", &curr_level, &curr_index, &next_index, &curr_allele, &curr_count);
      curr_allele--;  // fix to allele key later
      add_edge(tree, node[curr_level][curr_index], node[curr_level+1][next_index], curr_allele, curr_count);
    }
  }

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
    rand_u = (double)rand() / (double)RAND_MAX;
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
}

struct tree_node *new_node(int level, int id, int num_allele)
{
  int i;
  struct tree_node *new = malloc(sizeof(struct tree_node));

  new->level = level;
  new->id = id;
  new->count = new->num_back = new->num_haplo = new->marked = 0;
  new->num_allele = num_allele;
  new->count_child = calloc(num_allele, sizeof(int));
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
  new->num_allele = node->num_allele;
  new->num_back = node->num_back;
  new->num_haplo = node->num_haplo;
  new->marked = node->marked;

  new->count_child = malloc(node->num_allele*sizeof(int));
  memcpy(new->count_child, node->count_child, node->num_allele*sizeof(int));
  new->allele_miss = malloc(node->num_allele*sizeof(int));
  memcpy(new->allele_miss, node->allele_miss, node->num_allele*sizeof(int));

  new->p_node = node->p_node;

  new->p_child = malloc(node->num_allele*sizeof(double));
  memcpy(new->p_child, node->p_child, node->num_allele*sizeof(double));
  new->p_child_test = malloc(node->num_allele*sizeof(double));
  memcpy(new->p_child_test, node->p_child_test, node->num_allele*sizeof(double));
  new->back_weight = malloc(node->num_back*sizeof(double)); 
  memcpy(new->back_weight, node->back_weight, node->num_back*sizeof(double));

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
    if(reverse) fprintf(outfile, "Level %d: %d\n", tree->num_level-i-1, tree->num_node[i]);
    else fprintf(outfile, "Level %d: %d\n", i, tree->num_node[i]);
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      fprintf(outfile, "\tNode ID = %d, Node Weight = %8.8lf\n", curr_node->id+1, curr_node->p_node);
      for(k = 0; k < curr_node->num_allele; k++) {
	if(curr_node->allele_miss[k])
	  fprintf(outfile, "\t\tChild Node = %d, Child Allele = %d, Child Weight = %8.8lf\n", 0, k+1, curr_node->p_child[k]);
	else
	  fprintf(outfile, "\t\tChild Node = %d, Child Allele = %d, Child Weight = %8.8lf\n", curr_node->child[k]->id+1, k+1, curr_node->p_child[k]);
      }
    }
  }

  fclose(outfile);
}

double loss(struct tree *tree, double lambda)
{
  int i, j, k;
  double loss = 0;
  struct tree_node *curr_node;

  for(i = 0; i < tree->num_level; i++) {
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      for(k = 0; k < curr_node->num_allele; k++) 
        if(curr_node->count_child[k] > 0)
          loss -= curr_node->count_child[k] * log(curr_node->p_child[k]);
    }
    loss += lambda * tree->num_node[i] * (tree->num_allele[i]-1);
  }

  tree->loss = loss;
  return loss;
}

double merge_test(struct tree *tree, struct tree_node *node1, struct tree_node *node2, double lambda, double min_diff, int first)
{
  int i;
  double diff;

  if(first) tree->loss_test = tree->loss;

  if(node1 == node2) return min_diff;

  diff = 0;
  for(i = 0; i < node1->num_allele; i++) {
    node1->p_child_test[i] = node2->p_child_test[i] = (node1->p_child[i]*node1->p_node + node2->p_child[i]*node2->p_node) / (node1->p_node + node2->p_node);
    diff += (node1->count_child[i] + node2->count_child[i]) * log(node1->p_child_test[i]) - node1->count_child[i] * log(node1->p_child[i]) - node2->count_child[i] * log(node2->p_child[i]);
  }
  diff += lambda * (node1->num_allele-1);
  tree->loss_test -= diff;

  if(diff < min_diff) min_diff = diff;
  for(i = 0; i < node1->num_allele; i++)
    if(!node1->allele_miss[i] && !node2->allele_miss[i]) {
      diff = merge_test(tree, node1->child[i], node2->child[i], lambda, min_diff, 0);
      if(diff < min_diff) min_diff = diff;
    }
  
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
  int i, j, k;
  struct tree_node *curr_node, *new = malloc(sizeof(struct tree_node));

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
  new->num_allele = node1->num_allele;
  new->num_back = node1->num_back;
  new->num_haplo = node1->num_haplo + node2->num_haplo;
  new->marked = 0;
  new->p_node = node1->p_node + node2->p_node;

  new->count_child = malloc(new->num_allele*sizeof(int));
  new->allele_miss = malloc(new->num_allele*sizeof(int));
  new->p_child = malloc(new->num_allele*sizeof(double));
  new->p_child_test = malloc(new->num_allele*sizeof(double));
  new->child = malloc(new->num_allele*sizeof(struct tree_node*));
  for(i = 0; i < new->num_allele; i++) {
    new->count_child[i] = node1->count_child[i] + node2->count_child[i];
    new->allele_miss[i] = node1->allele_miss[i] && node2->allele_miss[i];
    new->p_child[i] = new->p_child_test[i] = (node1->p_child[i]*node1->p_node + node2->p_child[i]*node2->p_node) / new->p_node;
    tree->loss -= new->count_child[i] * log(new->p_child[i]) - node1->count_child[i] * log(node1->p_child[i]) - node2->count_child[i] * log(node2->p_child[i]);
  }
  tree->loss -= lambda * (new->num_allele-1);

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

  // haplo merging?

  if(first) {
    for(i = 0; i < tree->num_node[new->level-1]; i++) {
      curr_node = tree->node[new->level-1][i];
      for(j = 0; j < curr_node->num_allele; j++) 
	if(curr_node->child[j] == node1 || curr_node->child[j] == node2)
	  curr_node->child[j] = new;
    }
  }
  
  free_node(node1);
  free_node(node2);

  return new;
}

struct tree *merge_tree(struct tree *tree, double lambda)
{
  int i, j, k, max_j, max_k;
  double min_diff, max_diff;
  struct tree_node *curr_node1, *curr_node2;
  FILE *outfile;

  loss(tree, lambda);

  for(i = 0; i < tree->num_level; i++) {
    printf("%d\n", i);
    if(tree->num_node[i] > MAX_NODE) return NULL;
    do {
      max_diff = 0;
      //printf("curr_loss = %8.8lf\n", curr_loss);
      for(j = 0; j < tree->num_node[i]; j++) {
        //printf("Checking (%d, %d)...\n", i, j+1);
        for(k = j+1; k < tree->num_node[i]; k++) {
          curr_node1 = tree->node[i][j];
          curr_node2 = tree->node[i][k];
          min_diff = merge_test(tree, curr_node1, curr_node2, lambda, lambda, 1);
          //printf("\twith (%d, %d): Curr = %8.8lf", i, k+1, curr_loss - new_loss);
          //printf("\tnew_loss = %8.8lf\n", new_loss);
          if(min_diff > max_diff) {
              max_diff = min_diff;
              max_j = j;
              max_k = k;
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
      fprintf(outfile, "%d.%d %d %d %d\n", curr_node->level, curr_node->id+1, curr_node->level, curr_node->id+1, curr_node->count);
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
	  fprintf(outfile, "%d %d %d %d\n", node_list[i][j], node_list[curr_node->child[k]->level][curr_node->child[k]->id], k+1, curr_node->count_child[k]);
    }

  fclose(outfile);
  
  for(i = 0; i < tree->num_level; i++)
    free(node_list[i]);
  free(node_list);
}

double test_loss(char *filename, struct tree *tree)
{
  FILE *infile;
  char line[MAX_LINE], *line_ptr;
  int i, j, n, curr_level, curr_allele, curr_samp, num_samp;
  double loss, *temp, **node_list;
  struct tree_node *curr_node;

  /* open the input file */
  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return -1;
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

  print_tree("test", tree, 0);

  /* read input to calculate BIC */
  loss = 0;
  curr_level = 0;
  node_list = malloc(num_samp*sizeof(double*));
  for(curr_samp = 0; curr_samp < num_samp; curr_samp++) {
    node_list[curr_samp] = malloc(sizeof(double));
    node_list[curr_samp][0] = 1;
  }
  while(fgets(line, MAX_LINE, infile) != NULL) {
    if(line[0] != '\n') {
      line_ptr = line;
      sscanf(line_ptr, "%*s %*s %n", &n);
      line_ptr += n;

      curr_samp = 0;
      while(sscanf(line_ptr, "%d %n", &curr_allele, &n) == 1) {
	curr_allele--;  // fix to allele key later
	temp = calloc(tree->num_node[curr_level+1], sizeof(double));
	for(i = 0; i < tree->num_node[curr_level]; i++) {
	  curr_node = tree->node[curr_level][i];
	  if(node_list[curr_samp][i] > 0) 
	    loss -= log(curr_node->p_child[curr_allele]) * node_list[curr_samp][i];
	  if(curr_node->allele_miss[curr_allele])
	    for(j = 0; j < tree->num_node[curr_level+1]; j++)
	      temp[j] += node_list[curr_samp][i] * tree->node[curr_level+1][j]->p_node;
	  else
	    temp[curr_node->child[curr_allele]->id] += node_list[curr_samp][i];
	}
	free(node_list[curr_samp]);
	node_list[curr_samp] = temp;
	curr_samp++;
	line_ptr += n;
      }
      
      // loss += log(num_samp) * tree->num_node[curr_level] * (tree->num_allele[curr_level] - 1);
      curr_level++;
    }
  }
  
  fclose(infile);

  return(loss);
}

void calc_tree(struct tree *tree1, struct tree *tree2)
{
  int i, j, k, curr_level, curr_id, back_level;
  double child_sum, node_sum, back_mult;
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
    node_sum = 0;
    for(curr_id = 0; curr_id < tree1->num_node[curr_level]; curr_id++) 
      node_sum += tree1->node[curr_level][curr_id]->p_node;
    for(curr_id = 0; curr_id < tree1->num_node[curr_level]; curr_id++) {
      //printf("Level: %d, ID: %d\n", curr_level, curr_id+1);
      curr_node = tree1->node[curr_level][curr_id];
      curr_node->p_node /= node_sum;
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
	//printf("\tAllele %d = ", i+1);
	for(j = 0; j < tree2->num_node[back_level]; j++) {
	  back_node = tree2->node[back_level][j];
	  if(back_node->allele_miss[i]) {
	    back_mult = 0;
	    for(k = 0; k < tree2->num_node[back_level+1]; k++)
	      back_mult += curr_node->back_weight[k] * tree2->node[back_level+1][k]->p_node;
	  }
	  else back_mult = curr_node->back_weight[back_node->child[i]->id];
	  curr_node->p_child[i] += back_mult * back_node->p_child[i] * back_node->p_node / curr_node->p_node;
	  //printf("%4.4lf * %4.4lf * %4.4lf / %4.4lf", back_mult, back_node->p_child[i], back_node->p_node, curr_node->p_node); 
	  //if(j < tree2->num_node[back_level]-1) printf("\n\t         + ");
	  //else printf(" = %4.4lf\n", curr_node->p_child[i]);
	  if(!curr_node->allele_miss[i])
	    curr_node->child[i]->back_weight[back_node->id] += back_mult * back_node->p_child[i];
	}
	child_sum += curr_node->p_child[i];
      }
      for(i = 0; i < curr_node->num_allele; i++) {
	curr_node->p_child[i] /= child_sum;
	if(!curr_node->allele_miss[i])
	  curr_node->child[i]->p_node += curr_node->p_node * curr_node->p_child[i];
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
    curr_node->haplo[0]->freq = curr_node->p_node;
  }

  for(i = start-1; i < start+length; i++)
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      printf("%d - %d\n", curr_node->level, curr_node->id+1);
      for(k = 0; k < curr_node->num_allele; k++) {
	if(curr_node->allele_miss[k]) {
	  for(l = 0; l < tree->num_node[i+1]; l++) {
	    tree->node[i+1][l]->haplo = realloc(tree->node[i+1][l]->haplo, (tree->node[i+1][l]->num_haplo+curr_node->num_haplo)*sizeof(struct haplotype*));
	    for(m = 0; m < curr_node->num_haplo; m++) {
	      index = tree->node[i+1][l]->num_haplo + m;
	      tree->node[i+1][l]->haplo[index] = malloc(sizeof(struct haplotype));
	      tree->node[i+1][l]->haplo[index]->length = curr_node->haplo[m]->length + 1;
	      tree->node[i+1][l]->haplo[index]->allele = malloc(tree->node[i+1][l]->haplo[index]->length*sizeof(int));
	      memcpy(tree->node[i+1][l]->haplo[index]->allele, curr_node->haplo[m]->allele, curr_node->haplo[m]->length*sizeof(int));
	      tree->node[i+1][l]->haplo[index]->allele[curr_node->haplo[m]->length] = k+1;  /* allele fix */
	      tree->node[i+1][l]->haplo[index]->freq = curr_node->haplo[m]->freq * curr_node->p_child[k] * tree->node[i+1][l]->p_node;
	    }
	    tree->node[i+1][l]->num_haplo += curr_node->num_haplo;
	  }
	}
	else {
	  curr_node->child[k]->haplo = realloc(curr_node->child[k]->haplo, (curr_node->child[k]->num_haplo+curr_node->num_haplo)*sizeof(struct haplotype*));
	  for(l = 0; l < curr_node->num_haplo; l++) {
	    index = curr_node->child[k]->num_haplo + l;
	    curr_node->child[k]->haplo[index] = malloc(sizeof(struct haplotype));
	    curr_node->child[k]->haplo[index]->length = curr_node->haplo[l]->length + 1;
	    curr_node->child[k]->haplo[index]->allele = malloc(curr_node->child[k]->haplo[index]->length*sizeof(int));
	    memcpy(curr_node->child[k]->haplo[index]->allele, curr_node->haplo[l]->allele, curr_node->haplo[l]->length*sizeof(int));
	    curr_node->child[k]->haplo[index]->allele[curr_node->haplo[l]->length] = k+1;  /* allele fix */
	    curr_node->child[k]->haplo[index]->freq = curr_node->haplo[l]->freq * curr_node->p_child[k];
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
    if(compare((*results)+j, (*results)+i) == 0)
      (*results)[j]->freq += (*results)[i]->freq;
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

int compare(const void* a, const void* b)
{
  int i, n = (*(struct haplotype**)a)->length;

  for(i = 0; i < (*(struct haplotype**)a)->length; i++)
    if((*(struct haplotype**)a)->allele[i] != (*(struct haplotype**)b)->allele[i])
      return (*(struct haplotype**)a)->allele[i] - (*(struct haplotype**)b)->allele[i];

  return 0;
}


