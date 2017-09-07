#include "dualBEAGLE.h"
#include <sys/time.h>

struct tree *read_input(char *filename)
{
  FILE *infile;
  char line[MAX_LINE];
  int i, j, k, curr_level, curr_index, curr_allele, max_allele, num_samp, num_level, *num_node, *num_allele, *node_list;
  struct tree_node *curr_node, ***node;
  struct tree *tree;

  // open the input file
  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return NULL;
  }

  // reads opening line and calculates sample size
  fgets(line, MAX_LINE, infile);
  for(i = 0; line[i] > '0'; i++);
  num_samp = i;
  
  // read input to find number and size of levels
  num_level = 1;
  num_allele = malloc(num_level*sizeof(int));
  curr_level = 0;
  max_allele = 1; // this was 2? why?
  do {
    for(i = 0; line[i] > '0'; i++) {
      curr_allele = line[i] - '0';
      if(curr_allele > max_allele)
	max_allele = curr_allele;
    }
        
    num_allele[curr_level] = max_allele;
    max_allele = 1;
    
    curr_level++;
    if(curr_level == num_level) {
      num_level *= 2;
      num_allele = realloc(num_allele, num_level*sizeof(int));
    }
  } while(fgets(line, MAX_LINE, infile) != NULL);
  num_level = curr_level + 1;
  num_allele = realloc(num_allele, num_level*sizeof(int));
  num_allele[num_level-1] = 0;
  
  rewind(infile);

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
    if(curr_level+1 < num_level-1) {
      node[curr_level+1] = malloc(num_allele[curr_level]*num_node[curr_level]*sizeof(struct tree_node*));
      curr_index = 0;
      for(i = 0; line[i] > '0'; i++) {
	//printf("%d\n", curr_index);
	curr_allele = line[i] - '1';  // fix to allele key later
	curr_node = node[curr_level][node_list[curr_index]];
	curr_node->count++;
	curr_node->count_child[curr_allele]++;
	if(curr_node->allele_miss[curr_allele]) {
	  //printf("%d\n", curr_index);
	  node[curr_level+1][num_node[curr_level+1]] = new_node(curr_level+1, num_node[curr_level+1], num_allele[curr_level+1]);
	  add_edge(tree, curr_node, node[curr_level+1][num_node[curr_level+1]], curr_allele);
	  node_list[curr_index] = num_node[curr_level+1];
	  num_node[curr_level+1]++;
	  //printf("%d\n", curr_index);
	}
	else {
	  //printf("\t%d\n", curr_allele);
	  node_list[curr_index] = curr_node->child[curr_allele]->node_to->id;
	  //printf("\t%d\n", curr_allele);
	}
	curr_index++;
      }
    }
    
    else {
      node[curr_level+1] = malloc(sizeof(struct tree_node*));
      node[curr_level+1][0] = new_node(curr_level+1, 0, num_allele[curr_level+1]);
      num_node[curr_level+1] = 1;
      curr_index = 0;
      for(i = 0; line[i] > '0'; i++) {
        curr_allele = line[i] - '1';  // fix to allele key later
	curr_node = node[curr_level][node_list[curr_index]];
	curr_node->count++;
	curr_node->count_child[curr_allele]++;
	node[curr_level+1][0]->count++;
	if(curr_node->allele_miss[curr_allele]) 
	  add_edge(tree, curr_node, node[curr_level+1][0], curr_allele);
	curr_index++;
      }
    }
    
    node[curr_level+1] = realloc(node[curr_level+1], num_node[curr_level+1]*sizeof(struct tree_node*));
    curr_level++;
  } 

  tree->num_node = num_node;
  tree->node = node;
  free(node_list);
  
  // done with input
  fclose(infile);

  // calculate edge and node weights 
  for(i = 0; i < num_level; i++) {
    for(j = 0; j < num_node[i]; j++) {
      curr_node = node[i][j];
      for(k = 0; k < curr_node->num_allele; k++) {
        curr_node->p_child[k] = (double)curr_node->count_child[k] / (double)curr_node->count;
        if(!curr_node->allele_miss[k])
          curr_node->child[k]->node_to->p_node += curr_node->p_node * curr_node->p_child[k];
      }
    }
  }

  return(tree);
}

struct tree *read_bgl(char *filename)
{
  FILE *infile;
  char line[MAX_LINE];
  int i, j, k, prev_level, curr_level, curr_index, next_index, max_index, curr_allele, max_allele, num_level, curr_count, *num_node, *num_allele;
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
  num_level = 1;
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
      node[curr_level][curr_index]->count += curr_count;
      node[curr_level][curr_index]->count_child[curr_allele] += curr_count;
      add_edge(tree, node[curr_level][curr_index], node[curr_level+1][next_index], curr_allele);
    }
  }
  node[num_level-1][0]->count += node[0][0]->count;

  // done with input 
  fclose(infile);

  // calculate edge and node weights
  for(i = 0; i < num_level; i++) {
    for(j = 0; j < num_node[i]; j++) {
      curr_node = node[i][j];
      for(k = 0; k < curr_node->num_allele; k++) {
	curr_node->p_child[k] = (double)curr_node->count_child[k] / (double)curr_node->count;
	if(!curr_node->allele_miss[k])
	  curr_node->child[k]->node_to->p_node += curr_node->p_node * curr_node->p_child[k];
      }
    }
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

  while(fgets(line, MAX_LINE, infile) != NULL) num_line++;
  rewind(infile);

  for(i = 0; i < num_line; i++) {
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
  FILE *infile, *train_out, *test_out;
  char line[MAX_LINE], *temp_str;
  int i, j, temp, num_samp, *r;

  // open the input file
  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return;
  }

  // reads opening line and calculates sample size
  fgets(line, MAX_LINE, infile);
  for(i = 0; line[i] > '0'; i++);
  num_samp = i;
  
  srand(seed);
  r = malloc(num_samp*sizeof(int));
  for(i = 0; i < num_samp; i++)
    r[i] = i;
  for(i = num_samp-1; i >= 0; i--) {
    j = rand() % (i+1);
    temp = r[i];
    r[i] = r[j];
    r[j] = temp;
  }

  for(i = 0; i < num_split; i++) {
    temp_str = malloc((strlen(filename)+32)*sizeof(char));
    sprintf(temp_str, "%s%d.train", filename, i+1);
    train_out = fopen(temp_str, "w");
    if(train_out == NULL) {
      printf("Error: could not open file \"%s\".\n", temp_str);
      return;
    }
    free(temp_str);

    temp_str = malloc((strlen(filename)+32)*sizeof(char));
    sprintf(temp_str, "%s%d.test", filename, i+1);
    test_out = fopen(temp_str, "w");
    if(test_out == NULL) {
      printf("Error: could not open file \"%s\".\n", temp_str);
      return;
    }
    free(temp_str);
    
    rewind(infile);
    while(fgets(line, MAX_LINE, infile) != NULL) {
      for(j = 0; j < num_samp; j++) {
	if(j % num_split == i)
	  fprintf(test_out, "%c", line[r[j]]);
	else
	  fprintf(train_out, "%c", line[r[j]]);
      }
      fprintf(train_out, "\n");
      fprintf(test_out, "\n");
    }

    fclose(train_out);
    fclose(test_out);
  }

  free(r);
  fclose(infile);
}

struct tree_node *new_node(int level, int id, int num_allele)
{
  int i;
  struct tree_node *new = malloc(sizeof(struct tree_node));

  new->level = level;
  new->id = id;
  new->count = 0;
  new->num_allele = num_allele;
  new->num_edge = new->num_parent = 0;
  new->count_child = calloc(num_allele, sizeof(int));
  new->allele_miss = malloc(num_allele*sizeof(int));
  for(i = 0; i < num_allele; i++) new->allele_miss[i] = 1;
  new->p_node = (level==0);
  new->p_child = calloc(num_allele, sizeof(double));
  new->child = malloc(num_allele*sizeof(struct tree_edge*));
  for(i = 0; i < num_allele; i++) new->child[i] = NULL;
  new->parent = NULL;

  new->num_back = new->marked = 0;
  new->back_weight = NULL;

  return new;
}

void add_edge(struct tree *tree, struct tree_node *parent, struct tree_node *child, int allele)
{
  struct tree_edge *new = malloc(sizeof(struct tree_edge));

  new->node_from = parent;
  new->node_to = child;
  new->parent_id = child->num_parent;
  new->allele = allele;
  
  parent->num_edge++;
  parent->child[allele] = new;
  parent->allele_miss[allele] = 0;

  child->num_parent++;
  child->parent = realloc(child->parent, child->num_parent*sizeof(struct tree_edge*));  
  child->parent[child->num_parent-1] = new;
}

void free_tree(struct tree *tree)
{
  int i, j;

  if(tree == NULL) return;
  
  for(i = 0; i < tree->num_level; i++) {
    for(j = 0; j < tree->num_node[i]; j++)
      free_node(tree->node[i][j]);
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
  
  for(i = 0; i < node->num_allele; i++) free(node->child[i]);
  free(node->child);
  free(node->parent);

  if(node->num_back > 0) free(node->back_weight);
  
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
  new->loglik = tree->loglik;
  
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
	  add_edge(new, curr_node, new->node[i+1][tree->node[i][j]->child[k]->node_to->id], k);
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
  new->num_edge = new->num_parent = 0;

  new->count_child = malloc(node->num_allele*sizeof(int));
  memcpy(new->count_child, node->count_child, node->num_allele*sizeof(int));
  new->allele_miss = malloc(node->num_allele*sizeof(int));
  memcpy(new->allele_miss, node->allele_miss, node->num_allele*sizeof(int));

  new->p_node = node->p_node;

  new->p_child = malloc(node->num_allele*sizeof(double));
  memcpy(new->p_child, node->p_child, node->num_allele*sizeof(double));

  new->child = malloc(node->num_allele*sizeof(struct tree_node*));
  for(i = 0; i < node->num_allele; i++) 
    new->child[i] = NULL;
  new->parent = NULL;

  new->num_back = node->num_back;
  new->marked = node->marked;
  new->back_weight = malloc(node->num_back*sizeof(double));
  memcpy(new->back_weight, node->back_weight, node->num_back*sizeof(double));
  
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
      for(k = 0; k < curr_node->num_allele; k++) {
	if(!curr_node->allele_miss[k]) {
	  add_edge(new, curr_node, new->node[i+1][tree->node[i+start-1][j]->child[k]->node_to->id], k);
	}
      }
    }

  while(new->num_node[0] > 1) 
    merge_node(new, new->node[0][0], new->node[0][1], 0, 1);
  while(new->num_node[length] > 1) 
    merge_node(new, new->node[length][0], new->node[length][1], 0, 1);

  return(new);
}

void calc_loglik(struct tree *tree)
{
  int i, j, k;
  struct tree_node *curr_node;

  tree->loglik = 0;
  for(i = 0; i < tree->num_level; i++) {
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      for(k = 0; k < curr_node->num_allele; k++) {
	if(curr_node->p_child[k] > 0)
	  tree->loglik += curr_node->count_child[k] * log(curr_node->p_child[k]);
      }
    }
  }
}

void count_param(struct tree *tree)
{
  int i, j;
  struct tree_node *curr_node;

  tree->num_param_node = tree->num_param_edge = 0;
  for(i = 0; i < tree->num_level; i++) {
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      if(curr_node->num_edge > 0)
	tree->num_param_edge += curr_node->num_edge-1;
    }
    if(tree->num_allele[i] > 0)
      tree->num_param_node += tree->num_node[i] * (tree->num_allele[i]-1);
  }
}

double calc_loss(char *filename, struct tree *tree, double epsilon)
{
  FILE *infile;
  char line[MAX_LINE];
  int i, j, k, curr_level, curr_allele, num_samp, count, total, min_k;
  double add_loss, tot_loss, **loss, *curr_loss;
  struct tree_node *curr_node;

  // open the input file
  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return -1;
  }

  // reads opening line and calculates sample size
  fgets(line, MAX_LINE, infile);
  for(i = 0; line[i] != '\0'; i++);
  num_samp = i;

  loss = malloc(num_samp*sizeof(double*));
  for(i = 0; i < num_samp; i++)
    loss[i] = calloc(1, sizeof(double));
      
  curr_level = 0;
  count = total = 0; // test
  do {
    for(i = 0; line[i] > '0'; i++) {
      curr_allele = line[i] - '1';  // fix to allele key later
      curr_loss = malloc(tree->num_node[curr_level+1]*sizeof(double));
      for(j = 0; j < tree->num_node[curr_level+1]; j++)
	curr_loss[j] = INFINITY;
      for(j = 0; j < tree->num_node[curr_level]; j++) {
	curr_node = tree->node[curr_level][j];
	if(curr_node->allele_miss[curr_allele]) {
	  total++; // test
	  for(k = 0; k < tree->num_node[curr_level+1]; k++) {
	    add_loss = -log(epsilon / tree->num_node[curr_level+1]);
	    //add_loss = -log(EPSILON * tree->node[curr_level+1][k]->p_node);
	    if(loss[i][j] + add_loss < curr_loss[k]) {
	      min_k = k; // test
	      curr_loss[k] = loss[i][j] + add_loss;
	    }
	  }
	  if(curr_node->child[1-curr_allele]->node_to->id == min_k) count++; // test
	}
	else {
	  if(curr_node->allele_miss[1-curr_allele])
	    add_loss = -log(curr_node->p_child[curr_allele]-epsilon);
	  else
	    add_loss = -log(curr_node->p_child[curr_allele]);
	  if(loss[i][j] + add_loss < curr_loss[curr_node->child[curr_allele]->node_to->id])
	    curr_loss[curr_node->child[curr_allele]->node_to->id] = loss[i][j] + add_loss;
	}
      }
      free(loss[i]);
      loss[i] = curr_loss;
    }
    curr_level++;
  } while(fgets(line, MAX_LINE, infile) != NULL);

  tot_loss = 0;
  for(i = 0; i < num_samp; i++) {
    tot_loss += loss[i][0];
    free(loss[i]);
  }
  free(loss);
  
  // done with input
  fclose(infile);

  printf("epsilon = %8.8lf\tsame_prop = %8.8lf\n", epsilon, (double)count/(double)total); // test
  return(tot_loss);
}

/*
double calc_loss(char *filename, struct tree *tree, double *epsilon)
{
  FILE *infile;
  char line[MAX_LINE];
  int i, j, k, curr_level, curr_allele, num_samp, tot_count, **count, *curr_count, total;
  double tot_loss, **loss, *curr_loss;
  struct tree_node *curr_node;

  // open the input file
  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return -1;
  }

  // reads opening line and calculates sample size
  fgets(line, MAX_LINE, infile);
  for(i = 0; line[i] != '\0'; i++);
  num_samp = i;

  count = malloc(num_samp*sizeof(double*));
  loss = malloc(num_samp*sizeof(double*));
  for(i = 0; i < num_samp; i++) {
    count[i] = calloc(1, sizeof(int));
    loss[i] = calloc(1, sizeof(double));
  }
    
  curr_level = 0;
  total = 0;
  do {
    for(i = 0; line[i] > '0'; i++) {
      curr_allele = line[i] - '1';  // fix to allele key later
      curr_count = malloc(tree->num_node[curr_level+1]*sizeof(int));
      curr_loss = malloc(tree->num_node[curr_level+1]*sizeof(double));
      for(j = 0; j < tree->num_node[curr_level+1]; j++)
	curr_count[j] = tree->num_level;
      for(j = 0; j < tree->num_node[curr_level]; j++) {
	curr_node = tree->node[curr_level][j];
	if(curr_node->allele_miss[curr_allele]) {
	  total++;
	  for(k = 0; k < tree->num_node[curr_level+1]; k++) {
	    if(count[i][j] + 1 < curr_count[k]) {
	      curr_count[k] = count[i][j] + 1;
	      curr_loss[k] = loss[i][j] + log(tree->num_node[curr_level+1]);
	    }
	  }
	}
	else {
	  if(curr_node->allele_miss[1-curr_allele])
	    total++;
	  if(count[i][j] < curr_count[curr_node->child[curr_allele]->node_to->id]) {
	    curr_count[curr_node->child[curr_allele]->node_to->id] = count[i][j];
	    curr_loss[curr_node->child[curr_allele]->node_to->id] = loss[i][j] - log(curr_node->p_child[curr_allele]);
	  }
	}
      }
      free(count[i]);
      free(loss[i]);
      count[i] = curr_count;
      loss[i] = curr_loss;
    }
    curr_level++;
  } while(fgets(line, MAX_LINE, infile) != NULL);

  tot_count = 0;
  tot_loss = 0;
  for(i = 0; i < num_samp; i++) {
    tot_count += count[i][0];
    tot_loss += loss[i][0];
    free(count[i]);
    free(loss[i]);
  }
  free(count);
  free(loss);

  *epsilon = (double)tot_count / total;
  tot_loss -= tot_count * log(*epsilon) + (total - tot_count) * log(1 - *epsilon);
    
  // done with input
  fclose(infile);

  return(tot_loss);
}
*/

double search_epsilon(char *filename, struct tree *tree, double tol)
{
  double a, b, c, d, loss_c, loss_d, phi, step;

  b = 0;
  c = -2;
  d = -1;
  phi = (1+sqrt(5))/2;

  loss_c = calc_loss(filename, tree, pow(10, c));
  loss_d = calc_loss(filename, tree, pow(10, d));
  while(loss_c < loss_d) {
    step = d - c;
    d = c;
    loss_d = loss_c;
    c -= 2*step;
    loss_c = calc_loss(filename, tree, pow(10, c));
  }
  a = c;

  c = b - (b-a)/phi;
  d = a + (b-a)/phi;

  loss_c = calc_loss(filename, tree, pow(10, c));
  loss_d = calc_loss(filename, tree, pow(10, d));
  while((b-a)/(fabs(c)+fabs(d)) > tol) {
    if(loss_c < loss_d) {
      b = d;
      d = c;
      loss_d = loss_c;
      c = b - (b-a)/phi;
      loss_c = calc_loss(filename, tree, pow(10, c));
    }
    else {
      a = c;
      c = d;
      loss_c = loss_d;
      d = a + (b-a)/phi;
      loss_d = calc_loss(filename, tree, pow(10, d));
    }
  }

  return(pow(10, (a+b)/2));
}

double merge_test(struct tree_node *node1, struct tree_node *node2, double lambda, double diff, double min_diff, double *tot_diff)
{
  int i;
  double p_child_test, temp_diff = 0;

  if(node1 == node2) return min_diff;

  if(node1->num_allele == 0) return min_diff;
  
  if(min_diff < 0) return min_diff;

  //printf("Level %d: (%d,%d)\n", node1->level, node1->id, node2->id);
  for(i = 0; i < node1->num_allele; i++) {
    p_child_test = (node1->p_child[i]*node1->p_node + node2->p_child[i]*node2->p_node) / (node1->p_node + node2->p_node);
    if(p_child_test > 0)
      temp_diff += (node1->count_child[i] + node2->count_child[i]) * log(p_child_test);
    if(node1->p_child[i] > 0)
      temp_diff -= node1->count_child[i] * log(node1->p_child[i]);
    if(node2->p_child[i] > 0)
      temp_diff -= node2->count_child[i] * log(node2->p_child[i]);
    //printf("\t%d * log(%8.8lf) - %d * log(%8.8lf) - %d * log(%8.8lf)\n", (node1->count_child[i] + node2->count_child[i]), p_child_test, node1->count_child[i], node1->p_child[i], node2->count_child[i], node2->p_child[i]);
  }
  temp_diff += lambda * (node1->num_allele-1);
  //printf("\t%8.8lf * %d\n", lambda, node1->num_allele-1);

  //diff += lambda;
  //for(i = 0; i < node1->num_allele; i++)
  //  diff += lambda * (!node1->allele_miss[i] && !node2->allele_miss[i]);
  //temp_diff += lambda * (node1->num_edge + node2->num_edge - 1);
  //for(i = 0; i < node1->num_allele; i++)
  //  temp_diff -= lambda * (!node1->allele_miss[i] || !node2->allele_miss[i]);
  diff += temp_diff;
  *tot_diff += temp_diff;
  
  if(diff < min_diff) 
    min_diff = diff;

  for(i = 0; i < node1->num_allele; i++)
    if(!node1->allele_miss[i] && !node2->allele_miss[i]) {
      temp_diff = merge_test(node1->child[i]->node_to, node2->child[i]->node_to, lambda, diff, min_diff, tot_diff);
      if(temp_diff < min_diff)
	min_diff = temp_diff;
    }
  
  return min_diff;
} 

struct tree_node *merge_node(struct tree *tree, struct tree_node *node1, struct tree_node *node2, double lambda, int first)
{
  int i, j, k;
  double diff = 0;
  struct tree_node *curr_node, *new;

  if(node1->id == node2->id) return node1;

  if(node1->id > node2->id) {
    curr_node = node1;
    node1 = node2;
    node2 = curr_node;
  }

  new = malloc(sizeof(struct tree_node));
  new->level = node1->level;
  new->id = node1->id;
  new->count = node1->count + node2->count;
  new->num_allele = node1->num_allele;
  new->num_edge = 0;
  new->num_parent = node1->num_parent + node2->num_parent;
  new->p_node = node1->p_node + node2->p_node;
  
  new->count_child = malloc(new->num_allele*sizeof(int));
  new->allele_miss = malloc(new->num_allele*sizeof(int));
  new->p_child = malloc(new->num_allele*sizeof(double));
  new->child = malloc(new->num_allele*sizeof(struct tree_node*));
  for(i = 0; i < new->num_allele; i++) {
    new->count_child[i] = node1->count_child[i] + node2->count_child[i];
    new->allele_miss[i] = node1->allele_miss[i] && node2->allele_miss[i];
    new->p_child[i] = (node1->p_child[i]*node1->p_node + node2->p_child[i]*node2->p_node) / new->p_node;
    if(new->p_child[i] > 0)
      diff += new->count_child[i] * log(new->p_child[i]);
    if(node1->p_child[i] > 0)
      diff -= node1->count_child[i] * log(node1->p_child[i]);
    if(node2->p_child[i] > 0)
      diff -= node2->count_child[i] * log(node2->p_child[i]);
  }

  tree->loglik += diff;
  
  tree->node[node1->level][node1->id] = new;
  for(i = node2->id; i < tree->num_node[node1->level]-1; i++) {
    tree->node[node1->level][i] = tree->node[node1->level][i+1];
    tree->node[node1->level][i]->id--;
  }
  tree->num_node[node1->level]--;
  tree->node[node1->level] = realloc(tree->node[node1->level], tree->num_node[node1->level]*sizeof(struct tree_node*));
  
  new->parent = malloc(new->num_parent*sizeof(struct tree_edge*));
  for(i = 0; i < new->num_parent; i++) {
    if(i < node1->num_parent) 
      new->parent[i] = node1->parent[i];
    else 
      new->parent[i] = node2->parent[i-node1->num_parent];
    new->parent[i]->node_to = new;
    new->parent[i]->parent_id = i;
  }

  new->num_back = node1->num_back;
  new->back_weight = malloc(new->num_back*sizeof(double));
  for(i = 0; i < new->num_back; i++)
    new->back_weight[i] = (node1->back_weight[i]*node1->p_node + node2->back_weight[i]*node2->p_node) / new->p_node; 
  
  for(i = 0; i < new->num_allele; i++) {
    if(new->allele_miss[i]) 
      new->child[i] = NULL;
    else {
      if(!node1->allele_miss[i] && node2->allele_miss[i]) {
	new->num_edge++;
	new->child[i] = node1->child[i];
	new->child[i]->node_from = new;
	node1->child[i]= NULL;
      }
      else if(node1->allele_miss[i] && !node2->allele_miss[i]) {
	new->num_edge++;
	new->child[i] = node2->child[i];
	new->child[i]->node_from = new;
	node2->child[i] = NULL;
      }
      else {
	for(j = node1->child[i]->parent_id; j < node1->child[i]->node_to->num_parent-1; j++) {
	  node1->child[i]->node_to->parent[j] = node1->child[i]->node_to->parent[j+1];
	  node1->child[i]->node_to->parent[j]->parent_id--;
	}
	node1->child[i]->node_to->num_parent--;
	node1->child[i]->node_to->parent = realloc(node1->child[i]->node_to->parent, node1->child[i]->node_to->num_parent*sizeof(struct tree_edge*));
	for(j = node2->child[i]->parent_id; j < node2->child[i]->node_to->num_parent-1; j++) {
	  node2->child[i]->node_to->parent[j] = node2->child[i]->node_to->parent[j+1];
	  node2->child[i]->node_to->parent[j]->parent_id--;
	}
	node2->child[i]->node_to->num_parent--;
	node2->child[i]->node_to->parent = realloc(node2->child[i]->node_to->parent, node2->child[i]->node_to->num_parent*sizeof(struct tree_edge*));
	add_edge(tree, new, merge_node(tree, node1->child[i]->node_to, node2->child[i]->node_to, lambda, 0), i);
      }
    }
  }

  tree->num_param_node--;
  tree->num_param_edge += (new->num_edge > 0 ? new->num_edge-1 : 0) - (node1->num_edge > 0 ? node1->num_edge-1 : 0) - (node2->num_edge > 0 ? node2->num_edge-1 : 0);
    
  free_node(node1);
  free_node(node2);
  
  return new;
}
/*
struct tree *merge_tree(struct tree *tree, double lambda)
{
  int i, j, k, max_j, max_k, merge;
  double min_diff, max_diff, tot_diff, total = 0;
  struct tree_node *curr_node1, *curr_node2;
  FILE *outfile;
  struct timeval start, end;
  
  calc_loglik(tree);
  count_param(tree);
  for(i = 0; i < tree->num_level; i++) {
    //printf("Lambda %8.8lf: Level = %d, Nodes = %d\n", lambda, i, tree->num_node[i]);
    //printf("Merging Level %d\n", i+1);
    gettimeofday(&start, NULL);
    if(tree->num_node[i] > MAX_NODE) return NULL;
    do {
      max_diff = 0;
      merge = 0;
      for(j = 0; j < tree->num_node[i]-1; j++) {
	//if(i == 17) printf("Checking (%d, %d)...\n", i, j+1);
	for(k = j+1; k < tree->num_node[i]; k++) {
	  curr_node1 = tree->node[i][j];
	  curr_node2 = tree->node[i][k];
	  tot_diff = 0;
	  min_diff = merge_test(curr_node1, curr_node2, lambda, 0, lambda*(tree->num_allele[i]-1), &tot_diff); // should be times max of num_allele-1
	  //if(i == 17) printf("\twith (%d, %d): min_diff = %8.8lf, tot_diff = %8.8lf\n", i, k+1, min_diff, tot_diff);
	  if(min_diff >= 0 && tot_diff >= max_diff) {
	    max_diff = tot_diff;
	    max_j = j;
	    max_k = k;
	    merge = 1;
	  }
	}
      }
      if(merge) {
	//if(i == 17) printf("Merging (%d, %d) and (%d, %d).\n", i, max_j+1, i, max_k+1);
	curr_node1 = tree->node[i][max_j];
	curr_node2 = tree->node[i][max_k];
	merge_node(tree, curr_node1, curr_node2, lambda, 1);
      }
    } while(merge);
    gettimeofday(&end, NULL);
    //printf("\tTime = %8.8lf Seconds\n", (end.tv_sec + 0.000001 * end.tv_usec) - (start.tv_sec + 0.000001 * start.tv_usec));
    total += (end.tv_sec + 0.000001 * end.tv_usec) - (start.tv_sec + 0.000001 * start.tv_usec);
  }

  printf("Total Time = %8.8lf Seconds\n", total);
  return tree;
}
*/
struct tree *merge_tree(struct tree *tree, double lambda)
{
  int i, j, k, max_j, max_k, merge;
  double min_diff, max_diff, tot_diff, **diff_mat, total = 0;
  struct tree_node *curr_node1, *curr_node2;
  FILE *outfile;
  struct timeval start, end;
  
  calc_loglik(tree);
  count_param(tree);
  for(i = 0; i < tree->num_level; i++) {
    //printf("Lambda %8.8lf: Level = %d, Nodes = %d\n", lambda, i, tree->num_node[i]);
    printf("\tLevel %d\n", i);
    gettimeofday(&start, NULL);
    if(tree->num_node[i] > MAX_NODE) return NULL;
    diff_mat = malloc((tree->num_node[i]-1)*sizeof(double*));
    for(j = 0; j < tree->num_node[i]-1; j++) {
      diff_mat[j] = malloc((tree->num_node[i]-1-j)*sizeof(double));
      for(k = j+1; k < tree->num_node[i]; k++)
	diff_mat[j][k-1-j] = -1;
    }    
    for(j = 0; j < tree->num_node[i]-1; j++) {
      //if(i == 17) printf("Checking (%d, %d)...\n", i, j+1);
      for(k = j+1; k < tree->num_node[i]; k++) {
	curr_node1 = tree->node[i][j];
	curr_node2 = tree->node[i][k];
	tot_diff = 0;
	min_diff = merge_test(curr_node1, curr_node2, lambda, 0, lambda*(tree->num_allele[i]-1), &tot_diff); // should be times max of num_allele-1 
	//if(i == 17) printf("\twith (%d, %d): min_diff = %8.8lf, tot_diff = %8.8lf\n", i, k+1, min_diff, tot_diff);
	if(min_diff >= 0)
	  diff_mat[j][k-1-j] = tot_diff;
      }
    }
    do {
      max_diff = 0;
      merge = 0;
      for(j = 0; j < tree->num_node[i]-1; j++) {
	for(k = j+1; k < tree->num_node[i]; k++) {
	  if(diff_mat[j][k-1-j] >= max_diff) {
            max_diff = diff_mat[j][k-1-j];
	    max_j = j;
	    max_k = k;
	    merge = 1;
	  }
        }
      }
      if(merge) {
	//printf("nodes = %d. max_j = %d, max_k = %d\n", tree->num_node[i], max_j, max_k);
	//if(i == 17) printf("Merging (%d, %d) and (%d, %d).\n", lambda, i, max_j+1, i, max_k+1);
	curr_node1 = tree->node[i][max_j];
        curr_node2 = tree->node[i][max_k];
        merge_node(tree, curr_node1, curr_node2, lambda, 1);

	if(max_k < tree->num_node[i])
	  free(diff_mat[max_k]);
	for(j = 0; j < tree->num_node[i]-1; j++) {
	  if(j < max_k) {
	    for(k = max_k; k < tree->num_node[i]; k++)
	      diff_mat[j][k-1-j] = diff_mat[j][k-j];
	    //printf("\t\t%d\n", tree->num_node[i]-1-j);
	    diff_mat[j] = realloc(diff_mat[j], (tree->num_node[i]-1-j)*sizeof(double));
	  }
	  else
	    diff_mat[j] = diff_mat[j+1];

	  if(j < max_j) {
	    //if(i == 17) printf("Checking (%d, %d)...\n", i, j+1);
	    curr_node1 = tree->node[i][j];
	    curr_node2 = tree->node[i][max_j];
	    tot_diff = 0;
	    min_diff = merge_test(curr_node1, curr_node2, lambda, 0, lambda*(tree->num_allele[i]-1), &tot_diff); // should be times max of num_allele-1
	    //if(i == 17) printf("\twith (%d, %d): min_diff = %8.8lf, tot_diff = %8.8lf\n", i, max_j+1, min_diff, tot_diff);
	    if(min_diff >= 0)
	      diff_mat[j][max_j-1-j] = tot_diff;
	    else
	      diff_mat[j][max_j-1-j] = -1;
	  }
	  else if(j == max_j) {
	    //if(i == 17) printf("Checking (%d, %d)...\n", i, j+1);
	    for(k = j+1; k < tree->num_node[i]; k++) {
	      curr_node1 = tree->node[i][j];
	      curr_node2 = tree->node[i][k];
	      tot_diff = 0;
	      min_diff = merge_test(curr_node1, curr_node2, lambda, 0, lambda*(tree->num_allele[i]-1), &tot_diff); // should be times max of num_allele-1
	      //if(i == 17) printf("\twith (%d, %d): min_diff = %8.8lf, tot_diff = %8.8lf\n", i, k+1, min_diff, tot_diff);
	      if(min_diff >= 0)
		diff_mat[j][k-1-j] = tot_diff;
	      else
		diff_mat[j][k-1-j] = -1;
	    }
	  }
	}
	//printf("\t%d\n", tree->num_node[i]-1);
	diff_mat = realloc(diff_mat, (tree->num_node[i]-1)*sizeof(double*));
      }
    } while(merge);
    for(j = 0; j < tree->num_node[i]-1; j++)
      free(diff_mat[j]);
    free(diff_mat);
    
    gettimeofday(&end, NULL);
    //printf("\tTime = %8.8lf Seconds\n", (end.tv_sec + 0.000001 * end.tv_usec) - (start.tv_sec + 0.000001 * start.tv_usec));
    total += (end.tv_sec + 0.000001 * end.tv_usec) - (start.tv_sec + 0.000001 * start.tv_usec);
  }

  //printf("Total Time = %8.8lf Seconds\n", total);
  return tree;
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
      fprintf(outfile, "\tNode ID = %d, Node Weight = %8.8lf, Node Count = %d\n", curr_node->id+1, curr_node->p_node, curr_node->count);
      for(k = 0; k < curr_node->num_allele; k++) {
	if(curr_node->allele_miss[k])
	  fprintf(outfile, "\t\tChild Node = %d, Child Allele = %d, Child Weight = %8.8lf, Child Count = %d\n", 0, k+1, curr_node->p_child[k], curr_node->count_child[k]);
	else
	  fprintf(outfile, "\t\tChild Node = %d, Child Allele = %d, Child Weight = %8.8lf, Child Count = %d\n", curr_node->child[k]->node_to->id+1, k+1, curr_node->p_child[k], curr_node->count_child[k]);
      }
    }
  }

  fclose(outfile);
}

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
	  fprintf(outfile, "%d %d %d %3.3lf\n", node_list[i][j], node_list[curr_node->child[k]->node_to->level][curr_node->child[k]->node_to->id], k+1, curr_node->p_child[k]);
    }

  fclose(outfile);
  
  for(i = 0; i < tree->num_level; i++)
    free(node_list[i]);
  free(node_list);
}

void print_pop(char *outname, char *popname, char *chrname, struct tree *tree)
{
  int i, j, k, curr_level, curr_allele, sum, *pop_index, *node_list, ***pop_count;
  char *temp, pop[16], pop_line[1024], chr_line[MAX_LINE];
  struct tree_node *curr_node;
  FILE *pop_in, *chr_in, *outfile;
    
  pop_in = fopen(popname, "r");
  if(pop_in == NULL) {
    printf("Error: could not open file \"%s\".\n", popname);
    return;
  }

  i = 0;
  pop_index = malloc((tree->node[0][0]->count/2)*sizeof(int));
  fgets(pop_line, sizeof(pop_line), pop_in);
  while(fgets(pop_line, sizeof(pop_line), pop_in) != NULL) {
    sscanf(pop_line+12, "%s", pop);
    if(strcmp(pop, "AFR") == 0)
      pop_index[i] = 0;
    else if(strcmp(pop, "EAS") == 0)
      pop_index[i] = 1;
    else if(strcmp(pop, "EUR") == 0)
      pop_index[i] = 2;
    else if(strcmp(pop, "SAS") == 0)
      pop_index[i] = 3;
    else if(strcmp(pop, "AMR") == 0)
      pop_index[i] = 4;
    else
      pop_index[i] = 5;
    i++;
  }

  fclose(pop_in);

  chr_in = fopen(chrname, "r");
  if(chr_in == NULL) {
    printf("Error: could not open file \"%s\".\n", chrname);
    return;
  }

  curr_level = 0;
  pop_count = malloc(tree->num_level*sizeof(int**));
  node_list = calloc(tree->node[0][0]->count, sizeof(int));
  while(fgets(chr_line, MAX_LINE, chr_in) != NULL) {
    pop_count[curr_level] = malloc(tree->num_node[curr_level]*sizeof(int*));
    for(i = 0; i < tree->num_node[curr_level]; i++)
      pop_count[curr_level][i] = calloc(sizeof(int), 6);
    for(i = 0; chr_line[i] > '0'; i++) {
      curr_allele = chr_line[i] - '1';  // fix to allele key later
      pop_count[curr_level][node_list[i]][pop_index[i/2]]++;
      node_list[i] = tree->node[curr_level][node_list[i]]->child[curr_allele]->node_to->id;
    }
    curr_level++;
  }
  pop_count[curr_level] = malloc(tree->num_node[curr_level]*sizeof(int*));
  for(i = 0; i < tree->num_node[curr_level]; i++)
    pop_count[curr_level][i] = calloc(sizeof(int), 6);
  for(i = 0; i < tree->node[0][0]->count; i++) 
    pop_count[curr_level][node_list[i]][pop_index[i/2]]++;
  free(node_list);

  fclose(chr_in);

  temp = malloc((strlen(outname)+7)*sizeof(char));
  strcpy(temp, outname);
  strcat(temp, ".nodes");
  outfile = fopen(temp, "w");
  if(outfile == NULL) {
    printf("Error: could not open file \"%s\".\n", temp);
    return;
  }
  free(temp);

  fprintf(outfile, "ID LEVEL NUMBER COUNT TOT AFR EAS EUR SAS AMR NON\n");

  for(i = 0; i < tree->num_level; i++) {
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      fprintf(outfile, "%d.%d %d %d %3.3lf %d %d %d %d %d %d %d\n", curr_node->level, curr_node->id+1, curr_node->level, curr_node->id+1, curr_node->p_node, curr_node->count, pop_count[i][j][0], pop_count[i][j][1], pop_count[i][j][2], pop_count[i][j][3], pop_count[i][j][4], pop_count[i][j][5]);
    }
  }

  fclose(outfile);

  temp = malloc((strlen(outname)+7)*sizeof(char));
  strcpy(temp, outname);
  strcat(temp, ".edges");
  outfile = fopen(temp, "w");
  if(outfile == NULL) {
    printf("Error: could not open file \"%s\".\n", temp);
    return;
  }
  free(temp);

  fprintf(outfile, "NODE1 NODE2 ALLELE COUNT\n");

  sum = 0;
  for(i = 0; i < tree->num_level-1; i++) {
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      for(k = 0; k < curr_node->num_allele; k++)
	if(!curr_node->allele_miss[k])
	  fprintf(outfile, "%d %d %d %3.3lf\n", sum+curr_node->id, sum+tree->num_node[i]+curr_node->child[k]->node_to->id, k+1, curr_node->p_child[k]);
    }
    sum += tree->num_node[i];
  }

  for(i = 0; i < tree->num_level; i++) {
    for(j = 0; j < tree->num_node[i]; j++)
      free(pop_count[i][j]);
    free(pop_count[i]);
  }
  free(pop_count);
      
  fclose(outfile);
}

int calc_freq(struct haplotype ***results, struct tree* tree, int start, int length, int reverse)
{
  int i, j, k, l, m, index, tot_haplo = 0, **num_haplo, *temp;
  struct tree_node *curr_node;
  struct haplotype ****haplo_tree;

  if(reverse) start = tree->num_level - start - length + 1;

  haplo_tree = malloc((length+1)*sizeof(struct haplotype***));
  num_haplo = malloc((length+1)*sizeof(int*));

  haplo_tree[0] = malloc(tree->num_node[start-1]*sizeof(struct haplotype**));
  num_haplo[0] = malloc(tree->num_node[start-1]*sizeof(int));
  for(i = 0; i < tree->num_node[start-1]; i++) {
    curr_node = tree->node[start-1][i];
    num_haplo[0][i] = 1;
    haplo_tree[0][i] = malloc(sizeof(struct haplotype*));
    haplo_tree[0][i][0] = malloc(sizeof(struct haplotype));
    haplo_tree[0][i][0]->count = 0;
    haplo_tree[0][i][0]->length = 0;
    haplo_tree[0][i][0]->allele = NULL;
    haplo_tree[0][i][0]->freq = log(curr_node->p_node);
  }

  for(i = 0; i < length; i++) {
    printf("\tLevel %d\n", i+start-1);
    haplo_tree[i+1] = calloc(tree->num_node[i+start], sizeof(struct haplotype**));
    num_haplo[i+1] = calloc(tree->num_node[i+start], sizeof(int));
    for(j = 0; j < tree->num_node[i+start-1]; j++) {
      curr_node = tree->node[i+start-1][j];
      printf("\t\tNode %d: %d\n", j+1, num_haplo[i][j]); 
      for(k = 0; k < curr_node->num_allele; k++) {
	if(!curr_node->allele_miss[k]) {
	  haplo_tree[i+1][curr_node->child[k]->node_to->id] = realloc(haplo_tree[i+1][curr_node->child[k]->node_to->id], (num_haplo[i+1][curr_node->child[k]->node_to->id]+num_haplo[i][j])*sizeof(struct haplotype*));
	  for(l = 0; l < num_haplo[i][j]; l++) {
	    index = num_haplo[i+1][curr_node->child[k]->node_to->id] + l;
	    haplo_tree[i+1][curr_node->child[k]->node_to->id][index] = malloc(sizeof(struct haplotype));
	    if(i == length-1)
	      haplo_tree[i+1][curr_node->child[k]->node_to->id][index]->count = curr_node->count_child[k];
	    haplo_tree[i+1][curr_node->child[k]->node_to->id][index]->length = haplo_tree[i][j][l]->length + 1;
	    haplo_tree[i+1][curr_node->child[k]->node_to->id][index]->allele = malloc(haplo_tree[i+1][curr_node->child[k]->node_to->id][index]->length*sizeof(int));
	    memcpy(haplo_tree[i+1][curr_node->child[k]->node_to->id][index]->allele, haplo_tree[i][j][l]->allele, haplo_tree[i][j][l]->length*sizeof(int));
	    haplo_tree[i+1][curr_node->child[k]->node_to->id][index]->allele[haplo_tree[i][j][l]->length] = k+1;  /* allele fix */
	    haplo_tree[i+1][curr_node->child[k]->node_to->id][index]->freq = haplo_tree[i][j][l]->freq + log(curr_node->p_child[k]);
	  }
	  num_haplo[i+1][curr_node->child[k]->node_to->id] += num_haplo[i][j];
	}
      }
    }
  }
  /* free old haplos too */

  /* free old results? */
  *results = NULL;
  for(i = 0; i < tree->num_node[start+length-1]; i++) {
    *results = realloc(*results, (tot_haplo+num_haplo[length][i])*sizeof(struct haplotype*));
    for(j = 0; j < num_haplo[length][i]; j++)
      (*results)[tot_haplo+j] = haplo_tree[length][i][j];
    tot_haplo += num_haplo[length][i];
  }

  qsort(*results, tot_haplo, sizeof(struct haplotype*), compare_haplo);

  j = 0;
  for(i = 1; i < tot_haplo; i++) {
    if(compare_haplo((*results)+j, (*results)+i) == 0) {
      (*results)[j]->count += (*results)[i]->count;
      if((*results)[i]->freq == -1.0/0.0)
	continue;
      if((*results)[j]->freq == -1.0/0.0)
	(*results)[j]->freq = (*results)[i]->freq;
      else
	(*results)[j]->freq += log(1 + exp((*results)[i]->freq - (*results)[j]->freq));
    }
    else (*results)[++j] = (*results)[i];
  }
  tot_haplo = j + 1;
  *results = realloc(*results, tot_haplo*sizeof(struct haplotype*));
  
  if(reverse) {
    for(i = 0; i < tot_haplo; i++) {
      temp = malloc((*results)[i]->length*sizeof(int));
      memcpy(temp, (*results)[i]->allele, (*results)[i]->length*sizeof(int));
      for(j = 0; j < (*results)[i]->length; j++)
	(*results)[i]->allele[j] = temp[(*results)[i]->length-j-1];
      free(temp);
    }
    qsort(*results, tot_haplo, sizeof(struct haplotype*), compare_haplo);
  }

  return tot_haplo;
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

int compare_haplo(const void* a, const void* b)
{
  int i, n = (*(struct haplotype**)a)->length;
  
  for(i = 0; i < (*(struct haplotype**)a)->length; i++)
    if((*(struct haplotype**)a)->allele[i] != (*(struct haplotype**)b)->allele[i])
      return (*(struct haplotype**)a)->allele[i] - (*(struct haplotype**)b)->allele[i];
  
  return 0;
}

void calc_tree(struct tree *tree1, struct tree *tree2)
{
  int i, j, k, curr_level, curr_id, back_level, back_id, check;
  double child_sum, back_mult;
  struct tree_node *curr_node, *back_node;

  curr_node = tree1->node[0][0];
  if(curr_node->num_back != 1) {
    curr_node->num_back = 1;
    curr_node->back_weight = malloc(sizeof(double));
    curr_node->back_weight[0] = 1;
  }

  // Maybe back to tree1->num_level?
  for(curr_level = 0; curr_level < tree1->num_level-1; curr_level++) {
    back_level = tree2->num_level - curr_level - 2;
    for(curr_id = 0; curr_id < tree1->num_node[curr_level]; curr_id++) {
      curr_node = tree1->node[curr_level][curr_id];
      curr_node->marked = 0;
      child_sum = 0;
      for(i = 0; i < curr_node->num_allele; i++) {
	if(!curr_node->allele_miss[i] && !curr_node->child[i]->node_to->marked) {
	  curr_node->child[i]->node_to->marked = 1;
	  curr_node->child[i]->node_to->p_node = 0;
	  if(curr_node->child[i]->node_to->num_back != tree2->num_node[back_level]) {
	    curr_node->child[i]->node_to->num_back = tree2->num_node[back_level];
	    free(curr_node->child[i]->node_to->back_weight);
	    curr_node->child[i]->node_to->back_weight = calloc(curr_node->child[i]->node_to->num_back, sizeof(double));
	  }
	  else
	    for(j = 0; j < tree2->num_node[back_level]; j++)
	      curr_node->child[i]->node_to->back_weight[j] = 0;
	}
	curr_node->p_child[i] = 0;
	for(j = 0; j < tree2->num_node[back_level]; j++) {
	  back_node = tree2->node[back_level][j];
	  if(back_node->allele_miss[i])
	    back_mult = 0;
	  else
	    back_mult = curr_node->back_weight[back_node->child[i]->node_to->id];
	  if(curr_node->p_node > 0)
	    curr_node->p_child[i] += back_mult * back_node->p_child[i] * back_node->p_node / curr_node->p_node;
	  if(curr_node->allele_miss[i] && curr_node->p_child[i] > 0) {
	    if(curr_level+1 == tree1->num_level-1)
	      add_edge(tree1, curr_node, tree1->node[tree1->num_level-1][0], i);
	    else {
	      tree1->num_node[curr_level+1]++;
	      tree1->node[curr_level+1] = realloc(tree1->node[curr_level+1], tree1->num_node[curr_level+1]*sizeof(struct tree_node*));
	      tree1->node[curr_level+1][tree1->num_node[curr_level+1]-1] = new_node(curr_level+1, tree1->num_node[curr_level+1]-1, tree1->num_allele[curr_level+1]);
	      add_edge(tree1, curr_node, tree1->node[curr_level+1][tree1->num_node[curr_level+1]-1], i);
	      curr_node->child[i]->node_to->marked = 1;
	      curr_node->child[i]->node_to->p_node = 0;
	      curr_node->child[i]->node_to->num_back = tree2->num_node[back_level];
	      curr_node->child[i]->node_to->back_weight = calloc(curr_node->child[i]->node_to->num_back, sizeof(double));
	    }
	  }
	  if(!curr_node->allele_miss[i])
	    curr_node->child[i]->node_to->back_weight[back_node->id] += back_mult * back_node->p_child[i];
	}
	child_sum += curr_node->p_child[i];
      }
      for(i = 0; i < curr_node->num_allele; i++) {
	if(child_sum > 0)
	  curr_node->p_child[i] /= child_sum;
	if(!curr_node->allele_miss[i])
	  curr_node->child[i]->node_to->p_node += curr_node->p_node * curr_node->p_child[i];
      }
    }
  }
}

void random_phase(char *infilename, char *outfilename, int seed)
{
  FILE *infile, *outfile;
  char c, line[MAX_LINE];
  int i;

  // open the input file
  infile = fopen(infilename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", infilename);
    return;
  }

  // open the output file
  outfile = fopen(outfilename, "w");
  if(outfile == NULL) {
    printf("Error: could not open file \"%s\".\n", outfilename);
    return;
  }
  
  srand(seed);
  while(fgets(line, MAX_LINE, infile) != NULL) {
    for(i = 0; line[i] > '0' && line[i+1] > '0'; i += 2) {
      if(rand() % 2) {
	c = line[i];
	line[i] = line[i+1];
	line[i+1] = c;
      }
    }
    fprintf(outfile, "%s", line);
  }

  fclose(infile);
  fclose(outfile);
}

/*
double sample_phase(char *filename, struct tree *tree, int seed)
{
  FILE *infile, *outfile;
  char *curr_line;
  int i, j, k, n, max_line, curr_level, curr_allele1, curr_allele2, *num_edge, *edge_list1, *edge_list2;
  double log_sum, cum_sum, loglik, threshold, ****log_alpha;
  struct tree_edge *curr_edge1, *curr_edge2, *prev_edge1, *prev_edge2, ***edge;

  // open the input file
  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return 0;
  }
  
  // create edge structure
  edge = malloc((tree->num_level-1)*sizeof(struct tree_edge**));
  num_edge = calloc(tree->num_level-1, sizeof(int));
  for(curr_level = 0; curr_level < tree->num_level-1; curr_level++) {
    for(i = 0; i < tree->num_node[curr_level]; i++)
      num_edge[curr_level] += tree->node[curr_level][i]->num_edge;
    edge[curr_level] = malloc(num_edge[curr_level]*sizeof(struct tree_edge*));
    for(i = 0, k = 0; i < tree->num_node[curr_level]; i++)
      for(j = 0; j < tree->num_allele[curr_level]; j++)
	if(!tree->node[curr_level][i]->allele_miss[j]) {
	  tree->node[curr_level][i]->child[j]->id = k;
	  edge[curr_level][k++] = tree->node[curr_level][i]->child[j];
	}
  }

  // create log_alpha
  n = tree->node[0][0]->count / 2;
  log_alpha = malloc((tree->num_level-1)*sizeof(double***));
  for(curr_level = 0; curr_level < tree->num_level-1; curr_level++) {
    log_alpha[curr_level] = malloc(num_edge[curr_level]*sizeof(double**));
    for(i = 0; i < num_edge[curr_level]; i++) {
      log_alpha[curr_level][i] = malloc(num_edge[curr_level]*sizeof(double*));
      for(j = 0; j < num_edge[curr_level]; j++) {
	log_alpha[curr_level][i][j] = malloc(n*sizeof(double));
	for(k = 0; k < n; k++)
	  log_alpha[curr_level][i][j][k] = -INFINITY;
      }
    }
  }

  // initialize log_alpha
  curr_level = 0;
  max_line = 2*n + 5;
  curr_line = malloc(max_line*sizeof(char));
  fgets(curr_line, max_line, infile);
  for(k = 0; k < n; k++) {
    curr_allele1 = curr_line[2*k] - '1';
    curr_allele2 = curr_line[2*k+1] - '1';
    for(i = 0; i < num_edge[curr_level]; i++) {
      curr_edge1 = edge[curr_level][i];
      for(j = 0; j < num_edge[curr_level]; j++) {
	curr_edge2 = edge[curr_level][j];
	if((curr_allele1 == curr_edge1->allele && curr_allele2 == curr_edge2->allele) || (curr_allele2 == curr_edge1->allele && curr_allele1 == curr_edge2->allele))
	  log_alpha[0][i][j][k] = log(curr_edge1->node_from->p_child[curr_edge1->allele]) + log(curr_edge2->node_from->p_child[curr_edge2->allele]);
      }
    }
  }
  curr_level++;

  // recursively compute log_alpha (add forward instead?)
  while(fgets(curr_line, max_line, infile) != NULL) {
    for(k = 0; k < n; k++) {
      curr_allele1 = curr_line[2*k] - '1';
      curr_allele2 = curr_line[2*k+1] - '1';
      for(i = 0; i < num_edge[curr_level-1]; i++) {
	prev_edge1 = edge[curr_level-1][i];
	for(j = 0; j < num_edge[curr_level-1]; j++) {
	  prev_edge2 = edge[curr_level-1][j];
	  if(!prev_edge1->node_to->allele_miss[curr_allele1] && !prev_edge2->node_to->allele_miss[curr_allele2]) {
	    curr_edge1 = prev_edge1->node_to->child[curr_allele1];
	    curr_edge2 = prev_edge2->node_to->child[curr_allele2];
	    log_sum = log_alpha[curr_level-1][i][j][k] + log(curr_edge1->node_from->p_child[curr_edge1->allele]) + log(curr_edge2->node_from->p_child[curr_edge2->allele]);
	    if(log_alpha[curr_level][curr_edge1->id][curr_edge2->id][k] == -INFINITY)
	      log_alpha[curr_level][curr_edge1->id][curr_edge2->id][k] = log_sum;
	    else
	      log_alpha[curr_level][curr_edge1->id][curr_edge2->id][k] += log1p(exp(log_sum - log_alpha[curr_level][curr_edge1->id][curr_edge2->id][k]));
	  }	    
	  if(curr_allele1 != curr_allele2) {
	    if(!prev_edge1->node_to->allele_miss[curr_allele2] && !prev_edge2->node_to->allele_miss[curr_allele1]) {
	      curr_edge1 = prev_edge1->node_to->child[curr_allele2];
	      curr_edge2 = prev_edge2->node_to->child[curr_allele1];
	      log_sum = log_alpha[curr_level-1][i][j][k] + log(curr_edge1->node_from->p_child[curr_edge1->allele]) + log(curr_edge2->node_from->p_child[curr_edge2->allele]);
	      if(log_alpha[curr_level][curr_edge1->id][curr_edge2->id][k] == -INFINITY)
		log_alpha[curr_level][curr_edge1->id][curr_edge2->id][k] = log_sum;
	      else
		log_alpha[curr_level][curr_edge1->id][curr_edge2->id][k] += log1p(exp(log_sum - log_alpha[curr_level][curr_edge1->id][curr_edge2->id][k]));
	    }
	  }
	}
      }
    }
    curr_level++;
  }

  for(curr_level = 0; curr_level < tree->num_level-1; curr_level++)
    for(i = 0; i < num_edge[curr_level]; i++)
      for(j = 0; j < num_edge[curr_level]; j++)
	printf("(%d, %d, %d) = %8.8lf\n", curr_level, i, j, log_alpha[curr_level][i][j][0]);
    
  fclose(infile);
  infile = fopen(reverse_input(filename), "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s.rev\".\n", filename);
    return 0;
  }

  outfile = fopen(filename, "w");
  if(outfile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return 0;
  }

  srand(seed);
  curr_level = tree->num_level-2;
  edge_list1 = malloc(n*sizeof(int));
  edge_list2 = malloc(n*sizeof(int));

  // sample last marker
  fgets(curr_line, max_line, infile);
  loglik = 0;
  for(k = 0; k < n; k++) {
    curr_allele1 = curr_line[2*k] - '1';
    curr_allele2 = curr_line[2*k+1] - '1';
    cum_sum = -INFINITY;
    for(i = 0; i < num_edge[curr_level]; i++) 
      for(j = 0; j < num_edge[curr_level]; j++) {
	log_sum = log_alpha[curr_level][i][j][k];
	if(cum_sum == -INFINITY)
	  cum_sum = log_sum;
	else
	  cum_sum += log1p(exp(log_sum - cum_sum));
      }
    loglik += cum_sum;
    
    threshold = log((double)rand() / ((double)RAND_MAX + 1)) + cum_sum;
    cum_sum = -INFINITY;
    for(i = 0; i < num_edge[curr_level] && cum_sum <= threshold; i++)
      for(j = 0; j < num_edge[curr_level] && cum_sum <= threshold; j++) {
	log_sum = log_alpha[curr_level][i][j][k];
	if(cum_sum == -INFINITY)
	  cum_sum = log_sum;
	else
	  cum_sum += log1p(exp(log_sum - cum_sum));
      }
        
    curr_line[2*k] = edge[curr_level][--i]->allele + '1';
    curr_line[2*k+1] = edge[curr_level][--j]->allele + '1';
    edge_list1[k] = i;
    edge_list2[k] = j;
  }
  fprintf(outfile, "%s", curr_line);
  curr_level--;

  // recursively sample markers
  while(fgets(curr_line, max_line, infile) != NULL) {
    for(k = 0; k < n; k++) {
      curr_allele1 = curr_line[2*k] - '1';
      curr_allele2 = curr_line[2*k+1] - '1';
      prev_edge1 = edge[curr_level+1][edge_list1[k]];
      prev_edge2 = edge[curr_level+1][edge_list2[k]];
      
      threshold = log((double)rand() / ((double)RAND_MAX + 1));
      cum_sum = -INFINITY;
      for(i = 0; i < num_edge[curr_level] && cum_sum <= threshold; i++) {
	curr_edge1 = edge[curr_level][i];
	for(j = 0; j < num_edge[curr_level] && cum_sum <= threshold; j++) {
	  curr_edge2 = edge[curr_level][j];
	  if(prev_edge1 == curr_edge1->node_to->child[prev_edge1->allele] && prev_edge2 == curr_edge2->node_to->child[prev_edge2->allele]) {	  
	    log_sum = log(prev_edge1->node_from->p_child[prev_edge1->allele]) + log(prev_edge2->node_from->p_child[prev_edge2->allele]) + log_alpha[curr_level][i][j][k] - log_alpha[curr_level+1][edge_list1[k]][edge_list2[k]][k];
	    if(cum_sum == -INFINITY)
	      cum_sum = log_sum;
	    else
	      cum_sum += log1p(exp(log_sum - cum_sum));
	  }
	}
      }
      
      curr_line[2*k] = edge[curr_level][--i]->allele + '1';
      curr_line[2*k+1] = edge[curr_level][--j]->allele + '1';
      edge_list1[k] = i;
      edge_list2[k] = j;
    }
    fprintf(outfile, "%s", curr_line);
    curr_level--;
  } 

  free(curr_line);
  free(edge_list1);
  free(edge_list2);

  for(curr_level = 0; curr_level < tree->num_level-1; curr_level++) {
    for(i = 0; i < num_edge[curr_level]; i++) {
      for(j = 0; j < num_edge[curr_level]; j++)
	free(log_alpha[curr_level][i][j]);
      free(log_alpha[curr_level][i]);
    }
    free(log_alpha[curr_level]);
    free(edge[curr_level]);
  }
  free(num_edge);
  free(log_alpha);
  free(edge);

  fclose(infile);
  fclose(outfile);

  return loglik;
}
*/

double sample_phase(char *filename, struct tree *tree, int seed)
{
  FILE *infile, *outfile;
  char *curr_line, **line_list;
  int i, j, k, n, size, max_line, curr_level, curr_allele1, curr_allele2, *num_edge;
  long int *cum_edge;
  double log_sum, cum_sum, loglik, threshold;
  struct tree_edge *curr_edge1, *curr_edge2, *prev_edge1, *prev_edge2, ***edge;
  struct hashtable *log_alpha;
  
  // open the input file
  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return 0;
  }

  // read the input file
  n = tree->node[0][0]->count / 2;
  max_line = 2*n + 4;
  line_list = malloc((tree->num_level-1)*sizeof(char*));
  for(curr_level = 0; curr_level < tree->num_level-1; curr_level++) {
    line_list[curr_level] = malloc(max_line*sizeof(char));
    fgets(line_list[curr_level], max_line, infile);
  }

  // create edge structure
  edge = malloc((tree->num_level-1)*sizeof(struct tree_edge**));
  num_edge = calloc(tree->num_level-1, sizeof(int));
  for(curr_level = 0; curr_level < tree->num_level-1; curr_level++) {
    for(i = 0; i < tree->num_node[curr_level]; i++)
      num_edge[curr_level] += tree->node[curr_level][i]->num_edge;
    edge[curr_level] = malloc(num_edge[curr_level]*sizeof(struct tree_edge*));
    for(i = 0, k = 0; i < tree->num_node[curr_level]; i++)
      for(j = 0; j < tree->num_allele[curr_level]; j++)
	if(!tree->node[curr_level][i]->allele_miss[j]) {
	  tree->node[curr_level][i]->child[j]->id = k;
	  edge[curr_level][k++] = tree->node[curr_level][i]->child[j];
	}
  }
  cum_edge = malloc((tree->num_level-1)*sizeof(long int));
  cum_edge[0] = 0;
  for(curr_level = 1; curr_level < tree->num_level-1; curr_level++)
    cum_edge[curr_level] = cum_edge[curr_level-1] + (long int)(num_edge[curr_level-1]*num_edge[curr_level-1]);

  size = 0;
  for(curr_level = 0; curr_level < tree->num_level-1; curr_level++)
    size += 2*num_edge[curr_level];
  
  srand(seed);
  for(k = 0; k < n; k++) {
    // create log_alpha
    log_alpha = new_hashtable(size);

    // initialize log_alpha
    curr_level = 0;
    max_line = 2*n + 5;
    curr_line = line_list[curr_level];
    curr_allele1 = curr_line[2*k] - '1';
    curr_allele2 = curr_line[2*k+1] - '1';
    for(i = 0; i < num_edge[curr_level]; i++) {
      curr_edge1 = edge[curr_level][i];
      for(j = 0; j < num_edge[curr_level]; j++) {
	curr_edge2 = edge[curr_level][j];
	if((curr_allele1 == curr_edge1->allele && curr_allele2 == curr_edge2->allele) || (curr_allele2 == curr_edge1->allele && curr_allele1 == curr_edge2->allele))
	  insert_hashtable_logsum(log_alpha, edge_pair_key(0, i, j, num_edge, cum_edge), log(curr_edge1->node_from->p_child[curr_edge1->allele]) + log(curr_edge2->node_from->p_child[curr_edge2->allele]));
      }
    }
    curr_level++;

    // recursively compute log_alpha
    for(; curr_level < tree->num_level-1; curr_level++) {
      printf("\t%d - %d\n", k, curr_level);
      curr_line = line_list[curr_level];
      curr_allele1 = curr_line[2*k] - '1';
      curr_allele2 = curr_line[2*k+1] - '1';
      for(i = 0; i < num_edge[curr_level-1]; i++) {
	prev_edge1 = edge[curr_level-1][i];
	for(j = 0; j < num_edge[curr_level-1]; j++) {
	  prev_edge2 = edge[curr_level-1][j];
	  if(!prev_edge1->node_to->allele_miss[curr_allele1] && !prev_edge2->node_to->allele_miss[curr_allele2]) {
	    curr_edge1 = prev_edge1->node_to->child[curr_allele1];
	    curr_edge2 = prev_edge2->node_to->child[curr_allele2];
	    log_sum = search_hashtable_prob(log_alpha, edge_pair_key(curr_level-1, i, j, num_edge, cum_edge)) + log(curr_edge1->node_from->p_child[curr_edge1->allele]) + log(curr_edge2->node_from->p_child[curr_edge2->allele]);
	    insert_hashtable_logsum(log_alpha, edge_pair_key(curr_level, curr_edge1->id, curr_edge2->id, num_edge, cum_edge), log_sum);
	  }
	  if(curr_allele1 != curr_allele2) {
	    if(!prev_edge1->node_to->allele_miss[curr_allele2] && !prev_edge2->node_to->allele_miss[curr_allele1]) {
	      curr_edge1 = prev_edge1->node_to->child[curr_allele2];
	      curr_edge2 = prev_edge2->node_to->child[curr_allele1];
	      log_sum = search_hashtable_prob(log_alpha, edge_pair_key(curr_level-1, i, j, num_edge, cum_edge)) + log(curr_edge1->node_from->p_child[curr_edge1->allele]) + log(curr_edge2->node_from->p_child[curr_edge2->allele]);
	      insert_hashtable_logsum(log_alpha, edge_pair_key(curr_level, curr_edge1->id, curr_edge2->id, num_edge, cum_edge), log_sum);
	    }
	  }
	}
      }
    }

    // sample last marker
    curr_level = tree->num_level-2;
    curr_line = line_list[curr_level];
    curr_allele1 = curr_line[2*k] - '1';
    curr_allele2 = curr_line[2*k+1] - '1';
    cum_sum = -INFINITY;
    for(i = 0; i < num_edge[curr_level]; i++)
      for(j = 0; j < num_edge[curr_level]; j++) {
	log_sum = search_hashtable_prob(log_alpha, edge_pair_key(curr_level, i, j, num_edge, cum_edge));
	if(cum_sum == -INFINITY)
	  cum_sum = log_sum;
	else
	  cum_sum += log1p(exp(log_sum - cum_sum));
      }
    loglik += cum_sum;
    
    threshold = log((double)rand() / ((double)RAND_MAX + 1)) + cum_sum;
    cum_sum = -INFINITY;
    for(i = 0; i < num_edge[curr_level] && cum_sum <= threshold; i++)
      for(j = 0; j < num_edge[curr_level] && cum_sum <= threshold; j++) {
	log_sum = search_hashtable_prob(log_alpha, edge_pair_key(curr_level, i, j, num_edge, cum_edge));
	if(cum_sum == -INFINITY)
	  cum_sum = log_sum;
	else
	  cum_sum += log1p(exp(log_sum - cum_sum));
      }

    curr_line[2*k] = edge[curr_level][--i]->allele + '1';
    curr_line[2*k+1] = edge[curr_level][--j]->allele + '1';
    prev_edge1 = edge[curr_level][i];
    prev_edge2 = edge[curr_level][j];
    curr_level--;

    // recursively sample markers 
    for(; curr_level >= 0; curr_level--) {
      curr_line = line_list[curr_level];
      curr_allele1 = curr_line[2*k] - '1';
      curr_allele2 = curr_line[2*k+1] - '1';
      
      threshold = log((double)rand() / ((double)RAND_MAX + 1));
      cum_sum = -INFINITY;
      for(i = 0; i < num_edge[curr_level] && cum_sum <= threshold; i++) {
	curr_edge1 = edge[curr_level][i];
	for(j = 0; j < num_edge[curr_level] && cum_sum <= threshold; j++) {
	  curr_edge2 = edge[curr_level][j];
	  if(prev_edge1 == curr_edge1->node_to->child[prev_edge1->allele] && prev_edge2 == curr_edge2->node_to->child[prev_edge2->allele]) {
	    log_sum = log(prev_edge1->node_from->p_child[prev_edge1->allele]) + log(prev_edge2->node_from->p_child[prev_edge2->allele])
	      + search_hashtable_prob(log_alpha, edge_pair_key(curr_level, i, j, num_edge, cum_edge)) - search_hashtable_prob(log_alpha, edge_pair_key(curr_level+1, prev_edge1->id, prev_edge2->id, num_edge, cum_edge));
	    if(cum_sum == -INFINITY)
	      cum_sum = log_sum;
	    else
	      cum_sum += log1p(exp(log_sum - cum_sum));
	  }
	}
      }
      
      curr_line[2*k] = edge[curr_level][--i]->allele + '1';
      curr_line[2*k+1] = edge[curr_level][--j]->allele + '1';
      prev_edge1 = edge[curr_level][i];
      prev_edge2 = edge[curr_level][j];
    }

    // printf("k = %d: Used %d of %d.\n", k, log_alpha->filled, log_alpha->size);
    free_hashtable(log_alpha);
  }

  // open the output file
  outfile = fopen(filename, "w");
  if(outfile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return 0;
  }

  for(curr_level = tree->num_level-2; curr_level >= 0; curr_level--) {
    fprintf(outfile, "%s", line_list[curr_level]);
    free(line_list[curr_level]);
    free(edge[curr_level]);
  }
  free(line_list);
  free(edge);
  free(num_edge);
  fclose(outfile);
  
  return loglik;
}

void viterbi_phase(char *filename, struct tree *tree)
{
  FILE *infile, *outfile;
  char *curr_line, **line_list;
  int i, j, k, n, size, max_line, curr_level, curr_allele1, curr_allele2, max_i, max_j, *num_edge;
  long int *cum_edge;
  double log_sum;
  struct tree_edge *curr_edge1, *curr_edge2, *prev_edge1, *prev_edge2, ***edge;
  struct hashtable *log_alpha;
  struct hashtable_item *curr_item, *max_item;

  // open the input file
  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return;
  }

  // read the input file
  n = tree->node[0][0]->count / 2;
  max_line = 2*n + 4;
  line_list = malloc((tree->num_level-1)*sizeof(char*));
  for(curr_level = 0; curr_level < tree->num_level-1; curr_level++) {
    line_list[curr_level] = malloc(max_line*sizeof(char));
    fgets(line_list[curr_level], max_line, infile);
  }

  // create edge structure
  edge = malloc((tree->num_level-1)*sizeof(struct tree_edge**));
  num_edge = calloc(tree->num_level-1, sizeof(int));
  for(curr_level = 0; curr_level < tree->num_level-1; curr_level++) {
    for(i = 0; i < tree->num_node[curr_level]; i++)
      num_edge[curr_level] += tree->node[curr_level][i]->num_edge;
    edge[curr_level] = malloc(num_edge[curr_level]*sizeof(struct tree_edge*));
    for(i = 0, k = 0; i < tree->num_node[curr_level]; i++)
      for(j = 0; j < tree->num_allele[curr_level]; j++)
	if(!tree->node[curr_level][i]->allele_miss[j]) {
	  tree->node[curr_level][i]->child[j]->id = k;
	  edge[curr_level][k++] = tree->node[curr_level][i]->child[j];
	}
  }
  cum_edge = malloc((tree->num_level-1)*sizeof(long int));
  cum_edge[0] = 0;
  for(curr_level = 1; curr_level < tree->num_level-1; curr_level++)
    cum_edge[curr_level] = cum_edge[curr_level-1] + (long int)(num_edge[curr_level-1]*num_edge[curr_level-1]);

  size = 0;
  for(curr_level = 0; curr_level < tree->num_level-1; curr_level++)
    size += 2*num_edge[curr_level];

  for(k = 0; k < n; k++) {
    // create log_alpha
    log_alpha = new_hashtable(size);

    // initialize log_alpha
    curr_level = 0;
    max_line = 2*n + 5;
    curr_line = line_list[curr_level];
    curr_allele1 = curr_line[2*k] - '1';
    curr_allele2 = curr_line[2*k+1] - '1';
    for(i = 0; i < num_edge[curr_level]; i++) {
      curr_edge1 = edge[curr_level][i];
      for(j = 0; j < num_edge[curr_level]; j++) {
	curr_edge2 = edge[curr_level][j];
	if((curr_allele1 == curr_edge1->allele && curr_allele2 == curr_edge2->allele) || (curr_allele2 == curr_edge1->allele && curr_allele1 == curr_edge2->allele))
	  insert_hashtable_max(log_alpha, edge_pair_key(0, i, j, num_edge, cum_edge), log(curr_edge1->node_from->p_child[curr_edge1->allele]) + log(curr_edge2->node_from->p_child[curr_edge2->allele]), -1, -1);
      }
    }
    curr_level++;

    // recursively compute log_alpha
    for(; curr_level < tree->num_level-1; curr_level++) {
      curr_line = line_list[curr_level];
      curr_allele1 = curr_line[2*k] - '1';
      curr_allele2 = curr_line[2*k+1] - '1';
      for(i = 0; i < num_edge[curr_level-1]; i++) {
	prev_edge1 = edge[curr_level-1][i];
	for(j = 0; j < num_edge[curr_level-1]; j++) {
	  prev_edge2 = edge[curr_level-1][j];
	  if(!prev_edge1->node_to->allele_miss[curr_allele1] && !prev_edge2->node_to->allele_miss[curr_allele2]) {
	    curr_edge1 = prev_edge1->node_to->child[curr_allele1];
	    curr_edge2 = prev_edge2->node_to->child[curr_allele2];
	    log_sum = search_hashtable_prob(log_alpha, edge_pair_key(curr_level-1, i, j, num_edge, cum_edge)) + log(curr_edge1->node_from->p_child[curr_edge1->allele]) + log(curr_edge2->node_from->p_child[curr_edge2->allele]);
	    insert_hashtable_max(log_alpha, edge_pair_key(curr_level, curr_edge1->id, curr_edge2->id, num_edge, cum_edge), log_sum, i, j);
	  }
	  if(curr_allele1 != curr_allele2) {
	    if(!prev_edge1->node_to->allele_miss[curr_allele2] && !prev_edge2->node_to->allele_miss[curr_allele1]) {
	      curr_edge1 = prev_edge1->node_to->child[curr_allele2];
	      curr_edge2 = prev_edge2->node_to->child[curr_allele1];
	      log_sum = search_hashtable_prob(log_alpha, edge_pair_key(curr_level-1, i, j, num_edge, cum_edge)) + log(curr_edge1->node_from->p_child[curr_edge1->allele]) + log(curr_edge2->node_from->p_child[curr_edge2->allele]);
	      insert_hashtable_max(log_alpha, edge_pair_key(curr_level, curr_edge1->id, curr_edge2->id, num_edge, cum_edge), log_sum, i, j);
	    }
	  }
	}
      }
    }

    // backtrack to find max path
    curr_level = tree->num_level-2;
    max_item = NULL;
    for(i = 0; i < num_edge[curr_level]; i++)
      for(j = 0; j < num_edge[curr_level]; j++) {
	curr_item = search_hashtable(log_alpha, edge_pair_key(curr_level, i, j, num_edge, cum_edge));
	if(curr_item != NULL && (max_item == NULL || curr_item->prob > max_item->prob)) {
	  max_item = curr_item;
	  max_i = i;
	  max_j = j;
	}
      }

    for(; curr_level > 0; curr_level--) {
      curr_line = line_list[curr_level];
      curr_line[2*k] = edge[curr_level][max_i]->allele + '1';
      curr_line[2*k+1] = edge[curr_level][max_j]->allele + '1';
      max_i = max_item->back1;
      max_j = max_item->back2;
      max_item = search_hashtable(log_alpha, edge_pair_key(curr_level-1, max_i, max_j, num_edge, cum_edge));
    }
    curr_line = line_list[curr_level];
    curr_line[2*k] = edge[curr_level][max_i]->allele + '1';
    curr_line[2*k+1] = edge[curr_level][max_j]->allele + '1';

    free_hashtable(log_alpha);
  }

  // open the output file
  outfile = fopen(filename, "w");
  if(outfile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return;
  }

  for(curr_level = 0; curr_level < tree->num_level-1; curr_level++) {
    fprintf(outfile, "%s", line_list[curr_level]);
    free(line_list[curr_level]);
    free(edge[curr_level]);
  }
  free(line_list);
  free(edge);
  free(num_edge);
  fclose(outfile);
}

double switch_error(char *testfilename, char *targfilename) {
  FILE *testfile, *targfile;
  char curr_line[MAX_LINE], targ_line[MAX_LINE];
  int i, n, error, total, *phase;

  // open the target file
  targfile = fopen(targfilename, "r");
  if(targfile == NULL) {
    printf("Error: could not open file \"%s\".\n", targfilename);
    return;
  }

  // open the test file
  testfile = fopen(testfilename, "r");
  if(testfile == NULL) {
    printf("Error: could not open file \"%s\".\n", testfilename);
    return;
  }
    
  fgets(curr_line, MAX_LINE, testfile);
  fgets(targ_line, MAX_LINE, targfile);
  for(i = 0; curr_line[i] > '0' && curr_line[i+1] > '0'; i += 2);
  n = i / 2;

  phase = malloc(n*sizeof(int));
  for(i = 0; i < n; i++) {
    if(curr_line[2*i] == targ_line[2*i] && curr_line[2*i+1] == targ_line[2*i+1])
      phase[i] = 0;
    else if(curr_line[2*i] == targ_line[2*i+1] && curr_line[2*i+1] == targ_line[2*i])
      phase[i] = 1;
    else {
      printf("Error: Haplotype mismatch.\n");
      return -1;
    }	
  }
  
  error = total = 0;
  while(fgets(curr_line, MAX_LINE, testfile) != NULL) {
    fgets(targ_line, MAX_LINE, targfile);
    for(i = 0; i < n; i++) {
      if(curr_line[2*i] == targ_line[2*i+phase[i]] && curr_line[2*i+1] == targ_line[2*i+(phase[i]+1)%2]) {
	if(curr_line[2*i] != curr_line[2*i+1])
	  total++;
      }
      else if(curr_line[2*i] == targ_line[2*i+(phase[i]+1)%2] && curr_line[2*i+1] == targ_line[2*i+phase[i]]) {
	phase[i] = 1 - phase[i];
	error++;
	total++;
      }
      else {
	printf("Error: Haplotype mismatch.\n");
	return -1;
      }
    }
  }
  free(phase);
  
  fclose(testfile);
  fclose(targfile);

  return (double)error / (double)total;
}

struct hashtable *new_hashtable(int size)
{
  struct hashtable *new = malloc(sizeof(struct hashtable));
  new->array = calloc(size, sizeof(struct hashtable_item*));
  new->size = size;
  new->filled = 0;

  return new;
}

void free_hashtable(struct hashtable *hashtable)
{
  int i;

  for(i = 0; i < hashtable->size; i++)
    free(hashtable->array[i]);
  free(hashtable->array);
  free(hashtable);
}
 
long int edge_pair_key(int curr_level, int i, int j, int *num_edge, long int *cum_edge)
{
  return cum_edge[curr_level] + (long int)(i*num_edge[curr_level]) + (long int)j;
}

void insert_hashtable(struct hashtable *hashtable, long int key, double prob, int back1, int back2)
{
  int index;
  struct hashtable_item *item;

  if(prob == -INFINITY)
    return;

  index = key % hashtable->size;

  while(hashtable->array[index] != NULL) {
    if(hashtable->array[index]->key == key) {
      hashtable->array[index]->prob = prob;
      hashtable->array[index]->back1 = back1;
      hashtable->array[index]->back1 = back2;
      return;
    }
    index++;
    index %= hashtable->size;
  }

  item = malloc(sizeof(struct hashtable_item));
  item->prob = prob;
  item->back1 = back1;
  item->back2 = back2;
  item->key = key;
  hashtable->array[index] = item;
  hashtable->filled++;

  if(hashtable->filled > 0.75*(double)hashtable->size)
    resize_hashtable(hashtable, 2*hashtable->size);
}

void insert_hashtable_logsum(struct hashtable *hashtable, long int key, double prob)
{
  int index;
  struct hashtable_item *item;

  if(prob == -INFINITY)
    return;

  index = key % hashtable->size;

  while(hashtable->array[index] != NULL) {
    if(hashtable->array[index]->key == key) {
      hashtable->array[index]->prob += log1p(exp(prob - hashtable->array[index]->prob));
      return;
    }
    index++;
    index %= hashtable->size;
  }

  item = malloc(sizeof(struct hashtable_item));
  item->prob = prob;
  item->back1 = -1;
  item->back2 = -1;
  item->key = key;
  hashtable->array[index] = item;
  hashtable->filled++;

  if(hashtable->filled > 0.75*(double)hashtable->size)
    resize_hashtable(hashtable, 2*hashtable->size);
}

void insert_hashtable_max(struct hashtable *hashtable, long int key, double prob, int back1, int back2)
{
  int index;
  struct hashtable_item *item;

  if(prob == -INFINITY)
    return;

  index = key % hashtable->size;

  while(hashtable->array[index] != NULL) {
    if(hashtable->array[index]->key == key && prob > hashtable->array[index]->prob) {
      hashtable->array[index]->prob = prob;
      hashtable->array[index]->back1 = back1;
      hashtable->array[index]->back2 = back2;
      return;
    }
    index++;
    index %= hashtable->size;
  }

  item = malloc(sizeof(struct hashtable_item));
  item->prob = prob;
  item->back1 = back1;
  item->back2 = back2;
  item->key = key;
  hashtable->array[index] = item;
  hashtable->filled++;

  if(hashtable->filled > 0.75*(double)hashtable->size)
    resize_hashtable(hashtable, 2*hashtable->size);
}

struct hashtable_item *search_hashtable(struct hashtable *hashtable, long int key)
{
  int index;

  index = key % hashtable->size;

  while(hashtable->array[index] != NULL) {
    if(hashtable->array[index]->key == key)
      return hashtable->array[index];
    index++;
    index %= hashtable->size;
  }

  return NULL;
}

double search_hashtable_prob(struct hashtable *hashtable, long int key)
{
  struct hashtable_item *item = search_hashtable(hashtable, key);

  if(item != NULL)
    return item->prob;
  else
    return -INFINITY;
}

void resize_hashtable(struct hashtable *hashtable, int size)
{
  int i, old_size;
  struct hashtable_item **old_array;

  old_array = hashtable->array;
  old_size = hashtable->size;
  
  hashtable->array = calloc(size, sizeof(struct hashtable_item*));
  hashtable->size = size;
  hashtable->filled = 0;
  
  for(i = 0; i < old_size; i++) {
    if(old_array[i] != NULL)
      insert_hashtable(hashtable, old_array[i]->key, old_array[i]->prob, old_array[i]->back1, old_array[i]->back2);
    free(old_array[i]);
  }
  free(old_array);
}

