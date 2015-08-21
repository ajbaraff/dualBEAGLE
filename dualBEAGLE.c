#include "dualBEAGLE.h"

struct tree *read_input(char *filename)
{
  FILE *infile;
  char line[MAX_LINE], *line_ptr;
  int i, j, k, l, n, curr_level, curr_index, curr_allele, max_allele, num_samp, num_level = START_SIZE, *num_node, *num_allele, *node_list;
  struct tree_node *curr_node, ***node;
  struct tree *tree;
  struct queue *queue;

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
  curr_level = max_allele = 0;
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
      max_allele = 0;

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
	    node_list[curr_index] = curr_node->child[curr_allele][0]->id;
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

  /* calculate initial edge weights */
  for(i = 0; i < num_level; i++) {
    for(j = 0; j < num_node[i]; j++) {
      curr_node = node[i][j];
      for(k = 0; k < curr_node->num_allele; k++)
        curr_node->p_child[k] = (double)curr_node->count_child[k] / curr_node->count;
    }
  }

  /* calculate node weights */
  queue = new_queue();
  enqueue(queue, tree->node[0][0]);
  while(!empty(queue)) {
    curr_node = dequeue(queue);
    curr_node->marked = 0;
    for(i = 0; i < curr_node->num_allele; i++) {
      if(curr_node->allele_miss[i]) 
	continue;
      curr_node->child[i][0]->p_node += curr_node->p_node * curr_node->p_child[i];
      if(!curr_node->child[i][0]->marked) {
	enqueue(queue, curr_node->child[i][0]);
	curr_node->child[i][0]->marked = 1;
      }
    }
  }
  free_queue(queue);

  /* add missing edges */
  for(i = 0; i < num_level; i++)
    for(j = 0; j < num_node[i]; j++) {
      curr_node = node[i][j];
      curr_node->num_forw = num_node[i+1];
      curr_node->forw_weight = malloc(num_node[i+1]*sizeof(double));
      for(k = 0; k < num_node[i+1]; k++)
	curr_node->forw_weight[k] = node[i+1][k]->p_node;
      for(k = 0; k < num_allele[i]; k++)
	if(curr_node->allele_miss[k])
          for(l = 0; l < num_node[i+1]; l++)
            add_edge(tree, curr_node, node[i+1][l], k, 0);
    }
  
  return(tree);
}

char *reverse_input(char *filename) 
{
  FILE *infile, *outfile;
  char line[MAX_LINE];
  int i, j, num_line = 0;

  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return NULL;
  }

  /* open the output file */
  filename = strcat(filename, ".rev");
  outfile = fopen(filename, "w");
  if(outfile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
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

  return(filename);
}

struct tree_node *new_node(int level, int id, int num_allele)
{
  int i;
  struct tree_node *new = malloc(sizeof(struct tree_node));

  new->level = level;
  new->id = id;
  new->count = new->num_parent = new->num_forw = new->num_back = new->num_haplo = new->marked = 0;
  new->num_allele = num_allele;
  new->count_child = calloc(num_allele, sizeof(int));
  new->allele_parent = NULL;
  new->allele_miss = malloc(num_allele*sizeof(int));
  for(i = 0; i < num_allele; i++) new->allele_miss[i] = 1;
  new->p_node = (level==0);
  new->forw_weight = new->back_weight = NULL;
  new->p_child = malloc(num_allele*sizeof(double));
  new->haplo = NULL;
  new->parent = NULL;
  new->child = malloc(num_allele*sizeof(struct tree_node**));
  for(i = 0; i < num_allele; i++) new->child[i] = NULL;

  return new;
}

void add_edge(struct tree *tree, struct tree_node *parent, struct tree_node *child, int allele, int count)
{
  parent->count += count;
  parent->count_child[allele] += count;

  if(count > 0) {
    parent->allele_miss[allele] = 0;
    parent->child[allele] = malloc(sizeof(struct tree_node*));
    parent->child[allele][0] = child;
  }
  else {
    if(parent->child[allele] == NULL)
      parent->child[allele] = malloc(tree->num_node[child->level]*sizeof(struct tree_node*));
    parent->child[allele][child->id] = child; 
  }

  child->num_parent++;
  child->allele_parent = realloc(child->allele_parent, child->num_parent*sizeof(int));
  child->allele_parent[child->num_parent-1] = allele;
  child->parent = realloc(child->parent, child->num_parent*sizeof(struct tree_node*));
  child->parent[child->num_parent-1] = parent;
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
  free(node->allele_parent);
  free(node->allele_miss);
  
  free(node->p_child);
  free(node->forw_weight);
  // free(node->back_weight);

  /* for(i = 0; i < node->num_haplo; i++) {
    free(node->haplo[i]->allele);
    free(node->haplo[i]);
  }
  free(node->haplo); */
  
  free(node->parent);
  for(i = 0; i < node->num_allele; i++) free(node->child[i]);
  free(node->child);

  free(node);
}

struct tree *copy_tree(struct tree *tree)
{
  int i, j, k, l;
  struct tree_node *curr_node;
  struct tree *new;

  new = malloc(sizeof(struct tree));
  new->num_level = tree->num_level;
  new->num_node = malloc(tree->num_level*sizeof(int));
  memcpy(new->num_node, tree->num_node, tree->num_level*sizeof(int));
  new->num_allele = malloc(tree->num_level*sizeof(int));
  memcpy(new->num_allele, tree->num_allele, tree->num_level*sizeof(int));
 
  new->node = malloc(tree->num_level*sizeof(struct tree_node**));
  for(i = 0; i < tree->num_level; i++) {
    new->node[i] = malloc(tree->num_node[i]*sizeof(struct tree_node*));
    for(j = 0; j < tree->num_node[i]; j++)
      new->node[i][j] = copy_node(tree->node[i][j]);
  }

  for(i = 0; i < tree->num_level; i++)
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = new->node[i][j];
      for(k = 0; k < curr_node->num_parent; k++)
      	curr_node->parent[k] = new->node[i-1][tree->node[i][j]->parent[k]->id];
      for(k = 0; k < curr_node->num_allele; k++) {
	if(curr_node->allele_miss[k]) 
	  for(l = 0; l < curr_node->num_forw; l++)
	    curr_node->child[k][l] = new->node[i+1][tree->node[i][j]->child[k][l]->id];
	else 
	  curr_node->child[k][0] = new->node[i+1][tree->node[i][j]->child[k][0]->id];
      }
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
  new->num_parent = node->num_parent;
  new->num_allele = node->num_allele;
  new->num_forw = node->num_forw;
  new->num_haplo = node->num_haplo;
  new->marked = node->marked;

  new->count_child = malloc(node->num_allele*sizeof(int));
  memcpy(new->count_child, node->count_child, node->num_allele*sizeof(int));
  new->allele_parent = malloc(node->num_parent*sizeof(int));
  memcpy(new->allele_parent, node->allele_parent, node->num_parent*sizeof(int));
  new->allele_miss = malloc(node->num_allele*sizeof(int));
  memcpy(new->allele_miss, node->allele_miss, node->num_allele*sizeof(int));

  new->p_node = node->p_node;

  new->p_child = malloc(node->num_allele*sizeof(double));
  memcpy(new->p_child, node->p_child, node->num_allele*sizeof(double));
  new->forw_weight = malloc(node->num_forw*sizeof(double));
  memcpy(new->forw_weight, node->forw_weight, node->num_forw*sizeof(double));
  // new->back_weight = malloc(node->num_back*sizeof(double)); 
  // memcpy(new->back_weight, node->back_weight, node->num_back*sizeof(double));

  /* new->haplo = malloc(node->num_haplo*sizeof(struct haplotype*));
     for(i = 0; i < node->num_haplo; i++) {
     new->haplo[i] = malloc(sizeof(struct haplotype));
     new->haplo[i]->length = node->haplo[i]->length;
     new->haplo[i]->allele = malloc(node->haplo[i]->length*sizeof(int));
     memcpy(new->haplo[i]->allele, node->haplo[i]->allele, node->haplo[i]->length*sizeof(int));
     new->haplo[i]->freq = node->haplo[i]->freq;
     } */

  new->parent = malloc(node->num_parent*sizeof(struct tree_node*));
  new->child = malloc(node->num_allele*sizeof(struct tree_node**));
  for(i = 0; i < node->num_allele; i++) {
    if(node->allele_miss[i])
      new->child[i] = malloc(node->num_forw*sizeof(struct tree_node*));
    else 
      new->child[i] = malloc(sizeof(struct tree_node*));
  }

  return new;
}

void print_tree(char* filename, struct tree* tree, int reverse)
{
  int i, j, k, l;
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
          for(l = 0; l < curr_node->num_forw; l++)
            fprintf(outfile, "\t\tChild Node = %d, Child Allele = %d, Child Weight = %8.8lf\n", curr_node->child[k][l]->id+1, k+1, curr_node->p_child[k]*curr_node->forw_weight[l]);
        else
          fprintf(outfile, "\t\tChild Node = %d, Child Allele = %d, Child Weight = %8.8lf\n", curr_node->child[k][0]->id+1, k+1, curr_node->p_child[k]);
      }
    }
  }

  fclose(outfile);
}

struct tree_node *merge_node(struct tree *tree, struct tree_node *node1, struct tree_node *node2, int first)
{
  int i, j, k, l, index;
  struct tree_node *curr_node, *new = malloc(sizeof(struct tree_node));

  if(node1 == node2) return(node1);

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
  new->num_parent = 0;
  new->num_allele = node1->num_allele;
  new->num_forw = node1->num_forw;
  new->num_haplo = node1->num_haplo + node2->num_haplo;
  new->marked = 0;
  new->p_node = node1->p_node + node2->p_node;

  new->forw_weight = malloc(new->num_forw*sizeof(double));
  new->allele_parent = NULL;
  new->parent = NULL;

  /*
  for(i = 0; i < tree->num_node[new->level+1]; i++) {
    curr_node = tree->node[new->level+1][i];
    for(j = 0; j < curr_node->num_parent; j++)
      if(curr_node->parent[j] == node1 || curr_node->parent[j] == node2)
	curr_node->parent[j] = new;
  }

  new->allele_parent = malloc(new->num_parent*sizeof(int));
  new->parent = malloc(new->num_parent*sizeof(struct tree_node*));
  for(i = 0; i < node1->num_parent; i++) {
    new->allele_parent[i] = node1->allele_parent[i];
    new->parent[i] = node1->parent[i];
  }
  index = 0;
  for(i = 0; i < node2->num_parent; i++) {
    if(node2->parent[i]->allele_miss[node2->allele_parent[i]])
      new->num_parent--;
    else {
      new->allele_parent[node1->num_parent+index] = node2->allele_parent[i];
      new->parent[node1->num_parent+index] = node2->parent[i];
      index++;
    }
  }
  new->allele_parent = realloc(new->allele_parent, new->num_parent*sizeof(int));
  new->parent = realloc(new->parent, new->num_parent*sizeof(struct tree_node*));
  */

  new->count_child = malloc(new->num_allele*sizeof(int));
  new->allele_miss = malloc(new->num_allele*sizeof(int));
  new->p_child = malloc(new->num_allele*sizeof(double));
  new->child = malloc(new->num_allele*sizeof(struct tree_node**));
  for(i = 0; i < new->num_allele; i++) {
    new->count_child[i] = node1->count_child[i] + node2->count_child[i];
    new->allele_miss[i] = node1->allele_miss[i] && node2->allele_miss[i];
    new->p_child[i] = (node1->p_child[i]*node1->p_node + node2->p_child[i]*node2->p_node) / new->p_node;
    if(new->allele_miss[i])
      new->child[i] = malloc(new->num_forw*sizeof(struct tree_node*));
    else 
      new->child[i] = malloc(sizeof(struct tree_node*));
  }
  for(i = 0; i < new->num_allele; i++) {
    if(new->allele_miss[i]) 
      for(j = 0; j < new->num_forw; j++)
        new->child[i][j] = node1->child[i][j];
    else {
      if(!node1->allele_miss[i] && node2->allele_miss[i])
        new->child[i][0] = node1->child[i][0];
      else if(node1->allele_miss[i] && !node2->allele_miss[i])
        new->child[i][0] = node2->child[i][0];
      else
        new->child[i][0] = merge_node(tree, node1->child[i][0], node2->child[i][0], 0);
    }
  }

  // back_weight and haplo merging?

  if(first) {
    for(i = 0; i < tree->num_node[new->level-1]; i++) {
      curr_node = tree->node[new->level-1][i];
      for(j = 0; j < curr_node->num_allele; j++) 
	if(!curr_node->allele_miss[j]) 
	  if(curr_node->child[j][0] == node1 || curr_node->child[j][0] == node2)
	    curr_node->child[j][0] = new;
    }

    for(i = new->level; i < tree->num_level; i++)
      for(j = 0; j < tree->num_node[i]; j++)
	tree->node[i][j]->num_parent = 0;

    for(i = new->level-1; i < tree->num_level-1; i++)
      for(j = 0; j < tree->num_node[i]; j++) {
	curr_node = tree->node[i][j];
	curr_node->num_forw = tree->num_node[i+1];
      	curr_node->forw_weight = realloc(curr_node->forw_weight, curr_node->num_forw*sizeof(double));
	for(k = 0; k < curr_node->num_forw; k++)
	  curr_node->forw_weight[k] = tree->node[i+1][k]->p_node;
	for(k = 0; k < curr_node->num_allele; k++) {
	  if(curr_node->allele_miss[k]) {
	    curr_node->child[k] = realloc(curr_node->child[k], curr_node->num_forw*sizeof(struct tree_node*));
	    for(l = 0; l < curr_node->num_forw; l++) {
	      curr_node->child[k][l] = tree->node[i+1][l];
	      curr_node->child[k][l]->num_parent++;
	      curr_node->child[k][l]->allele_parent = realloc(curr_node->child[k][l]->allele_parent, curr_node->child[k][l]->num_parent*sizeof(int));
	      curr_node->child[k][l]->allele_parent[curr_node->child[k][l]->num_parent-1] = k;
	      curr_node->child[k][l]->parent = realloc(curr_node->child[k][l]->parent, curr_node->child[k][l]->num_parent*sizeof(struct tree_node*));
	      curr_node->child[k][l]->parent[curr_node->child[k][l]->num_parent-1] = curr_node;
	    }
	  }
	  else {
	    curr_node->child[k][0]->num_parent++;
	    curr_node->child[k][0]->allele_parent = realloc(curr_node->child[k][0]->allele_parent, curr_node->child[k][0]->num_parent*sizeof(int));
	    curr_node->child[k][0]->allele_parent[curr_node->child[k][0]->num_parent-1] = k;
	    curr_node->child[k][0]->parent = realloc(curr_node->child[k][0]->parent, curr_node->child[k][0]->num_parent*sizeof(struct tree_node*));
	    curr_node->child[k][0]->parent[curr_node->child[k][0]->num_parent-1] = curr_node;
	  }
	}
      }


    /*
    for(i = 0; i < new->num_parent; i++) {
      printf("Level = %d, Node = %d, Allele = %d\n", new->parent[i]->level, new->parent[i]->id, new->allele_parent[i]);
      if(new->parent[i]->allele_miss[new->allele_parent[i]])
        new->parent[i]->child[new->allele_parent[i]][new->id] = new;
      else
        new->parent[i]->child[new->allele_parent[i]][0] = new;
    }

    for(i = 0; i < node1->num_parent; i++)
      for(j = 0; j < node2->num_parent; j++)
        if((node1->parent[i]->id == node2->parent[j]->id) && (node1->allele_parent[i] == node2->allele_parent[j]))
	new->parent[i]->p_child[new->allele_parent[i]] += node2->parent[j]->p_child[new->allele_parent[i]]; */

  }

  free_node(node1);
  free_node(node2);

  return(new);
}

double loglik(struct tree* tree)
{
  int i, j, k;
  double loglik = 0;
  struct tree_node *curr_node;

  for(i = 0; i < tree->num_level; i++)
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->node[i][j];
      for(k = 0; k < curr_node->num_allele; k++)
        if(curr_node->count_child[k] > 0)
          loglik += curr_node->count_child[k] * log(curr_node->p_child[k]);
    }

  return loglik;
}

double loglik_pen(struct tree* tree, double lambda)
{
  int i, j, k, l;
  double loglik = 0;
  struct tree_node *curr_node1, *curr_node2;

  for(i = 0; i < tree->num_level; i++)
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node1 = tree->node[i][j];
      for(k = 0; k < curr_node1->num_allele; k++)
        if(curr_node1->count_child[k] > 0) 
          loglik += curr_node1->count_child[k] * log(curr_node1->p_child[k]);
      for(k = j+1; k < tree->num_node[i]; k++) {
	curr_node2 = tree->node[i][k];
	for(l = 0; l < tree->num_allele[i]; l++)
	  loglik -= lambda * fabs(curr_node1->p_child[l] - curr_node2->p_child[l]);
      }
    }
  
  return loglik;
}

struct tree *merge_tree(struct tree* tree, double lambda)
{
  int i, j, k, max_j, max_k;
  double curr_loglik, new_loglik, max_diff;
  struct tree_node *curr_node1, *curr_node2;
  struct tree *temp;

  for(i = 0; i < tree->num_level; i++) {
    printf("Level: %d\n", i);
    do {
      max_diff = 0;
      curr_loglik = loglik_pen(tree, lambda);
      for(j = 0; j < tree->num_node[i]; j++)
	for(k = j+1; k < tree->num_node[i]; k++) {
	  printf("  Try Merge: (%d, %d) and (%d, %d).\n", i, j+1, i, k+1);
	  temp = copy_tree(tree);
	  curr_node1 = temp->node[i][j];
	  curr_node2 = temp->node[i][k];
	  merge_node(temp, curr_node1, curr_node2, 1);
	  new_loglik = loglik_pen(temp, lambda);
	  printf("    curr_loglik = %8.8lf, new_loglik = %8.8lf\n", curr_loglik, new_loglik);
	  if(new_loglik - curr_loglik > max_diff) {
	    max_diff = new_loglik - curr_loglik;
	    max_j = j;
	    max_k = k;
	  }
	  free_tree(temp);
	}
      if(max_diff > 0) {
	printf("Merging (%d, %d) and (%d, %d).\n", i, max_j+1, i, max_k+1);
	curr_node1 = tree->node[i][max_j];
	curr_node2 = tree->node[i][max_k];
	merge_node(tree, curr_node1, curr_node2, 1);
      }
    } while(max_diff > 0);
  }
  return tree;
}
