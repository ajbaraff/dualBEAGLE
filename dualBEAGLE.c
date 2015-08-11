#include "dualBEAGLE.h"

struct tree *build_tree(char *filename)
{
  FILE *infile;
  char line[MAX_LINE], *line_ptr;
  int i, j, k, l, curr_level, curr_index, curr_allele, max_allele, num_samp, num_level = START_SIZE, *num_node, *num_allele, *node_list;
  struct tree_node *curr_node, ***level;
  struct tree *ret;

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
  while(sscanf(line_ptr, "%*s %n", &n) == 1) {
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
  level = malloc(num_level*sizeof(struct tree_node**));
  level[0] = malloc(sizeof(struct tree_node*));
  level[0][0] = new_node(0, 0, num_allele[0]);
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
	level[curr_level+1] = malloc(num_allele[curr_level]*num_node[curr_level]*sizeof(struct tree_node*));
	curr_index = 0;
	while(sscanf(line_ptr, "%d %n", &curr_allele, &n) == 1) {
	  curr_node <- level[curr_level][node_list[curr_index]];
	  if(curr_node->allele_miss[curr_allele-1]) {
	    level[curr_level+1][num_node[curr_level+1]] = new_node(curr_level+1, num_node[curr_level+1], num_allele[curr_level+1]);
	    add_edge(curr_node, level[curr_level+1][num_node[curr_level+1]], curr_allele, 1);
	    node_list[curr_index] = num_node[curr_level+1];
	    num_node[curr_level+1]++;
	  }
	  else {
	    curr_node->count++;
	    for(k = 0; k < curr_node->num_child && curr_node->allele_child[k] != curr_allele; k++);
	    curr_node->count_child[k]++;
	    node_list[curr_index] = curr_node->child[k]->id;
	  }
	  curr_index++;
	  line_ptr += n;
	}
      }

      else {
	level[curr_level+1] = malloc(sizeof(struct tree_node*));
	level[curr_level+1][0] = new_node(curr_level+1, 0, num_allele[curr_level+1]);
	num_node[curr_level+1] = 1;
        curr_index = 0;
        while(sscanf(line_ptr, "%d %n", &curr_allele, &n) == 1) {
          curr_node <- level[curr_level][node_list[curr_index]];
          if(curr_node->allele_miss[curr_allele-1]) add_edge(curr_node, level[curr_level+1][0], curr_allele, 1);
	  else {
            curr_node->count++;
            for(k = 0; k < curr_node->num_child && curr_node->allele_child[k] != curr_allele; k++);
            curr_node->count_child[k]++;
          }
          curr_index++;
          line_ptr += n;
        }
      }

      level[curr_level+1] = realloc(level[curr_level+1], num_node[curr_level+1]*sizeof(struct tree_node*));

      curr_level++;
      if(curr_level == num_level) {
        num_level *= 2;
        num_allele = realloc(num_allele, num_level*sizeof(int));
      }
    }
  }

  /* add missing edges */
  for(i = 0; i < num_level; i++) {
    for(j = 0; j < num_node[i]; j++) {
      curr_node = level[i][j];
      for(k = 0; k < num_allele[i]; k++)
        if(curr_node->allele_miss[k])
          for(l = 0; l < num_node[i+1]; l++)
            add_edge(level[i][j], level[i+1][l], k+1, 0, 0);
    }
  }

  /* calculate initial edge weights */
  for(i = 0; i < num_level; i++) {
    for(j = 0; j < num_node[i]; j++) {
      curr_node = level[i][j];
      for(k = 0; k < curr_node->num_child; k++)
        curr_node->p_child[k] = (double)curr_node->count_child[k] / curr_node->count;
      for(k = 0; k < curr_node->num_parent; k++)
        curr_node->p_parent[k] = curr_node->parent[k]->p_child[curr_node->child_id[k]];
    }
  }

  /* done with input */
  fclose(infile);

  ret = malloc(sizeof(struct tree));
  ret->level = level;
  ret->num_level = num_level;
  ret->num_node = num_node;
  ret->num_allele = num_allele;

  return(ret);
}

char *reverse_input(char *filename) 
{
  FILE *infile, *outfile;
  char line[MAX_LINE];
  int i, j, num_line = 0;

  infile = fopen(filename, "r");
  if(infile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return 1;
  }

  /* open the output file */
  filename = strcat(filename, ".rev");
  outfile = fopen(filename, "w");
  if(outfile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return 1;
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
  new->count = new->num_parent = new->num_child = new->num_haplo = new->marked = 0;
  new->count_child = new->allele_parent = new->allele_child = new->child_id = new->parent_id = NULL;
  new->allele_miss = malloc(num_allele*sizeof(int));
  for(i = 0; i < num_allele; i++) new->allele_miss[i] = 1;
  new->p_node = (level==0);
  new->p_parent = new->p_child = new->back_weight = NULL;
  new->haplo = NULL;
  new->parent = new->child = NULL;

  return new;
}

void add_edge(struct tree_node *parent, struct tree_node *child, int allele, int count)
{
  parent->count += count;
  parent->num_child += 1;
  parent->count_child = realloc(parent->count_child, parent->num_child*sizeof(int));
  parent->count_child[parent->num_child-1] = count;
  parent->allele_child = realloc(parent->allele_child, parent->num_child*sizeof(int));
  parent->allele_child[parent->num_child-1] = allele;
  if(count > 0) parent->allele_miss[allele-1] = 0; 
  parent->parent_id = realloc(parent->parent_id, parent->num_child*sizeof(int));
  parent->parent_id[parent->num_child-1] = child->num_parent;
  parent->p_child = realloc(parent->p_child, parent->num_child*sizeof(double));
  parent->child = realloc(parent->child, parent->num_child*sizeof(struct tree_node*));
  parent->child[parent->num_child-1] = child;

  child->num_parent += 1;
  child->allele_parent = realloc(child->allele_parent, child->num_parent*sizeof(int));
  child->allele_parent[child->num_parent-1] = allele;
  child->child_id = realloc(child->child_id, child->num_parent*sizeof(int));
  child->child_id[child->num_parent-1] = parent->num_child - 1;
  child->p_parent = realloc(child->p_parent, child->num_parent*sizeof(double));
  child->parent = realloc(child->parent, child->num_parent*sizeof(struct tree_node*));
  child->parent[child->num_parent-1] = parent;
}

void free_node(struct tree_node *node)

void prep_tree(struct tree *tree1, struct tree *tree2)
{
  int i, j, k;
  double sum;
  struct tree_node *curr_node;
  struct queue *queue; 
 
  /* calculate node weights */
  queue = new_queue();
  enqueue(queue, tree1->level[0][0]);
  while(!empty(queue)) {
    curr_node = dequeue(queue);
    curr_node->marked = 0;
    for(i = 0; i < curr_node->num_child; i++) {
      curr_node->child[i]->p_node += curr_node->p_node * curr_node->p_child[i];
      if(!curr_node->child[i]->marked) {
        enqueue(queue, curr_node->child[i]);
        curr_node->child[i]->marked = 1;
      }
    }
  }

  for(i = 0; i < tree1->num_level; i++) {
    sum = 0;
    for(j = 0; j < tree1->num_node[i]; j++)
      sum += tree1->level[i][j]->p_node;
    for(j = 0; j < tree1->num_node[i]; j++)
      tree1->level[i][j]->p_node /= sum;
  }

  /* allocate memory for back_weight arrays and set forw_weight arrays */
  for(i = 0; i < tree1->num_level; i++)
    for(j = 0; j < tree1->num_node[i]; j++) {
      curr_node = tree1->level[i][j];
      curr_node->back_weight = calloc(tree2->num_node[tree2->num_level-i-1], sizeof(double));
    }
}

void print_tree(char* filename, struct tree* tree, int reverse)
{
  int i, j, k;
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
      fprintf(outfile, "\tNode ID = %d, Node Weight = %8.8lf\n", tree->level[i][j]->id+1, tree->level[i][j]->p_node);
      for(k = 0; k < tree->level[i][j]->num_child; k++)
        fprintf(outfile, "\t\tChild Node = %d, Child Allele = %d, Child Weight = %8.8lf\n", tree->level[i][j]->child[k]->id+1, tree->level[i][j]->allele_child[k], tree->level[i][j]->p_child[k]);
    }
  }

  fclose(outfile);
}

void calc_tree(struct tree *tree1, struct tree *tree2)
{
  int i, j, k, curr_level, prev_level = 0;
  double sum;
  struct tree_node *curr_node;
  struct queue *queue;

  queue = new_queue();
  curr_node = tree1->level[0][0];
  curr_node->back_weight[0] = 1;
  enqueue(queue, curr_node);
  while(!empty(queue)) {
    curr_node = dequeue(queue);
    curr_node->marked = 0;
    curr_level = tree2->num_level - curr_node->level - 1;
    for(i = 0; i < curr_node->num_child; i++) {
      if(!curr_node->child[i]->marked) {
        enqueue(queue, curr_node->child[i]);
        curr_node->child[i]->marked = 1;
        curr_node->child[i]->p_node = 0;
        for(j = 0; j < tree2->num_node[curr_level-1]; j++)
          curr_node->child[i]->back_weight[j] = 0;
      }
      curr_node->p_child[i] = 0;
      for(j = 0; j < tree2->num_node[curr_level]; j++) {
        if(curr_node->back_weight[j] > 0) {
          for(k = 0; k < tree2->level[curr_level][j]->num_parent; k++) {
            if(tree2->level[curr_level][j]->allele_parent[k] == curr_node->allele_child[i]) {
              curr_node->p_child[i] += curr_node->back_weight[j] * tree2->level[curr_level][j]->p_parent[k] * tree2->level[curr_level][j]->parent[k]->p_node / curr_node->p_node;
              curr_node->child[i]->back_weight[tree2->level[curr_level][j]->parent[k]->id] += curr_node->back_weight[j] * tree2->level[curr_level][j]->p_parent[k];
            }
          }
        }
      }
      curr_node->child[i]->p_parent[curr_node->parent_id[i]] = curr_node->p_child[i];
      curr_node->child[i]->p_node += curr_node->p_node * curr_node->p_child[i];
    }
  }
}

double loglik(struct tree* tree) 
{
  int i, j, k;
  double loglik = 0;
  struct tree_node *curr_node;
  
  for(i = 0; i < tree->num_level; i++)
    for(j = 0; j < tree->num_node[i]; j++) {
      curr_node = tree->level[i][j];
      for(k = 0; k < curr_node->num_child; k++) 
	if(curr_node->count_child[k] > 0)
	  loglik += curr_node->count_child[k] * log(curr_node->p_child[k]);
    }
  
  return loglik;
}

int calc_freq(struct haplotype ***results, struct tree* tree, int start, int length, int reverse)
{
  int i, j, index, num_haplo = 0, *temp;
  struct tree_node *curr_node;
  struct queue *queue;

  if(reverse) start = tree->num_level - start - length + 1;

  queue = new_queue();
  for(i = 0; i < tree->num_node[start-1]; i++) {
    curr_node = tree->level[start-1][i];
    curr_node->num_haplo = 1;
    curr_node->haplo = malloc(sizeof(struct haplotype*));
    curr_node->haplo[0] = malloc(sizeof(struct haplotype));
    curr_node->haplo[0]->length = 0;
    curr_node->haplo[0]->allele = NULL;
    curr_node->haplo[0]->freq = curr_node->p_node;
    enqueue(queue, curr_node);
  }

  while(!empty(queue)) {
    curr_node = dequeue(queue);
    curr_node->marked = 0;
    for(i = 0; i < curr_node->num_child; i++) {
      if(!curr_node->child[i]->marked && curr_node->child[i]->level < start+length) {
        enqueue(queue, curr_node->child[i]);
        curr_node->child[i]->marked = 1;
      }
      curr_node->child[i]->haplo = realloc(curr_node->child[i]->haplo, (curr_node->child[i]->num_haplo+curr_node->num_haplo)*sizeof(struct haplotype*));
      for(j = 0; j < curr_node->num_haplo; j++) {
	index = curr_node->child[i]->num_haplo + j;
	curr_node->child[i]->haplo[index] = malloc(sizeof(struct haplotype));
	curr_node->child[i]->haplo[index]->length = curr_node->haplo[j]->length + 1;
	curr_node->child[i]->haplo[index]->allele = malloc(curr_node->child[i]->haplo[index]->length*sizeof(int));
	memcpy(curr_node->child[i]->haplo[index]->allele, curr_node->haplo[j]->allele, curr_node->haplo[j]->length*sizeof(int));
	curr_node->child[i]->haplo[index]->allele[curr_node->haplo[j]->length] = curr_node->allele_child[i];
	curr_node->child[i]->haplo[index]->freq = curr_node->haplo[j]->freq * curr_node->p_child[i];
      }
      curr_node->child[i]->num_haplo += curr_node->num_haplo;
    }
    if(curr_node->level < start+length-1) curr_node->num_haplo = 0;
  }
  /* free old haplos too */

  /* free old results? */
  *results = NULL;
  for(i = 0; i < tree->num_node[start+length-1]; i++) {
    curr_node = tree->level[start+length-1][i];
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

void sim_tree(char* filename, struct tree* tree, int n, int seed) 
{
  int i, j, k, max_index, *node_list;
  double sum_p, rand_p;
  struct tree_node *curr_node;
  FILE *outfile;

  outfile = fopen(filename, "w");
  if(outfile == NULL) {
    printf("Error: could not open file \"%s\".\n", filename);
    return;
  }

  fprintf(outfile, "I id");
  for(i = 0; i < n; i++)
    fprintf(outfile, " SIM%d", i+1);
  fprintf(outfile, "\n");

  node_list = calloc(n, sizeof(int));
  srand(seed);
  for(i = 0; i < tree->num_level-1; i++) {
    fprintf(outfile, "M m%d", i+1);
    for(j = 0; j < n; j++) {
      curr_node = tree->level[i][node_list[j]];
      sum_p = 0;
      rand_p = (double)rand() / (double)RAND_MAX;
      for(k = 0; k < curr_node->num_child && sum_p <= rand_p; k++)
	sum_p += curr_node->p_child[k];
      fprintf(outfile, " %d", curr_node->allele_child[k-1]);
      node_list[j] = curr_node->child[k-1]->id;
    }
    fprintf(outfile, "\n");
  }

  fclose(outfile);
  free(node_list);
}

