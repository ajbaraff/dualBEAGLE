#ifndef dualBEAGLE_H_  
#define dualBEAGLE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "queue.h"
#define MAX_LINE 65536
#define START_SIZE 16

struct haplotype
{
  int length, *allele;
  double freq;
};

struct tree_node
{
  int level, id, count, num_parent, num_allele, num_forw, num_back, num_haplo, marked, *count_child, *allele_parent, *allele_miss;
  double p_node, *p_child, *forw_weight, *back_weight;
  struct haplotype **haplo;
  struct tree_node **parent, ***child;
};

struct tree
{
  struct tree_node ***node;
  int num_level, *num_node, *num_allele;
};

struct tree *read_input(char*);
char *reverse_input(char*);
struct tree_node *new_node(int, int, int);
void add_edge(struct tree*, struct tree_node*, struct tree_node*, int, int);
void free_tree(struct tree*);
void free_node(struct tree_node*);
struct tree *copy_tree(struct tree*);
struct tree_node *copy_node(struct tree_node*);
void print_tree(char*, struct tree*, int);
struct tree_node *merge_node(struct tree*, struct tree_node*, struct tree_node*, int);
double loglik(struct tree*);
double loglik_pen(struct tree*, double);
struct tree *merge_tree(struct tree*, double);

void prep_tree(struct tree*, struct tree*);
void calc_tree(struct tree*, struct tree*);
int calc_freq(struct haplotype***, struct tree*, int, int, int);
int compare(const void*, const void*);
void sim_tree(char*, struct tree*, int, int);

#endif // dualBEAGLE_H_
