#ifndef dualBEAGLE_H_  
#define dualBEAGLE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "queue.h"
#define MAX_LINE 65536
#define MAX_ENTRY 64
#define START_SIZE 16
#define TOL 0.00001
#define OFFSET 0.0001
#define MAX_NODE 200

struct haplotype
{
  int length, *allele;
  double freq;
};

struct tree_node
{
  int level, id, count, num_allele, num_back, num_haplo, marked, *count_child, *allele_miss;
  double p_node, *p_child, *p_child_test, *back_weight;
  struct haplotype **haplo;
  struct tree_node **child;
};

struct tree
{
  struct tree_node ***node;
  int num_level, *num_node, *num_allele;
  double loss, loss_test;
};

struct tree *read_input(char*);
struct tree *read_bgl(char*);
char *reverse_input(char*);
void split_input(char*, int, int);
struct tree_node *new_node(int, int, int);
void add_edge(struct tree*, struct tree_node*, struct tree_node*, int, int);
void free_tree(struct tree*);
void free_node(struct tree_node*);
struct tree *copy_tree(struct tree*);
struct tree_node *copy_node(struct tree_node*);
void print_tree(char*, struct tree*, int);
double loss(struct tree*, double);
double merge_test(struct tree*, struct tree_node*, struct tree_node*, double, int);
struct tree_node *merge_node(struct tree*, struct tree_node*, struct tree_node*, double, int);
struct tree *merge_tree(struct tree*, double);
void print_view(char*, struct tree*);
double test_loss(char*, struct tree*);
void calc_tree(struct tree*, struct tree*);


void prep_tree(struct tree*, struct tree*);
int calc_freq(struct haplotype***, struct tree*, int, int, int);
int compare(const void*, const void*);
void sim_tree(char*, struct tree*, int, int);

#endif // dualBEAGLE_H_
