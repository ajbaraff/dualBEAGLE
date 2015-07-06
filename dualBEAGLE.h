#ifndef dualBEAGLE_H_  
#define dualBEAGLE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAX_LINE 65536
#define START_SIZE 16

struct haplotype
{
  int length, *allele;
  double freq;
};

struct tree_node
{
  int level, id, count, num_parent, num_child, num_haplo, marked, *count_child, *allele_parent, *allele_child, *allele_miss, *child_id, *parent_id;
  double p_node, *p_parent, *p_child, *back_weight;
  struct haplotype **haplo;
  struct tree_node **parent, **child;
};

struct tree
{
  struct tree_node ***level;
  int num_level, *num_node, *num_allele;
};

struct queue_node
{
  void *node;
  struct queue_node *next;
};

struct queue
{
  struct queue_node *front;
  struct queue_node *rear;
};

struct tree *read_input(char*);
struct tree_node *new_node(int, int, int);
void add_edge(struct tree_node*, struct tree_node*, int, int);
struct queue *new_queue(void);
int empty(struct queue*);
void enqueue(struct queue*, void*);
void *dequeue(struct queue*);
void prep_tree(struct tree*, struct tree*);
void print_tree(char*, struct tree*, int);
void calc_tree(struct tree*, struct tree*);
double loglik(struct tree*);
int calc_freq(struct haplotype***, struct tree*, int, int, int);
int compare(const void*, const void*);
void sim_tree(char*, struct tree*, int, int);

#endif // dualBEAGLE_H_
