#ifndef dualBEAGLE_H_  
#define dualBEAGLE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAX_LINE 262144
#define MAX_NODE 2000
#define EPSILON 1e-15

struct tree_node
{
  int level, id, count, num_allele, num_edge, num_parent, num_back, marked, *count_child, *allele_miss;
  double p_node, *p_child, *back_weight;
  struct tree_edge **child, **parent;
};

struct tree_edge
{
  int id, parent_id, allele;
  struct tree_node *node_from, *node_to;
};
 
struct tree
{
  struct tree_node ***node;
  int num_level, num_param_node, num_param_edge, *num_node, *num_allele;
  double loglik;
};

struct haplotype
{
  int count, length, *allele;
  double freq;
};

struct hashtable_item
{
  double prob;
  int back1, back2;
  long int key;
};

struct hashtable
{
  struct hashtable_item **array;
  int size, filled;
};
  
struct tree *read_input(char*);
struct tree *read_bgl(char*);
char *reverse_input(char*);
void split_input(char*, int, int);
struct tree_node *new_node(int, int, int);
void add_edge(struct tree*, struct tree_node*, struct tree_node*, int);
void free_tree(struct tree*);
void free_node(struct tree_node*);
struct tree *copy_tree(struct tree*);
struct tree_node *copy_node(struct tree_node*);
struct tree *cut_tree(struct tree*, int, int, int);
void calc_loglik(struct tree*);
void count_param(struct tree*);
double calc_loss(char*, struct tree*, double);
double search_epsilon(char*, struct tree*, double);
double merge_test(struct tree_node*, struct tree_node*, double, double, double, double*);
struct tree_node *merge_node(struct tree*, struct tree_node*, struct tree_node*, double, int);
struct tree *merge_tree(struct tree*, double);
void print_tree(char*, struct tree*, int);
void print_view(char*, struct tree*);
void print_pop(char*, char*, char*, struct tree*);
void print_haplo(char*, struct haplotype**, int);
int calc_freq(struct haplotype***, struct tree*, int, int, int);
int compare_haplo(const void*, const void*);
void calc_tree(struct tree*, struct tree*);
void random_phase(char*, char*, int);
double sample_phase(char*, struct tree*, int);
void viterbi_phase(char*, struct tree*);
double switch_error(char*, char*);  

struct hashtable *new_hashtable(int);
void free_hashtable(struct hashtable*);
long int edge_pair_key(int, int, int, int*, long int*);
void insert_hashtable(struct hashtable*, long int, double, int, int);
void insert_hashtable_logsum(struct hashtable*, long int, double);
void insert_hashtable_max(struct hashtable*, long int, double, int, int);
struct hashtable_item *search_hashtable(struct hashtable*, long int);
double search_hashtable_prob(struct hashtable*, long int);
void resize_hashtable(struct hashtable*, int);

#endif // dualBEAGLE_H_
