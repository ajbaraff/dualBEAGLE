#include <stdio.h>
#include "dualBEAGLE.h"

#define TOL 0.000001

int main(int argc, char* argv[])
{
  int num_haplo;
  double curr_loglik, prev_loglik, delta;
  struct tree *tree1, *tree2;
  struct haplotype **results;

  /* check for too few arguments */
  if(argc != 6) {
    printf("Error: wrong number of arguments - \"BEAGLEboth <input file 1> <input file 2> <output file for> <output file rev> <output file both>\".\n");
    return 1;
  }

  tree1 = read_input(argv[1]);
  tree2 = read_input(argv[2]);

  prep_tree(tree1, tree2);
  prep_tree(tree2, tree1);

  print_tree(argv[3], tree1, 0);
  print_tree(argv[4], tree2, 1);

  prev_loglik = loglik(tree1);
  do {
    calc_tree(tree1, tree2);
    calc_tree(tree2, tree1);
    curr_loglik = loglik(tree1);
    delta = (curr_loglik-prev_loglik)/curr_loglik;
    prev_loglik = curr_loglik;
    } while(delta > TOL);

  calc_tree(tree1, tree2);
  calc_tree(tree2, tree1);

  print_tree(argv[5], tree1, 0);

  return 0;
}
