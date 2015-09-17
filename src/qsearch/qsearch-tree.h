#ifndef __COMPLEARN_TREE_H
#define __COMPLEARN_TREE_H

#include <glib.h>
#include <gsl/gsl_matrix.h>
#include <complearn/complearn.h>
#include "qsearch-neighborlist.h"

#define NODE_TYPE_LEAF 0x01
#define NODE_TYPE_KERNEL 0x02
#define NODE_TYPE_ALL 0x03

/* A Tree has N leaves and N-2 nodes for a total of 2N - 2 nodes */

struct _MutationStatistics {
  int total_complex_mutations;
  int total_simple_mutations;
  int total_successful_mutations;
  int last_simple_mutations;
  int total_clonings;
  int total_order_simple_mutations;
  int total_order_complex_mutations;
  int last_order_simple_mutations;
};

typedef struct _MutationStatistics MutationStatistics;

struct _QSearchTree {
  int total_node_count;
  gboolean must_recalculate_paths;
  gboolean dist_calculated;
  gboolean f_score_good;
  double dist_min;
  double dist_max;
  MutationStatistics ms;
  double score;
  GPtrArray *n;
  GArray *spm, *nodeflags, *p1, *p2;
  GArray *leaf_placement;
};

typedef struct _QSearchTree QSearchTree;

struct _QLabeledTree {
  QSearchTree *qs;
  char **labels;
  gsl_matrix *mat;
};

typedef struct _QLabeledTree QLabeledTree;

QSearchTree *qsearch_tree_new(int how_many_leaves);
void qsearch_tree_free(QSearchTree *clt);
gboolean qsearch_tree_is_leaf_node(const QSearchTree *clt, int whichNode);
guint32 qsearch_tree_get_node_count(const QSearchTree *clt);
guint32 qsearch_tree_get_leaf_node_count(const QSearchTree *clt);
guint32 qsearch_tree_get_kernel_node_count(const QSearchTree *clt);
gsl_matrix *qsearch_tree_get_adjacency_matrix(const QSearchTree *clt);
gboolean qsearch_tree_is_connected(const QSearchTree *clt, guint32 a, guint32 b);
gboolean qsearch_tree_is_standard_tree(const QSearchTree *clt);
guint32 qsearch_tree_get_neighbor_count(const QSearchTree *clt, guint32 a);
void qsearch_tree_connect(QSearchTree *clt, guint32 a, guint32 b);
void qsearch_tree_disconnect(QSearchTree *clt, guint32 a, guint32 b);
GArray *qsearch_tree_find_path(QSearchTree *clt, guint32 a, guint32 b);
gboolean qsearch_tree_is_consistent_quartet(QSearchTree *clt, guint32 a, guint32 b, guint32 c, guint32 d);
double qsearch_tree_score_tree(QSearchTree *clt, const gsl_matrix *dm);
QSearchTree *qsearch_tree_clone(const QSearchTree *clt);
guint32 qsearch_tree_get_random_node(const QSearchTree *clt, guint32 what_kind);
guint32 qsearch_tree_get_random_node_but_not(const QSearchTree *clt, guint32 what_kind, guint32 but_not);
guint32 qsearch_tree_find_path_length(QSearchTree *clt,
                                        guint32 a, guint32 b);
guint32 qsearch_tree_get_random_neighbor(const QSearchTree *clt, guint32 who);
GArray *qsearch_tree_get_neighbors(const QSearchTree *clt, guint32 who);
gboolean qsearch_tree_is_valid_tree(QSearchTree *clt);
void qsearch_tree_complex_mutation(QSearchTree *clt);
void qsearch_tree_simple_mutation(QSearchTree *clt);
void qsearch_tree_simple_mutation_leaf_swap(QSearchTree *clt);
void qsearch_tree_simple_mutation_subtree_transfer(QSearchTree *clt);
void qsearch_tree_simple_mutation_subtree_interchange(QSearchTree *clt);
void qsearch_tree_simple_mutation(QSearchTree *clt);
gboolean qsearch_tree_can_subtree_transfer(const QSearchTree *clt);
gboolean qsearch_tree_can_subtree_interchange(const QSearchTree *clt);
QSearchTree *qsearch_tree_find_better_tree(QSearchTree *clt,
guint32 howManyTries, const gsl_matrix *dm);
GArray *qsearch_tree_walk_tree(const QSearchTree *clt, guint32 fromwhere, gboolean f_bfs);
GArray *qsearch_tree_walk_tree_bfs(const QSearchTree *clt, guint32 fromwhere);
GArray *qsearch_tree_walk_tree_dfs(const QSearchTree *clt, guint32 fromwhere);
double qsearch_tree_calculate_order_cost(const QSearchTree *clt,const gsl_matrix *dm);
guint32 qsearch_tree_get_column_number(const QSearchTree *clt, guint32 nodenum);
guint32 qsearch_tree_get_node_number(const QSearchTree *clt, guint32 columnnum);
void qsearch_tree_mutate_order_simple(QSearchTree *clt);
void qsearch_tree_mutate_order_complex(QSearchTree *clt);
char *qsearch_tree_to_s(const QSearchTree *clt);
GArray *qsearch_tree_flipped_node_order(const QSearchTree *clt);
QLabeledTree *qsearch_tree_add_labels(const QSearchTree *tree,
                                      const LabeledMatrix *lm);
GString *qsearch_tree_to_dot(const QLabeledTree *tree);
GArray *qsearch_tree_walk_filtered(const QSearchTree *clt, guint32 fromwhere,
                                   guint32 filterType);
char *qsearch_tree_to_nexus(QLabeledTree *tree);
char *qsearch_tree_to_nexus_full(QLabeledTree *tree, LabeledMatrix *lm);
void qsearch_labeledtree_free(QLabeledTree *ql);

/* Returns true if every node has exactly 1 or 3 neighbors */
gboolean qsearch_tree_is_tree_ternary(QSearchTree *clt);

/* Sets the "connectedness" state (TRUE or FALSE) between nodes a and b and
 * returns the old connectedness status that was overwritten.
 */
gboolean qsearch_tree_set_connected(QSearchTree *clt, guint32 a, guint32 b, gboolean newconstate);

void qsearch_tree_clear_all_connections(QSearchTree *clt);
gdouble qsearch_tree_read_from_dot(QSearchTree *clt, GString *treedot, LabeledMatrix *lm);

#endif
