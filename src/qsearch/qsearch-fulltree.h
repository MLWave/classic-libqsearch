#ifndef FULLTREE_H_
#define FULLTREE_H_

#include <glib.h>
#include "qsearch-tree.h"

struct _FullTree {
    void* data; 
    guint32      node_count;
    double       raw_score;
};

typedef struct _FullTree FullTree;

extern FullTree* qsearch_make_fulltree(QSearchTree *clt, const gsl_matrix *dm);
extern void qsearch_free_fulltree(FullTree *tree);

extern void qsearch_fulltree_to_searchtree(FullTree *src, QSearchTree *dest);

extern gboolean qsearch_fulltree_can_swap(FullTree *tree, int A, int B);
extern void qsearch_fulltree_swap_nodes(FullTree *tree, int A, int B, const gsl_matrix *dm);
                
extern guint32 qsearch_fulltree_move_to(FullTree *tree, guint32 from, guint32 to);
extern guint32 qsearch_fulltree_find_sibling(FullTree *tree, guint32 node, guint32 ancestor);
extern double  qsearch_fulltree_sum_distance(FullTree *tree, int a, int b);

// get children in the context of an ancestor 
extern void    qsearch_fulltree_get_children(FullTree *tree, int node, int ancestor, int *child1, int *child2);

#endif
