#ifndef __COMPLEARN_TREEMASTER_H
#define __COMPLEARN_TREEMASTER_H

#include <gsl/gsl_matrix.h>

struct _QSearchTreeObserver {
  void *user_data;
  void (*tree_search_started)(gsl_matrix *dm, void *user_data);
  void (*tried_to_improve)(gsl_matrix *dm, QSearchTree *old,
                                           QSearchTree *improved,
                           void *user_data);
  void (*tree_search_done)(gsl_matrix *dm, QSearchTree *final,
                           void *user_data);
};

typedef struct _QSearchTreeObserver QSearchTreeObserver;

struct _QSearchTreeMaster;
typedef struct _QSearchTreeMaster QSearchTreeMaster;

QSearchTreeMaster *qsearch_treemaster_new(gsl_matrix *dm);
void qsearch_treemaster_add_observer(
  QSearchTreeMaster *clt,
  void (*tree_search_started)(gsl_matrix *dm, void *user_data),
  void (*tried_to_improve)(gsl_matrix *dm, QSearchTree *old,
                                           QSearchTree *improved,
                           void *user_data),
  void (*tree_search_done)(gsl_matrix *dm, QSearchTree *final,
                           void *user_data),
void *user_data);
void qsearch_treemaster_free(QSearchTreeMaster *clt);
QSearchTree *qsearch_treemaster_find_best_tree(QSearchTreeMaster *clt);
gboolean qsearch_treemaster_was_search_stopped(QSearchTreeMaster *clt);
void qsearch_treemaster_stop_search(QSearchTreeMaster *clt);
double qsearch_treemaster_get_lmsd(QSearchTreeMaster *clt);

#endif
