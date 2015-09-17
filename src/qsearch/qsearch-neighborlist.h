#ifndef __COMPLEARN_NEIGHBOR_LIST_H
#define __COMPLEARN_NEIGHBOR_LIST_H

#include <glib.h>

struct _QSearchNeighborList {
  GArray *n; // list of guint32 neighbors
};

typedef struct _QSearchNeighborList QSearchNeighborList;

QSearchNeighborList *new_qsearch_neighborlist(void);
void free_qsearch_neighborlist(QSearchNeighborList *cl);
void qsearch_neighborlist_add_neighbor(QSearchNeighborList *cln, guint32 w);
void qsearch_neighborlist_remove_neighbor(QSearchNeighborList *cln, guint32 w);
gboolean qsearch_neighborlist_has_neighbor(QSearchNeighborList *cln, guint32 w);
gint32 qsearch_neighborlist_find_index(QSearchNeighborList *cln, guint32 w);
void qsearch_neighborlist_copy_from(QSearchNeighborList *d, QSearchNeighborList *s);
void qsearch_copy_intarray(GArray *d, GArray *s);
QSearchNeighborList *qsearch_neighborlist_clone(QSearchNeighborList *clnl);

#endif
