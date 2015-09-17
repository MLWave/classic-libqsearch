/*
Copyright (c) 2003-2008 Rudi Cilibrasi, Rulers of the RHouse
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of the University nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE RULERS AND CONTRIBUTORS ``AS IS'' AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED.  IN NO EVENT SHALL THE RULERS OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.
*/

#include <complearn.h>
#include "qsearch.h"
#include "qsearch/qsearch-maketree-private.h"

void qsearch_labeledtree_free(QLabeledTree *ql)
{
  qsearch_tree_free(ql->qs);
  g_strfreev(ql->labels);
  gsl_matrix_free(ql->mat);
  free(ql);
}

QLabeledTree *qsearch_tree_add_labels(const QSearchTree *tree, const LabeledMatrix *lm) {
  QLabeledTree *qlab = calloc(sizeof(*qlab), 1);
  qlab->qs = qsearch_tree_clone(tree);
  qlab->labels = calloc(qlab->qs->total_node_count + 1, sizeof(gpointer));
  int i;
  for (i = 0; i < qlab->qs->total_node_count; i += 1) {
    char *c;
    if (qsearch_tree_get_neighbor_count(qlab->qs, i) == 1) {
      //c = g_strdup(lm->labels1[qsearch_tree_get_column_number(qlab->qs, i)]);
      // TODO: careful
      c = g_strdup(lm->labels1[i]);
    }
    else
      c = g_strdup("");
    qlab->labels[i] = c;
  }
  qlab->mat = gsl_matrix_alloc(lm->mat->size1, lm->mat->size2);
  gsl_matrix_memcpy(qlab->mat, lm->mat);
  return qlab;
}

GString *qsearch_tree_to_dot(const QLabeledTree *tree) {
  char *str;
  int i, j;
  gboolean show_details =
          qsearch_maketree_get_dot_show_details(qsearch_maketree_top());
  gboolean show_ring =
          qsearch_maketree_get_dot_show_ring(qsearch_maketree_top());
  GString *r = g_string_new("");
  str = g_strdup_printf("graph \"%s\" {\n",
          qsearch_maketree_get_dot_title(qsearch_maketree_top()));
  g_string_append(r, str); g_free(str);
  if (show_details) {
    str = g_strdup_printf("label=\"S(T)=%f\";\n", qsearch_tree_score_tree(tree->qs, tree->mat));
    g_string_append(r, str); g_free(str);
  }
  for (i = 0; i < tree->qs->total_node_count; i += 1) {
    int nodenum = i;
    if (qsearch_tree_get_neighbor_count(tree->qs, nodenum) == 1)
      nodenum = qsearch_tree_get_column_number(tree->qs, nodenum);
    str = g_strdup_printf("%d [label=\"%s\"];\n", nodenum, tree->labels[i]);
    g_string_append(r, str); g_free(str);
  }
  for (i = 0; i < tree->qs->total_node_count; i += 1) {
    int nodenumi = i;
    if (qsearch_tree_get_neighbor_count(tree->qs, nodenumi) == 1)
      nodenumi = qsearch_tree_get_column_number(tree->qs, nodenumi);
    for (j = i+1; j < tree->qs->total_node_count; j += 1) {
      int nodenumj = j;
      if (qsearch_tree_get_neighbor_count(tree->qs, nodenumj) == 1)
        nodenumj = qsearch_tree_get_column_number(tree->qs, nodenumj);
      if (qsearch_tree_is_connected(tree->qs,i,j)) {
        str = g_strdup_printf("%d -- %d [weight=\"2\"];\n", nodenumi, nodenumj);
        g_string_append(r, str); g_free(str);
      }
    }
  }
  GArray *lcirc = qsearch_tree_walk_filtered(tree->qs, 0, NODE_TYPE_LEAF);
  int qq=tree->qs->total_node_count+10000;
  for (i = 0; i < lcirc->len; i += 1) {
    int n1 = g_array_index(lcirc, guint32, i);
    int n2 = g_array_index(lcirc, guint32, (i+1)%(lcirc->len));
    int c1, c2;
    double dist = gsl_matrix_get(tree->mat,
            c1=qsearch_tree_get_column_number(tree->qs, n1),
            c2=qsearch_tree_get_column_number(tree->qs, n2));
    str = g_strdup_printf("%d -- %d [style=\"dotted\"%s];\n",c1,qq,show_ring?"":",color=\"white\"");
    g_string_append(r, str); g_free(str);
    str = g_strdup_printf("%d -- %d [style=\"dotted\"%s];\n",c2,qq,show_ring?"":",color=\"white\"");
    g_string_append(r, str); g_free(str);
    if (show_ring)
      str = g_strdup_printf("%d [label=\"%03.3f\",color=\"white\"];\n",qq,dist);
    else
      str = g_strdup_printf("%d [label=\"\",color=\"white\"];\n",qq);
    g_string_append(r, str); g_free(str);
    qq += 1;
  }
  g_string_append(r, "}\n");
  return r;
}
#if 0
struct _QSearchTree {
  int total_node_count;
  int must_recalculate_paths;
  MutationStatistics ms;
  double score;
  gboolean f_score_good;
  GPtrArray *n;
  GArray *spm, *nodeflags, *p1, *p2;
  GArray *leaf_placement;
};
#endif
