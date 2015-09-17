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

#include <math.h>
#include "qsearch.h"

struct _QSearchTreeMaster
{
  GArray *forest;
  gsl_matrix *dm;
  GArray *obs;
  double lmsd;
  gboolean abort_search;
};

double qsearch_treemaster_get_lmsd(QSearchTreeMaster *clt)
{
  return clt->lmsd;
}

gboolean  qsearch_treemaster_was_search_stopped(QSearchTreeMaster *clt)
{
  return clt->abort_search;
}

void qsearch_treemaster_stop_search(QSearchTreeMaster *clt)
{
  clt->abort_search = TRUE;
}

void qsearch_treemaster_tree_observer_adaptor_search_started(gsl_matrix *a,
void *b) {
}

void qsearch_treemaster_tree_observer_adaptor_tried_to_improve(gsl_matrix *a,
QSearchTree *b, QSearchTree *c, void *d) {
}

void qsearch_treemaster_tree_observer_adaptor_search_done(gsl_matrix *a,
QSearchTree *b, void *c) {
}

void qsearch_treemaster_add_observer(
  QSearchTreeMaster *clt,
  void (*tree_search_started)(gsl_matrix *dm, void *user_data),
  void (*tried_to_improve)(gsl_matrix *dm, QSearchTree *old,
                                           QSearchTree *improved,
                           void *user_data),
  void (*tree_search_done)(gsl_matrix *dm, QSearchTree *final,
                           void *user_data),
void *user_data) {
  QSearchTreeObserver *cp = calloc(sizeof(*cp), 1);
  cp->tree_search_started = qsearch_treemaster_tree_observer_adaptor_search_started;
  cp->tried_to_improve = qsearch_treemaster_tree_observer_adaptor_tried_to_improve;
  cp->tree_search_done = qsearch_treemaster_tree_observer_adaptor_search_done;
  if (tree_search_started)
    cp->tree_search_started = tree_search_started;
  if (tried_to_improve)
    cp->tried_to_improve = tried_to_improve;
  if (tree_search_done)
    cp->tree_search_done = tree_search_done;
  cp->user_data = user_data;
  g_array_append_val(clt->obs, cp);
}

static int recommended_tree_duplicity(int how_many_leaves)
{
  double v = how_many_leaves;
  v = 26 * pow(v, -0.76) + 0.5;
  int i = (int) v;
  if (i < 2)
    i = 2;
  return i;
}

QSearchTreeMaster *qsearch_treemaster_new(gsl_matrix *dm)
{
  g_assert(dm != NULL && dm->size1 == dm->size2);
  QSearchTreeMaster *cltm = calloc(sizeof(*cltm), 1);
  cltm->forest = g_array_new(FALSE, TRUE, sizeof(gpointer));
  cltm->obs = g_array_new(FALSE, TRUE, sizeof(gpointer));
  int i, fs = recommended_tree_duplicity(dm->size1);
  for (i = 0; i < fs; i += 1) {
    QSearchTree *tree = qsearch_tree_new(dm->size1);
    qsearch_tree_complex_mutation(tree);
    qsearch_tree_complex_mutation(tree);
    g_array_append_val(cltm->forest, tree);
  }
  cltm->dm = gsl_matrix_calloc(dm->size1, dm->size2);
  gsl_matrix_memcpy(cltm->dm, dm);
  return cltm;
}

void qsearch_treemaster_free(QSearchTreeMaster *clt)
{
  int i;
  if (clt == NULL)
    return;
  for (i = 0; i < clt->obs->len; i += 1)
    free(g_array_index(clt->obs,gpointer,i));
  g_array_free(clt->obs, TRUE);
  for (i = 0; i < clt->forest->len; i += 1)
    qsearch_tree_free((QSearchTree *)g_array_index(clt->forest,gpointer,i));
  g_array_free(clt->forest, TRUE);
  gsl_matrix_free(clt->dm);
  free(clt);
}

gboolean qsearch_treemaster_is_done(QSearchTreeMaster *clt)
{
  const double MAXSCOREDIFF = 8e-14;
  if (clt->abort_search)
    return TRUE;
  clt->lmsd = -1.0;
  QSearchTree *orig = g_array_index(clt->forest, QSearchTree *, 0);
  double csco = qsearch_tree_score_tree(orig, clt->dm);
  int i;
  for (i = 1; i < clt->forest->len; i += 1) {
    double sco = qsearch_tree_score_tree(g_array_index(clt->forest, QSearchTree *, i), clt->dm);
    double deltasco = fabs(sco - csco);
    if (clt->lmsd == -1.0 || clt->lmsd < deltasco)
      clt->lmsd = deltasco;
    if (clt->lmsd > MAXSCOREDIFF)
      return FALSE;
  }
  return TRUE;
}

void qsearch_treemaster_try_to_improve_bucket(QSearchTreeMaster *clt, guint32 i)
{
  const int NUMTRIESPERBIGTRY = 24;
  int j;
  QSearchTree *cand, *old;
  old = g_array_index(clt->forest, QSearchTree *, i);
  cand = qsearch_tree_find_better_tree(old, NUMTRIESPERBIGTRY, clt->dm);
  if (!qsearch_treemaster_was_search_stopped(clt) && i == 0 && clt->obs->len > 0) {
    for (j = 0; j < clt->obs->len; j += 1) {
      QSearchTreeObserver *ob =
                   g_array_index(clt->obs,QSearchTreeObserver*,j);
                   ob->tried_to_improve(clt->dm, old, cand, ob->user_data);
    }
  }
  if (cand) {
    g_array_index(clt->forest, QSearchTree *, i) = cand;
    qsearch_tree_free(old);
  }
}

QSearchTree *qsearch_treemaster_find_best_tree(QSearchTreeMaster *clt)
{
  int i, j;
  double bestsco = -1.0;
  double osco, nsco;
  clt->abort_search = FALSE;
  if (clt->obs->len > 0) {
    for (j = 0; j < clt->obs->len; j += 1) {
      QSearchTreeObserver *ob =
                   g_array_index(clt->obs,QSearchTreeObserver*,j);
                   ob->tree_search_started(clt->dm, ob->user_data);
    }
  }
  do {
    for (i = 0; i < clt->forest->len; i += 1) {
      osco = qsearch_tree_score_tree(g_array_index(clt->forest, QSearchTree *, i), clt->dm);
      qsearch_treemaster_try_to_improve_bucket(clt, i);
      nsco = qsearch_tree_score_tree(g_array_index(clt->forest, QSearchTree *, i), clt->dm);
      if (nsco < osco) {
        fprintf(stderr, "Error, tree degraded: %f %f.\n", osco, nsco);
        exit(1);
      }
      if (nsco > bestsco)
        bestsco = nsco;
    }
  } while (!qsearch_treemaster_is_done(clt));
  QSearchTree *answer = qsearch_tree_clone(g_array_index(clt->forest, QSearchTree *, 0));
  double gsco;
  gsco = qsearch_tree_score_tree(answer, clt->dm);
  if (!qsearch_treemaster_was_search_stopped(clt)) {
  if (fabs(gsco - bestsco) > 1e-6) {
//    fprintf(stderr, "WARNING: Best score was %f but gsco was %f\n", bestsco, gsco);
    ;
  }
  }
  if (!qsearch_treemaster_was_search_stopped(clt) && clt->obs->len > 0) {
    for (j = 0; j < clt->obs->len; j += 1) {
      QSearchTreeObserver *ob =
                   g_array_index(clt->obs,QSearchTreeObserver*,j);
                   ob->tree_search_done(clt->dm, answer, ob->user_data);
    }
  }
  return answer;
}
