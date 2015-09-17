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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <check.h>
#include <math.h>
#include "../src/qsearch/qsearch.h"

#include "dm.c"

START_TEST (test_qsearch_tree_nexus)
{
  LabeledMatrix *bdm;
  CompLearnRealCompressor *rc = complearn_environment_load_compressor_named("bzlib");
  fail_unless(rc != NULL);
  GString *cdm = real_compressor_compress(rc, g_string_new(dmstr));
  fail_unless(cdm != NULL);
  bdm = complearn_load_any_matrix(cdm);
  QSearchTree *tree = qsearch_tree_new(bdm->mat->size1);
  QLabeledTree *t = qsearch_tree_add_labels(tree, bdm);
  char *n = qsearch_tree_to_nexus(t);
  char *n2 = qsearch_tree_to_nexus_full(t, bdm);
  fail_unless (n != NULL);
  fail_unless (n2 != NULL);
  fail_unless (strlen(n) > 10);
  fail_unless (strlen(n2) > 20);
}
END_TEST

START_TEST (test_qsearch_walk_tree)
{
  guint32 i, k, j;
  for (i = 0; i < 30; i += 1) {
    QSearchTree *clt = qsearch_tree_new(i < 4 ? 4 : i);
    gsl_matrix *m = qsearch_tree_get_adjacency_matrix(clt);
    gsl_matrix_free(m);
    guint32 s = qsearch_tree_get_node_count(clt);
    k = rand() % s;
    for (j = 0; j < 1; j += 1) {
      GArray *res = qsearch_tree_walk_tree(clt, k, (gboolean) j);
      fail_if (res == NULL);
      fail_unless (res->len == s);
      fail_unless (g_array_index(res,guint32,0) == k);
      int sum = 0;
      for (k = 0; k < res->len; k += 1)
        sum += g_array_index(res, guint32, k);
      fail_unless (sum == ((s*(s-1))/2)); // all nodes once
      g_array_free(res, TRUE);
    }
    qsearch_tree_free(clt);
  }
}

END_TEST
START_TEST (test_qsearch_tree_mutation)
{
  guint32 i, l, k, q, z;
  for (i = 0; i < 30; i += 1) {
    QSearchTree *clt = qsearch_tree_new(i < 4 ? 4 : i);
    l = qsearch_tree_get_random_node(clt, NODE_TYPE_LEAF);
    k = qsearch_tree_get_random_node(clt, NODE_TYPE_KERNEL);
    q = qsearch_tree_get_random_node_but_not(clt, NODE_TYPE_KERNEL, k);
    z = qsearch_tree_get_random_node_but_not(clt, NODE_TYPE_ALL, q);

    fail_if (l == k);
    fail_if (q == k);
    fail_if (q == z);
    fail_unless (qsearch_tree_get_neighbor_count(clt, l) == 1);
    fail_unless (qsearch_tree_get_neighbor_count(clt, k) == 3);
    fail_unless (qsearch_tree_get_neighbor_count(clt, q) == 3);
    qsearch_tree_free(clt);
  }
}
END_TEST
void printar(GArray *ga)
{
  int i;
  for (i = 0; i < ga->len; i += 1)
    printf("%d ", g_array_index(ga, guint32, i));
  printf("\n");
}

START_TEST (test_qsearch_tree_scoring)
{
  QSearchTree *clt = qsearch_tree_new(14), *cclt;
  double sco, sco2, osco, osco2;
  gsl_matrix *dm;
  dm = gsl_matrix_alloc(qsearch_tree_get_leaf_node_count(clt),qsearch_tree_get_leaf_node_count(clt));
  int i, j;
  for (i = 0; i < dm->size1; i += 1)
    for (j = 0; j < dm->size2; j += 1)
      gsl_matrix_set(dm, i,j,fabs(sin(i+j)));
  sco =  qsearch_tree_score_tree(clt, dm);
  fail_unless( sco >= 0.0 && sco <= 1.0 , "score is %f", sco);
  osco  = qsearch_tree_calculate_order_cost(clt, dm);
  qsearch_tree_mutate_order_complex(clt);
  qsearch_tree_mutate_order_complex(clt);
  qsearch_tree_mutate_order_complex(clt);
  osco2  = qsearch_tree_calculate_order_cost(clt, dm);
  fail_unless( osco > 0.0 && osco < dm->size1 * dm->size2 * 1.2 );
  fail_unless( osco2 > 0.0 && osco2 < dm->size1 * dm->size2 * 1.2 );
  fail_if( osco == osco2 );
  cclt = qsearch_tree_clone(clt);
  fail_unless( cclt != NULL );
  sco2 = qsearch_tree_score_tree(cclt, dm);
  g_assert(sco == sco2);
  for (i = 0; i < 10; i += 1) {
    qsearch_tree_simple_mutation_subtree_transfer(clt);
    qsearch_tree_simple_mutation_subtree_interchange(clt);
    sco =  qsearch_tree_score_tree(clt, dm);
    g_assert(sco != sco2);
    if (cclt) {
      qsearch_tree_free(cclt);
    }
    cclt = qsearch_tree_find_better_tree(clt, i+1, dm);
    if (cclt) {
      double sco3 = qsearch_tree_score_tree(cclt, dm);
      fail_unless(sco3 > sco);
    }
    fail_unless (qsearch_tree_is_valid_tree(clt));
  }
  qsearch_tree_free(clt);
  qsearch_tree_free(cclt);
}
END_TEST

START_TEST (test_qsearch_tree_quartets)
{
  guint32 i, j, k, l;
  QSearchTree *clt = qsearch_tree_new(14);
  for (i = 0; i < qsearch_tree_get_node_count(clt); i += 1) {
    if (qsearch_tree_get_neighbor_count(clt, i) != 1)
      continue;
    for (j = i+1; j < qsearch_tree_get_node_count(clt); j += 1) {
      if (qsearch_tree_get_neighbor_count(clt, j) != 1)
        continue;
    for (k = j+1; k < qsearch_tree_get_node_count(clt); k += 1) {
      if (qsearch_tree_get_neighbor_count(clt, k) != 1)
        continue;
    for (l = k+1; l < qsearch_tree_get_node_count(clt); l += 1) {
      if (qsearch_tree_get_neighbor_count(clt, l) != 1)
        continue;
      gboolean x1, x2, x3;
      x1 = qsearch_tree_is_consistent_quartet(clt, i, j, k, l);
      x2 = qsearch_tree_is_consistent_quartet(clt, i, k, j, l);
      x3 = qsearch_tree_is_consistent_quartet(clt, i, l, j, k);
      fail_unless((x1 + x2 + x3) == 1);
    }
    }
    }
  }
  qsearch_tree_free(clt);
}
END_TEST

START_TEST (test_qsearch_tree_paths)
{
  guint32 i, j;
  QSearchTree *clt = qsearch_tree_new(14);
  for (i = 0; i < qsearch_tree_get_node_count(clt); i += 1)
    for (j = 0; j < qsearch_tree_get_node_count(clt); j += 1) {
      GArray *path = qsearch_tree_find_path(clt, i, j);
      int sz = path->len;
      if (i == j)
        fail_unless( sz == 1 );
      else
        fail_unless( sz > 1 );
      g_array_free(path, TRUE);
    }
  qsearch_tree_free(clt);
}
END_TEST

START_TEST (test_qsearch_tree)
{
  QSearchTree *clt = qsearch_tree_new(4);
  fail_unless (clt != NULL,
    "Cannot have NULL QSearchNeighborList");
  fail_unless( qsearch_tree_get_kernel_node_count(clt) +
  qsearch_tree_get_leaf_node_count(clt) ==
  qsearch_tree_get_node_count(clt), "Node counts don't add up" );
  fail_unless(qsearch_tree_is_standard_tree(clt),
      "Incorrect neighbor count for new tree.");
  qsearch_tree_free(clt);
}
END_TEST

Suite *
qsearch_tree_suite (void)
{
  Suite *s = suite_create ("tree");

  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_qsearch_tree);
  tcase_add_test (tc_core, test_qsearch_tree_paths);
  tcase_add_test (tc_core, test_qsearch_tree_quartets);
  tcase_add_test (tc_core, test_qsearch_tree_scoring);
  tcase_add_test (tc_core, test_qsearch_tree_mutation);
  tcase_add_test (tc_core, test_qsearch_walk_tree);
  tcase_add_test (tc_core, test_qsearch_tree_nexus);
  suite_add_tcase(s, tc_core);

  return s;
}
