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
#include <stdlib.h>
#include <check.h>
#include <math.h>
#include "../src/qsearch/qsearch.h"

static double final_score;

static void endfunc(gsl_matrix *mat, QSearchTree *tree, void *udata)
{
  final_score = qsearch_tree_score_tree(tree, mat);
}

START_TEST (test_qsearch_treemaster_basic)
{
  gsl_matrix *dm;
  int k;
  for (k = 4; k < 6; k += 1) {
    final_score = -1;
    dm = gsl_matrix_alloc(k,k);
    int i, j;
    for (i = 0; i < dm->size1; i += 1)
      for (j = 0; j < dm->size2; j += 1)
        gsl_matrix_set(dm, i,j,fabs(sin(i+j)));
    QSearchTreeMaster *cltm = qsearch_treemaster_new(dm);
    qsearch_treemaster_add_observer(cltm, NULL, NULL, endfunc, NULL);
    QSearchTree *tree = qsearch_treemaster_find_best_tree(cltm);
    if (k == 4)
      fail_unless(qsearch_tree_score_tree(tree,dm) > 0.9999);
    fail_unless(fabs(qsearch_tree_score_tree(tree,dm) -  final_score) < 1e-14);
    qsearch_treemaster_free(cltm);
    qsearch_tree_free(tree);
    gsl_matrix_free(dm);
  }
}
END_TEST

Suite *
qsearch_treemaster_suite (void)
{
  Suite *s = suite_create ("treemaster");

  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_qsearch_treemaster_basic);
  suite_add_tcase(s, tc_core);

  return s;
}
