
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
#include <string.h>
#include <assert.h>
#include "qsearch/qsearch.h"
#include <math.h>
#include <libintl.h>
#include <stdio.h>

#define _(O) gettext(O)

int main(int argc, char **argv)
{
  g_type_init();
  GString *matfile = complearn_read_whole_file("distmatrix.clb");
  GString *treedot = complearn_read_whole_file("treefile.dot");
  QSearchTree *qst;
  assert(matfile && "must have distmatrix.clb");
  assert(treedot && "must have treefile.dot");
  LabeledMatrix *lm = complearn_load_any_matrix(matfile);
  assert(lm && "must have good distmatrix.clb");
  int l = lm->mat->size1;
  int tot = 2*l-2;
  printf("Got matrix of size %d yielding total nodes of %d\n", l, tot);
  qst = qsearch_tree_new(lm->mat->size1);
  printf("Initial score: %f\n", qsearch_tree_score_tree(qst, lm->mat));
  gdouble presumed_score = qsearch_tree_read_from_dot(qst, treedot, lm);
  gdouble real_score = qsearch_tree_score_tree(qst, lm->mat);
  gdouble dscore = presumed_score - real_score;
  printf("Real: %f   Presumed: %f   delta: %f\n", real_score, presumed_score, dscore);
  GString *newdot = qsearch_tree_to_dot(qsearch_tree_add_labels(qst, lm));
  printf("%s", newdot->str);
  //MakeTreeResult *mtr;
  //mtr = qsearch_maketree_process_options(argv);
/*  qsearch_maketree_write_tree_file(mtr); */
  /*qsearch_tree_free(tree);
  gsl_matrix_free(dm); */
  return 0;
}
