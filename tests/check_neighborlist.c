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
#include "../src/qsearch/qsearch.h"

START_TEST (test_qsearch_neighborlist)
{
  QSearchNeighborList *clnl = new_qsearch_neighborlist();
  QSearchNeighborList *cl2;
  fail_unless (clnl != NULL,
    "Cannot have NULL QSearchNeighborList");
  int n1 = 3, n2 = 7;
  fail_if ( qsearch_neighborlist_has_neighbor(clnl, n1) );
  fail_if ( qsearch_neighborlist_has_neighbor(clnl, n2) );
  qsearch_neighborlist_add_neighbor(clnl, n1);
  fail_unless ( qsearch_neighborlist_has_neighbor(clnl, n1) );
  fail_if ( qsearch_neighborlist_has_neighbor(clnl, n2) );
  qsearch_neighborlist_remove_neighbor(clnl, n1);
  fail_if ( qsearch_neighborlist_has_neighbor(clnl, n1) );
  fail_if ( qsearch_neighborlist_has_neighbor(clnl, n2) );
  qsearch_neighborlist_add_neighbor(clnl, n2);
  qsearch_neighborlist_add_neighbor(clnl, n1);
  cl2 = qsearch_neighborlist_clone(clnl);
  fail_unless ( qsearch_neighborlist_has_neighbor(clnl, n1) );
  fail_unless ( qsearch_neighborlist_has_neighbor(clnl, n2) );
  fail_unless ( qsearch_neighborlist_has_neighbor(cl2, n1) );
  fail_unless ( qsearch_neighborlist_has_neighbor(cl2, n2) );
  free_qsearch_neighborlist ( clnl );
  free_qsearch_neighborlist ( cl2 );
}
END_TEST

Suite *
qsearch_neighborlist_suite (void)
{
  Suite *s = suite_create ("neighborlist");

  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_qsearch_neighborlist);
  suite_add_tcase(s, tc_core);

  return s;
}
