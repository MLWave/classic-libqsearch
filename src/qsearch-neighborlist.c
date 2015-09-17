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

#include <stdlib.h>
#include <stdio.h>
#include <glib.h>

#include "qsearch.h"

void qsearch_init(void)
{
}

gint32 qsearch_neighborlist_find_index(QSearchNeighborList *cln, guint32 w)
{
  int i;
  for (i = 0; i < cln->n->len; i += 1)
    if (g_array_index(cln->n, guint32, i) == w)
      return i;
  return -1;
}

void qsearch_copy_intarray(GArray *d, GArray *s)
{
  g_array_set_size(d, 0);
  if (s->len > 0)
    g_array_append_vals(d, s->data, s->len);
}

void qsearch_neighborlist_copy_from(QSearchNeighborList *d, QSearchNeighborList *s) {
  qsearch_copy_intarray(d->n, s->n);
}

QSearchNeighborList *qsearch_neighborlist_clone(QSearchNeighborList *clnl)
{
  QSearchNeighborList *cl2 = new_qsearch_neighborlist();
  qsearch_neighborlist_copy_from(cl2, clnl);
  return cl2;
}

gboolean qsearch_neighborlist_has_neighbor(QSearchNeighborList *cln, guint32 w)
{
  return qsearch_neighborlist_find_index(cln, w) != -1;
}

void qsearch_neighborlist_remove_neighbor(QSearchNeighborList *cln, guint32 w)
{
  g_assert(qsearch_neighborlist_has_neighbor(cln, w) == TRUE);
  g_array_remove_index_fast(cln->n, qsearch_neighborlist_find_index(cln,w));
}

void qsearch_neighborlist_add_neighbor(QSearchNeighborList *cln, guint32 w)
{
  g_assert(qsearch_neighborlist_has_neighbor(cln, w) == FALSE);
  g_array_append_vals(cln->n, &w, 1);
}

void free_qsearch_neighborlist(QSearchNeighborList *cl)
{
  if (cl == NULL)
    return;
  g_array_free(cl->n, TRUE);
  g_free(cl);
}

QSearchNeighborList *new_qsearch_neighborlist(void)
{
  QSearchNeighborList *clnl = calloc(sizeof(*clnl), 1);
  clnl->n = g_array_new(FALSE, TRUE, sizeof(guint32));
  return clnl;
}

