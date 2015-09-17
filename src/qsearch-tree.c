#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <glib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include "qsearch.h"

#include <sys/time.h>

#define CHUNKSIZE 1

#define NODE_FLAG_FRINGE       0x01
#define NODE_FLAG_DONE         0x02
#define NODE_FLAG_NEXTFRINGE   0x04
#define NODE_FLAG_QUARTETINT   0x08
#define NODE_FLAG_ISWALKED     0x10
#define NODE_FLAG_ISFLIPPED    0x20

typedef struct _QSearchPath QSearchPath;

void qsearch_tree_freshen_spm(QSearchTree *clt);
void qsearch_tree_find_path_fast(QSearchTree *clt, guint32 a, guint32 b,
GArray *dest);

void qsearch_tree_mutate_order_simple(QSearchTree *clt)
{
  int k = qsearch_tree_get_random_node(clt, NODE_TYPE_KERNEL);
  g_array_index(clt->nodeflags, guint32, k) ^= NODE_FLAG_ISFLIPPED;
//  printf("made node %d flip with %d neighbors\n", k, qsearch_tree_get_neighbor_count(clt,k));
  clt->ms.last_order_simple_mutations += 1;
}

void qsearch_tree_mutate_order_complex(QSearchTree *clt)
{
  clt->ms.last_order_simple_mutations = 0;
  do {
    qsearch_tree_mutate_order_simple(clt);
  } while (rand() % 2 == 0);
  clt->ms.total_order_simple_mutations += clt->ms.last_order_simple_mutations;
  clt->ms.total_order_complex_mutations += 1;
}

double qsearch_tree_calculate_order_cost(const QSearchTree *clt,const gsl_matrix *dm);
void qsearch_copy_neighborlist(QSearchNeighborList *d, QSearchNeighborList *s)
{
  qsearch_neighborlist_copy_from(d, s);
}

GArray *qsearch_tree_walk_filtered(const QSearchTree *clt, guint32 fromwhere,
                                   guint32 filterType)
{
  int i;
  GArray *m = qsearch_tree_walk_tree_dfs(clt, fromwhere);
  GArray *r = g_array_new(FALSE, TRUE, sizeof(guint32));
  for (i = 0; i < m->len; i += 1) {
    int n = g_array_index(m, guint32, i);
    int nc = qsearch_tree_get_neighbor_count(clt, n);
    if (((nc == 1) && (filterType & NODE_TYPE_LEAF)) ||
        ((nc == 3) && (filterType & NODE_TYPE_KERNEL)))
      g_array_append_val(r, n);
  }
  g_array_free(m, TRUE);
  return r;
}

GArray *qsearch_tree_walk_tree_bfs(const QSearchTree *clt, guint32 fromwhere)
{
  return qsearch_tree_walk_tree(clt, fromwhere, TRUE);
}

GArray *qsearch_tree_walk_tree_dfs(const QSearchTree *clt, guint32 fromwhere)
{
  return qsearch_tree_walk_tree(clt, fromwhere, FALSE);
}

int qsearch_tree_get_mutation_distribution_sample(const QSearchTree *clt)
{
  const int MAXMUT = 80;
  static gsl_ran_discrete_t *d;
  static gsl_rng *r;
  if (d == NULL) {
    int i;
    double curp;
    double *p = calloc(sizeof(*d), MAXMUT);
    double k;
    for (i = 0; i < MAXMUT; i += 1) {
      k = i + 4; /* to make single-mutations somewhat less common */
      curp = 1.0 / (k * (log(k) / log(2.0)) * (log(k)/log(2.0)));
      p[i] = curp;
    }
    d = gsl_ran_discrete_preproc(MAXMUT, p);
    r = gsl_rng_alloc(gsl_rng_taus);
    free(p);
  }
  return gsl_ran_discrete(r, d) + 1;
}

void qsearch_tree_complex_mutation(QSearchTree *clt)
{
  clt->ms.last_simple_mutations = 0;
  int i, totmuts = qsearch_tree_get_mutation_distribution_sample(clt);
  for (i = 0; i < totmuts; i += 1)
    qsearch_tree_simple_mutation(clt);
  clt->ms.total_simple_mutations += clt->ms.last_simple_mutations;
  clt->ms.total_complex_mutations += 1;
}

void qsearch_tree_simple_mutation(QSearchTree *clt)
{
  gboolean hm = FALSE;
  int i;
  do {
    i = rand() % 3;
    switch (i) {
      case 0: qsearch_tree_simple_mutation_leaf_swap(clt); hm = TRUE; break;
      case 1: if (qsearch_tree_can_subtree_transfer(clt)) {
  qsearch_tree_simple_mutation_subtree_transfer(clt); hm = TRUE; }
  break;
      case 2: if (qsearch_tree_can_subtree_interchange(clt)) {
  qsearch_tree_simple_mutation_subtree_interchange(clt); hm = TRUE; }
  break;
    }
  } while (!hm);
}

gboolean qsearch_tree_is_valid_tree(QSearchTree *clt)
{
  int i, j;
  g_assert(clt != NULL);
  g_assert(clt->total_node_count > 3);
  for (i = 0; i < clt->leaf_placement->len; i += 1) {
    int nc = qsearch_tree_get_neighbor_count(clt, g_array_index(clt->leaf_placement, guint32, i));
    if (nc != 1)
      return FALSE;
  }
  for (i = 0; i < clt->total_node_count; i += 1) {
    int nc = qsearch_tree_get_neighbor_count(clt, i);
    if (nc != 1 && nc != 3)
      return FALSE;
  }
  for (i = 0; i < clt->total_node_count; i += 1)
    for (j = 0; j < clt->total_node_count; j += 1) {
      int pl1 = qsearch_tree_find_path_length(clt, i, j);
      int pl2 = qsearch_tree_find_path_length(clt, j, i);
      if (pl1 != pl2)
        return FALSE;
      if (pl1 < 1 || pl1 > clt->total_node_count)
        return FALSE;
    }
  return TRUE;
}

guint32 qsearch_tree_get_random_node(const QSearchTree *clt, guint32 what_kind)
{
  guint32 result;
  guint32 n;
  g_assert(what_kind == NODE_TYPE_LEAF || what_kind == NODE_TYPE_KERNEL ||
           what_kind == NODE_TYPE_ALL);
  do {
    result = rand() % (clt->total_node_count);
    n = qsearch_tree_get_neighbor_count(clt, result);
  } while ((what_kind & (n == 1 ? NODE_TYPE_LEAF : NODE_TYPE_KERNEL)) == 0);
  return result;
}

GArray *qsearch_tree_get_neighbors(const QSearchTree *clt, guint32 who)
{
  GArray *result = g_array_new(FALSE, TRUE, sizeof(guint32));
  int i;
  for (i = 0; i < clt->total_node_count; i += 1)
    if (qsearch_tree_is_connected(clt, i, who)) {
      guint32 v = i;
      g_array_append_val(result, v);
    }
  return result;
}

guint32 qsearch_tree_get_random_neighbor(const QSearchTree *clt, guint32 who)
{
  guint32 result;
  GArray *g = qsearch_tree_get_neighbors(clt, who);
  result = g_array_index(g, guint32, rand() % g->len);
  g_array_free(g, TRUE);
  return result;
}

guint32 qsearch_tree_get_random_node_but_not(const QSearchTree *clt, guint32 what_kind, guint32 but_not)
{
  guint32 result;
  do {
    result = qsearch_tree_get_random_node(clt, what_kind);
  } while (result == but_not);
  return result;
}

void qsearch_tree_simple_mutation_leaf_swap(QSearchTree *clt)
{
  int l1, l2, tmp;
  l1 = qsearch_tree_get_random_node(clt, NODE_TYPE_LEAF);
  l2 = qsearch_tree_get_random_node_but_not(clt, NODE_TYPE_LEAF, l1);
  tmp = g_array_index(clt->leaf_placement, guint32, l1);
  g_array_index(clt->leaf_placement, guint32, l1) =
  g_array_index(clt->leaf_placement, guint32, l2);
  g_array_index(clt->leaf_placement, guint32, l2) = tmp;
  clt->f_score_good = FALSE;
  clt->ms.last_simple_mutations += 1;
}

gboolean qsearch_tree_can_subtree_transfer(const QSearchTree *clt)
{
  return (clt->total_node_count >= 9);
}

gboolean qsearch_tree_can_subtree_interchange(const QSearchTree *clt)
{
  return (clt->total_node_count >= 11);
}

void qsearch_tree_simple_mutation_subtree_interchange(QSearchTree *clt)
{
  guint32 k1, k2, n1, n2;
  g_assert(qsearch_tree_can_subtree_interchange(clt));
  do {
    k1 = qsearch_tree_get_random_node(clt, NODE_TYPE_KERNEL);
    k2 = qsearch_tree_get_random_node_but_not(clt, NODE_TYPE_KERNEL, k1);
  } while (qsearch_tree_find_path_length(clt, k1, k2) <= 3);
  GArray *path = qsearch_tree_find_path(clt, k1, k2);
  n1 = g_array_index(path, guint32, 1);
  n2 = g_array_index(path, guint32, path->len-2);
  g_array_free(path, TRUE);
  qsearch_tree_disconnect(clt, n1, k1);
  qsearch_tree_disconnect(clt, n2, k2);
  qsearch_tree_connect(clt, n1, k2);
  qsearch_tree_connect(clt, n2, k1);
  clt->ms.last_simple_mutations += 1;
}

void qsearch_tree_simple_mutation_subtree_transfer(QSearchTree *clt)
{
  g_assert(qsearch_tree_can_subtree_transfer(clt));
  guint32 k1, k2, i1, m1, m2, m3;
  do {
    k1 = qsearch_tree_get_random_node(clt, NODE_TYPE_ALL);
    k2 = qsearch_tree_get_random_node_but_not(clt, NODE_TYPE_KERNEL, k1);
  } while (qsearch_tree_find_path_length(clt, k1, k2) <= 2);
  GArray *path = qsearch_tree_find_path(clt, k1, k2);
  i1 = g_array_index(path, guint32, 1);
  qsearch_tree_disconnect(clt, k1, i1);
  GArray *neighbors = qsearch_tree_get_neighbors(clt, i1);
  do {
    m3 = qsearch_tree_get_random_neighbor(clt, k2);
  } while (m3 == g_array_index(path, guint32, path->len-2));
  m1 = g_array_index(neighbors, guint32, 0);
  m2 = g_array_index(neighbors, guint32, 1);
  qsearch_tree_disconnect(clt, m1, i1);
  qsearch_tree_disconnect(clt, m2, i1);
  qsearch_tree_disconnect(clt, m3, k2);
  qsearch_tree_connect(clt, m1, m2);
  qsearch_tree_connect(clt, k2, i1);
  qsearch_tree_connect(clt, m3, i1);
  qsearch_tree_connect(clt, k1, i1);
  g_array_free(path, TRUE);
  g_array_free(neighbors, TRUE);
  clt->ms.last_simple_mutations += 1;
}

QSearchTree *qsearch_tree_clone(const QSearchTree *clt)
{
  int i;
  QSearchTree *nclt = qsearch_tree_new((clt->total_node_count + 2) / 2);
  nclt->must_recalculate_paths = TRUE;
  nclt->f_score_good = clt->f_score_good;
  nclt->score = clt->score;
  nclt->ms = clt->ms;
  nclt->ms.total_clonings += 1;
  nclt->dist_calculated = clt->dist_calculated;
  nclt->dist_min = clt->dist_min;
  nclt->dist_max = clt->dist_max;
  qsearch_copy_intarray(nclt->leaf_placement, clt->leaf_placement);
  for (i = 0; i < clt->total_node_count; i += 1) {
    qsearch_copy_neighborlist(g_ptr_array_index(nclt->n, i),
                  g_ptr_array_index( clt->n, i));
  }
  return nclt;
}

void calc_min_max(QSearchTree *clt,const gsl_matrix *dm) {

    int i,j,k,l;
    double minscore, maxscore;

    clt->dist_min = 0;
    clt->dist_max = 0;
  for (i = 0; i < clt->leaf_placement->len; i += 1)
    for (j = i+1; j < clt->leaf_placement->len; j += 1)
      for (k = j+1; k < clt->leaf_placement->len; k += 1)
        for (l = k+1; l < clt->leaf_placement->len; l += 1) {
          double c1, c2, c3;
          c1  = gsl_matrix_get(dm, i, j) + gsl_matrix_get(dm, k, l);
          c2  = gsl_matrix_get(dm, i, k) + gsl_matrix_get(dm, j, l);
          c3  = gsl_matrix_get(dm, i, l) + gsl_matrix_get(dm, j, k);
          minscore = c1; maxscore = c1;
          if (c2 < minscore) minscore = c2;
          if (c3 < minscore) minscore = c3;
          if (c2 > maxscore) maxscore = c2;
          if (c3 > maxscore) maxscore = c3;
          clt->dist_min += minscore; clt->dist_max += maxscore;
    }
    
}


struct _ConnectedNode {
    int done;
    int connections[3];
    int leaf_count[3];
    char *node_branch; // pointing in the direction where to find a node 
};

typedef struct _ConnectedNode ConnectedNode;

static int next_node(ConnectedNode *map, int from, int to) {
    return map[from].connections[ (int) map[from].node_branch[to] ];
}

inline static int find_branch(int connections[3], int to) {
    if (connections[0] == to) return 0;
    if (connections[1] == to) return 1;
    if (connections[2] == to) return 2;
    g_assert(0); // 
    return -1;
}

ConnectedNode *init_node_map(QSearchTree *clt) {
    int i,j; 
    int node_count = clt->total_node_count;
    int leaf_count = (node_count + 2)/2;

    ConnectedNode* map = malloc( node_count * sizeof(ConnectedNode));
       
    for (i = 0; i < node_count; ++i) {
        map[i].done = 0;
        for (j = 0; j < 3; ++j) {
            map[i].connections[j] = -1;
            map[i].leaf_count[j] = 0;
        }
        
        map[i].node_branch = malloc( node_count * sizeof(char));
        if (i < leaf_count) {
            for (j=0; j < node_count; ++j) map[i].node_branch[j] = 0;
            map[i].leaf_count[0] = leaf_count - 1;
            map[i].done = 1;
        } else {
            for (j = 0; j < node_count; ++j) map[i].node_branch[j] = -1;
        } 
    }
     
    // set connections (construct tree)
    for (i = 0; i < node_count; ++i) {

        QSearchNeighborList *cln = g_ptr_array_index(clt->n, i);
        // add to connected nodes
        for (j = 0; j < cln->n->len; ++j) {
            int node = g_array_index(cln->n, guint32, j);
        
            // find unfilled branch
            int branch = find_branch(map[node].connections, -1);
            map[node].connections[branch] = i; // set connection
            
            if (i < leaf_count) {
                map[node].leaf_count[branch] = 1; // set leaf
            }
            map[node].node_branch[i] = branch; // leaf can be found in branch

            // set connection back to this node
            branch = find_branch(map[i].connections, -1);
            map[i].connections[branch] = node;
            map[i].node_branch[node] = branch;
        }
    }

    // every iteration, this loop progresses by at least finishing one node. Worst case is n^3 (for 'linear' trees)
    for (;;) {
        for (i = leaf_count; i < node_count; ++i) {
            
            if (map[i].done) continue; // nothing to be done here
                        
            for (j = 0; j < 3; ++j) {
                int node = map[i].connections[j];
                if (map[node].done) continue; 

                int branch = find_branch(map[node].connections, i);
                if (map[node].leaf_count[branch] == 0) {
                    int first = (3 + j-1) % 3;
                    int second = (j + 1) % 3;
                    if (map[i].leaf_count[first] != 0 && map[i].leaf_count[second] != 0) {
                        map[node].leaf_count[branch] = map[i].leaf_count[first] + map[i].leaf_count[second];
                        
                        // set the nodes in the other
                        int k;
                        for (k = 0; k < node_count; ++k) {
                            // node present in one of the two branches pointing away from this
                            int node_present = k==i || map[i].node_branch[k] == first || map[i].node_branch[k] == second;
                            
                            if (node_present && map[node].node_branch[k] == -1) map[node].node_branch[k] = branch;
                        }
                    }
                }
            }
        }
        
        int all_done = 1;
        for (i = leaf_count; i < node_count; ++i) {
            if (map[i].done) continue;
            
            // we are done when all entries are set for this node, AND there are no pending assignments to neighbours

            // are all entries set for this node?
            int entry_completed = 1;
            for (j = 0; j < 3; ++j) {
                if (map[i].leaf_count[j] == 0) { 
                    entry_completed = 0;
                    break;
                }
            }

            if (entry_completed) {
                // are there pending assignments? (we could do the assignments here, but that would duplicate code)
                int done = 1;
                for (j = 0; j < 3; ++j) {
                    int node = map[i].connections[j];
                    if (node < leaf_count) continue;
                    // find connection
                    int branch = find_branch(map[node].connections, i);
                    if (map[node].leaf_count[branch] == 0) {
                        done = 0;
                        break;
                    }
                }
                
                map[i].done = done;
                if (!done) all_done = 0;
            } else {
                all_done = 0;
            }
        }
        
        // everything is set and accounted for
        if (all_done) break; 

    }
    
    
    return map;
}

void free_map(ConnectedNode *map, int node_count) { 
    int i;
    for (i = 0; i < node_count; ++i) if (map[i].node_branch) free(map[i].node_branch);
    free(map);
}
    
void print_map(ConnectedNode *map, int node_count) {
    
    int i,j;
    int leaf_count = (node_count+2)/2;

    printf("Nodes %d, leafs %d\n", node_count, leaf_count);
    for (i = 0; i < node_count; ++i) {
        printf("%d: connections: %d %d %d, counts: %d %d %d ", i, 
                map[i].connections[0],
                map[i].connections[1],
                map[i].connections[2],
                map[i].leaf_count[0],
                map[i].leaf_count[1],
                map[i].leaf_count[2]);
            printf("Leafs :");
            for (j = 0; j < node_count; ++j) {
                printf(" %d", i==j?5:map[i].node_branch[j]);
            }
        printf("\n");
    }
}


double qsearch_tree_score_map(ConnectedNode *map, int node_count, const gsl_matrix *dm) {
    
    // run over all pairs
    int leaf_count = (node_count + 2)/2;
    int i,j;

    /* loop over internal nodes. As we have know where leafs are, we can now, for each internal node
     * calculate the number of consistent pairs that participate in the sum
     *
     * Suppose we have 8 leafs, and an internal node structured as:
     *
     * branch1 points to leafs 0 3 6 7
     * branch2 points to leafs 2 5
     * branch3 points to leafs 1 4
     *
     * We can, for each branch, calculate the distances that needs to be added to the overall sum 
     *
     * for branch1, we can calculate that there are 6 pairs that are embedded in the tree that form quartets 
     * with all pairs from branch2 and branch3. Thus npairs = 6, and we compute:
     *
     * sum += 6 * ( d(2,1) + d(2,4) + d(5,1) + d(5, 4) )
     *
     * We do not calculate d(2,5), nor d(1,4), as these pairs are taken care of by the internal nodes that
     * split them in different branches
     *
     * By doing the same for branch2 and branch 4, we achieve an overall n^3 algorithm. 
     *
     */
    
    long long **tmpmat = calloc(node_count, sizeof(void *));

    int node;
    int branch, n, npairs, first, second;
//#pragma omp schedule(dynamic, CHUNKSIZE) private(branch,n,npairs,first,second,i,ni,j,nj)
{
    for (node = leaf_count; node < node_count; ++node) {
      tmpmat[node] = calloc( leaf_count * leaf_count, sizeof(long long));
        for (branch = 0; branch < 3; ++branch) {
            n = map[node].leaf_count[branch];
            
            if (n > 1) { // we need to accumulate the data for all pairs
                npairs = n * (n-1) / 2; // number of pairs
                
                first = (3 + branch - 1) % 3;
                second = (branch + 1) % 3;
                
                for (i = 0; i < leaf_count; ++i) {
                    //int ni = g_array_index(lp, guint32, i);                    
                    if (map[node].node_branch[i] != first) continue; // this leaf is in the wrong branch
                    
                    for (j = 0; j < leaf_count; ++j) {
                        //int nj = g_array_index(lp, guint32, j);
                        if ( map[node].node_branch[j] != second) continue; // this leaf is in the wrong branch
                        
                        tmpmat[node][i + leaf_count*j] += npairs;
                    }
                }
            }
        }
    }
}
    
    double sum = 0;
    for (node = leaf_count; node < node_count; ++node) {
    for (i = 0; i < leaf_count; ++i) {
        for (j = 0; j < leaf_count; ++j) {
            sum += tmpmat[node][i + j*leaf_count] * gsl_matrix_get(dm, i, j);
        }
      }
    }

    for (node = leaf_count; node < node_count; ++node)
      free(tmpmat[node]);
    free(tmpmat);
    return sum;
}

static double dist(ConnectedNode *cnm, int leaf_count, int A, int B, const gsl_matrix *dm) {
    int a2b = cnm[A].node_branch[B]; 
    int b2a = cnm[B].node_branch[A];
    
    double sum = 0;
    int i,j;
    for (i = 0; i < leaf_count; ++i) {
        
        if (cnm[A].node_branch[i] == a2b && i != A) continue;
        
        for (j = 0; j < leaf_count; ++j) {
            
            if ((j == i || cnm[B].node_branch[j] == b2a ) && j != B) continue;
            
            //printf("%d to %d, nodes %d %d, a2b %d b2a %d: cnma %d ncmb %d:  %f\n", A, B, i, j, a2b, b2a, cnm[A].node_branch[i], cnm[B].node_branch[j], gsl_matrix_get(dm, i, j));
            sum += gsl_matrix_get(dm, i, j);
        }
    }
    
    return sum;
}

// return all leafs that do are in branches *not* pointing to 'to'
static double count_leafs(ConnectedNode *cnm, int leaf_count, int from, int to) {
    int branch = cnm[from].node_branch[to];
    double nleafs = 0;
    
    if (from < leaf_count) nleafs++;

    int i;
    for (i = 0; i < leaf_count; ++i) {
        if (cnm[from].node_branch[i] != branch) nleafs++;
    }
    return nleafs;
}

inline static double npairs(double n) { return n * (n-1) / 2; }

static double predict_quartet_swap(ConnectedNode *cnm, int leaf_count, int A, int B, const gsl_matrix *dm) {
    
    int interior1 = next_node(cnm, A, B);
    int interior2 = next_node(cnm, B, A);
    
    if (interior1 == interior2 || interior1 == B) {
        printf("node mismatch\n");
        return 0.0; // swapping will have no effect
    }
    
    // check if interior1 and interior2 are connected
    if (cnm[interior1].connections[ (int) cnm[interior1].node_branch[interior2]] != interior2 ) {
        printf("node mismatch2");
        return 0.0;
    }
    
    // find the branch of interior1 which is not covered by A or B     
    int branchC = 3 - cnm[interior1].node_branch[A] - cnm[interior1].node_branch[B];
    int C = cnm[interior1].connections[branchC];

    // find the branch of interior2 which is not covered by A or B     
    int branchD = 3 - cnm[interior2].node_branch[A] - cnm[interior2].node_branch[B];
    int D = cnm[interior2].connections[branchD];
    
    // now we have nodes A, B, C, and D, swap cost can be determined by only looking
    // at the effect of changing <A,C>,<B,D> -> <B,C>,<A,D>
    // all nodes in the middle are irrelevant (will cancel out in the computation)
    
    double dac = dist(cnm, leaf_count, A, C, dm);
    double dad = dist(cnm, leaf_count, A, D, dm);
    double dbc = dist(cnm, leaf_count, B, C, dm);
    double dbd = dist(cnm, leaf_count, B, D, dm);
    
    double na = count_leafs(cnm, leaf_count, A, interior1);
    double nb = count_leafs(cnm, leaf_count, B, interior1);
    double nc = count_leafs(cnm, leaf_count, C, interior2);
    double nd = count_leafs(cnm, leaf_count, D, interior2);
    

    double d = 0;
    d += npairs(na) * (dbc - dbd);
    d += npairs(nb) * (dad - dac);
    d += npairs(nc) * (dad - dbd);
    d += npairs(nd) * (dbc - dac);
    
    d += npairs(na + nc) * dbd;
    d += npairs(nb + nd) * dac;

    d -= npairs(nb + nc) * dad;
    d -= npairs(na + nd) * dbc;

    //printf("%d %d %d %d : %d %d %d %d : %f %f %f %f %f %f %f %f \n", A, B, C, D, branchC, branchD, interior1, interior2, na, nb, nc, nd, dac, dad, dbc, dbd);
    /* 
    int i,j; 
    for (i = 0; i < leaf_count; ++i) {
        for (j = 0; j < leaf_count; ++j) {
            printf("%f ", gsl_matrix_get(dm, i, j));
        }
        printf("\n");
    }*/
     
    return -d;
}

static void swap_nodes(ConnectedNode *ctm, int node_count, int A, int B) {
   
   g_assert(A != B);
    
   int interiorA = ctm[A].connections[ (int) ctm[A].node_branch[B] ];
   int interiorB = ctm[B].connections[ (int) ctm[B].node_branch[A] ];
   
   g_assert(interiorA != interiorB);
   g_assert(interiorA != B);
    
   int unchangedABranch = ctm[A].node_branch[interiorA];
   int unchangedBBranch = ctm[B].node_branch[interiorB];
   
   int leaf_count = (node_count + 2)/2;
    
   int countA = leaf_count - ctm[A].leaf_count[unchangedABranch];
   int countB = leaf_count - ctm[B].leaf_count[unchangedBBranch];
    
   //printf("Counts: %d %d, branch counts: %d\n", countA, countB, ctm[A].leaf_count[unchangedABranch]);

   int interiorABranch = ctm[interiorA].node_branch[A];
   int interiorBBranch = ctm[interiorB].node_branch[B];
     
   // loop over internal nodes, swap node_branches from A <-> B
   int node = interiorA;
   while (node != B) {
        int aBranch = ctm[node].node_branch[A];
        int bBranch = ctm[node].node_branch[B];
        
        ctm[node].leaf_count[aBranch] += countB - countA;
        ctm[node].leaf_count[bBranch] += countA - countB;

        int i;
        for (i = 0; i < node_count; ++i) {
            
            if (i == A || ctm[A].node_branch[i] != unchangedABranch) { // this is a node that is a subtree of A
                g_assert(ctm[node].node_branch[i] == aBranch);
                ctm[node].node_branch[i] = bBranch;
        
            } 
            
            if (i == B || ctm[B].node_branch[i] != unchangedBBranch) { // this is a node that is a subtree of B
                    g_assert(ctm[node].node_branch[i] == bBranch);
                    ctm[node].node_branch[i] = aBranch;
            }

        }

        node = ctm[node].connections[ bBranch ];
   }

   // swap directions for A and B
   ctm[interiorA].connections[ interiorABranch ] = B;
   ctm[interiorB].connections[ interiorBBranch ] = A;
    
   // swap connections for A and b
   ctm[A].connections[ unchangedABranch ] = interiorB;
   ctm[B].connections[ unchangedBBranch ] = interiorA;

    // done

}

double qsearch_optimize_tree(QSearchTree *clt, const gsl_matrix *dm) {
    ConnectedNode *map = init_node_map(clt);
    
    // run over all pairs
    int node_count = clt->total_node_count;
    int leaf_count = (node_count + 2)/2;
    int i,j,k,l;
    
    double score = qsearch_tree_score_map(map, node_count, dm);
    printf("Init score %f\n", (clt->dist_max - score) / (clt->dist_max - clt->dist_min) );
    
    int done = 1;
    // now start swapping
    do {
RESTART:
        done = 1;
        for (i = leaf_count; i < node_count; ++i) {
            
            int interior1 = i;

            for (j = 0; j < 3; ++j) {
                
                int swap1 = map[i].connections[j];

                for (k = j+1; k < 3; ++k) {
                    
                    int interior2 = map[i].connections[k];

                    if (interior2 < leaf_count) continue;

                    for (l = 0; l < 3; ++l) {

                        int swap2 = map[interior2].connections[l];

                        if (swap2 == interior1) continue; // pointing back
                        //printf("Trying swap\n"); 
                        double swap_delta = predict_quartet_swap(map, leaf_count, swap1, swap2, dm);
                                     
                        if (swap_delta > 1e-6) {
                            //printf("Enter swap %d %d\n", swap1, swap2);
                            //print_map(map,node_count);
                            swap_nodes(map, node_count, swap1, swap2);
                            //print_map(map, node_count);
                            
                            //printf("Exit swap\n");
                            double newscore = qsearch_tree_score_map(map, node_count, dm);
                            printf("%f + %f = %f = %f\n", score, swap_delta, score+swap_delta, newscore);
                            
                            if (fabs(score + swap_delta - newscore) > 1e-6) {
                                printf("Internal error");
                                exit(0);
                            }
                            score = newscore;
                            goto RESTART;
                        }
                    }
                }
           }
        }
    } while (!done);
    
    printf("New score %f\n", (clt->dist_max - score) / (clt->dist_max - clt->dist_min) );
    /*
    // write out resulting tree in clt
    for (i = 0; i < leaf_count; ++i) {
        g_array_index(clt->leaf_placement, guint32, i) = i;
    }
    
    for (i = 0; i < node_count; ++i) {
        
       QSearchNeighborList *lst = g_ptr_array_index(clt->n, i);
       GArray *n = lst->n;
       // lst->n is GArray  
        g_array_set_size(n, 0);
           
        for (j = 0; j < 3; ++j) {
            int con = map[i].connections[j];
            if (con <= i) continue; // no need to write
            g_array_append_val(n, con);
        }
    }

    clt->must_recalculate_paths = TRUE;
    clt->f_score_good = FALSE;
    */
    free_map(map, node_count);
    return (clt->dist_max - score) / (clt->dist_max - clt->dist_min);
}

double qsearch_tree_score_tree_fast_v2(QSearchTree *clt, const gsl_matrix *dm) {
    
    ConnectedNode *map = init_node_map(clt);
    
    // run over all pairs
    int node_count = clt->total_node_count;
    int leaf_count = (node_count + 2)/2;
    int i,j;

    double sum = 0.0;
    
    /* loop over internal nodes. As we have know where leafs are, we can now, for each internal node
     * calculate the number of consistent pairs that participate in the sum
     *
     * Suppose we have 8 leafs, and an internal node structured as:
     *
     * branch1 points to leafs 0 3 6 7
     * branch2 points to leafs 2 5
     * branch3 points to leafs 1 4
     *
     * We can, for each branch, calculate the distances that needs to be added to the overall sum 
     *
     * for branch1, we can calculate that there are 6 pairs that are embedded in the tree that form quartets 
     * with all pairs from branch2 and branch3. Thus npairs = 6, and we compute:
     *
     * sum += 6 * ( d(2,1) + d(2,4) + d(5,1) + d(5, 4) )
     *
     * We do not calculate d(2,5), nor d(1,4), as these pairs are taken care of by the internal nodes that
     * split them in different branches
     *
     * By doing the same for branch2 and branch 4, we achieve an overall n^3 algorithm. 
     *
     */
    
    long long **tmpmat = calloc(node_count, sizeof(void *));

    int node;
    int branch, n, npairs, first, second, ni, nj;
{
#if QSOPENMP_ENABLED
#pragma omp parallel for shared(tmpmat) schedule(dynamic,CHUNKSIZE) private(branch,n,npairs,first,second,i,ni,j,nj)
#endif
    for (node = leaf_count; node < node_count; ++node) {
      tmpmat[node] = calloc( leaf_count * leaf_count, sizeof(long long));
        for (branch = 0; branch < 3; ++branch) {
            n = map[node].leaf_count[branch];
            
            if (n > 1) { // we need to accumulate the data for all pairs
                npairs = n * (n-1) / 2; // number of pairs
                
                first = (3 + branch - 1) % 3;
                second = (branch + 1) % 3;
                
//                double sumdistance = 0;
                 
                for (i = 0; i < leaf_count; ++i) {
                    
                    ni = g_array_index(clt->leaf_placement, guint32, i);
                    if (map[node].node_branch[ni] != first) continue; // this leaf is in the wrong branch
                    
                    for (j = 0; j < leaf_count; ++j) {
                        nj = g_array_index(clt->leaf_placement, guint32, j);
                        if ( map[node].node_branch[nj] != second) continue; // this leaf is in the wrong branch
                        
                        tmpmat[node][i + leaf_count*j] += npairs;
                        //double dist  = gsl_matrix_get(dm, i, j);
                        //sumdistance += dist;
                    }
                }
                
                //sum += npairs * sumdistance;
            }
        }
    }
}
    
    sum = 0;
    for (node = leaf_count; node < node_count; ++node) {
    for (i = 0; i < leaf_count; ++i) {
        for (j = 0; j < leaf_count; ++j) {
            sum += tmpmat[node][i + j*leaf_count] * gsl_matrix_get(dm, i, j);
        }
    }
    }

    for (node = leaf_count; node < node_count; ++node)
      free(tmpmat[node]);
    free(tmpmat);
    free_map(map, node_count);
    return sum;
}

double qsearch_tree_score_tree_fast(QSearchTree *clt, const gsl_matrix *dm) {
    ConnectedNode *map = init_node_map(clt);
    
    // run over all pairs
    int node_count = clt->total_node_count;
    int leaf_count = (node_count + 2)/2;
    int i,j,k;

    double sum = 0.0;

    // this is apparently very slow, why?
    qsearch_tree_freshen_spm(clt);
    
    for (i = 0; i < leaf_count; ++i) {
        
        int ni = g_array_index(clt->leaf_placement, guint32, i);
        GArray *ipath = g_array_index(clt->spm, GArray *, ni);

        for (j = i+1; j < leaf_count; ++j) {
          // walk from j to i

          double dist  = gsl_matrix_get(dm, i, j);
          int nj = g_array_index(clt->leaf_placement, guint32, j);
          
          guint32 prev = nj;
          guint32 current  = g_array_index(ipath, guint32, prev); 
          
          for (;;) {
            if (current == ni) break; // we're done
            guint32 next = g_array_index(ipath, guint32, current);
             int n = 0; 
             for (k = 0;k <3; ++k) {
                int con = map[current].connections[k];
                if (con != prev && con != next) {
                    n = map[current].leaf_count[k];
                    break;
                }
             }
             
             sum += dist * ( n * (n - 1) / 2 );
            
             prev = current;
             current = next;
          }
        }
    }
    free_map(map, node_count);
    return sum;
}


#define SKIP_ORIGINAL

double qsearch_tree_score_tree(QSearchTree *clt, const gsl_matrix *dm)
{
    if (!clt->dist_calculated) {
        calc_min_max(clt, dm);
        clt->dist_calculated = 1;
    }

  if (clt->f_score_good)
    return clt->score;
  guint32 i, j, k, l;
  double acc = 0.0;
  double amin=0, amax=0;
  g_assert(dm->size1 == dm->size2);
  for (i = 0; i < dm->size1; i += 1)
    for (j = 0; j < dm->size2; j += 1)
      g_assert(gsl_matrix_get(dm, i, j) >= 0.0);
    
    //qsearch_tree_freshen_spm(clt);
    
    //struct timespec start_time;
    //struct timespec end_time;
    //clockid_t clockid = CLOCK_REALTIME; 
    //clock_gettime(clockid, &start_time);
    

#ifdef SKIP_ORIGINAL
if(0) // will skip loop
#endif

  for (i = 0; i < clt->leaf_placement->len; i += 1)
    for (j = i+1; j < clt->leaf_placement->len; j += 1)
      for (k = j+1; k < clt->leaf_placement->len; k += 1)
        for (l = k+1; l < clt->leaf_placement->len; l += 1) {
          int ni = g_array_index(clt->leaf_placement, guint32, i);
          int nj = g_array_index(clt->leaf_placement, guint32, j);
          int nk = g_array_index(clt->leaf_placement, guint32, k);
          int nl = g_array_index(clt->leaf_placement, guint32, l);
          gboolean x1, x2;
          double c1, c2, c3;
          c1  = gsl_matrix_get(dm, i, j) + gsl_matrix_get(dm, k, l);
          c2  = gsl_matrix_get(dm, i, k) + gsl_matrix_get(dm, j, l);
          c3  = gsl_matrix_get(dm, i, l) + gsl_matrix_get(dm, j, k);
          /*minscore = c1; maxscore = c1;
          if (c2 < minscore) minscore = c2;
          if (c3 < minscore) minscore = c3;
          if (c2 > maxscore) maxscore = c2;
          if (c3 > maxscore) maxscore = c3;
          amin += minscore; amax += maxscore;*/
          x1 = qsearch_tree_is_consistent_quartet(clt,ni,nj,nk,nl);
          if (x1) { acc += c1; continue; }
          x2 = qsearch_tree_is_consistent_quartet(clt,ni,nk,nj,nl);
          if (x2) { acc += c2; continue; }
#ifdef G_DISABLE_ASSERT
          acc += c3; continue;
#else
          gboolean x3 = qsearch_tree_is_consistent_quartet(clt, ni,nl,nj,nk);
          if (x3)
            acc += c3;
          else {
            g_error("Error in program logic: no consistent quartets for "
          " %d,%d,%d,%d  yielded %d,%d,%d",ni,nj,nk,nl,x1,x2,x3);
          }
#endif
        }

  //long nanos_per_second = 1000000000L;
  //clock_gettime(clockid, &end_time);
  //int nanos1 = ((end_time.tv_sec - start_time.tv_sec) * nanos_per_second + end_time.tv_nsec - start_time.tv_nsec);
  //start_time = end_time;
    
   //qsearch_optimize_tree(clt, dm);
   static int nEvals = 0;
   if (++nEvals % 10 == 0) { printf("Evals: %d\r", nEvals); fflush(stdout); }
   
   double score2 = qsearch_tree_score_tree_fast_v2(clt, dm);
  //clock_gettime(clockid, &end_time);
  int inOrder = 0;

  if (inOrder) {  
      gsl_matrix *dm2 = gsl_matrix_calloc(dm->size1, dm->size2);
      for (i=0; i < dm->size1; ++i) {
          for (j=0;j<dm->size2;++j) {
            gsl_matrix_set(dm2, g_array_index(clt->leaf_placement, guint32, i), g_array_index(clt->leaf_placement, guint32, j), gsl_matrix_get(dm, i, j));
          }
      }
      FullTree *tree = qsearch_make_fulltree(clt, dm2);
      printf("Raw scores %f %f\n", score2, tree->raw_score);
      qsearch_free_fulltree(tree);
      gsl_matrix_free(dm2);
      exit(0);
  }
        
  //int nanos2 = ((end_time.tv_sec - start_time.tv_sec) * nanos_per_second + end_time.tv_nsec - start_time.tv_nsec);

#ifdef SKIP_ORIGINAL
    acc = score2;
#endif

    amin = clt->dist_min; amax=clt->dist_max;
  g_assert(amax >= amin - ERRTOL);
  g_assert(acc >= amin - ERRTOL);
  g_assert(acc <= amax + ERRTOL);
  clt->score = (amax-acc)/(amax-amin);
  clt->f_score_good = TRUE;
  g_assert(clt->score >= 0.0);
  g_assert(clt->score <= 1.0);

    /*static double best_score = 0;
    if (0 && clt->score > best_score) {
        printf("Old score: %f\n", clt->score); 
        double oscore = qsearch_optimize_tree(clt, dm);
        best_score = clt->score;
        if (oscore > best_score) best_score = oscore;
        
    }*/

   //printf("Optimized score = %f\n", clt->score);
    //printf("nanos org = %d nanos new = %d speedup factor = %f \n", nanos1, nanos2, (double)nanos1/(double)nanos2);
    //printf("Score: %f %f\n", acc, score2);

    if (fabs(score2-acc) > 1e-6) {
        
        fprintf(stderr, "Error, score should be %f, was %f\n", acc, score2);
            
        exit(EXIT_FAILURE);
    }

  return clt->score;
}

GArray *qsearch_tree_walk_tree(const QSearchTree *clt, guint32 fromwhere, gboolean f_bfs)
{
  GArray *result = g_array_new(FALSE, TRUE, sizeof(guint32));
  GArray *todo   = g_array_new(FALSE, TRUE, sizeof(guint32));
  guint32 i, d = 0, s = clt->total_node_count, v = fromwhere;
  g_array_append_val(todo, v);
  for (i = 0; i < s; i += 1)
    g_array_index(clt->nodeflags, guint32, i) &= ~NODE_FLAG_ISWALKED;
  while (d < s) {
    g_assert(todo->len > 0);
    int remind = (f_bfs ? 0 : (todo->len-1));
    guint32 nextguy = g_array_index(todo, guint32, remind);
    g_array_remove_index(todo, remind);
    g_array_append_val(result, nextguy);
    g_array_index(clt->nodeflags, guint32, nextguy) |= NODE_FLAG_ISWALKED;
    d += 1;
    GArray *nlist = qsearch_tree_get_neighbors(clt, nextguy);
    if (nlist->len == 3 && ((g_array_index(clt->nodeflags, guint32,nextguy) & NODE_FLAG_ISFLIPPED) != 0)) {
      GArray *n2 = g_array_new(FALSE, TRUE, sizeof(guint32));
      v = g_array_index(nlist, guint32, 2);
      g_array_append_val(n2, v);
      v = g_array_index(nlist, guint32, 1);
      g_array_append_val(n2, v);
      v = g_array_index(nlist, guint32, 0);
      g_array_append_val(n2, v);
      g_array_free(nlist, TRUE);
      nlist = n2;
    }
    for (i = 0; i < nlist->len; i += 1) {
      guint32 v = g_array_index(nlist, guint32, i);
      if (g_array_index(clt->nodeflags, guint32, v) & NODE_FLAG_ISWALKED)
        continue;
      g_array_append_val(todo, v);
    }
  g_array_free(nlist, TRUE);
  }
  g_array_free(todo, TRUE);
  return result;
}

guint32 qsearch_tree_get_column_number(const QSearchTree *clt, guint32 nodenum)
{
  guint32 i;
  for (i = 0; i < clt->leaf_placement->len; i += 1)
    if (g_array_index(clt->leaf_placement, guint32, i) == nodenum)
      return i;
  g_error("Bad column number %d", nodenum);
  return 0;
}

guint32 qsearch_tree_get_node_number(const QSearchTree *clt, guint32 columnnum)
{
  return g_array_index(clt->leaf_placement, guint32, columnnum);
}

double qsearch_tree_calculate_order_cost(const QSearchTree *clt, const gsl_matrix *dm)
{
  int i;
  double acc = 0.0;
  GArray *bres = qsearch_tree_flipped_node_order(clt);
  GArray *res = g_array_new(FALSE, TRUE, sizeof(guint32));
  for (i = 0; i < bres->len; i += 1) {
    guint32 a;
    a = g_array_index(bres, guint32, i);
    if (qsearch_tree_get_neighbor_count(clt, a) != 1)
      continue;
    a = qsearch_tree_get_column_number(clt, a);
    g_array_append_val(res, a);
  }
  for (i = 0; i < res->len; i += 1) {
    double c;
    guint32 a = g_array_index(res, guint32, i);
    guint32 b = g_array_index(res, guint32, (i+1)%(res->len));
    c = gsl_matrix_get(dm, a, b) + gsl_matrix_get(dm, b, a);
    acc += c;
  }
  g_array_free(bres, TRUE);
  g_array_free(res, TRUE);
  return acc;
}

gboolean qsearch_tree_is_consistent_quartet(QSearchTree *clt, guint32 a, guint32 b, guint32 c, guint32 d)
{
  g_assert(qsearch_tree_get_neighbor_count(clt, a) == 1);
  g_assert(qsearch_tree_get_neighbor_count(clt, b) == 1);
  g_assert(qsearch_tree_get_neighbor_count(clt, c) == 1);
  g_assert(qsearch_tree_get_neighbor_count(clt, d) == 1);
  guint32 i;
  for (i = 0; i < clt->total_node_count; i += 1)
    g_array_index(clt->nodeflags, guint32, i) &= ~NODE_FLAG_QUARTETINT;
  qsearch_tree_find_path_fast(clt, a, b, clt->p1);
  for (i = 0; i < clt->p1->len; i += 1)
    g_array_index(clt->nodeflags, guint32, g_array_index(clt->p1, guint32, i)) |= NODE_FLAG_QUARTETINT;
  qsearch_tree_find_path_fast(clt, c, d, clt->p2);
  for (i = 0; i < clt->p2->len; i += 1)
    if (g_array_index(clt->nodeflags, guint32, g_array_index(clt->p2, guint32, i)) & NODE_FLAG_QUARTETINT)
      return FALSE;
  return TRUE;
}

int qsearch_count_nodes_with(const QSearchTree *clt, GArray *flags, guint32 whichflag)
{
  guint32 i, acc = 0;
  for (i = 0; i < clt->total_node_count; i += 1)
    if (g_array_index(flags, guint32, i) & whichflag)
      acc += 1;
  return acc;
}

gboolean qsearch_tree_is_standard_tree(const QSearchTree *clt)
{
  int i;
  for (i = 0; i < qsearch_tree_get_node_count(clt); i += 1) {
    int nc = qsearch_tree_get_neighbor_count(clt, i);
    if (nc != 1 && nc != 3)
      return FALSE;
  }
  return TRUE;
}

void qsearch_calculate_spm_for(QSearchTree *clt, GArray *flags, guint32 target)
{
  GArray *spm = g_array_index(clt->spm, GArray *, target);
  guint32 i, j;
  g_assert(clt->total_node_count > 1);
  g_assert(target >= 0 && target < clt->total_node_count);
  for (i = 0; i < clt->total_node_count; i += 1) {
    guint32 aval = NODE_FLAG_NEXTFRINGE;
    guint32 tval = NODE_FLAG_FRINGE | NODE_FLAG_DONE;
    guint32 q = g_array_index(flags, guint32, i) & (~aval);
    guint32 v = ((i == target) ? (q | tval) : (q & ~tval));
    g_array_index(flags, guint32, i) = v;
  }
  while (qsearch_count_nodes_with(clt, flags, NODE_FLAG_DONE) < clt->total_node_count) {
    for (i = 0; i < clt->total_node_count; i += 1) {
      guint32 v = g_array_index(flags, guint32, i) & (~NODE_FLAG_NEXTFRINGE);
      g_array_index(flags, guint32, i) = v;
    }
    for (i = 0; i < clt->total_node_count; i += 1) {
      if (g_array_index(flags, guint32, i) & (NODE_FLAG_DONE | NODE_FLAG_FRINGE | NODE_FLAG_NEXTFRINGE))
        continue;
      for (j = 0; j < clt->total_node_count; j += 1) {
        if ((g_array_index(flags, guint32, j) & NODE_FLAG_FRINGE) == 0)
          continue;
        if (qsearch_tree_is_connected(clt, i, j)) {
          guint32 v = g_array_index(flags, guint32, i) | (NODE_FLAG_NEXTFRINGE);
          g_array_index(flags, guint32, i) = v;
          g_array_index(spm, guint32, i) = j;
        }
      }
    }
    for (i = 0; i < clt->total_node_count; i += 1) {
      guint32 v = g_array_index(flags, guint32, i);
      if (v & NODE_FLAG_DONE) {
        v = v & (~NODE_FLAG_NEXTFRINGE);
        g_array_index(flags, guint32, i) = v;
      }
    }
    for (i = 0; i < clt->total_node_count; i += 1) {
      guint32 v = g_array_index(flags, guint32, i);
      if (v & NODE_FLAG_FRINGE) {
        v = (v & (~NODE_FLAG_FRINGE)) | NODE_FLAG_DONE;
        g_array_index(flags, guint32, i) = v;
      }
    }
    for (i = 0; i < clt->total_node_count; i += 1) {
      guint32 v = g_array_index(flags, guint32, i);
      if (v & NODE_FLAG_NEXTFRINGE) {
        v = (v & (~NODE_FLAG_NEXTFRINGE)) | NODE_FLAG_FRINGE;
        g_array_index(flags, guint32, i) = v;
      }
    }
  }
}

guint32 qsearch_tree_find_path_length(QSearchTree *clt,
                                        guint32 a, guint32 b)
{
  qsearch_tree_find_path_fast(clt, a, b, clt->p1);
  return clt->p1->len;
}

GArray *qsearch_tree_find_path(QSearchTree *clt, guint32 a, guint32 b) {
  GArray *result = g_array_new(FALSE, TRUE, sizeof(guint32));
  qsearch_tree_find_path_fast(clt, a, b, result);
  return result;
}

void qsearch_tree_find_path_fast(QSearchTree *clt, guint32 a, guint32 b,
GArray *dest)
{
  g_assert(a >= 0 && b >= 0 && a < clt->total_node_count && b < clt->total_node_count);
  g_array_set_size(dest, 0);
  qsearch_tree_freshen_spm(clt);
  
  int step_counter = -1;
  for (;;) {
    g_array_append_val(dest, a);
    if (step_counter > clt->total_node_count)
      break;
    step_counter += 1;
    if (a == b)
      break;
    a = g_array_index(g_array_index(clt->spm, GArray *, b), guint32, a);
  }
  if (a != b)
    g_error("Error, broken path from %d to %d for tree.\n", a, b);
}

GArray *qsearch_tree_flipped_node_order(const QSearchTree *clt)
{
  return qsearch_tree_walk_tree_dfs(clt, 0);
}

char *qsearch_tree_to_s(const QSearchTree *clt)
{
  int s = clt->total_node_count;
  char *buf = calloc(s * 1000, 1);
  char *ptr = buf;
  int i;
  GArray *res = qsearch_tree_flipped_node_order(clt);
  ptr += sprintf(ptr, "size %d\n", clt->total_node_count);
  for (i = 0; i < res->len; i += 1) {
    int j = g_array_index(res, guint32, i);
    ptr += sprintf(ptr,"%2d(%02x) ", j, g_array_index(clt->nodeflags, guint32, j));
  }
  ptr += sprintf(ptr,"\n");
  return buf;
}

void qsearch_tree_freshen_spm(QSearchTree *clt)
{
  //guint32 target;
  if (!clt->must_recalculate_paths)
    return;
  clt->must_recalculate_paths = 0;
  g_assert(clt->total_node_count > 1);
  //for (target = 0; target < clt->total_node_count; target += 1)
  //  qsearch_calculate_spm_for(clt, clt->nodeflags, target);

  ConnectedNode* map = init_node_map(clt);
    
  int i,j;
  for (i=0;i<clt->total_node_count;++i) {
      
      GArray *spm = g_array_index(clt->spm, GArray *, i);
      for (j=0;j<clt->total_node_count; ++j) {
         if (j==i) continue;
         // path from j to i
         g_array_index(spm, guint32, j) = map[j].connections[ (int) map[j].node_branch[i] ]; 
      }
  }

  free_map(map, clt->total_node_count);
}

QSearchTree *qsearch_tree_new(int how_many_leaves)
{
  QSearchTree *clt = calloc(sizeof(*clt), 1);
  g_assert(how_many_leaves >= 4);
  clt->total_node_count = how_many_leaves * 2 - 2;
  int s = clt->total_node_count;
  clt->n = g_ptr_array_sized_new(s);
  clt->dist_calculated = FALSE;
  clt->p1 = g_array_new(FALSE, TRUE, sizeof(guint32));
  clt->p2 = g_array_new(FALSE, TRUE, sizeof(guint32));
  clt->spm = g_array_sized_new(FALSE, TRUE, sizeof(gpointer), s);
  clt->nodeflags = g_array_new(FALSE, TRUE, sizeof(guint32));
  int i;
  for (i = 0; i < s; i += 1) {
    guint32 v = 0;
    g_array_append_val(clt->nodeflags, v);
  }
  clt->leaf_placement = g_array_sized_new(FALSE, TRUE, sizeof(guint32),
                                     how_many_leaves);
  clt->must_recalculate_paths = TRUE;
  for (i = 0; i < s; i += 1) {
    g_ptr_array_add(clt->n, new_qsearch_neighborlist());
    g_array_index(clt->spm, GArray *, i) = g_array_sized_new(FALSE, TRUE,
                  sizeof(guint32), s);
  }
  for (i = 0; i < how_many_leaves - 2; i += 1) {
    qsearch_tree_connect(clt, i, how_many_leaves + i);
    if (i > 0)
      qsearch_tree_connect(clt, i + how_many_leaves-1, i + how_many_leaves);
  }
  qsearch_tree_connect(clt, how_many_leaves - 2, how_many_leaves);
  qsearch_tree_connect(clt, how_many_leaves-1, s-1);
  g_assert(s > 0);
//  j = 0;
  for (i = 0; i < s; i += 1) {
    if (qsearch_tree_get_neighbor_count(clt, i) == 1) {
      guint32 v = i;
      g_array_append_val(clt->leaf_placement, v);
    }
  }
  return clt;
}

guint32 qsearch_tree_get_neighbor_count(const QSearchTree *clt, guint32 a) {
  int i, acc = 0;
  for (i = 0; i < clt->total_node_count; i += 1)
    if (qsearch_tree_is_connected(clt, i, a))
      acc += 1;
  return acc;
}

gboolean qsearch_tree_is_connected(const QSearchTree *clt, guint32 a, guint32 b) {
  g_assert(a >= 0 && b >= 0 && a < clt->total_node_count && b < clt->total_node_count);
  if (a == b)
    return FALSE;
  return a > b ?
    qsearch_neighborlist_has_neighbor(g_ptr_array_index(clt->n, b), a) :
    qsearch_neighborlist_has_neighbor(g_ptr_array_index(clt->n, a), b);
}

guint32 qsearch_tree_get_node_count(const QSearchTree *clt)
{
  return clt->total_node_count;
}

guint32 qsearch_tree_get_leaf_node_count(const QSearchTree *clt)
{
  return (clt->total_node_count+2)/2;
}

guint32 qsearch_tree_get_kernel_node_count(const QSearchTree *clt)
{
  return (clt->total_node_count-2)/2;
}

void qsearch_tree_free(QSearchTree *clt)
{
  int i;
  if (clt == NULL)
    return;
  for (i = 0; i < clt->total_node_count; i += 1) {
    free_qsearch_neighborlist(g_ptr_array_index(clt->n, i));
    g_array_free(g_array_index(clt->spm, GArray *, i), TRUE);
  }
  g_array_free(clt->spm, TRUE);
  g_array_free(clt->nodeflags, TRUE);
  g_array_free(clt->leaf_placement, TRUE);
  g_array_free(clt->p1, TRUE);
  g_array_free(clt->p2, TRUE);
  g_ptr_array_free(clt->n, TRUE);
  g_free(clt);
}

void qsearch_tree_connect(QSearchTree *clt, guint32 a, guint32 b)
{
  g_assert(qsearch_tree_is_connected(clt, a,b) == FALSE);
  g_assert(a != b);
  if (a < b)
    qsearch_neighborlist_add_neighbor(g_ptr_array_index(clt->n, a), b);
  else
    qsearch_neighborlist_add_neighbor(g_ptr_array_index(clt->n, b), a);
  clt->must_recalculate_paths = TRUE;
  clt->f_score_good = FALSE;
}

void qsearch_tree_disconnect(QSearchTree *clt, guint32 a, guint32 b)
{
  g_assert(qsearch_tree_is_connected(clt, a,b) == TRUE);
  g_assert(a != b);
  if (a < b)
    qsearch_neighborlist_remove_neighbor(g_ptr_array_index(clt->n, a), b);
  else
    qsearch_neighborlist_remove_neighbor(g_ptr_array_index(clt->n, b), a);
  clt->must_recalculate_paths = TRUE;
  clt->f_score_good = FALSE;
}

static void random_pair(FullTree *tree, int *A, int *B) {
    
    int p1 = (int)(rand() % tree->node_count);
    int p2 = p1;
    
    while (p2 == p1 || qsearch_fulltree_move_to(tree, p1, p2) == p2 ) p2 = (int)(rand() % tree->node_count);
    
    *A = p1;
    *B = p2;
}

/*
static void random_internal_pair(FullTree *tree, int *A, int *B) {
    int leaf_count = (tree->node_count + 2)/2;
    int p1 = (int)(rand() % (tree->node_count - leaf_count)) + leaf_count;
    int p2 = p1;
    
    while (p2 == p1 || qsearch_fulltree_move_to(tree, p1, p2) == p2 ) p2 = (int)(rand() % (tree->node_count - leaf_count)) + leaf_count;
    
    *A = p1;
    *B = p2;
}
*/

QSearchTree *qsearch_tree_find_better_tree(QSearchTree *clt,
guint32 howManyTries, const gsl_matrix *dm)
{
    if (!clt->dist_calculated) {
        calc_min_max(clt, dm);
        clt->dist_calculated = 1;
    }
  
  QSearchTree *cand, *result = NULL;
  int i, totmuts;
  double candscore;
  double best_score;
  double curscore = qsearch_tree_score_tree(clt, dm);
#if QSOPENMP_ENABLED
  if (!g_thread_supported ()) g_thread_init (NULL);
#pragma omp parallel shared(result, curscore) private(i, totmuts, cand, candscore, best_score) firstprivate(howManyTries, dm, clt)
 {
#pragma omp for schedule (dynamic, 1)
#else
  howManyTries = 1;
  do {
#endif
  for (i = 0; i < howManyTries; i += 1) {
    cand = qsearch_tree_clone(clt);
    cand->dist_min = clt->dist_min;
    cand->dist_max = clt->dist_max;
     
    //qsearch_tree_complex_mutation(cand);
    FullTree *tree = qsearch_make_fulltree(cand, dm);

    // perform node_count swaps, keep track of best
    best_score = tree->raw_score;
    totmuts = tree->node_count;//qsearch_tree_get_mutation_distribution_sample(clt);
    
    int j;
     
    for (j = 0; j < totmuts; ++j) {
        double beta = 1.0; // set to 0.0 to mimick random behaviour. This behaviour is a metropolis markov chain
         
        int p1, p2;
        random_pair(tree, &p1, &p2);

        double cur = tree->raw_score;

        if ( rand() % 3 < 2) { 
            
            qsearch_fulltree_swap_nodes(tree, p1, p2, dm);
            
            if (tree->raw_score <= best_score || fabs(tree->raw_score - best_score) < 1e-6) { 
                qsearch_fulltree_to_searchtree(tree, cand);
                best_score = tree->raw_score;
                //printf("Score improved from %f to %f, raw: %f \n", curscore, cand->score, tree->raw_score);
            }

            // calculate acceptance
            double now = tree->raw_score;
            if (drand48() >= exp(beta * (cur-now) )) { // reject
                qsearch_fulltree_swap_nodes(tree, p1, p2, dm);
            } 

        } else { // transfer tree
            
            int interior = qsearch_fulltree_move_to(tree, p1, p2);
            g_assert(interior != p2);
             
            int sibling = qsearch_fulltree_find_sibling(tree, p1, p2);
           
            // move entire subtree containing p1 and sibling in the place of p2 
            qsearch_fulltree_swap_nodes(tree, interior, p2, dm);
            
            if (tree->raw_score <= best_score || fabs(tree->raw_score - best_score) < 1e-6) { 
                qsearch_fulltree_to_searchtree(tree, cand);
                best_score = tree->raw_score;
                //printf("Score improved from %f to %f, raw: %f \n", curscore, cand->score, tree->raw_score);
            }
            
            // swap the sibling back in its original place (making 'node' a sibling of 'p2'
            qsearch_fulltree_swap_nodes(tree, sibling, p2, dm);
           
            // postcondition: 
            g_assert( qsearch_fulltree_find_sibling(tree, p1, sibling) == p2);

            if (tree->raw_score <= best_score || fabs(tree->raw_score - best_score) < 1e-6) { 
                qsearch_fulltree_to_searchtree(tree, cand);
                best_score = tree->raw_score;
                //printf("Score improved from %f to %f, raw: %f \n", curscore, cand->score, tree->raw_score);
            }
            
            // calculate acceptance
            double now = tree->raw_score;
            if (drand48() >= exp(beta*(cur-now) )) { // reject
                qsearch_fulltree_swap_nodes(tree, sibling, p2, dm);
                qsearch_fulltree_swap_nodes(tree, interior, p2, dm);
            } 
        }
    }
    
    //qsearch_fulltree_to_searchtree(tree, cand);
    qsearch_free_fulltree(tree);
    
    //cand->f_score_good = 0;
    //double org = cand->score;
    candscore  = qsearch_tree_score_tree(cand, dm);
    
    //if (fabs(org-cand->score) > 1e-6) {
    //    printf("Error, score should be %f, was %f\n", cand->score, org);
    //    exit(0);
    //}

    if (candscore <= curscore) {
      qsearch_tree_free(cand);
      continue;
    }
#if QSOPENMP_ENABLED
#pragma omp critical
#endif
    if (candscore > curscore) {
      if (result)
        qsearch_tree_free(result);
      result = cand;
      curscore = candscore;
    }
    else
      qsearch_tree_free(cand);
    }
#if QSOPENMP_ENABLED
}
#else
  } while (0);
#endif
  return result;
}
  
  /*
  QSearchTree *cand, *result = NULL;
  int i;
  double candscore;
  double curscore = qsearch_tree_score_tree(clt, dm);
#if QSOPENMP_ENABLED
  if (!g_thread_supported ()) g_thread_init (NULL);
#pragma omp parallel shared(result, curscore) private(i, cand, candscore) firstprivate(howManyTries, dm, clt)
 {
#pragma omp for schedule (dynamic, 1)
#else
  howManyTries = 1;
#endif
  for (i = 0; i < howManyTries; i += 1) {
    cand = qsearch_tree_clone(clt);
    qsearch_tree_complex_mutation(cand);
    candscore  = qsearch_tree_score_tree(cand, dm);
    if (candscore <= curscore) {
      qsearch_tree_free(cand);
      continue;
    }
#if QSOPENMP_ENABLED
#pragma omp critical
#endif
    if (candscore > curscore) {
      if (result)
        qsearch_tree_free(result);
      result = cand;
      curscore = candscore;
    }
    else
      qsearch_tree_free(cand);
  }
#if QSOPENMP_ENABLED
}
#endif
  return result;
}
*/

static char *nodeToNewick(QLabeledTree *tree, int me, int fromWhere) {
  GArray *a;
  int i, todo;
  char *result[3];
  char *bigr;
  int done = 0;
  if (tree->labels[me] == NULL)
    g_error("Tree must be labeled to convert to Newick format.");
  if (qsearch_tree_get_neighbor_count(tree->qs, me) == 1)
    return g_strdup(tree->labels[me]);
  todo = sizeof(a)/sizeof(a[0]);
  a = qsearch_tree_get_neighbors(tree->qs, me);
  for (i = 0; i < a->len; i += 1) {
    int curnode = g_array_index(a, guint32, i);
    if (curnode != fromWhere) {
      result[done++] = nodeToNewick(tree, curnode, me);
    }
  }
  assert(done == 2 || done == 3);
  if (done == 2)
    bigr = g_strdup_printf("(%s,%s)", result[0], result[1]);
  else
    bigr = g_strdup_printf("(%s,%s,%s)",result[0],result[1],result[2]);
  for (i = 0; i < done; i += 1)
    g_free(result[i]);
  g_array_free(a, TRUE);
  return bigr;
}

char *qsearch_tree_to_nexus(QLabeledTree *tree) {
   int root;
  char *nstr, *rstr;
  for (root = 0; qsearch_tree_get_neighbor_count(tree->qs, root) == 1;)
    root += 1;
  nstr = nodeToNewick(tree, root, -1);
  rstr = g_strdup_printf("BEGIN Trees;\n[1] tree 'CLQCS'= %s;\nEND;",nstr);
  g_free(nstr);
  return rstr;
}


char *qsearch_tree_to_nexus_full(QLabeledTree *tree, LabeledMatrix *lm)
{
  char *cres;
  GString *res;
  char *ctb = qsearch_tree_to_nexus(tree);
  GString *tb = g_string_new(ctb);
  res = complearn_matrix_prettyprint_nex(lm, tb);
  g_string_free(tb, TRUE);
  g_free(ctb);
  cres = res->str;
  g_string_free(res, FALSE);
  return cres;
}

/* Returns true if every node has exactly 1 or 3 neighbors */
gboolean qsearch_tree_is_tree_ternary(QSearchTree *clt)
{
  int i;
  for (i = 0; i < clt->total_node_count; i += 1) {
    int nc = qsearch_tree_get_neighbor_count(clt, i);
    if (nc != 1 && nc != 3)
      return FALSE;
  }
  return TRUE;
}

/* Sets the "connectedness" state (TRUE or FALSE) between nodes a and b and
 * returns the old connectedness status that was overwritten.
 */
gboolean qsearch_tree_set_connected(QSearchTree *clt, guint32 a, guint32 b, gboolean newconstate)
{
  g_assert(a >= 0 && b >= 0 && a < clt->total_node_count && b < clt->total_node_count);
  if (a == b)
    return FALSE;
  gboolean oldconstate = qsearch_tree_is_connected(clt, b, a);
  if (oldconstate != newconstate) {
    if (newconstate)
      qsearch_tree_connect(clt, a, b);
    else
      qsearch_tree_disconnect(clt, a, b);
  }
  return oldconstate;
}

void qsearch_tree_clear_all_connections(QSearchTree *clt)
{
  guint32 i, j;
  for (i = 0; i < clt->total_node_count; i += 1)
    for (j = i+1; j < clt->total_node_count; j += 1)
      qsearch_tree_set_connected(clt, j, i, FALSE);
}

static void qsearch_tree_flatten_leafperm(QSearchTree *clt)
{
  int i;
  for (i = 0; i < clt->leaf_placement->len; i += 1)
    g_array_index(clt->leaf_placement, guint32, i) = i;
}

/*
static gint32 qsearch_tree_find_index(LabeledMatrix *lm, const gchar *lab)
{
  int i;
  for (i = 0; lm->labels1[i]; i += 1) {
    gchar *cur = lm->labels1[i];
    if (strcmp(cur, lab) == 0)
      return i;
  }
  return -1;
}
*/

/*
label="S(T)=0.994961";
*/
gdouble qsearch_tree_read_from_dot(QSearchTree *clt, GString *treedot, LabeledMatrix *lm)
{
  gchar **lines, *cur, *str;
  const gchar *scoreprefix = "label=\"S(T)=";
  lines = g_strsplit(treedot->str, "\n", 0);
  int i;
  double score;
  gboolean gotscore = FALSE;
  qsearch_tree_clear_all_connections(clt);
  for (i = 0; (cur = lines[i]) != NULL; i += 1) {
    guint32 a, b;
    if (g_str_has_prefix(cur, scoreprefix)) {
      str = cur + strlen(scoreprefix);
      strtok(str, "\"");
      score = atof(str);
      gotscore = TRUE;
      continue;
    }
    if (strlen(cur) < 3)
      continue;
    if (strstr(cur, "label") || strstr(cur, "dotted") || strstr(cur, "graph"))
      continue;
    sscanf(cur, "%d -- %d", &a, &b);
    qsearch_tree_connect(clt, b, a);
  }
  for (i = 0; (cur = lines[i]) != NULL; i += 1)
    free(cur);
  free(lines);
  qsearch_tree_flatten_leafperm(clt);
  if (qsearch_tree_is_tree_ternary(clt))
    return gotscore ? score : -2.0;
  else
    return -1.0;
}

gsl_matrix *qsearch_tree_get_adjacency_matrix(const QSearchTree *clt)
{
  int s = qsearch_tree_get_node_count(clt);
  gsl_matrix *m = gsl_matrix_calloc(s,s);
  int i, j;
  for (i = 0; i < s; i += 1) {
    for (j = i+1; j < s; j += 1) {
      int b = qsearch_tree_is_connected(clt, i, j) ? 1 : 0;
      gsl_matrix_set(m, i, j, b);
      gsl_matrix_set(m, j, i, b);
    }
  }
  return m;
}

