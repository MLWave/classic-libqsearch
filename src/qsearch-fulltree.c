
#include "qsearch.h"

// All data is statically allocated, so there's no need to resize things

struct _FullNode {
    int connections[3];
    
    // cached counts and distance sums
    int         leaf_count[3];
    double      dist[3];

    guint8 *node_branch; // pointing in the direction where to find a node
};

typedef struct _FullNode FullNode;

struct _Misc {
    FullNode *nodes;
    
    // some temporary storage
    GArray *tmpA;
    GArray *tmpB;
};

typedef struct _Misc Misc;

static FullNode *get_nodes(FullTree *tree) { return (FullNode*) ((Misc *)tree->data)->nodes; }
static GArray   *get_tmpA(FullTree *tree)  { return ((Misc *)tree->data)->tmpA; }
static GArray   *get_tmpB(FullTree *tree)  { return ((Misc *)tree->data)->tmpB; }

int qsearch_next_node(FullTree *tree, int from, int to) {
    FullNode *map = get_nodes(tree);
    return map[from].connections[ (int) map[from].node_branch[to] ];
}

static inline double npairs(double n) { return n * (n-1)/2; }

/*
static void print_tree(FullTree *tree) {
    
    FullNode *map = get_nodes(tree);
    int node_count = tree->node_count;

    int i,j;
    int leaf_count = (node_count+2)/2;

    printf("Nodes %d, leafs %d\n", node_count, leaf_count);
    for (i = 0; i < node_count; ++i) {
        printf("%d: connections: %d %d %d, counts: %d %d %d, distances %f %f %f ", i, 
                map[i].connections[0],
                map[i].connections[1],
                map[i].connections[2],
                map[i].leaf_count[0],
                map[i].leaf_count[1],
                map[i].leaf_count[2],
                map[i].dist[0],
                map[i].dist[1],
                map[i].dist[2]
                );
            printf("Leafs :");
            for (j = 0; j < node_count; ++j) {
                printf(" %c", i==j?'x':'0' + map[i].node_branch[j]);
            }
        printf("\n");
    }
}
*/

inline static int find_branch(int connections[3], int to) {
    if (connections[0] == to) return 0;
    if (connections[1] == to) return 1;
    if (connections[2] == to) return 2;
    g_assert(0); // 
    return -1;
}

static void set_score(FullTree *tree) {
    // calculate the score
    tree->raw_score = 0;
    FullNode *map = get_nodes(tree);
    int i;
    for (i = (tree->node_count + 2)/2; i < tree->node_count; ++i) {
        tree->raw_score += npairs(map[i].leaf_count[0]) * map[i].dist[0];
        tree->raw_score += npairs(map[i].leaf_count[1]) * map[i].dist[1];
        tree->raw_score += npairs(map[i].leaf_count[2]) * map[i].dist[2];
    }
}


FullTree* qsearch_make_fulltree(QSearchTree *clt, const gsl_matrix *dm) {
    
    int i,j; 
    int node_count = clt->total_node_count;
    int leaf_count = (node_count + 2)/2;
    
    FullTree *tree = malloc(sizeof(FullTree));
    tree->node_count = node_count;
    
    tree->data = malloc(sizeof(Misc));
    
    ((Misc *)tree->data)->nodes = malloc(sizeof(FullNode) * node_count);
    ((Misc *)tree->data)->tmpA  = g_array_sized_new(FALSE, FALSE, sizeof(guint32), node_count);
    ((Misc *)tree->data)->tmpB  = g_array_sized_new(FALSE, FALSE, sizeof(guint32), node_count);

    FullNode *map = get_nodes(tree);
    
    GArray *todo = g_array_sized_new(FALSE, FALSE, sizeof(guint32), node_count - leaf_count);
    g_array_set_size(todo, node_count-leaf_count);

    for (i = 0; i < node_count; ++i) {
        FullNode *node = (map + i);
        
        for (j = 0; j < 3; ++j) {
            map[i].connections[j] = -1;
            map[i].leaf_count[j] = 0;
            map[i].dist[j] = 0;
        }
        
        guint8 *node_branch = malloc(sizeof(guint8) * node_count);
        node->node_branch = node_branch;

        if (i < leaf_count) {
            for (j=0; j < node_count; ++j) node_branch[j] = 0;
            node->leaf_count[0] = leaf_count - 1;
        } else {
            for (j = 0; j < node_count; ++j) node_branch[j] = -1;
            g_array_index(todo, guint32, i - leaf_count) = i;
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
    while (todo->len > 0) {
        for (i = 0; i < todo->len; ++i) {
            
            int this_node = g_array_index(todo, guint32, i); 
                        
            for (j = 0; j < 3; ++j) {
                int connected_node = map[this_node].connections[j];
                int branch = find_branch(map[connected_node].connections, this_node);

                if (map[connected_node].leaf_count[branch] == 0) {
                    int first = (3 + j-1) % 3;
                    int second = (j + 1)  % 3;
                    if (map[this_node].leaf_count[first] != 0 && map[this_node].leaf_count[second] != 0) {
                        map[connected_node].leaf_count[branch] = map[this_node].leaf_count[first] + map[this_node].leaf_count[second];
                        
                        // set the node_branch information
                        int k;
                        for (k = 0; k < node_count; ++k) {
                            // node present in one of the two branches pointing away from this
                            int node_present = k==this_node || map[this_node].node_branch[k] == first || map[this_node].node_branch[k] == second;
                            
                            if (node_present) map[connected_node].node_branch[k] = branch;
                        }
                    }
                }
            }
        }
        
        for (i = 0; i < todo->len; ++i) {
            
            int this_node = g_array_index(todo, guint32, i);
             
            // we are done when all entries are set for this node, AND there are no pending assignments to neighbours

            // are all entries set for this node?
            for (j = 0; j < 3; ++j) {
                if (map[this_node].leaf_count[j] == 0) { 
                    break;
                }
            }

            // are there pending assignments? (we could do the assignments here, but that would duplicate code)
            int done = 1;
            for (j = 0; j < 3; ++j) {
                if (this_node < leaf_count) continue;
                int connected_node = map[this_node].connections[j];
                // find connection
                int branch = find_branch(map[connected_node].connections, this_node);
                if (map[connected_node].leaf_count[branch] == 0) {
                    done = 0;
                    break;
                }
            }
            
            if (done) {
                //printf("Removing node %d, counts %d %d %d\n", i, map[this_node].leaf_count[0], map[this_node].leaf_count[1], map[this_node].leaf_count[2]);
                g_array_remove_index_fast(todo, i);
                --i;
            }
        }
    }
    
    g_array_free(todo, TRUE); 
    
    int k; 
    // now fill in the distances
    for (i=leaf_count; i < node_count; ++i) {
        for (j = 0; j < leaf_count; ++j) {
            for (k = j+1; k < leaf_count; ++k) {

                int b1 = map[i].node_branch[j];
                int b2 = map[i].node_branch[k];
                
                if (b1 == b2) continue;

                int b3 = 3 - b1 - b2;

                map[i].dist[b3] += gsl_matrix_get(dm, j, k);
            }
        }
    }
    
    set_score(tree);
    
    return tree;
}

void qsearch_free_fulltree(FullTree *tree) {
    
    g_array_free(get_tmpA(tree), TRUE);
    g_array_free(get_tmpB(tree), TRUE);

    FullNode *map = get_nodes(tree);
     
    int i;
    for (i = 0; i < tree->node_count; ++i) free(map[i].node_branch);
    free(map);
    free(tree->data);
    free(tree);
}
    
gboolean qsearch_fulltree_can_swap(FullTree *tree, int A, int B) {
   FullNode *ctm = get_nodes(tree);

   if (A == B) return FALSE; // no point in doing anything
    
   int interiorA = ctm[A].connections[ ctm[A].node_branch[B] ];
   int interiorB = ctm[B].connections[ (int) ctm[B].node_branch[A] ];
   
   if (interiorA == interiorB || interiorA == B) return FALSE; // swap does not change score
    
   return TRUE;
}

void qsearch_fulltree_swap_nodes(FullTree *tree, int A, int B, const gsl_matrix *dm) {
   FullNode *ctm = get_nodes(tree);

   if (A == B) return; // no point in doing anything
    
   int interiorA = ctm[A].connections[ ctm[A].node_branch[B] ];
   int interiorB = ctm[B].connections[ (int) ctm[B].node_branch[A] ];
   
   if (interiorA == interiorB || interiorA == B) return; // swap does not change score
   
   int node_count = tree->node_count;
   int leaf_count = (node_count + 2)/2;

   int aToInteriorBranch = ctm[A].node_branch[interiorA];
   int bToInteriorBranch = ctm[B].node_branch[interiorB];
   
   int countA = leaf_count - ctm[A].leaf_count[aToInteriorBranch];
   int countB = leaf_count - ctm[B].leaf_count[bToInteriorBranch];
   
   int interiorToABranch = ctm[interiorA].node_branch[A];
   int interiorToBBranch = ctm[interiorB].node_branch[B];
     
   // loop over internal nodes, swap node_branches from A <-> B
   int node = interiorA;
   
   //printf("Score before %f\n", tree->raw_score);
    
   // store the nodes that need to be updated
   GArray *aNodes = get_tmpA(tree); //g_array_sized_new(FALSE, FALSE, sizeof(guint32), countA*2 - 1);
   GArray *bNodes = get_tmpB(tree); //g_array_sized_new(FALSE, FALSE, sizeof(guint32), countB*2 - 1);
    
   g_array_set_size(aNodes, 0);
   g_array_set_size(bNodes, 0);

   int i,j;
   for (i = 0; i < node_count; ++i) {
        if (i == A || ctm[A].node_branch[i] != aToInteriorBranch) {
            g_array_append_val(aNodes, i);
        }
        if (i == B || ctm[B].node_branch[i] != bToInteriorBranch) {
            g_array_append_val(bNodes, i);
        }
   }

   // move towards B 
   while (node != B) {

        int aBranch = ctm[node].node_branch[A];
        int bBranch = ctm[node].node_branch[B];
        int cBranch = 3 - aBranch - bBranch;
        
        tree->raw_score -= npairs(ctm[node].leaf_count[0]) * ctm[node].dist[0];
        tree->raw_score -= npairs(ctm[node].leaf_count[1]) * ctm[node].dist[1];
        tree->raw_score -= npairs(ctm[node].leaf_count[2]) * ctm[node].dist[2];
        
        //printf("Node %d\n", node);
        //printf("Score now %f, distances %f %f %f\n", tree->raw_score, ctm[node].dist[0], ctm[node].dist[1], ctm[node].dist[2]);
        
        ctm[node].leaf_count[aBranch] += countB - countA;
        ctm[node].leaf_count[bBranch] += countA - countB;
       
        // update the branches that point to elements from A 
        for (i = 0; i < aNodes->len; ++i) {
            guint32 aNode = g_array_index(aNodes, guint32, i);
            g_assert(ctm[node].node_branch[aNode] == aBranch);
            ctm[node].node_branch[aNode] = bBranch;
           
            // update distances 
            if (aNode < leaf_count) { // it's  a leaf
                for (j = 0; j < leaf_count; ++j) {
                    if (aNode==j) continue;
                    
                    double d = gsl_matrix_get(dm, aNode, j);

                    if (ctm[node].node_branch[j] == cBranch) {
                        ctm[node].dist[aBranch] += d;
                        ctm[node].dist[bBranch] -= d;
                    } else if (ctm[node].node_branch[j] == aBranch) {
                        ctm[node].dist[cBranch] += d;
                    } else {
                        ctm[node].dist[cBranch] -= d;
                    }
                }
            }
        } 
        
            
        // update the branches that point to elements from B
        for (i = 0; i < bNodes->len; ++i) {
            guint32 bNode = g_array_index(bNodes, guint32, i);

            g_assert(ctm[node].node_branch[bNode] == bBranch);
            ctm[node].node_branch[bNode] = aBranch;
                
            if (bNode < leaf_count) { // it's  a leaf
                for (j = 0; j < leaf_count; ++j) {
                    if (bNode==j) continue;
                    
                    double d = gsl_matrix_get(dm, bNode, j);

                    if (ctm[node].node_branch[j] == cBranch) {
                        ctm[node].dist[bBranch] += d;
                        ctm[node].dist[aBranch] -= d;
                    } else if (ctm[node].node_branch[j] == bBranch) {
                        ctm[node].dist[cBranch] += d;
                    } else {
                        ctm[node].dist[cBranch] -= d;
                    }
                }
            }
        }

        tree->raw_score += npairs(ctm[node].leaf_count[0]) * ctm[node].dist[0];
        tree->raw_score += npairs(ctm[node].leaf_count[1]) * ctm[node].dist[1];
        tree->raw_score += npairs(ctm[node].leaf_count[2]) * ctm[node].dist[2];
        
        //printf("Score now2 %f, distances %f %f %f\n", tree->raw_score, ctm[node].dist[0], ctm[node].dist[1], ctm[node].dist[2]);
        
        node = ctm[node].connections[ bBranch ];
   }
     
   // swap directions for A and B
   ctm[interiorA].connections[ interiorToABranch ] = B;
   ctm[interiorB].connections[ interiorToBBranch ] = A;
    
   // swap connections for A and b
   ctm[A].connections[ aToInteriorBranch ] = interiorB;
   ctm[B].connections[ bToInteriorBranch ] = interiorA;
}

void qsearch_fulltree_to_searchtree(FullTree *src, QSearchTree *clt) {
    
    FullNode *map = get_nodes(src);
    int node_count = src->node_count;
    int leaf_count = (node_count + 2)/2;
    int i,j;
    
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
    
    clt->score = (clt->dist_max - src->raw_score) / (clt->dist_max - clt->dist_min); 

    clt->must_recalculate_paths = TRUE;
    clt->f_score_good = TRUE;
}

guint32 qsearch_fulltree_move_to(FullTree *tree, guint32 from, guint32 to) {
    FullNode *map = get_nodes(tree);
    return map[from].connections[ map[from].node_branch[to] ];
}

guint32 qsearch_fulltree_find_sibling(FullTree *tree, guint32 node, guint32 ancestor) {
    FullNode *map = get_nodes(tree);
    
    g_assert(node != ancestor);

    guint32 parent = qsearch_fulltree_move_to(tree, node, ancestor);
    
    g_assert(parent != ancestor);

    int branch2node = map[parent].node_branch[node];
    int branch2ancestor = map[parent].node_branch[ancestor];

    int branch2sibling = 3 - branch2node - branch2ancestor;

    return map[parent].connections[ branch2sibling ];

}

double  qsearch_fulltree_sum_distance(FullTree *tree, int a, int b) {
    
    FullNode *map = get_nodes(tree);

    double sum = 0.0;

    int node = map[a].connections[ map[a].node_branch[b] ];
    int n = 0;
    while (node != b) {
        int toa = map[node].node_branch[a];
        int tob = map[node].node_branch[b];
        
        sum += npairs( map[node].leaf_count[ toa ] ) * map[node].dist[toa];
        sum += npairs( map[node].leaf_count[ tob ] ) * map[node].dist[tob];
        n++;
        node = map[node].connections[ tob ];
    }

    return sum / n;
}

double  qsearch_fulltree_sum_distance_org(FullTree *tree, int a, int b) {
    
    FullNode *map = get_nodes(tree);
    int branch2b = map[a].node_branch[b];
    int branch2a = map[b].node_branch[a];
     
    return npairs( map[a].leaf_count[branch2b] ) * map[a].dist[branch2b] + npairs( map[b].leaf_count[branch2a] ) * map[b].dist[branch2a];

}

void qsearch_fulltree_get_children(FullTree *tree, int node, int ancestor, int *child1, int *child2) {
    
    FullNode *map = get_nodes(tree);
    int branch = map[node].node_branch[ancestor];
    
    *child1 = map[node].connections[ (3 + branch-1) % 3];
    *child2 = map[node].connections[ (branch + 1) % 3];

}


