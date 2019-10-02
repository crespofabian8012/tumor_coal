//
//  tree.h
//  TumorCellCoal
//
//  Created by Fausto Fabian Crespo Fernandez on 7/8/19.
//  Copyright © 2019 Fausto Fabian Crespo Fernandez. All rights reserved.
//

#ifndef tree_h
#define tree_h

static void BuildTree(Population **populations,Population *CurrentPop,
               long int *seed,
               ProgramOptions *programOptions,
               //               int numNodes,
               //               int TotalNumSequences,
               //               int    numClones,
               //               long int userSeed,
               //               double  mutationRate,
               //               int noisy,
               //               int outgroupSelection,
               //               int thereisOutgroup,
               //               double outgroupBranchLength_Root1Root2,
               //               double outgroupBranchLength_RootSample,
               TreeNode    **nodes,
               TreeNode   **treeTips,
               TreeNode    **treeRootInit,
               int *nextAvailable,
               int *newInd,
               double *currentTime,
               int *labelNodes
                      );

static void connectNodelets(TreeNode *node );
static void connectNodes(TreeNode *left, TreeNode *right, TreeNode *ancester  );
static void RelabelNodes (TreeNode *p, TreeNode **treeRootInit, int *intLabel);
inline static int   Index (TreeNode *p);
inline static int    Label (TreeNode *p);
static void PrintTrees(int replicate, TreeNode **treeRootInit, FILE  *fpTrees, double mutationRate, int doUseObservedCellNames);
static void     PrintTrees2(int replicate, TreeNode **treeRootInit, FILE *fpTrees2, double mutationRate, char * cellNames[], int doUseObservedCellNames);


static void WriteTree2 (TreeNode *p, double mutationRate, FILE    *fpTrees2, char *ObservedCellNames[], int *indexCurrentCell, int doUseObservedCellNames);
static void WriteTree (TreeNode *p, double mutationRate, FILE    *fpTrees, int doUseObservedCellNames);
double SumBranches (TreeNode *p,  double mutationRate);

void toNewickString ( TreeNode *p, double mutationRate, char NewickString[],    int doUseObservedCellNames);
char * toNewickString2 ( TreeNode *p, double mutationRate,    int doUseObservedCellNames);
char * toNewickString4 ( TreeNode *p, double mutationRate,     int doUseObservedCellNames);

void AssignTreeNodestoObservedCellNamesRecursive(TreeNode *p, char* ObservedCellNames[], int *currentCellName);
void AssignTreeNodestoObservedCellNames(TreeNode **treeRootInit, char* ObservedCellNames[]);


pll_utree_t * load_tree(const char *fname);

static void MakeCoalescenceTree2 (long int *seed, Population **populations,
                                  int *numNodes,
                                  int numClones,
                                  ProgramOptions *programOptions,
                                  double      cumNumCA,
                                  double meanNumCA,
                                  double cumNumMIG,
                                  double meanNumMIG,
                                  int  *numMIG,
                                  int  *numCA,
                                  double *numEventsTot,
                                  TreeNode** nodes,
                                  TreeNode** treeTips,
                                  TreeNode    **treeRootInit,
                                  char* ObservedCellNames[],
                                  int *sampleSizes
                                  ) ;
static void MakeCoalescenceEvent(Population **populations, Population *popI, TreeNode **nodes, int numClones, long int* seed, int noisy, int *numActiveGametes,
                                 int* nextAvailable,
                                 int*labelNodes, double *currentTime, int *numNodes);

void  destroyTree(pll_utree_t * tree, void (*cb_destroy)(void *));

bool isLeaf(TreeNode *node);
TreeNode* getInternalNodebyIndexSubTree(TreeNode *currentNode,int indexNode, Population *pop);
TreeNode *getHealthyTip(TreeNode *treeRootInit);
void dealloc_data(pll_unode_t * node, void (*cb_destroy)(void *));

int CountTrueVariants (TreeNode *nodes, int numSites, int numCells, TreeNode *HEALTHY_ROOT, SiteStr* allSites, int* variantSites,  int *SNVsites);
double ComputeAvgHeterocigocity (int numSites, int numCells, TreeNode **treeTips, TreeNode *treeRoot, SiteStr* allSites, int* SNVsites);

void setLength(TreeNode *node );

PLL_EXPORT int pllmod_utree_spr1(pll_unode_t * p_edge,
                                 pll_unode_t * r_edge,
                                 pll_tree_rollback_t * rollback_info,
                                 pll_utree_rb_t * rb,
                                 TreeNode **nodes,
                                 double time_from_r_to_p,
                                 double old_father_pop_effect_pop_size,
                                 double new_father_pop_effect_pop_size,
                                 TreeNode *container_of_p,
                                 TreeNode *container_of_r,
                                 TreeNode *container_of_p_prime,//mrca
                                 TreeNode *container_of_r_prime,//the other end ppoint of the edge to reconnect
                                 TreeNode * container_of_u,
                                 TreeNode * container_of_v
                                 );

PLL_EXPORT int pll_utree_spr1(pll_unode_t * p,
                              pll_unode_t * r,
                              pll_utree_rb_t * rb,
                              double * branch_lengths,
                              unsigned int * matrix_indices,
                              TreeNode **nodes,
                              double time_from_r_to_p,
                              double old_father_pop_effect_pop_size,
                              double new_father_pop_effect_pop_size,
                              TreeNode *container_of_p,
                              TreeNode *container_of_r,
                              TreeNode *container_of_p_prime,//mrca
                              TreeNode *container_of_r_prime,//the other end ppoint of the edge to reconnect
                              TreeNode * container_of_u,
                              TreeNode * container_of_v
                              );

static void utree_link(pll_unode_t * a,pll_unode_t * b,double length,unsigned pmatrix_index);
double * pllmod_msa_empirical_frequencies1(pll_partition_t * partition);


static int set_tipclv1(pll_partition_t * partition,
                       unsigned int tip_index,
                       const pll_state_t * map,
                       const char * sequence,
                       double _seq_error_rate,
                       double _dropout_rate
                       );

int pll_set_tip_states1(pll_partition_t * partition,
                        unsigned int tip_index,
                        const pll_state_t * map,
                        const char * sequence);
#endif /* tree_h */
