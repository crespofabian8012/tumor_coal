//
//  output_functions.hpp
//  simul
//
//  Created by Seong-Hwan Jun on 2019-10-08.
//

#ifndef output_functions_hpp
#define output_functions_hpp

#include <stdio.h>

#include "data_types.hpp"
#include "tree_node.hpp"

TreeNode *getHealthyTip(TreeNode *treeRootInit);

void WriteTree2 ( TreeNode *p, double mutationRate, FILE    *fpTrees2, char *cellNames[], int *indexCurrentCell, int doUseObservedCellNames);
void WriteTree (TreeNode *p, double mutationRate, FILE    *fpTrees, int doUseObservedCellNames);
void PrintTrees2(int replicate, TreeNode *treeRootInit,   FILE   *fpTrees2 , double mutationRate,char * ObservedCellNames[],int doUseObservedCellNames);
void PrintTrees(int replicate, TreeNode *treeRootInit,   FILE  *fpTrees, double mutationRate, int doUseObservedCellNames);
void PrintTimes(int replicate, FILE   *fpTimes, double mutationRate, vector<TreeNode *> &nodes,  int thereisOutgroup);
void PrintTimes2(int replicate, FILE  *fpTimes2, double mutationRate,  vector<TreeNode *> &nodes,  int thereisOutgroup);

void PrintTrueFullHaplotypes (FILE *fp, vector<TreeNode *> &nodes, TreeNode* treeRoot, int numNodes, int doPrintIUPAChaplotypes, int doPrintAncestors, int numSites, int numCells, int alphabet, int doUserTree , int doNGS,   char **cellNames, CellStr            *cell, int        HEALTHY_ROOT, int TUMOR_ROOT , char *cellnames[], int doUseObservedCellName);

void ListTimes (int j, double mutationRate, vector<TreeNode *> &nodes, FILE *fpTimes, int thereisOutgroup);
void ListTimes2 (int j,  double mutationRate, vector<TreeNode *> &nodes,  FILE *fpTimes2, int thereisOutgroup);
int Label (TreeNode *p);
int Index (TreeNode *p);


#endif /* output_functions_hpp */
