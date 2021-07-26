/*################################################################################
 ##
 ##   Copyright (C) 2018-2020 Fausto Fabian Crespo Fernandez
 ##
 ##   This file is part of the tumor_coal C++ library.
 ##
 ##   Licensed under the Apache License, Version 2.0 (the "License");
 ##   you may not use this file except in compliance with the License.
 ##   You may obtain a copy of the License at
 ##
 ##       http://www.apache.org/licenses/LICENSE-2.0
 ##
 ##   Unless required by applicable law or agreed to in writing, software
 ##   distributed under the License is distributed on an "AS IS" BASIS,
 ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 ##   See the License for the specific language governing permissions and
 ##   limitations under the License.
 ##
 ################################################################################*/
/*
 * output functions
 */
#ifndef output_functions_hpp
#define output_functions_hpp


#include "data_types.hpp"
#include "data_utils.hpp"
#include "tree_node.hpp"

#include <stdio.h>
namespace Output
{
 TreeNode *getHealthyTip(TreeNode *treeRootInit);

 void PrintUsage();
 void WriteTree2 ( TreeNode *p, double mutationRate, FILE    *fpTrees2, char *cellNames[], int *indexCurrentCell, int doUseObservedCellNames);
 void WriteTree (TreeNode *p, double mutationRate, FILE    *fpTrees, int doUseObservedCellNames);
 void PrintTrees2(int replicate, TreeNode *treeRootInit,   FILE   *fpTrees2 , double mutationRate,char * ObservedCellNames[],int doUseObservedCellNames);
 void PrintTrees(int replicate, TreeNode *treeRootInit,   FILE  *fpTrees, double mutationRate, int doUseObservedCellNames);
 void PrintTimes(int replicate, FILE   *fpTimes, double mutationRate, std::vector<TreeNode *> &nodes,  int thereisOutgroup);
 void PrintTimes2(int replicate, FILE  *fpTimes2, double mutationRate,  std::vector<TreeNode *> &nodes,  int thereisOutgroup);

 void PrintTrueFullHaplotypes (FILE *fp, std::vector<TreeNode *> &nodes, TreeNode* treeRoot, int numNodes, int doPrintIUPAChaplotypes, int doPrintAncestors, int numSites, int numCells, int alphabet, int doUserTree , int doNGS,   char **cellNames, CellStr            *cell, int        HEALTHY_ROOT, int TUMOR_ROOT , char *cellnames[], int doUseObservedCellName);

void PrintSNVGenotypes (FILE *fp, std::vector<TreeNode *> &nodes, TreeNode* treeRoot, int numNodes, int doPrintIUPAChaplotypes, int doPrintAncestors, int numSites, int numCells, int alphabet, int doUserTree , int doNGS,   char **cellNames, CellStr    *cell, int        HEALTHY_ROOT, int TUMOR_ROOT , char *cellnames[], int doUseObservedCellName,
                                int numSNVs, std::vector<int> &SNVsites);

void PrintFullGenotypes(FILE *fp, std::vector<TreeNode *> &nodes, TreeNode* treeRoot, int numNodes, int doPrintIUPAChaplotypes, int doPrintAncestors, int numSites, int numCells, int alphabet, int doUserTree , int doNGS,   char **cellNames, CellStr    *cell, int        HEALTHY_ROOT, int TUMOR_ROOT , char *cellnames[], int doUseObservedCellName,
int numSNVs, std::vector<int> &SNVsites);

 void ListTimes (int j, double mutationRate, std::vector<TreeNode *> &nodes, FILE *fpTimes, int thereisOutgroup);
 void ListTimes2 (int j,  double mutationRate, std::vector<TreeNode *> &nodes,  FILE *fpTimes2, int thereisOutgroup);
 int Label (TreeNode *p);
 int Index (TreeNode *p);
 int Index (pll_rnode_t *p);
 void PrintTrees(int replicate, pll_rnode_t *treeRootInit,   FILE  *fpTrees, double mutationRate, int doUseObservedCellNames);
 void PrintTrees2(int replicate, pll_rnode_t *treeRootInit,   FILE   *fpTrees2 , double mutationRate,char * ObservedCellNames[],int doUseObservedCellNames);
 void WriteTree (pll_rnode_t *p, double mutationRate, FILE    *fpTrees, int doUseObservedCellNames);
 void WriteTree2 ( pll_rnode_t *p, double mutationRate, FILE    *fpTrees2, char *cellNames[], int *indexCurrentCell, int doUseObservedCellNames);
 void PrintTimes(int replicate, FILE   *fpTimes, double mutationRate, std::vector<pll_rnode_t *> &nodes,  int thereisOutgroup);
 void PrintTimes2(int replicate, FILE  *fpTimes2, double mutationRate,  std::vector<pll_rnode_t *> &nodes,  int thereisOutgroup);
 void ListTimes (int j, double mutationRate, std::vector<pll_rnode_t *> &nodes, FILE *fpTimes, int thereisOutgroup);
 void ListTimes2 (int j,  double mutationRate, std::vector<pll_rnode_t *> &nodes,  FILE *fpTimes2, int thereisOutgroup);
//static   void writetextToFile(FilePath &filePath,  char * text);
    
};
#endif /* output_functions_hpp */
