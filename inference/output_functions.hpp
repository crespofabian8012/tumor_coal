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
#include <stdio.h>

#include "data_types.hpp"
#include "data_utils.hpp"
#include "tree_node.hpp"


//struct CellStr;
//struct SiteStr;
//class  TreeNode;


struct Output
{
static TreeNode *getHealthyTip(TreeNode *treeRootInit);

static void PrintUsage();
static void WriteTree2 ( TreeNode *p, double mutationRate, FILE    *fpTrees2, char *cellNames[], int *indexCurrentCell, int doUseObservedCellNames);
static void WriteTree (TreeNode *p, double mutationRate, FILE    *fpTrees, int doUseObservedCellNames);
static void PrintTrees2(int replicate, TreeNode *treeRootInit,   FILE   *fpTrees2 , double mutationRate,char * ObservedCellNames[],int doUseObservedCellNames);
static void PrintTrees(int replicate, TreeNode *treeRootInit,   FILE  *fpTrees, double mutationRate, int doUseObservedCellNames);
static void PrintTimes(int replicate, FILE   *fpTimes, double mutationRate, std::vector<TreeNode *> &nodes,  int thereisOutgroup);
static void PrintTimes2(int replicate, FILE  *fpTimes2, double mutationRate,  std::vector<TreeNode *> &nodes,  int thereisOutgroup);

static void PrintTrueFullHaplotypes (FILE *fp, std::vector<TreeNode *> &nodes, TreeNode* treeRoot, int numNodes, int doPrintIUPAChaplotypes, int doPrintAncestors, int numSites, int numCells, int alphabet, int doUserTree , int doNGS,   char **cellNames, CellStr            *cell, int        HEALTHY_ROOT, int TUMOR_ROOT , char *cellnames[], int doUseObservedCellName);

static void ListTimes (int j, double mutationRate, std::vector<TreeNode *> &nodes, FILE *fpTimes, int thereisOutgroup);
static void ListTimes2 (int j,  double mutationRate, std::vector<TreeNode *> &nodes,  FILE *fpTimes2, int thereisOutgroup);
static int Label (TreeNode *p);
static int Index (TreeNode *p);
static int Index (pll_rnode_t *p);
static void PrintTrees(int replicate, pll_rnode_t *treeRootInit,   FILE  *fpTrees, double mutationRate, int doUseObservedCellNames);
static void PrintTrees2(int replicate, pll_rnode_t *treeRootInit,   FILE   *fpTrees2 , double mutationRate,char * ObservedCellNames[],int doUseObservedCellNames);
static void WriteTree (pll_rnode_t *p, double mutationRate, FILE    *fpTrees, int doUseObservedCellNames);
static void WriteTree2 ( pll_rnode_t *p, double mutationRate, FILE    *fpTrees2, char *cellNames[], int *indexCurrentCell, int doUseObservedCellNames);
static void PrintTimes(int replicate, FILE   *fpTimes, double mutationRate, std::vector<pll_rnode_t *> &nodes,  int thereisOutgroup);
static void PrintTimes2(int replicate, FILE  *fpTimes2, double mutationRate,  std::vector<pll_rnode_t *> &nodes,  int thereisOutgroup);
static void ListTimes (int j, double mutationRate, std::vector<pll_rnode_t *> &nodes, FILE *fpTimes, int thereisOutgroup);
static void ListTimes2 (int j,  double mutationRate, std::vector<pll_rnode_t *> &nodes,  FILE *fpTimes2, int thereisOutgroup);
//static   void writetextToFile(FilePath &filePath,  char * text);
    
};
#endif /* output_functions_hpp */
