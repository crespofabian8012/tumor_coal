//
//  data_utils.hpp
//  tumor_coal
//
//  Created by Fausto Fabian Crespo Fernandez on 1/10/19.
//  Copyright © 2019 Fausto Fabian Crespo Fernandez. All rights reserved.
//

#ifndef data_utils_hpp
#define data_utils_hpp

#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <stdarg.h>
#include <search.h>
#include <time.h>

#include "data_types.hpp"
#include "definitions.hpp"
//#include "pllmod_common.h"
#include "eigen.hpp"
#include "population.hpp"
#include "tree_node.hpp"

extern "C"
{
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
#include "libpll/pll_optimize.h"
#include "libpll/pllmod_algorithm.h"
#include "libpll/pllmod_common.h"
}

using namespace std;

//pll_msa_t * pll_phylip_load(const char * fname, pll_bool_t interleaved);

void Initialize( double (*Eij)[4], double (*Mij)[4], double *freq,  ProgramOptions &programOptions );

// reading from configuration file
void ReadParametersFromFile(ProgramOptions &programOptions, FilePaths &filePaths,
                            vector<int> &CloneNameBegin,
                            vector<int> &CloneSampleSizeBegin,
                            vector<int> &ClonePopSizeBegin,
                            vector<double> &CloneBirthRateBegin,
                            vector<double> &CloneDeathRateBegin,
                            vector<double> &CloneTimeOriginInput,
                            double Mij[4][4],
                            double freq[4]);
void ReadUntil(FILE *fv, char stopChar, string what);
void PrintUsage();
void ValidateParameters(ProgramOptions &programOptions,
                        vector<int> CloneNameBegin , vector<int> CloneSampleSizeBegin, vector<int> ClonePopSizeBegin) ;

// initializing the clones
void InitListClones(vector<Population *> &populations, int numClones, int noisy, const vector<int> &CloneNameBegin, const vector<int> &CloneSampleSizeBegin, const vector<double> &CloneBirthRateBegin,  const vector<double> &CloneDeathRateBegin, const vector<int> &ClonePopSizeBegin, const vector<double> &CloneTimeOriginInput, int TotalNumSequences  );
void InitNumberNodes(double &TotalBirthRate, double &TotalDeathRate, int &TotalN,  vector<Population *> &populations, ProgramOptions &programOptions);

bool comparePopulationsByTimeOrigin(const void *s1, const void *s2);
void ListClonesAccordingTimeToOrigin(vector<Population *> &populations, int numClones);

void InitFilesPathsOptions( FilePaths &filePaths, ProgramOptions &programOptions);

int SimulateData(ProgramOptions &programOptions, vector<int> &CloneNameBegin, vector<int> &CloneSampleSizeBegin,
                 vector<int> &ClonePopSizeBegin,
                 vector<Population *> &populations,
                 FilePaths &filePaths,
                 Files &files,
                 double freq[4],
                 double Mij[4][4]
                 );

// helper function for simulating the data -- called by SimulateData
void InitListPossibleMigrations(vector<Population *> &populations, int numClones);
void resetMigrationsList(vector<Population*> &populations, int numClones);
void UpdateListMigrants(Population **populations, int numClones, Population *PopChild, Population *PopFather );
Population* ChooseFatherPopulation(Population **populations, int numClones, Population  *PopChild,  long int *seed, int noisy);
void AssignCurrentSequencesToPopulation(vector<Population *> &populations, vector<TreeNode*> &nodes,
                                        ProgramOptions &programOptions,
                                        int numClones, int numNodes, int noisy,  int TotalNumSequences,
                                        int &numActiveGametes, int &nextAvailable,
                                        int &labelNodes, char* ObservedCellNames[], int doUseObservedCellNames);

int bbinClones (double dat, double *v, int n);
void ChooseRandomIndividual(int *firstInd,   int numClones, Population *popI,  int *secondInd, long *seed, int choosePairIndividuals);
void MakeCoalescenceEvent(vector<Population*> &populations, Population *popI, vector<TreeNode *> &nodes, int numClones, long int* seed, int noisy,   int &numActiveGametes, int &nextAvailable,
                          int &labelNodes, double &currentTime, int &numNodes);
TreeNode *BuildTree(vector<Population* > &populations,
                    Population *CurrentPop,
                    long int *seed,
                    ProgramOptions &programOptions,
                    vector<TreeNode *> &nodes,
                    vector<TreeNode *> &treeTips,
                    TreeNode *treeRootInit,
                    int &nextAvailable,
                    int &newInd,
                    double &currentTime,
                    int &labelNodes);

void PrepareSeparateFiles(int ChainNumber, int paramSetNumber, int replicate,const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files);

void PrepareSeparateFilesGenotypes(int paramSetNumber, int TreeNum,int MutationAssignNum,
                                   const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files);

void InitPopulationsCoalescentEvents( int numClones,  vector<Population *> &populations) ;

void InitializeGenomes (TreeNode *p, long int *seed,  int alphabet, int doUserGenome, int numSites, vector<SiteStr> &allSites, int doGeneticSignatures, double cumfreq[4], double *triNucFreq, char **cellNames);

double SumBranches (TreeNode *p, double mutationRate);

TreeNode *MakeCoalescenceTree2 (long int *seed, vector<Population *> &populations,
                                int &numNodes,
                                int numClones,
                                ProgramOptions &programOptions,
                                double cumNumCA,
                                double meanNumCA,
                                double cumNumMIG,
                                double meanNumMIG,
                                int  &numMIG,
                                int  &numCA,
                                double &numEventsTot,
                                vector<TreeNode *> &nodes,
                                vector<TreeNode *> &treeTips,
                                TreeNode    *treeRootInit
                                ) ;
void SimulatePopulation( Population *popI, vector<Population*> &populations,
                        ProgramOptions &programOptions,
                        long int *seed,
                        int &numNodes,
                        int numClones,
                        double      cumNumCA,
                        double meanNumCA,
                        double cumNumMIG,
                        double meanNumMIG,
                        int  &numMIG,
                        int  &numCA,
                        double &numEventsTot,
                        vector<TreeNode *> &nodes,
                        int &nextAvailable,
                        int &numActiveGametes,
                        int &labelNodes,
                        double &currentTime,
                        int &eventNum);
char * toNewickString2 ( TreeNode *p, double mutationRate,     int doUseObservedCellNames);
void RelabelNodes(TreeNode *p, int &intLabel);

double ProbabilityCloneiFromClonej2 (Population *PopI, Population* PopJ, vector<Population*> &populations, int numClones);

//void connectNodes(TreeNode *left, TreeNode *right, TreeNode *ancester  );
//void setLength(TreeNode *node );

TreeNode *getHealthyTip(TreeNode *treeRootInit);
int openFile(FILE **file, char path[MAX_NAME] );
double  LogConditionalLikelihoodSequences(pll_msa_t * msa, char* NewickString, ProgramOptions &programOptions, double seqError,double dropoutError);

void set_partition_tips_costum( pll_partition_t * partition, pll_msa_t * msa, ProgramOptions &programOptions, double seqError, double dropoutError);
void  destroyTree(pll_utree_t * tree, void (*cb_destroy)(void *));
void dealloc_data_costum(pll_unode_t * node, void (*cb_destroy)(void *));
int set_tipclv_custom_error_model(pll_partition_t * partition,
                                  unsigned int tip_index,
                                  const pll_state_t * map,
                                  const char * sequence,
                                  double _seq_error_rate,
                                  double _dropout_rate);
double LogUniformDensity(double value, double from, double to);
void ReadMCMCParametersFromFile(ProgramOptions &programOptions, FilePaths &filePaths, MCMCoptions &mcmcOptions);
void computeUnfoldedISMSFS(int numSites,vector<SiteStr> &allSites,int numSNVs, vector<int> &SNVsites, vector<int> &SFS, vector<int> &numberDifferences);
int countTrueVariants (vector<TreeNode *> &nodes,  int numSites, int numCells, TreeNode *HEALTHY_ROOT, vector<SiteStr> &allSites, vector<int> &variantSites, vector<int> &SNVsites );
#endif /* data_utils_hpp */
