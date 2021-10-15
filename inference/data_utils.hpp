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
 * data utils functions
 */


#ifndef data_utils_hpp
#define data_utils_hpp


#include "data_types.hpp"
#include "eigen.hpp"
#include "population.hpp"
#include "tree_node.hpp"
#include "mutationModel.h"

#include <gsl/gsl_rng.h>
#include <boost/random.hpp>

#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <stdarg.h>
#include <search.h>

extern "C"
{
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
#include "libpll/pll_optimize.h"
#include "libpll/pllmod_algorithm.h"
#include "libpll/pllmod_common.h"
}

void Initialize( double (*Eij)[4], double (*Mij)[4], double *freq,  ProgramOptions &programOptions );

// reading from configuration file
void ReadParametersFromFile(ProgramOptions &programOptions, FilePaths &filePaths,
                            std::vector<int> &CloneNameBegin,
                            std::vector<int> &CloneSampleSizeBegin,
                            std::vector<int> &ClonePopSizeBegin,
                            std::vector<double> &CloneBirthRateBegin,
                            std::vector<double> &CloneDeathRateBegin,
                            std::vector<double> &CloneTimeOriginInput,
                            double Mij[4][4],
                            double freq[4]);
void ReadUntil(FILE *fv, char stopChar, std::string what);
void PrintUsage();
void ValidateParameters(ProgramOptions &programOptions,
                        std::vector<int> CloneNameBegin , std::vector<int> CloneSampleSizeBegin, std::vector<int> ClonePopSizeBegin) ;

// initializing the clones
void InitListClones(std::vector<Population *> &populations, int numClones, int noisy, const std::vector<int> &CloneNameBegin, const std::vector<int> &CloneSampleSizeBegin, const std::vector<double> &CloneBirthRateBegin,  const std::vector<double> &CloneDeathRateBegin, const std::vector<int> &ClonePopSizeBegin, const std::vector<double> &CloneTimeOriginInput, int &TotalNumSequences,int  doEstimateTimesOriginClones , ProgramOptions &programOptions,  std::vector<gsl_rng *> rngGslvector);
void InitNumberNodes(double &TotalBirthRate, double &TotalDeathRate, int &TotalN,  std::vector<Population *> &populations, ProgramOptions &programOptions);

bool comparePopulationsByTimeOrigin(const void *s1, const void *s2);
void ListClonesAccordingTieToOrigin(std::vector<Population *> &populations, int numClones);

void InitFilesPathsOptions( FilePaths &filePaths, ProgramOptions &programOptions);

int SimulateData(ProgramOptions &programOptions, std::vector<int> &CloneNameBegin, std::vector<int> &CloneSampleSizeBegin,
                 std::vector<int> &ClonePopSizeBegin,
                 std::vector<Population *> &populations,
                 FilePaths &filePaths,
                 Files &files,
                  double freq[4],
                  double Mij[4][4],
                 double Eij[4][4],
                 std::vector<gsl_rng *> &rngGslvector,
                 std::vector<boost::random::mt19937 *> &rngBoostvector
                 );

// helper function for simulating the data -- called by SimulateData
void InitListPossibleMigrations(std::vector<Population *> &populations, int numClones);
void resetMigrationsList(std::vector<Population*> &populations, int numClones);
void UpdateListMigrants(Population **populations, int numClones, Population *PopChild, Population *PopFather );
Population* ChooseFatherPopulation(Population **populations, int numClones, Population  *PopChild,  long int *seed, int noisy);
void AssignCurrentSequencesToPopulation(std::vector<Population *> &populations, std::vector<TreeNode*> &nodes,
                                        ProgramOptions &programOptions,
                                        int numClones, int numNodes, int noisy,  int TotalNumSequences,
                                        int &numActiveGametes, int &nextAvailable,
                                        int &labelNodes, char* ObservedCellNames[], int doUseObservedCellNames);

int bbinClones (long double dat, long double *v, int n);
void ChooseRandomIndividual(int *firstInd,   int numClones, Population *popI,  int *secondInd, long *seed, int choosePairIndividuals, const gsl_rng * rngGsl);
void MakeCoalescenceEvent(std::vector<Population*> &populations, Population *popI, std::vector<TreeNode *> &nodes, int numClones, long int* seed, int noisy,   int &numActiveGametes, int &nextAvailable,
                          int &labelNodes, double &currentTime, int &numNodes,const gsl_rng * rngGsl);
TreeNode *BuildTree(std::vector<Population* > &populations,
                    Population *CurrentPop,
                    long int *seed,
                    ProgramOptions &programOptions,
                    std::vector<TreeNode *> &nodes,
                    std::vector<TreeNode *> &treeTips,
                    TreeNode *treeRootInit,
                    int &nextAvailable,
                    int &newInd,
                    double &currentTime,
                    int &labelNodes,
                    const gsl_rng *rngGsl);

void PrepareSeparateFiles(int ChainNumber, int paramSetNumber, int replicate,const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files, std::vector<Population*> &populations);

void PrepareSeparateFilesGenotypes(int paramSetNumber, int TreeNum,int MutationAssignNum,
                                   const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files, std::vector<Population*> populations, double numMU);

void PrepareLikelihoodOutputFile(const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files);

void InitPopulationsCoalescentEvents( int numClones,  std::vector<Population *> &populations) ;

void writeHeaderLikelihoodFile(const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files, int numClones );

void writeLineLikelihoodFile( int  simulationNumber, const FilePaths &filePaths, const ProgramOptions    &programOptions,Files &files ,  std::vector<Population *> &populations,
     long double logLikCoalTree, long double logLikTrueSequences, long double logLikErrorSequences, double treeLength);

void InitializeGenomes (TreeNode *p, long int *seed,  int alphabet, int doUserGenome, int numSites, std::vector<SiteStr> &allSites, int doGeneticSignatures, double cumfreq[4], long double *triNucFreq, char **cellNames,
                       const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost);

double SumBranches (TreeNode *p, double mutationRate, std::string &healthyTipLabel);

TreeNode *MakeCoalescenceTree2 (long int *seed, std::vector<Population *> &populations,
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
                                std::vector<TreeNode *> &nodes,
                                std::vector<TreeNode *> &treeTips,
                                TreeNode    *treeRootInit,
                                long double K,
                                const gsl_rng *rngGsl,
                                boost::random::mt19937 * rngBoost
                                ) ;
void SimulatePopulation( Population *popI, std::vector<Population*> &populations,
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
                        std::vector<TreeNode *> &nodes,
                        int &nextAvailable,
                        int &numActiveGametes,
                        int &labelNodes,
                        double &currentTime,
                        int &eventNum,
                         long double K,
                        const gsl_rng * rngGsl,
                        boost::mt19937* rngBoost);
char * toNewickString2 ( TreeNode *p, double mutationRate,     int doUseObservedCellNames);
void RelabelNodes(TreeNode *p, int &intLabel);
void ListClonesAccordingTimeToOrigin(std::vector<Population *> &populations, int numClones);

double ProbabilityCloneiFromClonej2 (Population *PopI, Population* PopJ, std::vector<Population*> &populations, int numClones, long double K);

//void connectNodes(TreeNode *left, TreeNode *right, TreeNode *ancester  );
//void setLength(TreeNode *node );

TreeNode *getHealthyTip(TreeNode *treeRootInit);
int openFile(FILE **file, char path[MAX_NAME] );


void ReadMCMCParametersFromFile(ProgramOptions &programOptions, FilePaths &filePaths, MCMCOptions &mcmcOptions);
void computeUnfoldedISMSFS(int numSites,std::vector<SiteStr> &allSites,int numSNVs, std::vector<int> &SNVsites, std::vector<int> &SFS, std::vector<int> &numberDifferences);
int countTrueVariants (std::vector<TreeNode *> &nodes,  int numSites, int numCells, TreeNode *HEALTHY_ROOT, std::vector<SiteStr> &allSites, std::vector<int> &variantSites, std::vector<int> &SNVsites );


long double computeParamPowerDistribQuantileUntil(long double areaUntilb, long double b, long double from);
void setDefaultOptions(ProgramOptions &programOptions, MCMCOptions &mcmcOptions );
void printProgramHeader();
void printSimulatorProgramHeader();

long double  initMutationRate( MCMCOptions &mcmcOptions, ProgramOptions &programOptions, const gsl_rng * randomGenerator);
long double  sampleMutationRateSimulation( MCMCOptions &mcmcOptions, ProgramOptions &programOptions,const gsl_rng *randomGsl, boost::mt19937* rngBoost) ;
void simulateTrees(int numberTrees,std::vector<StructuredCoalescentTree *> &structuredCoalTrees,  std::vector<pll_rtree_t *> &trees,          std::vector<long double> &realThetas,
                   std::vector<std::vector<long double>> &realDeltaTs,
                   std::vector<std::vector<long double>> &realTs,
                   std::vector<int> & sampleSizes, ProgramOptions &programOptions, MCMCOptions & mcmcOptions, std::vector<gsl_rng * > rngGsl,std::vector<boost::mt19937* > rngBoost, std::vector<std::vector<int> > &ObservedData,char* ObservedCellNames[], pll_msa_t *msa, std::string& healthyTipLabel,
                   long double seqErrorRate,
                   long double dropoutRate );
void InitFiles(Files &files);
void SetPopulationTimeOriginSTD(std::vector<Population *> &populations, int numClones, const gsl_rng* rngGsl, bool doEstimateTorigins);
void InitNumberNodes( std::vector<Population *> &populations, ProgramOptions &programOptions);
void computeStatisticsNumberMutations(std::vector<SiteStr> allSites, long double &meanMaternalMutationPerSite, long double &meanPaternalMutationPerSite);
void SetPopulationParametersFromPriors(std::vector<Population *> &populations, int numClones,const gsl_rng* rngGsl, ProgramOptions &programOptions);
#endif /* data_utils_hpp */
