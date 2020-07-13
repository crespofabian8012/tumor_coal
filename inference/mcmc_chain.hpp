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
 * MCMC chain class
 */
#ifndef mcmc_chain_hpp
#define mcmc_chain_hpp


#include <map>
#include <unordered_map>
#include "data_utils.hpp"

//class MCMCmove;
//class Population;
//class ProgramOptions;
//class MCMCoptions;


class MCMCmove;
class Chain{
public:
    int chainNumber;
    int numClones;
    int oldnumClones;
    int gammaParam;
    int currentNumberIerations;
    int totalAccepted;
    int totalRejected;
    pll_unode_t *root;
    pll_rnode_t *rootRootedTree;
    std::vector<pll_unode_t*> nodes;
    std::vector<pll_unode_t*> treeTips;
    
    std::vector<pll_rnode_t*> rnodes;
    std::vector<pll_rnode_t*> rtreeTips;
    
    pll_unode_t *oldroot;
    std::vector<pll_unode_t *>oldnodes;
    std::vector<pll_unode_t*>  oldtreeTips;
    
    int numNodes;
    int oldnumNodes;
    
    long double mutationRate;
    long double  seqErrorRate;
    long double  dropoutRate;
    long double oldmutationRate;
    long double theta;
    long double oldtheta;
    
    std::vector<Population*> populations;
    
    std::vector<long double  > proportionsVector;
    std::vector<long double  > oldproportionsVector;
    unsigned long int totalEffectPopSize;//total effect population size
    unsigned long int oldTotalEffectPopSize;
    //long double  lambda ;
    long double  currentlogConditionalLikelihoodTree;
    long double  currentlogConditionalLikelihoodSequences;
    
    pll_utree_t *initialUnrootedTree;
    pll_rtree_t *initialRootedTree;
    
    std::unordered_map<pll_unode_t*, Population*> tipsAssign;
    std::unordered_map<std::string, pll_unode_t*> labelsAssign;
    
    std::map<pll_rnode_t*, std::vector<Population*>> rMRCAPopulation;
    std::map<pll_rnode_t*, std::vector<Population*>> proposedrMRCAPopulation;
    std::map<pll_rnode_t*, std::vector<Population*>> currentrMRCAPopulation;
    
    std::vector<pll_tree_edge_t *> edges;
    
    std::vector<pair<double, pll_tree_edge_t *> > edgeLengths;
    
    std::vector<double> sampledTotalEffectPopSize;
    std::vector<std::vector<double> > sampledPoportionVector;
    Files files;
    vector<MCMCmove*>  moves;
    //std::vector<pll_edge_node_t*> edges;
public:
    
    Chain(int chainNumber,
          int numClones,
          int gammaParam,
          long double  mutationRate,
          long double  seqErrorRate,
          long double  dropoutRate
          );
    
    void MakeCoalescenceEvent( Population *Population, std::vector<pll_unode_t *> &nodes, int numClones, long int* seed, int noisy,   int &numActiveGametes, int &nextAvailable,
                              int &labelNodes, long double  &currentTime, int &numNodes);
    
    
    void SimulatePopulation( Population *popI, std::vector<Population*> &populations,
                                   ProgramOptions &programOptions,
                                   long int *seed,
                                   int &numNodes,
                                   int numClones,
                                   long double       cumNumCA,
                                   long double  meanNumCA,
                                   long double  cumNumMIG,
                                   long double  meanNumMIG,
                                   int  &numMIG,
                                   int  &numCA,
                                   long double  &numEventsTot,
                                   std::vector<pll_unode_t *> &nodes,
                                   int &nextAvailable,
                                   int &numActiveGametes,
                                   int &labelNodes,
                                   long double  &currentTime,
                            int &eventNum);
    
    int setInitialTreeFromNewick(char * NewickString);
    
    int setInitialTreeUnrootedTree(pll_utree_t *unrootedTree);
    
    pll_unode_t* BuildTree(std::vector<Population* > &populations,
                           Population *CurrentPop,
                           long int *seed,
                           ProgramOptions &programOptions,
                           std::vector<pll_unode_t *> &nodes,
                           std::vector<pll_unode_t *> &treeTips,
                           pll_unode_t *tumour_mrca,
                           int &nextAvailable,
                           int &newInd,
                           long double  &currentTime,
                           int &labelNodes);
    
    pll_unode_t * MakeCoalescenceTree (long int *seed,
                                       int &numNodes,
                                       int numClones,
                                       ProgramOptions &programOptions,
                                       long double       cumNumCA,
                                       long double  meanNumCA,
                                       long double  cumNumMIG,
                                       long double  meanNumMIG,
                                       int  &numMIG,
                                       int  &numCA,
                                       long double  &numEventsTot,
                                       char* ObservedCellNames[],
                                       std::vector<int> &sampleSizes
                                       );
    
    Population* ChooseFatherPopulation( int numClones, Population  *PopChild,  long int *seed, int noisy);
    
    void AssignSequencesToPopulations( ProgramOptions& programOptions,
                                      int numNodes, int noisy,  int TotalNumSequences, int &numActiveGametes, int &nextAvailable,
                                      int &labelNodes, char* SimulatePopulationObservedCellNames[], int doUseObservedCellNames, std::vector<int> &sampleSizes);
    
    void SetPopulationsBirthRate( long double  lambda);
    void GenerateEffectPopSizesFromPriors2( int noisy,     int doGenerateProportionsVector);
    
    void FillChainPopulationsFromPriors( ProgramOptions &programOptions,  MCMCoptions &mcmcOptions, std::vector<int> &sampleSizes, long int *seed);
    
    void setChainPopulationSampleSizes(std::vector<int > &sampleSizes,  ProgramOptions &programOptions);
    void ListClonesAccordingTimeToOrigin(std::vector<Population *> &populations);
    void GenerateTimesFromPriorsOriginal(int noisy,  long int *seed);
    void InitChainPopulations( int noisy,  int TotalNumSequences  ) ;
    void RelabelNodes(pll_unode_t *p, int &intLabel);
    void InitPopulationSampleSizes(std::vector<Population*> populations, int TotalSampleSize, int numClones, std::vector<long double> &proportionsVector, long int *seed);
    long double  SumBranches(pll_unode_t *root, long double  mutationRate);
    char * toNewickString ( pll_unode_t *p, long double  mutationRate,     int doUseObservedCellNames);
    
    long double  LogDensityCoalescentTimesForPopulation();
    long double  LogConditionalLikelihoodTree( ProgramOptions &programOptions, MCMCoptions &mcmcOptions  );
    long double   LogConditionalLikelihoodSequences(pll_msa_t * msa, char* NewickString, ProgramOptions &programOptions, long double  seqError,long double  dropoutError);
    int set_tipclv_custom_error_model(pll_partition_t * partition,
                                      unsigned int tip_index,
                                      const pll_state_t * map,
                                      const char * sequence,
                                      long double  _seq_error_rate,
                                      long double  _dropout_rate);
    void set_partition_tips_costum( pll_partition_t * partition, pll_msa_t * msa, ProgramOptions &programOptions, long double  seqError, long double  dropoutError);
    void dealloc_data_costum(pll_unode_t * node, void (*cb_destroy)(void *));
    void  destroyTree(pll_utree_t * tree, void (*cb_destroy)(void *));
    static void extracted(Chain *&chain, MCMCoptions &mcmcOptions, ProgramOptions &programOptions);
    
    static Chain *initializeChain( ProgramOptions &programOptions,  MCMCoptions &mcmcOptions, std::vector<int> &sampleSizes, long int *seed, char* ObservedCellNames[], pll_msa_t *msa, char* treefileName, string& healthyTipLabel);
    
    
    vector<MCMCmove *> extracted();
    
    void runChain(   MCMCoptions &opt,  long int *seed,  FilePaths &filePaths, Files &files,  ProgramOptions &programOptions,
                  char* ObservedCellNames[], pll_msa_t * msa, std::vector<int> &sampleSizes, int currentIteration
                  );
    void newScaledGrowthRateMoveforPopulation( Population *popI, long int *seed,  ProgramOptions &programOptions, char *ObservedCellNames[], pll_msa_t * msa, MCMCoptions & mcmcOptions, std::vector<int> &sampleSizes);
    void initializeCoalescentEventTimesFormSampleSizes(pll_utree_t *utree, std::vector<int > &sampleSizes, string &healthyCellLabel );
    void initializeMapPopulationAssignFromTree();
    Population * getPopulationbyIndex(int indexPopulation);
    int getPopulationIndex(char * label);
    std::map<pll_unode_t*, Population*> chooseTimeOfOriginsOnTree( long int *seed);
    std::map<pll_rnode_t*, std::vector<Population*>>  initTimeOfOriginsOnRootedTree(  std::vector<pair<double, pll_tree_edge_t *> > edgeLengths, int numberPoints, string &healthyCellLabel, MCMCoptions &mcmcOptions);
    void initNodeDataFromTree();
    
    void initPopulationSampleSizesFromRootNodeOnTree(pll_unode_t *p, Population *population );
    bool initPopulationsSampleSizes(std::map<pll_rnode_t*, std::vector<Population*>>  rmrcaOfPopulation, string &healthyTipLabel );
    void initPopulationsCoalescentAndMigrationEvents(std::map<pll_unode_t*, Population*> mrcaOfPopulation );
    void initPopulationCoalescentAndMigrationEventsFromRootNodeOnTree(pll_unode_t *p, Population *population, std::map<pll_unode_t*, Population*> mrcaOfPopulation );
    void initProportionsVector();
    void generateProportionsVectorFromDirichlet(long double alpha[]);
    void initPopulationMigration();
    void initTotalEffectivePopulationSize(MCMCoptions &mcmcOptions, long int *seed);
    void initPopulationsCoalTimes();
    void initEffectPopulationSizesFromProportionsVector();
    void initPopulationsTipsFromTree(pll_utree_t *utree, bool assignationKnown, string &healthyCellLabel );
    static void computeNumberTipsSubTree(pll_unode_t *node, void *data);
    void initPopulationsTipsFromRootedTree(pll_rtree_t *rtree, bool assignationKnown, string &healthyCellLabel  );
    void initNodeDataFromRootedTree();
    void initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(pll_rnode_t *p, Population *currentPopulation, std::map<pll_rnode_t*, std::vector<Population*>> &rmrcaOfPopulation,  string& healthyTipLabel);
    void initPopulationsCoalescentAndMigrationEventsFromRootedTree(std::map<pll_rnode_t*, std::vector<Population*>> rmrcaOfPopulation, string& healthyTipLabel);
    void initPopulationSampleSizesFromNodeOnRootedTree(pll_rnode_t *p, Population *population, std::map<pll_rnode_t*, std::vector<Population*> > rmrcaOfPopulation, string &healthyTipLabel);
    void  initNumberTipsSubTree(pll_rnode_t *node);
    void initBranches(string& healthyCellLabel,std::vector<pair<double, pll_tree_edge_t *> > &edgeLengths, std::vector<pll_tree_edge_t *> &edges);
    void initTimeOriginSTD(MCMCoptions &mcmcOptions);
    long double  SumBranches2(pll_rnode_t *p, long double  mutationRate);
    void filterSortPopulationsCoalescentEvents();
    void samplePopulationDeltaFromPriors(MCMCoptions &mcmcOptions, long int *seed );
    void rescaleRootedTreeBranchLengths(long double mutationRate);
    void newTotalEffectivePopulationSizeMove( ProgramOptions &programOptions, char *ObservedCellNames[], pll_msa_t * msa, MCMCoptions &mcmcOptions, std::vector<int> &sampleSizes);
    void newProportionsVectorMove(ProgramOptions &programOptions, char *ObservedCellNames[], pll_msa_t * msa, MCMCoptions &mcmcOptions, std::vector<int> &sampleSizes);
    
    long double   proposalSlidingWindow( long double  oldvalue,  long double  windowSize);
    void proposalProportionsVector(std::vector<long double  > &newProportionvector, long double  tuningParameter );
    long double  DirichletDensity(std::vector<long double> &proportionsVector,  std::vector<long double> &concentrationVector, int sizeVector);
    void updateEffectPopSizesCurrentProportionsVector();
    int totalSampleSize();
    std::map<pll_rnode_t*, std::vector<Population*>> chooseAvailableEdgeOnRootedTreeForPopulation(Population *pop, std::map<pll_rnode_t*, std::vector<Population*>> &mrcaOfPopulation, string &healthyCellLabel);
    long double  sumAvailableBranchLengths(std::map<pll_rnode_t*, std::vector<Population*>> currentMRCAPopulation);
    std::map<pll_rnode_t*, std::vector<Population*>> chooseNewTimeofOriginOnEdge(Population *pop, MCMCoptions &mcmcOptions);
    void PrepareFiles(const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files, int chainNumber);
    void writeMCMCState( int  currentIteration, const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files,  MCMCoptions &mcmcOptions );
    void writeHeaderOutputChain(  const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files );
    long double  autoCorrelation(int lag, std::vector<double> values);
    long double  ESS(int lag, std::vector<double> values);
    void copyProportionsVector(long double  alpha[]);
    void computeAvailableEdges( std::vector<pair<double, pll_tree_edge_t *> > &availableEdges, std::map<pll_rnode_t *, Population *> &currentMrcaOfPopulation, std::string &healthyCellLabel);
    bool checkMigrationsOrder();
    long double  LogDensityCoalescentTimesForPopulation2();
    void rescaleNodeDataFromRootedTree(long double scale);
    long double   computeAdjacentEdges( std::vector<pair<double, pll_tree_edge_t *> > &availableEdges, std::map<pll_rnode_t *, std::vector<Population *>> &currentMrcaOfPopulation, std::string &healthyCellLabel, Population *pop, pll_rnode_t * mrca);
    void safeTreeNodeCurrentTimePUnits();
    void rollbackTreeNodeCurrentTimePUnits();
    void updateNodeScaledTimeForRootedTree(long double  newScale);
    Population* getYoungestPopulationOnEdge(pll_rnode_t* p, std::map<pll_rnode_t*, std::vector<Population*> >  rmrcaOfPopulation);
    void samplePopulationGrowthRateFromPriors(MCMCoptions &mcmcOptions, long int *seed );
    void updateFatherPopOnSameEdge(std::map<pll_rnode_t*, std::vector<Population*> > &rmrcaOfPopulation, pll_rnode_t *p, Population * populationOfCurrentNode);
    long double  sumAdjacentEdges(std::map<pll_rnode_t*, std::vector<Population*>> currentMRCAPopulation, pll_rnode_t*MRCA, Population *pop);
    bool isOldestPopulation(Population *pop, std::map<pll_rnode_t*, std::vector<Population*>> &currentMrcaOfPopulation);
    void closeFiles(const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files, int chainNumber);
    void updateNodeInfoOldestPopulation(Population * oldestPop, double newMRCAtimeInputTreeUnits);
    void initMutationRate( MCMCoptions &mcmcOptions, ProgramOptions &programOptions);
    void drawModelTimeOriginFromConditionalDensity(Population * oldestPop, MCMCoptions &mcmcOptions);
    void addOldestPopulation(std::map<pll_rnode_t*, std::vector<Population*> > &mrcaOfPopulation, string &healthyCellLabel, MCMCoptions &mcmcOptions);
    void initOriginTimeOldestPopulation( string &healthyCellLabel, MCMCoptions &mcmcOptions);
    void  initListMoves();
    void printMovesSummary();
    
private:
    static double * expand_uniq_rates(int states, const double * uniq_rates, const int * rate_sym);
};

#endif /* mcmc_chain_hpp */
