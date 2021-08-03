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

#include "mcmc_move.hpp"
#include "genotype_error_model.hpp"
#include "treeLikelihood.hpp"

//#include <Eigen/Dense>

#include <map>
#include <unordered_map>


class MCMCmove;
class Chain{
public:
    int  chainNumber;           /* id Chain                                      */
    int  sampleFreq;            /* frequency to sample chain                     */
    int  printFreq;             /* frequency to print chain                      */
    int  swapFreq;              /* frequency to attempt swap of states           */
    int  numRuns;               /* number of runs                                */
       
    bool  converged;
    int   numClones;
    int   oldnumClones;
    bool  endReached;
    bool  samplingFromPosteriorAfterConvergence;
    
    int   gammaParam;
    int   thinning;
    int   currentNumberIerations;
    int   currentNumberEdgesRootedTree;
    int   totalAccepted;
    int   totalRejected;
    int   numberIndependentLongUpdates;
    pll_unode_t *root;
    pll_rnode_t *rootRootedTree;
    pll_rnode_t *healthyTip;
    
    std::string initialValuesString;
    
    std::vector<pll_unode_t*> nodes;
    std::vector<pll_unode_t*> treeTips;
    
    std::vector<pll_rnode_t*> rnodes;
    std::vector<pll_rnode_t*> rtreeTips;
    
    pll_unode_t *oldroot;
    
    std::vector<pll_unode_t *>oldnodes;
    std::vector<pll_unode_t*>  oldtreeTips;
    
 
    std::vector<std::string> paramNames;
    int numNodes;
    int oldnumNodes;
    
    long double mutationRate;
    long double  seqErrorRate;
    long double  dropoutRate;
    long double oldmutationRate;
    long double theta;
    long double oldtheta;
    
    MCMCParameterWithKernel *thetaPar;
    MCMCParameterWithKernel *mutationRatePar;
    MCMCParameterWithKernel *seqErrorRatePar;
    MCMCParameterWithKernel *dropoutRatePar;
    
    std::vector<std::pair<MCMCParameterWithKernel *, int>> indexMCMCParameters;
    int currentIdxMCMCParameters =0;
       
    //Eigen::MatrixXd stored;

    std::vector<Population*> populations;
    
    std::vector<long double  > proportionsVector;
    std::vector<long double  > oldproportionsVector;
    //unsigned long int totalEffectPopSize;//total effect population size
    //unsigned long int oldTotalEffectPopSize;
 
    long double  currentlogConditionalLikelihoodTree;
    long double  currentlogSumDensitiesTimeOriginSTDPopulations;
    long double  currentlogProbFatherPopulations;
    long double  currentlogDensityCoalescentTimesForPopulation;
    long double  currentlogConditionalLikelihoodSequences;
    
    pll_utree_t *initialUnrootedTree;
    pll_rtree_t *initialRootedTree;
    
    std::unordered_map<pll_unode_t*, Population*> tipsAssign;
    std::unordered_map<std::string, pll_unode_t*> labelsAssign;
    
    std::map<pll_rnode_t*, std::vector<Population*>> rMRCAPopulation;
    std::map<pll_rnode_t*, std::vector<Population*>> proposedrMRCAPopulation;
    std::map<pll_rnode_t*, std::vector<Population*>> currentrMRCAPopulation;
    
    std::vector<pll_tree_edge_t *> edges;
    
    std::vector<std::pair<double, pll_tree_edge_t *> > edgeLengths;

    //std::vector<std::vector<double> > sampledPoportionVector;
    Files files;
    std::vector<MCMCmove*>  moves;
    
   ;
    
    std::vector<std::vector<long double>>  posteriorDeltaT;
    std::vector<std::vector<long double>>  posteriorT;
    std::vector<std::vector<long double>>  posteriorProportions;
    std::vector<long double>  posteriorTheta;
    
    
    GenotypeErrorModel *gtErrorModel;
    TreeLikelihood *treeLik;
    //std::ofstream fout;
    //std::vector<pll_edge_node_t*> edges;

    Chain(int chainNumber,
          int numClones,
          int gammaParam,
          long double  mutationRate,
          long double  seqErrorRate,
          long double  dropoutRate,
          MCMCoptions &mcmcOptions
          );
    
    void MakeCoalescenceEvent( Population *Population, std::vector<pll_unode_t *> &nodes, int numClones,const  gsl_rng *randomGenerator, int noisy,   int &numActiveGametes, int &nextAvailable,
                              int &labelNodes, long double  &currentTime, int &numNodes);
    
    
    void  SimulatePopulation( Population *popI,
                                                ProgramOptions &programOptions,
                                                long int *seed,
                                                const gsl_rng *randomGenerator,
                                                int &numNodes,
                                                int numClones,
                                                int &nextAvailable,
                                                int &numActiveGametes,
                                                int &labelNodes,
                                                long double  &currentTime,
                                                int &eventNum,
                                                long double K);
    
    int setInitialTreeFromNewick(char * NewickString);
    
    int setInitialTreeUnrootedTree(pll_utree_t *unrootedTree);
    
    pll_rnode_t* BuildTree(Population *CurrentPop,
                                  const gsl_rng *randomGenerator,
                                  ProgramOptions &programOptions,
                                  pll_rnode_t *tumour_mrca,
                                  int &nextAvailable,
                                  int &newInd,
                                  long double   &currentTime,
                                  int &labelNodes);
    
    pll_rnode_t * MakeCoalescenceTree (long int *seed,
                                       const gsl_rng * randomGenerator,
                                       pll_msa_t *msa,
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
                                       std::vector<std::vector<int> > &ObservedData,
                                       char* ObservedCellNames[],
                                       std::vector<int> &sampleSizes
                                       );
    
    Population* ChooseFatherPopulation( int numClones, Population  *PopChild, const gsl_rng *randomGenerator, int noisy, long double K);
    
    void AssignSequencesToPopulations( ProgramOptions& programOptions,
                                      int numNodes, int noisy,  int TotalNumSequences, int &numActiveGametes, int &nextAvailable,
                                      int &labelNodes, char* SimulatePopulationObservedCellNames[], int doUseObservedCellNames, std::vector<int> &sampleSizes);
    
    void SetPopulationsBirthRate( long double  lambda);
    void GenerateEffectPopSizesFromPriors2( int noisy,     int doGenerateProportionsVector,const gsl_rng* randomGenerator);
    
    void FillChainPopulationsFromPriors( ProgramOptions &programOptions,  MCMCoptions &mcmcOptions, std::vector<int> &sampleSizes, const gsl_rng* randomGenerator, long double K);
    
    void setChainPopulationSampleSizes(std::vector<int > &sampleSizes,  ProgramOptions &programOptions);
    void ListClonesAccordingTimeToOrigin(std::vector<Population *> &populations);
    void GenerateTimesFromPriorsOriginal(int noisy, const gsl_rng * randomGenerator, long double K);
    void InitChainPopulations( int noisy,  int TotalNumSequences, MCMCoptions& mcmcOptions  ) ;
    void RelabelNodes(pll_unode_t *p, int &intLabel);
    void RelabelNodes2(pll_rnode_t *p, int &intLabel);
    void InitPopulationSampleSizes(std::vector<Population*> &populations, int TotalSampleSize, int numClones, std::vector<long double> &proportionsVector,const gsl_rng* randomGenerator);
    long double  SumBranches(pll_unode_t *root, long double  mutationRate);
    char * toNewickString ( pll_unode_t *p, long double  mutationRate,     int doUseObservedCellNames);
    char * toNewickString2 ( pll_rnode_t *p,  std::string healthyTipLabel,    int doUseObservedCellNames);
    long double  LogDensityCoalescentTimesForPopulation(long double K);
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
    
 
    
    static Chain *initializeChain(int chainNumber,  ProgramOptions &programOptions,  MCMCoptions &mcmcOptions, std::vector<int> &sampleSizes,const gsl_rng* randomGenerator, std::vector<std::vector<int>> &ObservedData,
                                  char* ObservedCellNames[], pll_msa_t *msa, pll_rtree_t * initialRootedTree, StructuredCoalescentTree *structCoalTree, std::string& healthyTipLabel, FilePaths &filePaths, Partition *partition);
    
    
   
    
    void stepAllMoves(   MCMCoptions &mcmcOptions, const gsl_rng *randomGenerator,    ProgramOptions &programOptions);
    //void newScaledGrowthRateMoveforPopulation( Population *popI, long int *seed,  ProgramOptions &programOptions, char *ObservedCellNames[], pll_msa_t * msa, MCMCoptions & mcmcOptions, std::vector<int> &sampleSizes);
    void initializeCoalescentEventTimesFormSampleSizes(pll_utree_t *utree, std::vector<int > &sampleSizes, std::string &healthyCellLabel );
    void initializeMapPopulationAssignFromTree();
    Population * getPopulationbyIndex(int indexPopulation);
    int getPopulationIndex(char * label);
    std::map<pll_unode_t*, Population*> chooseTimeOfOriginsOnTree( const gsl_rng * randomGenerator);
    std::map<pll_rnode_t*, std::vector<Population*>>  initTimeOfOriginsOnRootedTree(  std::vector<std::pair<double, pll_tree_edge_t *> > &edgeLengths, int numberPoints, std::string &healthyCellLabel, MCMCoptions &mcmcOptions, const gsl_rng * randomGenerator);
    void initNodeDataFromTree();
    
    void initPopulationSampleSizesFromRootNodeOnTree(pll_unode_t *p, Population *population );
    bool initPopulationsSampleSizes(std::map<pll_rnode_t*, std::vector<Population*>>  &rmrcaOfPopulation, std::string &healthyTipLabel );
    void initPopulationsCoalescentAndMigrationEvents(std::map<pll_unode_t*, Population*> mrcaOfPopulation );
    void initPopulationCoalescentAndMigrationEventsFromRootNodeOnTree(pll_unode_t *p, Population *population, std::map<pll_unode_t*, Population*> &mrcaOfPopulation );
    void initProportionsVector();
    void initProportionsVectorFromSampleSizes(std::vector<int> &sampleSizes);
    void generateProportionsVectorFromDirichlet(long double alpha[], const gsl_rng* randomGenerator);
   
    
    void initPopulationsCoalTimes();
  
    void initPopulationsTipsFromTree(pll_utree_t *utree, bool assignationKnown, std::string &healthyCellLabel );
    static void computeNumberTipsSubTree(pll_unode_t *node, void *data);
    void initPopulationsTipsFromRootedTree(pll_rtree_t *rtree, bool assignationKnown, std::string &healthyCellLabel  );
    void initNodeDataFromRootedTree();
    void initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(pll_rnode_t *p, Population *currentPopulation, std::map<pll_rnode_t*, std::vector<Population*>> &rmrcaOfPopulation,  std::string& healthyTipLabel);
    void initPopulationsCoalescentAndMigrationEventsFromRootedTree(std::map<pll_rnode_t*, std::vector<Population*>> &rmrcaOfPopulation, std::string& healthyTipLabel);
    void initPopulationSampleSizesFromNodeOnRootedTree(pll_rnode_t *p, Population *population, std::map<pll_rnode_t*, std::vector<Population*> > &rmrcaOfPopulation, std::string &healthyTipLabel);
    void  initNumberTipsSubTree(pll_rnode_t *node);
    void initBranches(std::string& healthyCellLabel,std::vector<std::pair<double, pll_tree_edge_t *> > &edgeLengths, std::vector<pll_tree_edge_t *> &edges);
    void initTimeOriginSTDYoungerPopulations(MCMCoptions &mcmcOptions);
    long double  SumBranches2(pll_rnode_t *p, long double  mutationRate);
    void filterSortPopulationsCoalescentEvents();
    void samplePopulationDeltaFromPriors(MCMCoptions &mcmcOptions, long int *seed, const  gsl_rng * rngGsl);
    void rescaleRootedTreeBranchLengths(long double mutationRate);
   
    void newProportionsVectorMove(ProgramOptions &programOptions, char *ObservedCellNames[], pll_msa_t * msa, MCMCoptions &mcmcOptions, std::vector<int> &sampleSizes);
    
    long double   proposalSlidingWindow( long double  oldvalue,  long double  windowSize);
    void proposalProportionsVector(std::vector<long double  > &newProportionvector, long double  tuningParameter );
    long double  DirichletDensity(std::vector<long double> &proportionsVector,  std::vector<long double> &concentrationVector, int sizeVector);
    void updateEffectPopSizesCurrentProportionsVector();
    int totalSampleSize();
    std::map<pll_rnode_t*, std::vector<Population*>> chooseAvailableEdgeOnRootedTreeForPopulation(Population *pop, std::map<pll_rnode_t*, std::vector<Population*>> &mrcaOfPopulation, std::string &healthyCellLabel, const gsl_rng* rngGsl);
    long double  sumAvailableBranchLengths(std::map<pll_rnode_t*, std::vector<Population*>> &currentMRCAPopulation);
    std::map<pll_rnode_t*, std::vector<Population*>> chooseNewTimeofOriginOnEdge(Population *pop, MCMCoptions &mcmcOptions, const gsl_rng* rngGsl, boost::mt19937* rngBoost);
    void PrepareFiles(const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files, int chainNumber);
    void writeMCMCState( int  currentIteration, const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files,  MCMCoptions &mcmcOptions );
   void writeHeaderOutputChain(const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files, const pll_rtree_t* trueTree, const long double trueTheta,
    const std::vector<long double> trueDeltaTs,
    const std::vector<long double> trueTs,
                                      const std::vector<int> trueSampleSizes, StructuredCoalescentTree* tree);
    long double  autoCorrelation(int lag, std::vector<long double> &values, long double meanP, long double varianceP);
    long double  autoCorrelation(std::vector<int> &lags, std::vector<long double> &values, long double correlationThreshold,  int &indexLag, long double meanP, long double varianceP) ;
    long double  ESS(int lag, std::vector<long double> &values);
    void copyProportionsVector(long double  alpha[]);
    void computeAvailableEdges( std::vector<std::pair<double, pll_tree_edge_t *> > &availableEdges, std::map<pll_rnode_t *, Population *> &currentMrcaOfPopulation, std::string &healthyCellLabel);
    bool checkMigrationsOrder();
    long double  LogDensityCoalescentTimesForPopulation2(long double K);
    void rescaleNodeDataFromRootedTree(long double scale);
    long double   computeAdjacentEdges( std::vector<std::pair<double, pll_tree_edge_t *> > &availableEdges, std::map<pll_rnode_t *, std::vector<Population *>> &currentMrcaOfPopulation, std::string &healthyCellLabel, Population *pop, pll_rnode_t * mrca);
    void saveTreeNodeCurrentTimePUnits();
    void rollbackTreeNodeCurrentTimePUnits();
    void updateNodeScaledTimeForRootedTree(long double  newScale);
    Population* getYoungestPopulationOnEdge(pll_rnode_t* p, std::map<pll_rnode_t*, std::vector<Population*> >  &rmrcaOfPopulation);
    void samplePopulationGrowthRateFromPriors(MCMCoptions &mcmcOptions, const gsl_rng* randomGenerator );
    void updateFatherPopOnSameEdge(std::map<pll_rnode_t*, std::vector<Population*> > &rmrcaOfPopulation, pll_rnode_t *p, Population * populationOfCurrentNode);
    long double  sumAdjacentEdges(std::map<pll_rnode_t*, std::vector<Population*>> &currentMRCAPopulation, pll_rnode_t*MRCA, Population *pop);
    bool isOldestPopulation(Population *pop, std::map<pll_rnode_t*, std::vector<Population*>> &currentMrcaOfPopulation);
    void closeFiles(const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files, int chainNumber);
    void updateNodeInfoOldestPopulation(Population * oldestPop, double newMRCAtimeInputTreeUnits);
    void initMutationRate( MCMCoptions &mcmcOptions, ProgramOptions &programOptions, const gsl_rng * randomGenerator);
    void drawModelTimeOriginFromConditionalDensity(Population * oldestPop, MCMCoptions &mcmcOptions, const   gsl_rng * randomGenerator);
    void addOldestPopulation(std::map<pll_rnode_t*, std::vector<Population*> > &mrcaOfPopulation, std::string &healthyCellLabel, MCMCoptions &mcmcOptions);
    void initOriginTimeOldestPopulation( std::string &healthyCellLabel, MCMCoptions &mcmcOptions);
    void  initListMoves(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void printMovesSummary(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void initTreeTips(std::string &healthyTipLabel, MCMCoptions &mcmcOptions,ProgramOptions &programOptions );
  
    void initTreeBranchLengths(std::string &healthyTipLabel);
    
    void initChainTree( std::string &healthyTipLabel, MCMCoptions &mcmcOptions, pll_msa_t *msa, ProgramOptions &programOptions );
    void initLogLikelihoods(MCMCoptions &mcmcOptions, pll_msa_t *msa, ProgramOptions &programOptions, Partition *partition) ;
    long double SumLogDensitiesTimeOriginSTDPopulations() ;
    long double SumLogProbFatherPopulations(long double K) ;
    void initMRCAOldestPopulation(std::string& healthyTipLabel);
    std::vector<long double> initVectorSampleSizes(std::string &healthyTipLabel, MCMCoptions &mcmcOptions, ProgramOptions &programOptions);
    void initPopulationsThetaDelta();
    void  InitListPossibleMigrations();
    void InitPopulationsCoalescentEvents( );
    void InitPopulationSampleSizes(std::vector<int> &sampleSizes);
    void deleteRNodes();
    void addEdgeFromNode(pll_rnode_t *node );
    int setInitialRootedTreeFromNewick(char * NewickString);
    int setInitialTreeRootedTree(pll_rtree_t *rootedTree);
    void InitPopulationGametes();
    void InitPopulationRTips();
    void resetMRCAMap();
    void initEdgesRootedTree(int numberTips);
    void insertMRCAMap(pll_rnode_t   *r, Population *population);
    void initMRCAMap();
    void closeFiles(ProgramOptions &programOptions);
    void saveMCMCState( int  currentIteration, ProgramOptions &programOptions,  MCMCoptions &mcmcOptions);
    void computeThinnig(MCMCoptions & mcmcOptions);
    void resetPosteriorValues();
    void saveTrueTreeInfo(pll_rtree_t * trueTree, char *  rootedNewick3, ProgramOptions& programOptions);
    void closeTrueTreeFiles( ProgramOptions &programOptions);
    void printLastMovesSummary();
    std::vector<int> vectorLags(int valuesSize);
    static long double LogDensityCoalTimes(std::vector<long double> allEventsSorted, std::vector<long double> coalEventTimes, std::vector<long double> immigrantsTimes, long double timeOriginSTD, long double delta, int sampleSize, long double K);
//    void initEdgesRootedTree(int numberTips);
    int numParams();
    int index(const std::string& name) const;
    std::vector<std::vector<long double>> getSamples();
    std::vector<long double> getSampleByName(const std::string& name);
    int getIndexMCMCParameter(MCMCParameterWithKernel *parameter);
    
    //Eigen::MatrixXd getStoredSamples() const;
    void resizeStoredMCMCparameters(int rows, int cols);
    void resizeStoredMCMCparametersAfterWarmUp(MCMCoptions &mcmcOptions);
    int  getLagWithCorrelationCloseToZero(std::vector<long double>& correlationVector, long double correlationThreshold);
    std::vector<long double> getSampleByIndex(const int idx) const;
    std::string getParameterNameByIndex(const int idx) const;
    int getSampleSizeByIndex(const int idx) const;
   ~Chain()
    {
        if (initialRootedTree!=NULL)
        {
            pll_rtree_destroy(initialRootedTree, free );
            //free(initialRootedTree);
        }
        
        for (auto ptr1 : populations)
        {
            delete ptr1;
        }
        populations.clear();
//        for (auto ptr2 : rnodes)
//        {
//            delete (TreeNode*)(ptr2->data);
//            delete ptr2;
//        }
         rnodes.clear();
        for (auto ptr3 : edges)
        {
            delete ptr3;
        }
        edges.clear();
        moves.clear();
        
        //ptr->moves.clear();
       // delete ptr;
      
//        for (std::vector<pll_rnode_t *>::iterator i = rnodes.begin(); i != rnodes.end(); ++i)
//        {
//            delete *i;
//        }
//        //if (healthyTip)
//       //   delete healthyTip;
//        rnodes.clear();
//        rtreeTips.clear();
//
//        for (std::vector<Population *>::iterator i = populations.begin(); i != populations.end(); ++i)
//        {
//            delete *i;
//        }
//        populations.clear();
        
   }
    
private:
    static double * expand_uniq_rates(int states, const double * uniq_rates, const int * rate_sym);
};

#endif /* mcmc_chain_hpp */
