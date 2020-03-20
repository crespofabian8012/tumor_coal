//
//  mcmc_chain.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 7/10/19.
//

#ifndef mcmc_chain_hpp
#define mcmc_chain_hpp

#include <map>
#include <unordered_map>
#include "data_utils.hpp"

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
    vector<pll_unode_t*> nodes;
    vector<pll_unode_t*> treeTips;
    
    vector<pll_rnode_t*> rnodes;
    vector<pll_rnode_t*> rtreeTips;
    
    pll_unode_t *oldroot;
    vector<pll_unode_t *>oldnodes;
    vector<pll_unode_t*>  oldtreeTips;
    
    int numNodes;
    int oldnumNodes;
    
    long double mutationRate;
    long double  seqErrorRate;
    long double  dropoutRate;
    long double oldmutationRate;
    long double theta;
    long double oldtheta;
    
    vector<Population*> populations;
    
    vector<long double  > proportionsVector;
    vector<long double  > oldproportionsVector;
    unsigned long int totalEffectPopSize;//total effect population size
    unsigned long int oldTotalEffectPopSize;
    //long double  lambda ;
    long double  currentlogConditionalLikelihoodTree;
    long double  currentlogConditionalLikelihoodSequences;
    
    pll_utree_t *initialUnrootedTree;
    pll_rtree_t *initialRootedTree;
    
    std::unordered_map<pll_unode_t*, Population*> tipsAssign;
    std::unordered_map<std::string, pll_unode_t*> labelsAssign;
    
    std::map<pll_rnode_t*, vector<Population*>> rMRCAPopulation;
    std::map<pll_rnode_t*, vector<Population*>> proposedrMRCAPopulation;
    std::map<pll_rnode_t*, vector<Population*>> currentrMRCAPopulation;
    
    vector<pll_tree_edge_t *> edges;
    
    vector<pair<double, pll_tree_edge_t *> > edgeLengths;
    
    vector<double> sampledTotalEffectPopSize;
    vector<vector<double> > sampledPoportionVector;
    Files files;
    //vector<MCMCmove*>  moves
    //vector<pll_edge_node_t*> edges;
public:
    Chain(int chainNumber,
          int numClones,
          int gammaParam,
          long double  mutationRate,
          long double  seqErrorRate,
          long double  dropoutRate
          );
    
    void MakeCoalescenceEvent( Population *Population, vector<pll_unode_t *> &nodes, int numClones, long int* seed, int noisy,   int &numActiveGametes, int &nextAvailable,
                              int &labelNodes, long double  &currentTime, int &numNodes);
    
    
    void SimulatePopulation( Population *popI, vector<Population*> &populations,
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
                                   vector<pll_unode_t *> &nodes,
                                   int &nextAvailable,
                                   int &numActiveGametes,
                                   int &labelNodes,
                                   long double  &currentTime,
                            int &eventNum);
    
    int setInitialTreeFromNewick(char * NewickString);
    
    int setInitialTreeUnrootedTree(pll_utree_t *unrootedTree);
    
    pll_unode_t* BuildTree(vector<Population* > &populations,
                           Population *CurrentPop,
                           long int *seed,
                           ProgramOptions &programOptions,
                           vector<pll_unode_t *> &nodes,
                           vector<pll_unode_t *> &treeTips,
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
                                       vector<int> &sampleSizes
                                       );
    
    Population* ChooseFatherPopulation( int numClones, Population  *PopChild,  long int *seed, int noisy);
    
    void AssignSequencesToPopulations( ProgramOptions& programOptions,
                                      int numNodes, int noisy,  int TotalNumSequences, int &numActiveGametes, int &nextAvailable,
                                      int &labelNodes, char* SimulatePopulationObservedCellNames[], int doUseObservedCellNames, vector<int> &sampleSizes);
    
    void SetPopulationsBirthRate( long double  lambda);
    void GenerateEffectPopSizesFromPriors2( int noisy,     int doGenerateProportionsVector);
    
    void FillChainPopulationsFromPriors( ProgramOptions &programOptions,  MCMCoptions &mcmcOptions, vector<int> &sampleSizes, long int *seed);
    
    void setChainPopulationSampleSizes(vector<int > &sampleSizes,  ProgramOptions &programOptions);
    void ListClonesAccordingTimeToOrigin(vector<Population *> &populations);
    void GenerateTimesFromPriorsOriginal(int noisy,  long int *seed);
    void InitChainPopulations( int noisy,  int TotalNumSequences  ) ;
    void RelabelNodes(pll_unode_t *p, int &intLabel);
    void InitPopulationSampleSizes(vector<Population*> populations, int TotalSampleSize, int numClones, vector<long double> &proportionsVector, long int *seed);
    long double  SumBranches(pll_unode_t *root, long double  mutationRate);
    char * toNewickString ( pll_unode_t *p, long double  mutationRate,     int doUseObservedCellNames);
    
    long double  LogDensityCoalescentTimesForPopulation();
    long double  LogConditionalLikelihoodTree( ProgramOptions &programOptions  );
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
    static Chain *initializeChain( ProgramOptions &programOptions,  MCMCoptions &mcmcOptions, vector<int> &sampleSizes, long int *seed, char* ObservedCellNames[], pll_msa_t *msa, char* treefileName, string& healthyTipLabel);
    void runChain(   MCMCoptions &opt,  long int *seed,  FilePaths &filePaths, Files &files,  ProgramOptions &programOptions,
                  char* ObservedCellNames[], pll_msa_t * msa, vector<int> &sampleSizes, int currentIteration
                  );
    void newScaledGrowthRateMoveforPopulation( Population *popI, long int *seed,  ProgramOptions &programOptions, char *ObservedCellNames[], pll_msa_t * msa, MCMCoptions & mcmcOptions, vector<int> &sampleSizes);
    void initializeCoalescentEventTimesFormSampleSizes(pll_utree_t *utree, vector<int > &sampleSizes, string &healthyCellLabel );
    void initializeMapPopulationAssignFromTree();
    Population * getPopulationbyIndex(int indexPopulation);
    int getPopulationIndex(char * label);
    std::map<pll_unode_t*, Population*> chooseTimeOfOriginsOnTree( long int *seed);
    std::map<pll_rnode_t*, vector<Population*>>  initTimeOfOriginsOnRootedTree(  vector<pair<double, pll_tree_edge_t *> > edgeLengths, int numberPoints, string &healthyCellLabel);
    void initNodeDataFromTree();
    
    void initPopulationSampleSizesFromRootNodeOnTree(pll_unode_t *p, Population *population );
    bool initPopulationsSampleSizes(std::map<pll_rnode_t*, vector<Population*>>  rmrcaOfPopulation, string &healthyTipLabel );
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
    void initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(pll_rnode_t *p, Population *currentPopulation, std::map<pll_rnode_t*, vector<Population*>> &rmrcaOfPopulation,  string& healthyTipLabel);
    void initPopulationsCoalescentAndMigrationEventsFromRootedTree(std::map<pll_rnode_t*, vector<Population*>> rmrcaOfPopulation, string& healthyTipLabel);
    void initPopulationSampleSizesFromNodeOnRootedTree(pll_rnode_t *p, Population *population, std::map<pll_rnode_t*, vector<Population*> > rmrcaOfPopulation, string &healthyTipLabel);
    void  initNumberTipsSubTree(pll_rnode_t *node);
    void initBranches(string& healthyCellLabel,vector<pair<double, pll_tree_edge_t *> > &edgeLengths, vector<pll_tree_edge_t *> &edges);
    void initTimeOriginSTD();
    long double  SumBranches2(pll_rnode_t *p, long double  mutationRate);
    void filterSortPopulationsCoalescentEvents();
    void samplePopulationDeltaFromPriors(MCMCoptions &mcmcOptions, long int *seed );
    void rescaleRootedTreeBranchLengths(long double mutationRate);
    void newTotalEffectivePopulationSizeMove( ProgramOptions &programOptions, char *ObservedCellNames[], pll_msa_t * msa, MCMCoptions &mcmcOptions, vector<int> &sampleSizes);
    void newProportionsVectorMove(ProgramOptions &programOptions, char *ObservedCellNames[], pll_msa_t * msa, MCMCoptions &mcmcOptions, vector<int> &sampleSizes);
    
    long double   proposalSlidingWindow( long double  oldvalue,  long double  windowSize);
    void proposalProportionsVector(vector<long double  > &newProportionsvector, long double  tuningParameter );
    long double  DirichletDensity(vector<long double> &proportionsVector,  vector<long double> &concentrationVector, int sizeVector);
    void updateEffectPopSizesCurrentProportionsVector();
    int totalSampleSize();
    std::map<pll_rnode_t*, vector<Population*>> chooseAvailableEdgeOnRootedTreeForPopulation(Population *pop, std::map<pll_rnode_t*, vector<Population*>> &mrcaOfPopulation, string &healthyCellLabel);
    long double  sumAvailableBranchLengths(std::map<pll_rnode_t*, vector<Population*>> currentMRCAPopulation);
    std::map<pll_rnode_t*, vector<Population*>> chooseNewTimeofOriginOnEdge(Population *pop);
    void PrepareFiles(const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files, int chainNumber);
    void writeMCMCState( int  currentIteration, const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files,  MCMCoptions &mcmcOptions );
    void writeHeaderOutputChain(  const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files );
    long double  autoCorrelation(int lag, vector<double> values);
    long double  ESS(int lag, vector<double> values);
    void copyProportionsVector(long double  alpha[]);
    void computeAvailableEdges( vector<pair<double, pll_tree_edge_t *> > &availableEdges, std::map<pll_rnode_t *, Population *> &currentMrcaOfPopulation, std::string &healthyCellLabel);
    bool checkMigrationsOrder();
    long double  LogDensityCoalescentTimesForPopulation2();
    void rescaleNodeDataFromRootedTree(long double scale);
    long double   computeAdjacentEdges( vector<pair<double, pll_tree_edge_t *> > &availableEdges, std::map<pll_rnode_t *, vector<Population *>> &currentMrcaOfPopulation, std::string &healthyCellLabel, Population *pop, pll_rnode_t * mrca);
    void safeTreeNodeCurrentTimePUnits();
    void rollbackTreeNodeCurrentTimePUnits();
    void updateNodeScaledTimeForRootedTree(long double  newScale);
    Population* getYoungestPopulationOnEdge(pll_rnode_t* p, std::map<pll_rnode_t*, vector<Population*> >  rmrcaOfPopulation);
    void samplePopulationGrowthRateFromPriors(MCMCoptions &mcmcOptions, long int *seed );
    void updateFatherPopOnSameEdge(std::map<pll_rnode_t*, vector<Population*> > &rmrcaOfPopulation, pll_rnode_t *p, Population * populationOfCurrentNode);
    long double  sumAdjacentEdges(std::map<pll_rnode_t*, vector<Population*>> currentMRCAPopulation, pll_rnode_t*MRCA, Population *pop);
    bool isOldestPopulation(Population *pop, std::map<pll_rnode_t*, vector<Population*>> &currentMrcaOfPopulation);
    void closeFiles(const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files, int chainNumber);
private:
    static double * expand_uniq_rates(int states, const double * uniq_rates, const int * rate_sym);
};

#endif /* Chain_hpp */
