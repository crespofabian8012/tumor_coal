//
//  Chain.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 7/10/19.
//

#ifndef Chain_hpp
#define Chain_hpp

#include <stdio.h>
#include "data_utils.hpp"
#include "chain.hpp"
//#include "population.hpp"


#include "utils.hpp"
#include "definitions.hpp"
class Chain{
public:
    int chainNumber;
    int numClones;
    int oldnumClones;
    int gammaParam;
    int currentNumberIerations;

    pll_unode_t *root;
    vector<pll_unode_t*> nodes;
    vector<pll_unode_t*> treeTips;
    
    pll_unode_t *oldroot;
    vector<pll_unode_t *>oldnodes;
    vector<pll_unode_t*>  oldtreeTips;
    
    int numNodes;
    int oldnumNodes;
    
    double mutationRate;
    double seqErrorRate;
    double dropoutRate;
    double oldmutationRate;
    vector<Population*> populations;
    
    vector<double > proportionsVector;
    vector<double > oldproportionsVector;
    int totalPopSize;
    int oldtotalPopSize;
    double lambda ;
    double currentlogConditionalLikelihoodTree;
    double currentlogConditionalLikelihoodSequences;
    
    pll_utree_t *initialTree;
    
public:
    Chain( int chainNumber,
          int numClones,
          int gammaParam,
          double mutationRate,
          double seqErrorRate,
          double dropoutRate
          );
    
    void MakeCoalescenceEvent( Population *Population, vector<pll_unode_t *> &nodes, int numClones, long int* seed, int noisy,   int &numActiveGametes, int &nextAvailable,
                              int &labelNodes, double &currentTime, int &numNodes);
    

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
                            vector<pll_unode_t *> &nodes,
                            int &nextAvailable,
                            int &numActiveGametes,
                            int &labelNodes,
                            double &currentTime,
                            int &eventNum);
    
    
    
    
    int setIntitialTree(char * NewickString);
    
    pll_unode_t* BuildTree(vector<Population* > &populations,
                                  Population *CurrentPop,
                                  long int *seed,
                                  ProgramOptions &programOptions,
                                  vector<pll_unode_t *> &nodes,
                                  vector<pll_unode_t *> &treeTips,
                                  pll_unode_t *tumour_mrca,
                                  int &nextAvailable,
                                  int &newInd,
                                  double &currentTime,
                                  int &labelNodes);
    
  


    pll_unode_t * MakeCoalescenceTree (long int *seed,
                                              int &numNodes,
                                              int numClones,
                                              ProgramOptions &programOptions,
                                              double      cumNumCA,
                                              double meanNumCA,
                                              double cumNumMIG,
                                              double meanNumMIG,
                                              int  &numMIG,
                                              int  &numCA,
                                              double &numEventsTot,
                                              char* ObservedCellNames[],
                                              vector<int> &sampleSizes
                                              );
    
    Population* ChooseFatherPopulation( int numClones, Population  *PopChild,  long int *seed, int noisy);
    
    void AssignSequencesToPopulations( ProgramOptions& programOptions,
                                             int numNodes, int noisy,  int TotalNumSequences, int &numActiveGametes, int &nextAvailable,
                                             int &labelNodes, char* SimulatePopulationObservedCellNames[], int doUseObservedCellNames, vector<int> &sampleSizes);
    
    void SetPopulationsBirthRate( double lambda);
    void GenerateEffectPopSizesFromPriors2( int noisy,   long int *seed,  int doGenerateProportionsVector);
  
    void FillChainPopulationsFromPriors( ProgramOptions &programOptions,  MCMCoptions &mcmcOptions, vector<int> &sampleSizes, long int *seed);
 
    void setChainPopulationSampleSizes(vector<int > &sampleSizes,  ProgramOptions &programOptions);
    void ListClonesAccordingTimeToOrigin(vector<Population *> &populations);
    void GenerateTimesFromPriorsOriginal(int noisy,  long int *seed);
    void InitChainPopulations( int noisy,  int TotalNumSequences  ) ;
    void RelabelNodes(pll_unode_t *p, int &intLabel);
    void InitPopulationSampleSizes(vector<Population*> populations, int TotalSampleSize, int numClones, vector<double> &proportionsVector, long int *seed);
    double SumBranches(pll_unode_t *root, double mutationRate);
    char * toNewickString ( pll_unode_t *p, double mutationRate,     int doUseObservedCellNames);
    
    double LogDensityCoalescentTimesForPopulation(pll_unode_t  *tree);
    double LogConditionalLikelihoodTree(pll_unode_t  *tree, ProgramOptions &programOptions  );
    double  LogConditionalLikelihoodSequences(pll_msa_t * msa, char* NewickString, ProgramOptions &programOptions, double seqError,double dropoutError);
    int set_tipclv_custom_error_model(pll_partition_t * partition,
                           unsigned int tip_index,
                           const pll_state_t * map,
                           const char * sequence,
                           double _seq_error_rate,
                    double _dropout_rate);
    void set_partition_tips_costum( pll_partition_t * partition, pll_msa_t * msa, ProgramOptions &programOptions, double seqError, double dropoutError);
    void dealloc_data_costum(pll_unode_t * node, void (*cb_destroy)(void *));
    void  destroyTree(pll_utree_t * tree, void (*cb_destroy)(void *));
   static void InitializeChains(vector<Chain*> &chains,   ProgramOptions &programOptions,  MCMCoptions &mcmcOptions, vector<int> &sampleSizes, long int *seed, char* ObservedCellNames[], pll_msa_t *msa, pll_utree_t * initialTree);
};

#endif /* Chain_hpp */
