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
#include "population.hpp"
#include "libpll/pll_optimize.h"
#include "libpll/pll_tree.h"
#include "libpll/pllmod_algorithm.h"
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
    pll_unode_t *nodes;
    pll_unode_t  **treeTips;
    
    pll_unode_t *oldroot;
    pll_unode_t *oldnodes;
    pll_unode_t  **oldtreeTips;
    
    int numNodes;
    int oldnumNodes;
    
    double mutationRate;
    double seqErrorRate;
    double dropoutRate;
    double oldmutationRate;
    Population **populations;
    
    double *proportionsVector;
    double *oldproportionsVector;
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
          int totalPopSize,
          double mutationRate,
          double seqErrorRate,
          double dropoutRate
          );
    
    void MakeCoalescenceEvent(Population *population,pll_unode_t **nodes, int numClones, long int* seed, int noisy, int *numActiveGametes,   int* nextAvailable,
                              int*labelNodes, double *currentTime, int *numNodes);
    
    void SimulatePopulation(Population *population, ProgramOptions *programOptions,
                            long int *seed,
                            int *numNodes,
                            int numClones,
                            double      cumNumCA,
                            double meanNumCA,
                            double cumNumMIG,
                            double meanNumMIG,
                            int  *numMIG,
                            int  *numCA,
                            double *numEventsTot,
                            pll_unode_t    **nodes,
                            int *nextAvailable,
                            int*  numActiveGametes,
                            int* labelNodes,
                            double *currentTime,
                            int* eventNum
                            );
    
    
    int setIntitialTree(char * NewickString);

    void BuildTree(Population *olderPopulation, long int *seed,
                          ProgramOptions *programOptions,
                          pll_unode_t    **nodes,
                          pll_unode_t   **treeTips,
                          pll_unode_t    **treeRootInit,
                          int *nextAvailable,
                          int *newInd,
                          double *currentTime,
                          int *labelNodes
                          );

    void MakeCoalescenceTree (long int *seed,
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
                                     char* ObservedCellNames[],
                                     int *sampleSizes
                              ) ;
    
    Population* ChooseFatherPopulation( int numClones, Population  *PopChild,  long int *seed, int noisy);
    
    void AssignSequencesToPopulations( ProgramOptions* programOptions,
                                             int numNodes, int noisy,  int TotalNumSequences, int *numActiveGametes, int* nextAvailable,
                                             int *labelNodes, char* ObservedCellNames[], int doUseObservedCellNames, int *sampleSizes);
    
    void SetPopulationsBirthRate( double lambda);
    void GenerateEffectPopSizesFromPriors2( int noisy,   long int *seed,  int doGenerateProportionsVector);
    void FillChainPopulationsFromPriors( ProgramOptions *programOptions,  MCMCoptions *mcmcOptions, int *sampleSizes, long int *seed);
    void setChainPopulationSampleSizes(int *sampleSizes,  ProgramOptions *programOptions);
    void ListClonesAccordingTimeToOrigin2();
    void GenerateTimesFromPriorsOriginal(int noisy,  long int *seed);
    void InitChainPopulations( int noisy,  int TotalNumSequences  ) ;
};

#endif /* Chain_hpp */
