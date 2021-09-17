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
 * tumor population class
 */

#ifndef population_hpp
#define population_hpp


#include "tree_node.hpp"
//#include <vector>

#include <boost/random.hpp>

extern "C"
{
#include <gsl/gsl_rng.h>
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
#include "libpll/pll_optimize.h"
#include "libpll/pllmod_algorithm.h"
#include "libpll/pllmod_common.h"
}


//#include "tree_node.hpp"

class MCMCoptions;
class ProgramOptions;
class MCMCParameterWithKernel;
template<typename T>
class MCMCParameter;

class Population
{
public:
    int index; // label of the population {1, 2, 3, ...}
    int order; // population order by time of origin
    int oldOrder; // population order by time of origin
    long double timeOriginSTD; // standardized time of origin
    long double oldTimeOriginSTD; // standardized time of origin
    long double timeOriginInput; // time of origin
    long double oldTimeOriginInput; // time of origin
    long double lowerBoundTimeOriginInput;
    long double scaledtimeOriginInput; // time of origin
    long double oldScaledTimeOriginInput; // time of origin
    long double delta; // growth rate * effectPopSize = r * x
    long double olddelta;
    long double deltaT;
    long double olddeltaT;
    long double x;
    long double oldx;
    long double theta;
    long double oldTheta;
    long double currentModelTime;

    long double birthRate, deathRate, growthRate;
    long double oldDeathRate, oldGrowthRate;
    int sampleSize; // number of cells in the population
    int oldSampleSize; // number of cells in the population

    int numActiveGametes; // number of (currently) living cells
    int numGametes; // ???
    int numCompletedCoalescences;
    int nextAvailableIdInmigrant; // TODO: check later
    int numIncomingMigrations, numPossibleMigrations;
    int indexNextInmigrant;
    bool doEstimateTimeOrigin;
    bool isAlive, CellAssignationCompleted;
    //    double timeMigrationSTDCurrentPop;
    std::vector<std::pair<long double, Population *>> immigrantsPopOrderedByModelTime; // migrationTime, Population
    
    std::vector<std::pair<long double, Population *>> oldimmigrantsPopOrderedByModelTime; // migrationTime, Population
    
    std::vector<long double>  posteriorDeltaT;
    std::vector<long double>  posteriorTimeOriginSTD;
    std::vector<long double>  posteriorProportion;
    
    MCMCParameterWithKernel *DeltaT;
    MCMCParameterWithKernel *TimeOriginSTD;
    MCMCParameterWithKernel *Theta;
    MCMCParameterWithKernel *TimeOriginInput;
    MCMCParameter<long double> *X;
    MCMCParameter<long double> * sampleSizePar;
    
    std::vector<int> idsActiveGametes;
    std::vector<int> idsGametes;
    int indexFirstObservedCellName;
    int nodeIdAncestorMRCA;
    //pll_unode_t *MRCA;
    TreeNode *MRCA;
    std::vector<pll_unode_t *> tips;
    std::vector<pll_rnode_t *> rtips;
    std::vector<pll_rnode_t *> oldrtips;
    
    pll_unode_t * nodeletMRCA;
    pll_rnode_t * rMRCA;
    pll_rnode_t * oldrMRCA;
    
    Population *FatherPop;
    Population *oldFatherPop;
    std::vector<long double> CoalescentEventTimes;
    std::vector<long double> oldCoalescentEventTimes;
    
public:
    Population(int ind, int ord, long double timeOriginInput,
               int sampleSize,int totalSampleSize, int popSize,int totalPopSize, long double birthRate,
               long double deathRate, bool estimateTOR);
    Population(int ind, int ord, long double timeOriginInput,
               int sampleSize, int popSize, long double birthRate,
               long double deathRate, bool estimateTOR, MCMCoptions &mcmcOptions);
    Population(int ind, int ord, int sampleSize, long double delta, long double theta,
               ProgramOptions &programOptions);
    
    Population(const Population &original);
    
    long double ProbabilityComeFromPopulation(Population *PopJ, std::vector<Population*> &populations, int numClones, long double K);
    static long double FmodelTstandard (long double t, long  double TOrigin, long double delta, long double K);
    static long double GstandardTmodel (long double V, long double TOrigin, long double delta, long double K);
    static long double LogProbNoCoalescentEventBetweenTimes(long double from, long double to, int numberActiveInd, long double TOrigin, long double delta,  long double K);
    static long double CalculateH (long double t, long double TOrigin, long double delta, long double K);
    static bool comparePopulationsPairByTimeOrigin(const std::pair<long double, Population *> s1, const std::pair<long double, Population *> s2);
    static long double DensityTimeSTD(long double u, long double deltaPar, long double from);
    static int compare (const void * a, const void * b);
    void InitListPossibleMigrations(int order);
    int resetMigrationsList();
    static void  UpdateListMigrants( int numClones, Population *PopChild, Population *PopFather  );
    
    
    void ChooseRandomIndividual(int *firstInd,   int numClones,   int *secondInd, long *seed, int choosePairIndividuals);
    void ChooseRandomIndividual(int *firstInd,   int numClones,   int *secondInd, const gsl_rng *randomGenerator, int choosePairIndividuals);
    void InitCoalescentEvents(int numClones);
    void resetGametesCounters();
    
    void resetActiveGametes();
    static int bbinClones (long double dat, long double *v, int n);
    long double DensityTime( long double u);
    long double LogProbNoCoalescentEventBetweenTimes(long double from, long double to, int numberActiveInd,  long double K);
    void filterAndSortCoalescentEvents();
    void multiplyCoalescentsEventByFactor(long double factor);
    void multiplyMigrationsTimesByFactor(long double factor);
    long double LogDensityTime(long double u);
    static long double LogCalculateH (long double t, long double TOrigin, long double delta, long double K);
   static  long double LogLambda (long double t, long double TOrigin, long double delta, long double K);
    long double LogDensityTime2(long double u);
    long double DensityTimeSTD(long double u, long double from);
    long double LogDensityTimeSTDFrom(long double u, long double from);
    void setLowerBoundTimeOriginInput(long double from);
    void restoreOldCoalescentTimes();
    void restoreOldImmigrationTimes();
    void InitIdsActiveGametes();
    void InitIdsGametes(int numClones);
    void InitRTips();
    void savePosteriorValues();
    void resetPosteriorValues();
    long double getCurrentModelTime() const{return currentModelTime;} ;
    void setCurrentModelTime(long double newModelTime );
  
    
    void setTimeOriginInputTree(long double timeOriginInputTree);
    void setTimeOriginSTD(long double timeOriginSTD);
    void setScaledTimeOriginInputTree(long double scaledTimeOriginInputTree);
    void setTheta(long double theta);
    void setProportion(long double x);
    void setDelta(long double parDelta);
    void setPopulationToriginConditionalDelta( const gsl_rng *rngGsl );
    long double proposeTimeNextCoalEvent(gsl_rng* rngGsl, int numActiveLineages,  double K);
    long double proposeWaitingKingmanTimeNextCoalEvent(gsl_rng* rngGsl );
    long double proposeWaitingKingmanTimeNextCoalEvent(gsl_rng* rngGsl, size_t numActiveGametesPar   );
    long double logConditionalLikelihoodNextCoalescentTime(long double timeNextEvent,  double K);
    long double logConditionalDensityWaitingTimeScaledByTheta(long double waitingTimeScaledByTheta,  double K);
  
    long double logConditionalDensityWaitingTimeScaledByTheta(long double waitingTimeScaledbyTheta, size_t numActiveGametesP, double currentModelTimeP,  double K);
    void sampleEventTimesScaledByProportion(gsl_rng *random, double K);
    double  nextCoalEventTime(int idxNextCoal, int indexNextMigration,  double currentTime, bool& isThereInmigration, Population* inmigrantPop );
    void printCoalEventTimes(std::ostream &stream);
    void printInmigrationEventTimes(std::ostream &stream);
    
private:
    static bool isNotPositive(long double d);
};
class PopulationSet{
public:
    int numClones;

    std::vector<Population *> populations;
    std::vector<long double  > proportionsVector;
    std::vector<long double  > oldproportionsVector;
public:
    PopulationSet(int numClones);
    
    PopulationSet(const PopulationSet &original);
    //PopulationSet &operator=(const PopulationSet &original);
    
    void initPopulation();
    
    std::vector<long double> samplePopulationGrowthRateFromPriors(MCMCoptions &mcmcOptions, const gsl_rng * randomGenerator );
    void  setPopulationsBirthRate( long double  lambda);
    void  setPopulationsToriginConditionalDelta( const gsl_rng *rngGsl );
    void initPopulationSampleSizes(std::vector<int> &sampleSizes);
    Population * getPopulationbyIndex(int indexPopulation);
    void initPopulationGametes();
    void initPopulationRTips();
    void initProportionsVectorFromSampleSizes(std::vector<int> &sampleSizes);
    void initPopulationsThetaDelta(long double theta);
    void initListPossibleMigrations();
    void initPopulationsCoalescentEvents();
    void initProportionsVector();
    
    void initPopulationDeltas(std::vector<double> &deltas);
    void initPopulationTOriginSTD(std::vector<double> &TOriginSTD);
    void initTheta( long double &theta);
    std::vector <long double> getDeltaTs();
    std::vector <long double> getDeltas();
    std::vector <long double> getTs();
    std::vector <long double> getSampleSizes();
    void initDeltaThetaFromPriors( const gsl_rng *rngGsl, long double &theta);
    Population* ChooseFatherPopulation( Population  *PopChild,  const gsl_rng *randomGenerator, int noisy, long double K);
    void AssignSequencesToPopulations(std::vector<pll_rnode_t*> rnodes,
                                                     ProgramOptions &programOptions,
                                                     int noisy,  int TotalTumorSequences,
                                                     int &numActiveGametes, int &nextAvailable,
                                                     int &labelNodes, char* ObservedCellNames[], int doUseObservedCellNames, std::vector<int> &sampleSizes);
    std::vector<Population *>& getPopulations();
    void sortPopulationsByTorigin();
    void sampleEventTimesScaledByProportion(gsl_rng * random, double K, int noisy);
    void resetNumActiveGametesCounter();
    double proposeNextCoalEventTime( gsl_rng *random, int& idxLeftNodePop, int& idxRightNodePop, int& idxIncomingNode, double &logLik, double K);
    void  acceptNextCoalEventTime(gsl_rng *random, double receiverModelTime, int idxReceiverPop, int idxIncomingPop, double K);
    Population* getCurrentPopulation();

};
class StructuredCoalescentTree{
public:
    int numClones;
    long double theta;
    long double oldTheta;
    int currentNumberEdgesRootedTree;
    PopulationSet * populationSet;
    std::vector<int> sampleSizes;
    pll_rnode_t *root;
    pll_rtree_t *rtree;
    pll_rnode_t *healthyTip;
    std::vector<pll_rnode_t *> rnodes;
    std::vector<pll_rnode_t *> rtreeTips;
    std::vector<pll_tree_edge_t *> edges;
    std::vector<std::pair<double, pll_tree_edge_t *> > edgeLengths;
    long double logLikelihoodTree;
    long double logLikelihoodSequences;
    long double seqErrorRate;
    long double dropoutRate;
public:
    StructuredCoalescentTree(int numClones, std::vector<int> &sampleSizes, long double theta, MCMCoptions &mcmcOptions,ProgramOptions &programOptions,const gsl_rng *rngGsl,
    boost::mt19937 *rngBoost, std::vector<std::vector<int> > &ObservedData,char* ObservedCellNames[], pll_msa_t *msa, std::string& healthyTipLabel,  long double seqErrorRate,
       long double dropoutRate);
    StructuredCoalescentTree(PopulationSet *populationSet, std::vector<int> &sampleSizes, long double theta, MCMCoptions &mcmcOptions,ProgramOptions &programOptions, const gsl_rng *rngGsl,
       boost::mt19937 *rngBoost, std::vector<std::vector<int> > &ObservedData,char* ObservedCellNames[], pll_msa_t *msa, std::string& healthyTipLabel,  long double seqErrorRate,
       long double dropoutRate );
    PopulationSet& getPopulationSet();
    pll_rtree_t* getTree();
    pll_rnode_t* getRoot();
    void addEdgeFromNode(pll_rnode_t *node );
    void RelabelNodes2(pll_rnode_t *p, int &intLabel);
    void initEdgesRootedTree(int numberTips);
    std::vector <long double> getDeltaTs();
    std::vector <long double> getDeltas();
    std::vector <long double> getTs();
    std::vector <long double> getSampleSizes();
    void SimulatePopulation(Population &popI,
                                     ProgramOptions &programOptions,
                                     const gsl_rng *rngGsl,
                                     boost::mt19937 *rngBoost,
                                     int &numNodes,
                                     int numClones,
                                     int &nextAvailable,
                                     int &numActiveGametes,
                                     int &labelNodes,
                                     long double  &currentTime,
                                     int &eventNum,
                                     long double K);
    pll_rnode_t * MakeCoalescenceTree (const gsl_rng *rngGsl,
                                       boost::mt19937 *rngBoost,
                                       pll_msa_t *msa,
                                       int &numNodes,
                                       ProgramOptions &programOptions,
                                        std::vector<std::vector<int> > &ObservedData,
                                        char* ObservedCellNames[],
                                        std::vector<int> &sampleSizes) ;

    void MakeCoalescenceEvent( Population &population, const gsl_rng *randomGenerator, int noisy,   int &numActiveGametes, int &nextAvailable,
                              int &labelNodes, long double  &currentTime, int &numNodes);
    pll_rnode_t* BuildTree(Population *CurrentPop,
                           const gsl_rng *randomGenerator,
                           ProgramOptions &programOptions,
                           pll_rnode_t *tumour_mrca,
                           int &nextAvailable,
                           int &newInd,
                           long double   &currentTime,
                           int &labelNodes);
    
    long double SumLogProbFatherPopulations(long double K);
    static long double SumLogProbFatherPopulations(std::vector<Population *> populations, int numClones, long double K);
    long double SumLogDensitiesTimeOriginSTDPopulations();
    static long double SumLogDensitiesTimeOriginSTDPopulations(std::vector<Population *> populations, int numClones, long double K);
    void initLogLikelihoodSequences(const ProgramOptions &programOptions, const pll_msa_t *msa);
    long double  LogDensityCoalescentTimesForPopulation(long double K);
    static long double  LogDensityCoalescentTimesForPopulation(std::vector<Population *> populations, int numClones, long double K);
 
    long double getLogLikelihoodTree() const{return logLikelihoodTree;}
    long double getLogLikelihoodSequences() const{return logLikelihoodSequences;}
    
};
#endif /* Population_hpp */
