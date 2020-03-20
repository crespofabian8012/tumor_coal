//
//  Population.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 4/10/19.
//

#ifndef Population_hpp
#define Population_hpp

#include <stdio.h>
#include <vector>


extern "C"
{
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
#include "libpll/pll_optimize.h"
#include "libpll/pllmod_algorithm.h"
#include "libpll/pllmod_common.h"
}

#include "tree_node.hpp"

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
    long double scaledtimeOriginInput; // time of origin
    long double oldScaledTimeOriginInput; // time of origin
    long double delta; // growth rate * effectPopSize = r * x
    long double olddelta;
    long double r;
    long double oldr;
    long double x;
    long double oldx;
    long double theta;
    long double oldTheta;
    long double effectPopSize;
    long double oldeffectPopSize;
    long double birthRate, deathRate, growthRate;
    long double oldDeathRate, oldGrowthRate;
    int sampleSize; // number of cells in the population
    int oldSampleSize; // number of cells in the population
    unsigned long popSize; // total number of cells (unknown at inference)
    unsigned long oldPopSize;
    int numActiveGametes; // number of (currently) living cells
    int numGametes; // ???
    int numCompletedCoalescences;
    int nextAvailableIdInmigrant; // TODO: check later
    int numIncomingMigrations, numPossibleMigrations;
    bool doEstimateTimeOrigin;
    bool isAlive, CellAssignationCompleted;
    //    double timeMigrationSTDCurrentPop; // TODO: check
    vector<pair<long double, Population *> > immigrantsPopOrderedByModelTime; // migrationTime, Population
    
    vector<pair<long double, Population *> > oldimmigrantsPopOrderedByModelTime; // migrationTime, Population
    vector<int> idsActiveGametes;
    vector<int> idsGametes;
    int indexFirstObservedCellName;
    int nodeIdAncestorMRCA;
    //pll_unode_t *MRCA;
    TreeNode *MRCA;
    vector<pll_unode_t *> tips;
    vector<pll_rnode_t *> rtips;
    vector<pll_rnode_t *> oldrtips;
    
    pll_unode_t * nodeletMRCA;
    pll_rnode_t * rMRCA;
    pll_rnode_t * oldrMRCA;
    
    Population *FatherPop;
    Population *oldFatherPop;
    vector<long double> CoalescentEventTimes;
    vector<long double> oldCoalescentEventTimes;
    
public:
    Population(int ind, int ord, long double timeOriginInput,
               int sampleSize, int popSize, long double birthRate,
               long double deathRate, bool estimateTOR);
    
    long double ProbabilityComeFromPopulation(Population *PopJ, vector<Population*> &populations, int numClones);
    static long double FmodelTstandard (long double t, long  double TOrigin, long double delta);
    static long double GstandardTmodel (long double V, long double TOrigin, long double delta);
    static long double CalculateH (long double t, long double TOrigin, long double delta);
    static bool comparePopulationsPairByTimeOrigin(const pair<long double, Population *> s1, const pair<long double, Population *> s2);
    static int compare (const void * a, const void * b);
    void InitListPossibleMigrations(int order);
    int resetMigrationsList();
    static void  UpdateListMigrants( int numClones, Population *PopChild, Population *PopFather  );
    
    
    void ChooseRandomIndividual(int *firstInd,   int numClones,   int *secondInd, long *seed, int choosePairIndividuals);
    void InitCoalescentEvents(int numClones);
    void resetActiveGametes();
    static int bbinClones (long double dat, long double *v, int n);
    long double DensityTime( long double u);
    long double LogProbNoCoalescentEventBetweenTimes(long double from, long double to, int numberActiveInd);
    void filterAndSortCoalescentEvents();
    void multiplyCoalescentsEventByFactor(long double factor);
    void multiplyMigrationsTimesByFactor(long double factor);
    long double LogDensityTime(long double u);
    static long double LogCalculateH (long double t, long double TOrigin, long double delta);
    long double LogDensityTime2(long double u);
private:
    static bool isNotPositive(long double d);
};
#endif /* Population_hpp */
