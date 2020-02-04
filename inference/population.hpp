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
    double timeOriginSTD; // standardized time of origin
    double oldTimeOriginSTD; // standardized time of origin
    double timeOriginInput; // time of origin
    double oldTimeOriginInput; // time of origin
    double scaledtimeOriginInput; // time of origin
    double oldScaledTimeOriginInput; // time of origin
    double delta; // growth rate * effectPopSize = r * x
    double olddelta;
    double r;
    double oldr;
    double x;
    double oldx;
    double theta;
    double oldTheta;
    double effectPopSize;
    double oldeffectPopSize;
    double birthRate, deathRate, growthRate;
    double oldDeathRate, oldGrowthRate;
    int sampleSize; // number of cells in the population
    int oldSampleSize; // number of cells in the population
    int popSize; // total number of cells (unknown at inference)
    int oldPopSize;
    int numActiveGametes; // number of (currently) living cells
    int numGametes; // ???
    int numCompletedCoalescences;
    int nextAvailableIdInmigrant; // TODO: check later
    int numIncomingMigrations, numPossibleMigrations;
    bool doEstimateTimeOrigin;
    bool isAlive, CellAssignationCompleted;
//    double timeMigrationSTDCurrentPop; // TODO: check
    vector<pair<double, Population *> > immigrantsPopOrderedByModelTime; // migrationTime, Population
    
    vector<pair<double, Population *> > oldimmigrantsPopOrderedByModelTime; // migrationTime, Population
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
    vector<double> CoalescentEventTimes;
    vector<double> oldCoalescentEventTimes;
    
public:
    Population(int ind, int ord, double timeOriginInput,
               int sampleSize, int popSize, double birthRate,
               double deathRate, bool estimateTOR);

    double ProbabilityComeFromPopulation(Population *PopJ, vector<Population*> &populations, int numClones);
    static double FmodelTstandard (double t, double TOrigin, double delta);
    static double GstandardTmodel (double V, double TOrigin, double delta);
    static double CalculateH (double t, double TOrigin, double delta);
    static bool comparePopulationsPairByTimeOrigin(const pair<double, Population *> s1, const pair<double, Population *> s2);
    static int compare (const void * a, const void * b);
    void InitListPossibleMigrations(int order);
    int resetMigrationsList();
    static void  UpdateListMigrants( int numClones, Population *PopChild, Population *PopFather  );
    
  
    void ChooseRandomIndividual(int *firstInd,   int numClones,   int *secondInd, long *seed, int choosePairIndividuals);
    void InitCoalescentEvents(int numClones);
    void resetActiveGametes();
    static int bbinClones (double dat, double *v, int n);
    double DensityTime( double u);
    double LogProbNoCoalescentEventBetweenTimes(double from, double to, int numberActiveInd);
    void filterAndSortCoalescentEvents();
    void multiplyCoalescentsEventByFactor(double factor);
    void multiplyMigrationsTimesByFactor(double factor);
    double LogDensityTime(double u);
    static double LogCalculateH (double t, double TOrigin, double delta);
private:
    static bool isNotPositive(double d);
};
#endif /* Population_hpp */
