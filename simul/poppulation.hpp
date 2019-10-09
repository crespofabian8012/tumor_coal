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

#include <libpll/pll_tree.h>

#include "tree_node.hpp"

using namespace std;

class Population
{
public:
    int index; // label of the population {1, 2, 3, ...}
    int order; // population order by time of origin
    double timeOriginSTD; // standardized time of origin
    double timeOriginInput; // time of origin
    double delta; // growth rate * effectPopSize
    double effectPopSize;
    double birthRate, deathRate, growthRate;
    int sampleSize; // number of cells in the population
    int popSize; // total number of cells (unknown at inference)
    int numActiveGametes; // number of (currently) living cells
    int numGametes; // ???
    int numCompletedCoalescences;
    int nextAvailableIdInmigrant; // TODO: check later
    int numIncomingMigrations, numPossibleMigrations;
    bool doEstimateTimeOrigin;
    bool isAlive, CellAssignationCompleted;
//    double timeMigrationSTDCurrentPop; // TODO: check
//    vector<double> migrationTimes;
//    vector<Population *> immigrantsPopOrderedModelTime;
    vector<pair<double, Population *> > immigrantsPopOrderedByModelTime; // migrationTime, Population
    vector<int> idsActiveGametes;
    vector<int> idsGametes;
    int indexFirstObservedCellName;
    int nodeIdAncestorMRCA;
    //pll_unode_t *MRCA;
    TreeNode *MRCA;
    Population *FatherPop;
    vector<double> CoalescentEventTimes;

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
};
#endif /* Population_hpp */
