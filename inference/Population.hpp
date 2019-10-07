//
//  Population.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 4/10/19.
//

#ifndef Population_hpp
#define Population_hpp

#include <stdio.h>

#include "data_utils.hpp"
#include "Population.hpp"

using namespace std;

class population
{
public:
    int                 index, order;
    double          timeOriginSTD, delta, effectPopSize, oldeffectPopSize, olddelta;
    double              birthRate, deathRate, growthRate, oldbirthRate, olddeathRate, oldgrowthRate;
    double              timeOriginInput, oldtimeOriginInput;
    int           sampleSize, popSize, numActiveGametes, oldsampleSize, oldpopSize, numGametes;
    int                 numCompletedCoalescences;
    int           nextAvailableIdInmigrant;
    int                 numIncomingMigrations, numPossibleMigrations;
    int                 doEstimateTimeOrigin;
    int                 isAlive, CellAssignationCompleted;
    double              timeMigrationSTDCurrentPop;
    double             *migrationTimes;
    int                *idsActiveGametes;
    int                *idsGametes;
    int                 indexFirstObservedCellName;
    int                 nodeIdAncesterMRCA;
    pll_unode_t         *MRCA;
    population                *FatherPop;
    population                *oldFatherPop;
    double            *CoalescentEventTimes;
    double            *oldCoalescentEventTimes;
    population** immigrantsPopOrderedModelTime;
    

    population(int ind, int ord,   double  timeOriginInputp,
               int           sampleSizep, int popSizep, double birthRatep,
               double deathRatep);
public:
     double    ProbabilityComeFromPopulation ( population* PopJ, population **populations, int numClones);
    static double FmodelTstandard (double t, double TOrigin, double delta);
    static double CalculateH (double t, double TOrigin, double delta);
    static int comparePopulationsByTimeOrigin(const void *s1, const void *s2);
    static int compare (const void * a, const void * b);
    int   InitListPossibleMigrations( );
    int  resetMigrationsList();
    void  UpdateListMigrants( int numClones, population *PopChild, population *PopFather  );
    
    void SimulatePopulation(  Population** populations,
                                        ProgramOptions *programOptions,
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
    void ChooseRandomIndividual(int *firstInd,   int numClones,   int *secondInd, long *seed, int choosePairIndividuals);
    void MakeCoalescenceEvent(pll_unode_t **nodes, int numClones, long int* seed, int noisy, int *numActiveGametes,   int* nextAvailable,
                              int*labelNodes, double *currentTime, int *numNodes);
};
#endif /* Population_hpp */
