//
//  Population.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 4/10/19.
//

#ifndef Population_hpp
#define Population_hpp

#include <stdio.h>
#include "definitions.h"

using namespace std;
class Population
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
    TreeNode            *MRCA;
    struct Population                *FatherPop;
    struct Population                *oldFatherPop;
    double            *CoalescentEventTimes;
    double            *oldCoalescentEventTimes;
    struct Population  **immigrantsPopOrderedModelTime;
    
    //Default Constructor
    Population();
    //Parametrized Constructor
    Population(int ind, int ord,   double  timeOriginInputp,
               int           sampleSizep, int popSizep, double birthRatep
               double deathRatep);
    ~Population()
    {
        
    }
    
private:
   
}
#endif /* Population_hpp */
