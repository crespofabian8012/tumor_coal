//
//  population.h
//  TumorCellCoal
//
//  Created by Fausto Fabian Crespo Fernandez on 7/8/19.
//  Copyright © 2019 Fausto Fabian Crespo Fernandez. All rights reserved.
//

#ifndef population_h
#define population_h

#include <vector>

#include "data_utils.hpp"

static void InitListClones(Population **populations, int numClones, int noisy, int *CloneNameBegin, int *CloneSampleSizeBegin, double *CloneBirthRateBegin,  double *CloneDeathRateBegin,  int *ClonePopSizeBegin, int TotalNumSequences  );
static void InitListClonesTimes(Population  **populations, int numClones, int *doEstimateTimesOriginClones,   double *CloneTimeOriginInput);

static void CheckNumberOfClonesOrigin(int numClones, double  *CloneTimeOriginInputSTDPOP, double  *CloneTimeOriginInput );
void ListClonesAccordingTimeToOrigin(int numClones, int noisy, int      **CloneNameBeginOrderByTimesOnlyControl, int **CloneNameBeginOrderByTimes, double *CloneTimeOriginInput);

static void InitPopulationSampleSizes(Population **populations, int TotalSampleSize, int numClones, double *proportionsVector, long int *seed);

static void ValidateParameters(ProgramOptions *programOptions,
                               // int noisy, int numDataSets, int Nscaling, int numClones,
                               int *CloneNameBegin , int *CloneSampleSizeBegin, int *ClonePopSizeBegin);

static void SimulatePopulation( Population *popI, Population** populations,
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
                               TreeNode    **nodes,
                               int *nextAvailable,
                               int*  numActiveGametes,
                               int* labelNodes,
                               double *currentTime,
                               int* eventNum
                               );
static int comparePopulationsByTimeOrigin(const void *s1, const void *s2);
static int comparePopulationsBySTDPopTime(const void *s1, const void *s2); static void ListClonesAccordingTimeToOrigin2(Population **populations, int numClones);
static void  InitListPossibleMigrations(Population **populations,  int numClones);
static void  UpdateListMigrants(Population **populations, int numClones, Population *PopChild, Population *PopFather  );
static Population* ChooseFatherPopulation(Population **populations, int numClones, Population*  PopChild, long int *seed, int noisy);
static  double ProbabilityCloneiFromClonej2 (Population *PopI, Population* PopJ, Population **populations, int numClones);

static  void AssignCurrentSequencesToPopulation(Population **populations, TreeNode **nodes,  ProgramOptions* programOptions, int numClones, int numNodes, int noisy,  int TotalNumSequences, int *numActiveGametes, int* nextAvailable,
                                                int *labelNodes, char* ObservedCellNames[], int doUseObservedCellNames, int *sampleSizes);
void resetMigrationsList(Population **populations, int numClones);
void SetPopulationsBirthRate(Population** populations, double lambda, int numClones);
void ChooseRandomIndividual(int *firstInd,   int numClones, Population *popI,  int *secondInd, long *seed, int choosePairIndividuals);
void InitPopulationsCoalescentEvents( int numClones,  Population **populations) ;

void AssignObservedCellNamestoTips(TreeNode **nodes, Population ** populations, int *sampleSizes,char* ObservedCellNames[],   ProgramOptions *programOptions);

void AssignObservedCellNamestoTips2(TreeNode **nodes, Population ** populations, int *sampleSizes, char* ObservedCellNames[],   ProgramOptions *programOptions);

#endif /* population_h */
