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
 * data types
 */

#ifndef data_types_hpp
#define data_types_hpp

#include "definitions.hpp"


#include <iostream>
#include <unordered_map>
#include <unordered_set>

//#include "definitions.hpp"


// SiteStr: information about a site
typedef struct
{
    int        isSNV, isSNP, isVariant;
    int        numMutations, numMutationsMaternal, numMutationsPaternal;
    int        numDeletions, numDeletionsMaternal, numDeletionsPaternal;
    int        numCNLOH, numCNLOHmaternal, numCNLOHpaternal;
    int        hasADO;
    int        hasGenotypeError;
    int        countA, countC, countG, countT, countACGT, countCellswithData, countDropped;
    int        referenceAllele;
    int        *alternateAlleles; // ??? Ask David?
    int        numAltAlleles;
    double  rateMultiplier;
    int    numberDiffReference;
}
SiteStr;

typedef struct
{
    int     tempLength;
    int     numAvailablePositions;
    int     *position;
}
TriNucStr;

typedef struct csite
{
    int        numReads;
    int        *readCount;
    int        trueMaternalAllele, truePaternalAllele;
    int        maternalAllele, paternalAllele;
    int        thereIsMaternalAllele, thereIsPaternalAllele;
    double  maternalSiteAmplificationError, paternalSiteAmplificationError;
    double    **genLike;
    double    **scaledGenLike;
    int        MLmatAllele;
    int        MLpatAllele;
    int        numReadsDoublet;
    int        *readCountDoublet;
    int        maternalAlleleDoublet, paternalAlleleDoublet;
    double    **genLikeDoublet;
    double    **scaledGenLikeDoublet;
    int        MLmatAlleleDoublet;
    int        MLpatAlleleDoublet;
}
CellSiteStr;

typedef struct
{
    struct csite    *site;
    int                hasDoublet;
}
CellStr;

class ProgramOptions 
{
public:
    int doUsefixedMutationRate= 1;
    int doUseGenotypes = 0;
    int doUseFixedTree = 0;
    int doPrintTrees = 1;
    int doPrintSNVgenotypes = 0;
    int doPrintSNVhaplotypes = 0;
    int doPrintTrueHaplotypes = 1;
    int doPrintFullGenotypes = 0;
    int doPrintFullHaplotypes = 0;
    int doNGS = 0;
    int doPrintCATG = 0;
    int doPrintAncestors = 0;
    int doPrintMLhaplotypes = 0;
    int doSimulateReadCounts = 0;
    int doPrintTimes = 1;
    int doSimulateData = 1;
    int doPrintIUPAChaplotypes = 0;
    int doGeneticSignatures = 0;
    int doJC = 1;
    int doHKY = 0;
    int doGTR = 0;
    int doGTnR = 0;
    int thereIsMij = 0;
    int doSimulateFixedNumMutations = 0;
    int    thereisOutgroup = 0;
    int outgroupSelection = 0;
    int doEstimateTimesOriginClones = 0;
    int doAcceptTimes = 0;
    int doPrintSeparateReplicates = 1;
    double outgroupBranchLength_RootSample = 0;
    double outgroupBranchLength_Root1Root2 = 0;
    int doUseObservedCellNames = 0;
    double ADOrate = 0;
    double Eij[4];
    double Mij[4];
    double SNPrate = 0;
    double allelicImbalance = 0;
    double alphaBranches = 0;
    double alphaCoverage = 0;
    double alphaSites = 0;
    int alphabet = DNA;
    int altModel = 0;
    int coverage = 0;
    int doDemographics = 0;
    int doExponential = 0;
    int doUserGenome = 0;
    int doUserTree = 0;
    double doubletRate = 0;
    int equalBaseFreq = 0;
    double freq;
    double genotypingError = 0;
    double growthRate = 0;
    double haploidCoverageReduction = 0;
    double healthyTipBranchLength = 0;
    double meanAmplificationError = 0;
    double mutationRate = 0;
    double CNLOHrate =0;
    double meanADOcell  = 0;
    double varADOcell= 0;
    double meanADOsite= 0;
    double varADOsite=0;
    double deletionRate =0;
    int noisy = 1;
    double fixedADOrate =0 ;
    int doADOcell = 0;
    int doADOsite = 0;
    double nonISMRelMutRate = 0;
    int numCells = 0;
    int numDataSets = 0;
    int numNodes = 0;
    int Nscaling = 0;
    int numPeriods = 0;
    int numSites = 0;
    int numUserSignatures = 0;
    int ploidy = 0;
    long int seed = 0;
    long int userSeed = 0;
    int simulateOnlyTwoTemplates = 0;
    double varAmplificationError = 0;
    double sequencingError = 0;
    int rateVarCoverage = 0;
    double transformingBranchLength = 0;
    int rateVarAmongSites = 0;
    int rateVarAmongLineages = 0;
    double propAltModelSites = 0;
    double titv = 0;
    double seqErrorRate = 0;
    double dropoutRate = 0;
    int thereIsEij = 0;
    int numClones = 0;
    long double altModelMutationRate = 0;
    int numFixedMutations = 0;
    int TotalNumSequences = 0;
    int TotalTumorSequences =0;
    int doSimulateFromPriors = 0;
    int populationSampleSizesKnown = 0;
    int assignationKnown = 0;
    int numberClonesKnown = 0;
    int MutationAssignNum = 0;
    std::string healthyTipLabel;
    long double K;
    long double K_inference;
    int totalParameterValuesFromPriors;
    int minSampleSize;
    int maxSampleSize;
    double meanGenotypingError= 0;
    double varGenotypingError=0;
    //gsl_rng * r;
};

class MCMCoptions {
public:
     double numChains;
     double Niterations;
     double thinning;
     bool useGSLRandomGenerator;
     bool splitThetaDeltaTmoves;
     bool doThinning;
     int numberTumorCells;
     long double Deltafrom;
     long double Deltato;
     long double GrowthRatefrom;
     long double GrowthRateto;
     long  double MutRatefrom;
     long double MutRateto;
     long double totalEffectPopSizefrom;
     long double totalEffectPopSizeto;
     long double lambdafrom;
     long double lambdato;
     double tuningParameterDirichlet;
     int    startChainsRandomTree;
     int    doTopologicalMoves;
     int    printChainStateEvery;
     int    priorsType;//0: logUniform, 1: exponential, 2: powerLaw
     int    kernelType;//0: multiplier move, 1: normal centered current value,
     double slidingWindowSizeTotalEffectPopSize;
     double slidingWindowSizeGrowtRate;
     double meanGrowthRate;
     double dispersion;
     double tuningParameter;
     double fixedLambda;
     int maxNumberProposalAttempts;
     int numberWarmUpIterations;
     long double paramMultiplierMoveTheta;
     long double paramMultiplierMutationRate;
     long double paramMultiplierEffectPopSize;
     long double paramMultiplierGrowthRate;
     long double paramMultiplierTheta;
     long double paramMultiplierTimeOriginOldestPop;
     bool noData;//1: without data(just priors), 0: with data
     int useSequencesLikelihood;//1: true, 0 false
     int verbose;// 0 nothing, 1, 2 more verbose
     int fixTimeOriginInputTreeUnits;
     long double lambdaExponentialPriorTime;
     long double lambdaExponentialPriorMutationRate;
     long double lambdaExponentialPriorTotalEffectivePopSize;
     long double lambdaExponentialPriorGrowthRate;
     long double lambdaExponentialPriorSeqError;
     long double lambdaExponentialPriorDropoutError;
     long double sigmaNormalKernelTimeofOrigin;
     long double parameterPowerLawDistributionTotalEffectPopSize;
     long double parameterPowerLawDistributionMutationRate;
     long double parameterPowerLawDistributionGrowthRate;
     long double parameterPowerLawDistributionTimeOriginInputOldestPop;
     long double sigmaNormalKernelTotalEffectivePopulationSize;
     long double sigmaNormalKernelMutationRate;
     long double sigmaNormalKernelGrowthRate;
     double lengthIntervalMultiplier ;
     double lengthIntervalMultiplierDeltaT;
     double lengthIntervalMultiplierTheta;
     double lengthIntervalMultiplierTimeOriginOldestPop;
     double upperBoundTimeOriginInputOldestPop;
     double percentIterationsToComputeThinnig;
     int iterationToComputeThinnig;
     double thresholdAutoCorrelation;
     int maxNumberIndependentPosteriorValues;
     double thresholdAccceptanteRate;
     int iterationsToMonitorChain;
     int numberChainsPerTree;
     bool doMCMCMoveTimeOriginInputOldestPop;
     long double updateLengthMultiplierMCMCMove;
     bool fixedValuesForSimulation;
     long double lambdaExponentialMutationRateSimulation;
     long double lambdaExponentialGrowthRateSimulation;
     int burnInIterations;
     int numberIterationsAfterConvergence;
};

typedef struct {
    char  SNVgenotypesFile[MAX_NAME];
    char  SNVhaplotypesFile[MAX_NAME];
    char  trueHaplotypesFile[MAX_NAME];
    char  MLhaplotypesFile[MAX_NAME];
    char  fullGenotypesFile[MAX_NAME];
    char  fullHaplotypesFile[MAX_NAME];
    char  CATGfile[MAX_NAME];
    char  VCFfile[MAX_NAME];
    char  logFile[MAX_NAME];
    char  settingsFile[MAX_NAME];
    char  userTreeFile[MAX_NAME];
    char  treeFile[MAX_NAME];
    char  timesFile[MAX_NAME];
    char  userGenomeFile[MAX_NAME];
    char  SNVgenotypesDir[MAX_NAME];
    char  SNVhaplotypesDir[MAX_NAME];
    char  trueHaplotypesDir[MAX_NAME];
    char  MLhaplotypesDir[MAX_NAME];
    char  fullGenotypesDir[MAX_NAME];
    char  fullHaplotypesDir[MAX_NAME];
    char  treeDir[MAX_NAME];
    char  timesDir[MAX_NAME];
    char  CATGdir[MAX_NAME];
    char  VCFdir[MAX_NAME];
    char  resultsDir[MAX_NAME];
    char  dir[MAX_NAME];
    char  File[MAX_NAME];
    char  inputTreeFile[500];
    char  inputGenotypeFileFasta[500];
    char  inputGenotypeFilePhylip[500];
    char  likelihoodOuput[500];
}FilePaths;


typedef struct  {
    FILE *f;
    char  path[MAX_NAME];
}FilePath;

typedef struct {
    FilePath *fpTrees;
    FilePath *fpTrees2;
    FilePath *fpSNVgenotypes;
    FilePath *fpTimes;
    FilePath *fpTimes2;
    FilePath *fpSNVhaplotypes;
    FilePath *fpTrueHaplotypes;
    FilePath *fpFullGenotypes;
    FilePath *fpFullHaplotypes;
    FilePath *fpVCF;
    FilePath *fpCATG;
    FilePath *fpMLhaplotypes;
    FilePath *fplog;
    FilePath *fpTreeOutput;
    FilePath *fpLikelihood;
}Files;

typedef std::unordered_map<std::string,size_t> NameIdMap;

#endif /* data_types_h */
