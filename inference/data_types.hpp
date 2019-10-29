//
//  data_types.h
//  tumor_coal
//
//  Created by Fausto Fabian Crespo Fernandez on 2019-10-08.
//

#ifndef data_types_h
#define data_types_h

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include "definitions.hpp"

using namespace std;

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
    double freq = 0;
    double genotypingError = 0;
    double growthRate = 0;
    double haploidCoverageReduction = 0;
    double healthyTipBranchLength = 0;
    double meanAmplificationError = 0;
    double mutationRate = 0;
    int noisy = 1;
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
    double altModelMutationRate = 0;
    int numFixedMutations = 0;
    int TotalNumSequences = 0;
    int doSimulateFromPriors = 0;
    int populationSampleSizesKnown = 0;
    int assignationKnown = 0;
    int numberClonesKnown = 0;
    int MutationAssignNum = 0;
    string healthyTipLabel;
};

class MCMCoptions {
public:
    double numChains;
    double Niterations;
    double thinning;
    double Deltafrom;
    double Deltato;
    double MutRatefrom;
    double MutRateto;
    double totalEffectPopSizefrom;
    double totalEffectPopSizeto;
    double lambdafrom;
    double lambdato;
    double tuningParameterDirichlet;
    int    startChainsRandomTree;
    int    doTopologicalMoves;
    double slidingWindowSizeTotalEffectPopSize;
    double meanGrowthRate;
    double dispersion;
    double tuningParameter;
    double fixedLambda;
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
} FilePaths;

typedef struct {
    FILE *fpTrees;
    FILE *fpTrees2;
    FILE  *fpSNVgenotypes;
    FILE  *fpTimes;
    FILE  *fpTimes2;
    FILE  *fpSNVhaplotypes;
    FILE  *fpTrueHaplotypes;
    FILE  *fpFullGenotypes;
    FILE *fpFullHaplotypes;
    FILE *fpVCF;
    FILE *fpCATG;
    FILE *fpMLhaplotypes;
}Files;

#endif /* data_types_h */
