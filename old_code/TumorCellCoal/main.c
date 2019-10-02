//
//  CoalCC.c
//  CoalCC
//  TumorCellCoal
//
//  Created by Fausto Fabian Crespo Fernandez, Miguel Arenas, David Posada on 11/12/18.

//  Copyright © 2018 Fausto Fabian Crespo Fernandez, Miguel Arenas, Carsten Wiuf & David Posada. All rights reserved.
//
//   is a coalescent-based program oriented to simulate the evolution of somatic cells under:
//
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <sys/timeb.h>
#include <float.h>
#include "pll_msa.h"
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include "definitions.h"

#include "kseq.h"
#include "tumorcellcoal.h"
#include <zlib.h>
#include <unistd.h>
//#include "pll.h"
//#include "pll_binary.h"
//#include "pll_msa.h"
#include "tree.h"
#include "likelihood.h"
#include "model.h"
#include "random.h"
#include "population.h"
#include "prior.h"
#include "mutationModel.h"
#include "util.h"
#include "fileFuntions.h"
#include "proposal.h"
#include "MCMCchain.h"

#include "pll_optimize.h"
#include "pll_tree.h"
#include "pllmod_algorithm.h"
#include "pllmod_common.h"
//#include "pllmod_util.h"
#include <stdlib.h>

#include <stdarg.h>
#include <search.h>
#include <time.h>

//KSEQ_INIT(FILE*, read);
KSEQ_INIT(int, read);

//#include "signatures.h"
//#include <mpi.h>

#ifndef macintosh
#include <sys/types.h>
#include <sys/stat.h>
#endif
#ifdef MAC
#include <sioux.h>
#include <console.h>
#include <unix.h>
#endif

/* functions */


void InitializeGenomes (TreeNode *p, long int *seed,  int alphabet, int doUserGenome, int numSites, SiteStr* allSites, int doGeneticSignatures, double cumfreq[4],double *triNucFreq, char **cellNames);





/***************************** Main *******************************/
/* Main function */
int main(int argc, char **argv)
{
    
    //variables
    
    static int    dataSetNum, i, j, k, z, m, w;
    static float     start, secs;
    long int  seed, seedFirst, originalSeed;
    char    File[80];
    char    dir[80];
    double    thistimeNi, TotalProbability, LocalProbability, aboveTerm, belowTerm, ranHere;
    struct    timeb tmb;
    
    ProgramOptions programOptions;
    Files files;
    FilePaths filePaths;
    MCMCoptions mcmcOptions;
    
    int         AttemptsAcceptation;
    FILE    *fp;
    /* File pointers */
    FILE            *fpSNVgenotypes, *fpFullGenotypes, *fpSNVhaplotypes, *fpTrueHaplotypes, *fpFullHaplotypes, *fpMLhaplotypes, *fpCATG, *fpVCF, *fpLog, *fpUserTree, *fpUserGenome;
    
    TreeNode    *nodes;
    //TreeNode **treeRootInit;
    long int    userSeed;
    static int             numClones;
    static int       TotalNumSequences, TotalN, numDataSets, noisy, numCA, numMIG;
    static double          TotalBirthRate, TotalDeathRate;
    static int      *CloneNameBegin, *CloneSampleSizeBegin, *ClonePopSizeBegin;
    static double   *CloneBirthRateBegin, *CloneDeathRateBegin, *CloneTimeOriginInput;
    static double   *CloneGrowthRateBegin;
    static double   *CloneTimeOriginInputSTDPOP, *CloneDelta, *ClonePopSizeMeffectBegin;
    static double   *ListMigrationTimesInitial, *ListMigrationTimesUpdated;
    static int      *CloneNameBeginOrderByTimes;
    static int      *CloneNameBeginOrderByTimesOnlyControl;
    static int             control1;
    static int             numSimClones;
    static double          expectedPopSize;
    double      cumNumCA, meanNumCA, cumNumMIG, meanNumMIG, numEventsTot;
    static double      mutationRate;
    char      treeFile[20], timesFile[20];
    static double outgroupBranchLength_RootSample, expTMRCA, expVarTMRCA, outgroupBranchLength_Root1Root2;
    static int    thereisOutgroup, doPrintTrees, doPrintTimes, outgroupSelection, doEstimateTimesOriginClones, doAcceptTimes;
    double      counterTime, counterTimeInit, countTMRCA;
    static double *varTimeGMRCA;
    static int    *varEvent;
    int      numNodes;
    int      Nscaling;
    int      nextAvailable, intLabel;
    int numSites;
    int alphabet;
    double  titv;/* transition/transversion rate ratio */
    FILE            *fpTrees, *fpTimes, *fpTrees2, *fpTimes2;
    double freq[4];
    double cumfreq[4];
    double Mij[4][4];
    double cumMij[4][4];
    
    double Eij[4][4];
    double cumEij[4][4];
    int thereIsMij, thereIsEij;
    int *** data;
    
    data=NULL;
    
     Population **populations;

    //////////////////////////////////////////////////////////
    
     TreeNode    *treeNodes, *coalTreeMRCA, *healthyRoot, *healthyTip;
   // TreeNode  *treeTips;
     SiteStr     *allSites;
    TriNucStr   *triNucleotideMaternal, *triNucleotidePaternal;
    CellStr     *cell;
  
    static int      *SNVsites, *DefaultModelSites, *AltModelSites, *variantSites, *SFS, *numberDifferences;
    static int    tipLabel;
    static int    ploidy, numCells, *Nbegin, *Nend,  *cumDuration, numPeriods;
    static int     numAltModelSites, numDefaultModelSites, numISMmutations, altModel;
    static int    numMU, numDEL, numCNLOH, numProposedMU, numSNVs, numFixedMutations, numSNVmaternal, zeroSNVs;
    static int    numISMdeletions, numISMCNLOH;
    static double meanNumSNVs, meanNumMU, meanNumDEL, meanNumCNLOH;
    static double   cumNumSNVs, cumNumMU, cumNumDEL, cumNumCNLOH, cumCountMLgenotypeErrors;
    static double cumNumMUSq, cumNumSNVsSq, cumNumDELSq, cumNumCNLOHSq;
    static double varNumMU, varNumSNVs, varNumDEL, varNumCNLOH;
    static double expNumMU, expVarNumMU;
    static double theta, healthyTipBranchLength, transformingBranchLength, totalTreeLength;
    static double nonISMRelMutRate, propAltModelSites, altModelMutationRate, deletionRate, CNLOHrate;
    static char   SNVgenotypesFile[MAX_NAME], SNVhaplotypesFile[MAX_NAME], trueHaplotypesFile[MAX_NAME],MLhaplotypesFile[MAX_NAME], fullGenotypesFile[MAX_NAME], fullHaplotypesFile[MAX_NAME];
    static char  CATGfile[MAX_NAME], VCFfile[MAX_NAME], logFile[MAX_NAME], settingsFile[MAX_NAME], userTreeFile[MAX_NAME], userGenomeFile[MAX_NAME];
    static char   SNVgenotypesDir[MAX_NAME], SNVhaplotypesDir[MAX_NAME], trueHaplotypesDir[MAX_NAME], MLhaplotypesDir[MAX_NAME], fullGenotypesDir[MAX_NAME], fullHaplotypesDir[MAX_NAME];
    static char   treeDir[MAX_NAME], timesDir[MAX_NAME], CATGdir[MAX_NAME], VCFdir[MAX_NAME];
    static char   resultsDir[MAX_NAME],  *CommandLine, *treeString, *taxonName, **cellNames;
    static int    doPrintSNVgenotypes, doPrintSNVhaplotypes, doPrintTrueHaplotypes, doPrintMLhaplotypes, doPrintFullHaplotypes, doPrintFullGenotypes, doPrintTree, doUserTree, doUserGenome, doInference;
    static int    doPrintAncestors, doPrintCATG, doPrintSeparateReplicates, doPrintIUPAChaplotypes, doNGS;
    static int    doExponential, doDemographics, doSimulateData, doSimulateFixedNumMutations,doSimulateReadCounts, taxonNamesAreChars;
    static int    doJC, doHKY, doGTR, doGTnR, doGeneticSignatures;
    //GTR parameters
     static double    Rmat[6], NRmat[12], Cijk[256], Root[4];
    //////////////
    static int      rateVarAmongSites, rateVarAmongLineages, rateVarCoverage, equalBaseFreq,  coverage, countMLgenotypeErrors;
    static double *periodGrowth, growthRate, sequencingError, ADOrate, allelicImbalance, haploidCoverageReduction, genotypingError, meanAmplificationError, varAmplificationError, doubletRate;
    static double TMRCA, cumTMRCA, cumTMRCASq, meanTMRCA,  varTMRCA;
    static double    kappa, beta, freqR, freqY, freqAG, freqCT,  alphaSites, alphaBranches;
   
    static double SNPrate, alphaCoverage;
    static int    HEALTHY_ROOT, TUMOR_ROOT;
    static int    readingParameterFile, simulateOnlyTwoTemplates;
    static int    TipNodeNum, IntNodeNum;
    static char   *maternalUserGenome, *paternalUserGenome;
    static int    complementBase[4] = {3,2,1,0};
    static int    targetTriChange[6] = {0,2,3,0,1,2};
    static int    *triMutationsCounter;
    static int    numUserSignatures;
    static double *signatureWeight;
    static int    *signatureID;
    
    
    // read real data from fasta file
    int** ObservedData ;
    char TreeFile[MAX_NAME];
    char SeqFile[MAX_NAME];
  
    char *ObservedCellNames[programOptions.numCells];
    
    char newickString[20000];
    char *newickString2;
    
//    extern double   Qij[16], mr;
//   extern double   ***selectedSignature;
//   extern double ****geneticSignature;
//
//    extern double **signatureProbs;
//    extern double *triNucFreq;
    
#ifdef CHECK_MUT_DISTRIBUTION
    static int    *MutCount, *SiteMut, dataSetsWithSNVs;
    static double sumPos, meansumPos;
  #endif
    
    /* Default settings */
    Initialize( Eij, Mij, freq,  &programOptions );

    //////////////////////////////////////////////
    
  
//    /* Pointer initializacion */
    
    CloneNameBegin = NULL;
    CloneSampleSizeBegin = NULL;
    ClonePopSizeBegin = NULL;
    
    CloneBirthRateBegin = NULL;
    CloneDeathRateBegin = NULL;
    CloneTimeOriginInput = NULL;
    CloneGrowthRateBegin = NULL;
    CloneTimeOriginInputSTDPOP = NULL;
    CloneDelta = NULL;
    ClonePopSizeMeffectBegin = NULL;
    ListMigrationTimesInitial = NULL;
    ListMigrationTimesUpdated = NULL;
    CloneNameBeginOrderByTimes = NULL;
    CloneNameBeginOrderByTimesOnlyControl = NULL;
    
    /* Read input settings */
    /* arguments from external place */
#ifdef MAC
    if ((fp = freopen("parameters", "r", stdin)) == NULL)       /* input from parameters file and it reads */
    {
        fprintf (stderr, "\nERROR: Can't read parameters file.");
        PrintUsage();
    }
    //ReadParametersFromFile();
    
          ReadParametersFromFile(&programOptions, &filePaths, &CloneNameBegin, &CloneSampleSizeBegin, &ClonePopSizeBegin, &CloneBirthRateBegin, &CloneDeathRateBegin, &CloneTimeOriginInput,Mij,freq);
    Sioux();
    fclose(fp);
#else
    // ReadParametersFromCommandLine (argc, argv);               /* input from external arguments */
    if (argc <= 2)
    {
         if ((fp = freopen("parameters", "r", stdin)) != NULL)
        {
            ReadParametersFromFile(&programOptions, &filePaths, &CloneNameBegin, &CloneSampleSizeBegin, &ClonePopSizeBegin, &CloneBirthRateBegin, &CloneDeathRateBegin, &CloneTimeOriginInput,Mij,freq);
        }
        else
        {
            fprintf (stderr, "\nERROR: No parameters specified (use command line or parameter file)");
            PrintUsage();
        }
    }
#endif
    //////////////////////////////////////////////////////////////////////////
    
    
    int cumNumMUperTree;
    /* initialize a few variables */
    
    cumNumCA = cumNumMU = cumNumMUSq = cumNumSNVs= cumNumSNVsSq = cumCountMLgenotypeErrors = 0;
    zeroSNVs = cumTMRCA = cumTMRCASq = 0;
    cumNumMUperTree=0;
    programOptions.altModelMutationRate = programOptions.mutationRate * programOptions.nonISMRelMutRate;
    
    if ( programOptions.alphabet == DNA)
    {
        /* initialize cumfreq */
        cumfreq[0] = freq[0];
        for (i=1; i<4; i++)
            cumfreq[i] = cumfreq[i-1] + freq[i];
        
        /* initialize cumMij */
        for (i=0; i<4; i++)
        {
            cumMij[i][0] = Mij[i][0];
            for (j=1; j<4; j++)
                cumMij[i][j] = cumMij[i][j-1] + Mij[i][j];
        }
        
        /* initialize cumEij */
        for (i=0; i<4; i++)
        {
            cumEij[i][0] = Eij[i][0];
            for (j=1; j<4; j++)
                cumEij[i][j] = cumEij[i][j-1] + Eij[i][j];
        }
    }
    ftime(&tmb);
    //printf("tmb.time     = %ld (seconds), %ld \n", tmb.time, fabs(tmb.time));
    //printf("tmb.millitm  = %ld (mlliseconds), %ld \n", tmb.millitm, fabs(tmb.millitm));
  
    programOptions.seed =  programOptions.seed + tmb.millitm;
    if ( programOptions.userSeed > 0)
        //seed = userSeed;
        originalSeed =  programOptions.seed;
    for (i = 0; i < 10; i++)
        RandomUniform(& (programOptions.seed));   /* function that generates random seed */
    seedFirst =  programOptions.seed;
    
    /* Presentation */
    start = clock();
    
    //allocate memory for the population structs
    populations = (Population**)malloc (sizeof(struct Population*)  * programOptions.numClones);
    if (!populations)
    {
        fprintf (stderr, "Could not allocate populations (%lu bytes)\n", (programOptions.numClones)  * (long) sizeof(Population*));
        exit (1);
    }
     InitListClones(populations, programOptions.numClones, programOptions.noisy, CloneNameBegin, CloneSampleSizeBegin, CloneBirthRateBegin,  CloneDeathRateBegin, ClonePopSizeBegin, programOptions.TotalNumSequences);
     InitListClonesTimes(populations, programOptions.numClones,  &programOptions.doEstimateTimesOriginClones,  CloneTimeOriginInput  );
    
    InitNumberNodes(&TotalBirthRate, &TotalDeathRate, &TotalN, populations, &programOptions);
    ListClonesAccordingTimeToOrigin2(populations, programOptions.numClones);
    
    /* set file dirs and names */
    InitFilesPathsOptions(&filePaths, &programOptions);
    if(programOptions.doPrintSeparateReplicates==NO)
    {
        CreateFolder(filePaths,programOptions, &files) ;
    }
    /* more defaults and memories */
    cumNumCA = 0.0;
    cumNumMIG = 0.0;
    counterTime = 0.0;
    numEventsTot = 0.0;
    i = j = 0;
    
    //programOptions.numCells = programOptions.numCells -1;//not counting the outgroup / healthycell
    programOptions.numCells = programOptions.numCells;
    programOptions.TotalNumSequences =programOptions.numCells;
    
    HEALTHY_ROOT = 2 * programOptions.numCells;
    TUMOR_ROOT = (2 * programOptions.numCells) - 1;
    
    MCMCoptions opt;
    if  (programOptions.doSimulateData==YES)
    {
        programOptions.doSimulateFromPriors =NO;
       programOptions.doUseObservedCellNames=NO;
    }
    else {
        programOptions.doSimulateFromPriors =YES;
        programOptions.doUseObservedCellNames=YES;
        
    }
    int parameterNum = 10;
    z=0;
    double randomDelta;
    double randomT;
    mcmcOptions.meanGrowthRate=0.6;
    mcmcOptions.dispersion=0.4;
    double randomTimeOfOrigin=0;
    cumNumMUperTree=0;
    double shape=9;
    programOptions.MutationAssignNum=1;
    double  oldestPopSize;
    double totalPopSize;
    double gammaParam;
    double lambdaPoissonOldestPopSize= 10000;
    //parameters for the log uniform prior distribution for deltas
    mcmcOptions.Deltafrom = -4;
    mcmcOptions.Deltato = 1;
    
    mcmcOptions.numChains =2;
    mcmcOptions.Niterations = 200000;
    mcmcOptions.thinning =  500;
    
    //parameters for the log uniform prior distribution for mutation rate
    mcmcOptions.MutRatefrom = -13;
    mcmcOptions.MutRateto = -8;
    mcmcOptions.totalEffectPopSizefrom=7;
    mcmcOptions.totalEffectPopSizeto=13;
    double logConditionalLikelihoodTree;
    
    int paramNum=0;
    int TajimaD=0;
    double AvgHeterocigocity =0.0;
    double logConditionalLikelihoodSequences;
    int currentSizeNewick;
    
    /////////// memory allocation ///////
    
    /* allocate memory for site information (equal for maternal and paternal) */
    allSites = (SiteStr*) calloc (programOptions.numSites, sizeof(SiteStr));
    if (!allSites)
    {
        fprintf (stderr, "Could not allocate the allSites structure\n");
        exit (-1);
    }
    for (i=0; i< programOptions.numSites; i++)
    {
        allSites[i].alternateAlleles = (int *) calloc (4, sizeof(int));
        if (!allSites[i].alternateAlleles)
        {
            fprintf (stderr, "Could not allocate the allSites[i].alternateAlleles[] structure\n");
            exit (-1);
        }
    }
    
    /* the arrays below keep the index for different types of sites */
    SNVsites = (int*) malloc (programOptions.numSites* sizeof(int));
    if (!SNVsites)
    {
        fprintf (stderr, "Could not allocate the SNVsites structure\n");
        exit (-1);
    }
    SFS = (int*) malloc (programOptions.numSites* sizeof(int));
    if (!SFS)
    {
        fprintf (stderr, "Could not allocate the SNVsites structure\n");
        exit (-1);
    }
    /* the arrays below keep the index for different types of sites */
    variantSites = (int*) malloc (programOptions.numSites* sizeof(int));
    if (!variantSites)
    {
        fprintf (stderr, "Could not allocate the variantSites structure\n");
        exit (-1);
    }
    
    DefaultModelSites = (int*) malloc (programOptions.numSites* sizeof(int));
    if (!DefaultModelSites)
    {
        fprintf (stderr, "Could not allocate the DefaultModelSites structure\n");
        exit (-1);
    }
    
    AltModelSites = (int*) malloc (programOptions.numSites* sizeof(int));
    if (!AltModelSites)
    {
        fprintf (stderr, "Could not allocate the AltModelSites structure\n");
        exit (-1);
    }
    TreeNode  **treeTips;
    treeTips =  malloc (programOptions.TotalNumSequences* sizeof(TreeNode*));
    if (!treeTips)
    {
        fprintf (stderr, "Could not allocate the treeTips array\n");
        exit (-1);
    }
    double    *proportionsVector = (double *) calloc((programOptions.numClones), (long) sizeof(double));
    if (!proportionsVector)
    {
        fprintf (stderr, "Could not allocate proportions vector (%lu bytes)\n", (programOptions.numClones ) * (long) sizeof(double));
        exit (-1);
    }
    /* Variance memories */
    varEvent = (int *) malloc(programOptions.numDataSets * (long) sizeof(int));
    if (!varEvent)
    {
        fprintf (stderr, "Could not allocate varEvent (%lu bytes)\n", programOptions.numDataSets * (long) sizeof(int));
        exit (1);
    }
    varTimeGMRCA = (double *) malloc(programOptions.numDataSets* (long) sizeof(double));
    if (!varTimeGMRCA)
    {
        fprintf (stderr, "Could not allocate varTimeGMRCA (%lu bytes)\n", programOptions.numDataSets  * (long) sizeof (double));
        exit (1);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////
    TreeNode *root;
    if (programOptions.doSimulateData == YES)
    {
        //do simulations
        programOptions.doUseObservedCellNames=NO;
        PrintTitle();   /* write title */
        PrintDate();    /* write date */
        ValidateParameters(&programOptions,
                           // programOptions.noisy, programOptions.numDataSets, programOptions.Nscaling, programOptions.numClones,
                           CloneNameBegin , CloneSampleSizeBegin, ClonePopSizeBegin);
        CloneGrowthRateBegin = (double *) malloc((programOptions.numClones + 1)* (long) sizeof(double));
        ClonePopSizeMeffectBegin = (double *) malloc((programOptions.numClones + 1)* (long) sizeof(double));
        if (programOptions.noisy > 1)
        {
            printf("Demographic parameters per clone\n");
            printf("\tClone number\t|\tSample size\t|\tPopulation size\t|\tEfective population size\t|\tBirth rate\t|\tDeath rate\t|\tGrowth rate (br - dr)\n");
        }
        InitListPossibleMigrations(populations,programOptions.numClones);
        InitPopulationsCoalescentEvents( programOptions.numClones,  populations);
        for (dataSetNum = 0; dataSetNum < programOptions.numDataSets; dataSetNum++)// dataSetNum refers to a simulated tree number
        {
            if (programOptions.doPrintSeparateReplicates == YES)
                PrepareSeparateFiles(0, 1, dataSetNum,  &filePaths, &programOptions, &files);
            /* Reorganize seed per replicate */
            //seed = seedFirst+dataSetNum+10;
            numCA = numMU = numDEL = numProposedMU = TMRCA = numSNVmaternal = 0;
            seed = seedFirst + dataSetNum * 1000;
            //fprintf(stderr, "\n seed = %lu \n", seed);
            /* reset variables per replicate */
            if (programOptions.noisy > 0)
            {
                fprintf (stderr, "\rReplicate #%3d/%d", dataSetNum + 1, programOptions.numDataSets);
                fflush (stdout);
            }
            varEvent[dataSetNum] = 0;
            varTimeGMRCA[dataSetNum] = 0.0;
            counterTimeInit = 0.0;
            numCA = numMIG = 0;
            countTMRCA = 0.0;
            /* coalescent tree */
           
            if (programOptions.noisy > 1)
                fprintf (stderr, "\n>> Start coalescent tree .. \n");
            MakeCoalescenceTree2 (&programOptions.seed, populations,
                                  &programOptions.numNodes,
                                  programOptions.numClones,
                                  &programOptions,
                                  cumNumCA,
                                  meanNumCA,
                                  cumNumMIG,
                                  meanNumMIG,
                                  &numMIG,
                                  &numCA,
                                  &numEventsTot,
                                  &nodes,
                                  treeTips,
                                  &root,
                                  ObservedCellNames, NULL
                                  ) ;
            if (programOptions.noisy > 1)
                fprintf (stderr, "\n>> Finishing coalescent tree ... DONE");
            /* counters */
            
            //AssignTreeNodestoObservedCellNames(treeRootInit[0],  ObservedCellNames);
             //logConditionalLikelihoodTree= LogConditionalLikelihoodTree(root, nodes, populations,  programOptions.numClones);
            
            cumNumCA += numCA;
            cumNumMIG += numMIG;
           // countTMRCA = treeRootInit[0]->timePUnits;
            countTMRCA = root->timePUnits;
            
            //fprintf ( stderr, "\n countTMRCA = %lf\n", countTMRCA);
            varTimeGMRCA[dataSetNum] = countTMRCA;
            varEvent[dataSetNum] = numCA + numMIG;
            // counterTime = counterTime + counterTimeInit;
            /*************** output files *************/

            newickString2=NULL;
            newickString2 = toNewickString2 ( root, programOptions.mutationRate,     programOptions.doUseObservedCellNames);
            printf("\n newick = %s  \n", newickString2);
            
            if (programOptions.doPrintTrees == YES)
            {
//                PrintTrees(dataSetNum, treeRootInit, files.fpTrees, programOptions.mutationRate);
//                PrintTrees2(dataSetNum, treeRootInit, files.fpTrees2, programOptions.mutationRate,  ObservedCellNames);
                PrintTrees(dataSetNum, &root, files.fpTrees, programOptions.mutationRate, programOptions.doUseObservedCellNames);
                PrintTrees2(dataSetNum, &root, files.fpTrees2, programOptions.mutationRate,  ObservedCellNames, programOptions.doUseObservedCellNames);
            }
            if (programOptions.doPrintTimes == YES)
            {
                PrintTimes(dataSetNum, files.fpTimes, programOptions.mutationRate, nodes, programOptions.thereisOutgroup);
                PrintTimes2(dataSetNum, files.fpTimes2, programOptions.mutationRate, nodes, programOptions.thereisOutgroup);
            }
            if (programOptions.noisy > 1)
            {
                fprintf (stderr, "\nData set %d", dataSetNum + 1);
                fprintf (stderr, "\n\tNumber of coalescence events   =   %d", numCA);
                fprintf (stderr, "\n\tNumber of migration events     =   %d", numMIG);
            }
           // totalTreeLength = SumBranches(treeRootInit[0], programOptions.mutationRate);
            totalTreeLength = SumBranches(root, programOptions.mutationRate);
            cumNumMUperTree=0;
            
            free(newickString2);
            newickString2=NULL;
            
            if (programOptions.doSimulateData ==YES)
            {
                for (z = 0; z < programOptions.MutationAssignNum; z++)
                {
                    numMU=0;
                    if (programOptions.doPrintSeparateReplicates == YES)
                        PrepareSeparateFilesGenotypes(1, dataSetNum, z,
                                                      &filePaths, &programOptions,&files);
                    
                    //here there was the code
                    if (programOptions.doSimulateData == YES)
                    {
                        numISMmutations = 0;
                        numISMdeletions = 0;
                        numISMCNLOH = 0;
                    }
                    for (i=0; i< programOptions.numSites; i++)
                    {   allSites[i].numMutations =0;
                        allSites[i].numMutationsMaternal =0;
                        allSites[i].numMutationsPaternal =0;
                    }
                    //InitializeGenomes (treeRootInit[0], &seed, programOptions.alphabet, programOptions.doUserGenome,programOptions.numSites,  allSites, programOptions.doGeneticSignatures,cumfreq, triNucFreq, cellNames);
                    InitializeGenomes (root, &programOptions.seed, programOptions.alphabet, programOptions.doUserGenome,programOptions.numSites,  allSites, programOptions.doGeneticSignatures,cumfreq, triNucFreq, cellNames);
                    // SNPrate=0.01;
                   // HEALTHY_ROOT=treeRootInit[0]->label;
                    //TUMOR_ROOT=treeRootInit[0]->label;
                    HEALTHY_ROOT=root->label;
                    TUMOR_ROOT=root->label;
                    //        if (SNPrate > 0)
                    //            AddGermlineVariation (treeRootInit[0], &seed,  numSites, SNPrate, allSites, alphabet,  data,   HEALTHY_ROOT, cumMij );
                    EvolveSitesOnTree (root, MATERNAL, & programOptions.seed, programOptions.rateVarAmongSites,  programOptions.numSites,  allSites, programOptions.doGeneticSignatures, programOptions.alphaSites, programOptions.propAltModelSites ,  numDefaultModelSites, numAltModelSites, DefaultModelSites, AltModelSites,  totalTreeLength , &numISMmutations, programOptions.numFixedMutations, numSNVmaternal,  programOptions.doSimulateFixedNumMutations,  programOptions.alphabet,  data,  &numMU, cumMij,  programOptions.altModel, programOptions.altModelMutationRate, programOptions.doUserTree,  programOptions.doJC,  programOptions.doHKY,  programOptions.doGTR,
                                       programOptions.doGTnR,  freqR,  freqY,
                                       freqAG, freqCT, programOptions.titv, freq, Mij ,   Root,  Cijk);
                    EvolveSitesOnTree (root, PATERNAL, & programOptions.seed, programOptions.rateVarAmongSites,  programOptions.numSites,  allSites, programOptions.doGeneticSignatures, programOptions.alphaSites, programOptions.propAltModelSites ,  numDefaultModelSites, numAltModelSites, DefaultModelSites, AltModelSites,  totalTreeLength , &numISMmutations, numFixedMutations, numSNVmaternal,  programOptions.doSimulateFixedNumMutations,  programOptions.alphabet,  data,  &numMU, cumMij,  programOptions.altModel, programOptions.altModelMutationRate, doUserTree,  programOptions.doJC,  programOptions.doHKY,  programOptions.doGTR,
                                       programOptions.doGTnR,  freqR,  freqY,
                                       freqAG, freqCT, programOptions.titv, freq, Mij,   Root,  Cijk );
                    cumNumMU += numMU;
                    cumNumMUSq += pow(numMU,2);
                    numSNVs = CountTrueVariants (nodes, programOptions.numSites,  programOptions.numCells,   root, allSites, variantSites, SNVsites);
                    TajimaD =  computeTajimaD(treeTips, programOptions.numSites,  programOptions.numCells);
                    AvgHeterocigocity= ComputeAvgHeterocigocity( programOptions.numSites,  programOptions.numCells, treeTips, root, allSites, SNVsites);
                    numberDifferences = (int*) calloc (numSNVs, sizeof(int));
                    if (!numberDifferences)
                    {
                        fprintf (stderr, "Could not allocate the variantSites structure\n");
                        exit (-1);
                    }
                    if (programOptions.doPrintTrueHaplotypes == YES)
                    {
                        if (programOptions.doPrintSeparateReplicates == NO)
                            fprintf (files.fpTrueHaplotypes, "[#%d]\n", z+1);
                        //
                        PrintTrueFullHaplotypes (files.fpTrueHaplotypes,  nodes, root , programOptions.numNodes, programOptions.doPrintIUPAChaplotypes, programOptions.doPrintAncestors, programOptions.numSites,  programOptions.numCells, programOptions.alphabet, programOptions.doUserTree , data,    programOptions.doNGS,   cellNames, cell, HEALTHY_ROOT, TUMOR_ROOT, ObservedCellNames, programOptions.doUseObservedCellNames);
                    }
                    
                    /// here there was a code snippet
                    computeUnfoldedISMSFS(programOptions.numSites, allSites, numSNVs, SNVsites,  SFS, numberDifferences);
                    
                    if (programOptions.doPrintTrees ==YES && programOptions.doPrintSeparateReplicates == YES){
                        
                        fclose(files.fpTrees);
                        fclose(files.fpTrees2);
                    }
                    if (programOptions.doPrintTimes ==YES && programOptions.doPrintSeparateReplicates == YES){
                        
                        fclose(files.fpTimes);
                        fclose(files.fpTimes2);
                    }
                    if (programOptions.doPrintTrueHaplotypes ==YES){
                        if (programOptions.doPrintSeparateReplicates == YES)
                            fclose(files.fpTrueHaplotypes);
                        
                    }
                    sprintf(TreeFile,"%s/%s/%s_2_%04d_%04d.tre", filePaths.resultsDir, filePaths.treeDir, filePaths.treeFile, paramNum+1, dataSetNum+1);
                    
                    sprintf(SeqFile,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.trueHaplotypesDir, filePaths.trueHaplotypesFile, paramNum+1, dataSetNum+1, z +1);
                    
                    //        cumNumSNVs += numSNVs;
                    //        cumNumSNVsSq += pow(numSNVs,2);
                    
                    if (programOptions.doPrintSNVgenotypes == YES && numSNVs > 0) /* we only print replicates with variation */
                    {
                        if (programOptions.doPrintSeparateReplicates == NO)
                        {
                            fprintf (files.fpSNVgenotypes, "[#%d]\n", dataSetNum+1);
                        }
                        PrintSNVGenotypes (files.fpSNVgenotypes,programOptions.numClones, populations,  nodes, root, programOptions.doPrintAncestors,  programOptions.doNGS, numCells,  numSNVs,  SNVsites,  programOptions.alphabet,  doUserTree,data, cellNames,  TUMOR_ROOT,  HEALTHY_ROOT,  cell, SFS, numberDifferences, TajimaD, programOptions.mutationRate);
                    }
                    //here another code
                    free(numberDifferences);
                    numberDifferences=NULL;
                
                }/* end of mutation simulation process */
            }//end of if(doSimulateData)
            cumNumMUperTree=numMU;
            meanNumMU  = cumNumMUperTree /  (double) programOptions.MutationAssignNum;
            if (programOptions.noisy >1)
                fprintf(stderr, "avg number of mutations for tree %d is  %lf \n",dataSetNum, meanNumMU );
  
            TreeNode *p;
            
            for(i=0; i<programOptions.numNodes; i++){
                p= nodes + i;
                //             free( p->cellName);
                //    p->cellName=NULL;
                free( p->maternalSequence);
                p->maternalSequence=NULL;
                free( p->paternalSequence);
                p->paternalSequence=NULL;
                free(p->numbersMutationsUnderSubtreePerSite);
                p->numbersMutationsUnderSubtreePerSite=NULL;
                free(p->numbersMaternalMutationsPerSite);
                p->numbersMaternalMutationsPerSite=NULL;
                free(p->numbersPaternalMutationsPerSite);
                p->numbersPaternalMutationsPerSite=NULL;
            }
            free (nodes);
            nodes=NULL;
            //free (treeRootInit);
            //treeRootInit=NULL;
            
//            free (treeTips);
//            treeTips=NULL;
        } /* end of replicates */
        Population *pop;
        for( i = 0 ; i < programOptions.numClones; i++)
        {
            pop=*(populations + i);
            free(pop->idsActiveGametes);
            pop->idsActiveGametes =NULL;
            free(pop->CoalescentEventTimes);
            pop->CoalescentEventTimes = NULL;
            for (j = 1; j < ( pop->order); j++)
            {
                if (pop->order >0){
                    //free(pop->immigrantsPopOrderedModelTime[j]);
                    pop->immigrantsPopOrderedModelTime[j]=NULL;
                }
            }
            free( pop->migrationTimes);
            pop->migrationTimes=NULL;
            if (pop->order >0)
            {
                free(pop->immigrantsPopOrderedModelTime);
                pop->immigrantsPopOrderedModelTime=NULL;
            }
         }
    }//end simulation
    /////////////////////////////////////////////////////////////////////////////////////////////
    else { //do MCMC inference with known number of clones and known assignation
        programOptions.numberClonesKnown=YES;
        programOptions.populationSampleSizesKnown = YES;
        mcmcOptions.slidingWindowSizeTotalEffectPopSize = 20000;
        programOptions.doUseGenotypes = YES;
        programOptions.doUseFixedTree =NO;
        
        mcmcOptions.tuningParameter = 1;
        mcmcOptions.thinning  = 1000;
        
        programOptions.seqErrorRate=programOptions.sequencingError;
        programOptions.dropoutRate=programOptions.ADOrate;
        
        pll_msa_t * msa;
        
        char *fileName; //="Houetal_HS.fasta";
        char *fileNamePhylip ;//="Houetal_HS1.phylips";
        
        
        fileName ="SimulatedDataJC.fasta";
        fileNamePhylip ="SimulatedDataJC.phylip";
        //char * fileName ="Simulated_data.phylip";
        //sprintf(SeqFile,"Houetal_HS1.phylips" );
        sprintf(SeqFile,"Simulated_dataJC.phylip" );
        
        int *sampleSizes =(int *) malloc(programOptions.numClones* (long) sizeof(int));
        if (!sampleSizes)
        {
            fprintf (stderr, "Could not allocate samplesSizes");
            exit (1);
        }
        
//        sampleSizes[0]=57;
//        sampleSizes[1]=5;
//        sampleSizes[2]=3;
//        programOptions.numClones = 3;
//        programOptions.numCells = 65;
//        programOptions.TotalNumSequences= 65;
        
        sampleSizes[0]=5;
        sampleSizes[1]=3;
        sampleSizes[2]=4;
        programOptions.numClones = 3;
     
        
        programOptions.numNodes = 2 * programOptions.TotalNumSequences + programOptions.numClones+ 10;
        programOptions.numCells = programOptions.TotalNumSequences;
        
        if (programOptions.numberClonesKnown==YES)
        {
            mcmcOptions.totalEffectPopSizefrom = round(log(programOptions.numCells +2));
            mcmcOptions.totalEffectPopSizefrom = round(log(100* programOptions.numCells +200 ));
            
        }
        
        //char *fileName ="Nietal_HS.fasta";
        ReadParametersFromFastaFile(fileName,  &programOptions);
        ObservedData = (int**) malloc (programOptions.numCells * sizeof(int*));
        if (!ObservedData)
        {
            fprintf (stderr, "Could not allocate the ObservedData memory\n");
            exit (-1);
        }
        for (i=0; i< programOptions.numCells; i++)
        {
            ObservedData[i] = (int *) calloc (programOptions.numSites +1, sizeof(int));
            if (!ObservedData[i])
            {
                fprintf (stderr, "Could not allocate the observed data structure\n");
                exit (-1);
            }
        }
        ReadFastaFile(fileName, ObservedData,  ObservedCellNames, &programOptions);
        
         for( i = 0 ; i < programOptions.numCells; i++)
             fprintf (stderr, "observed cell name %s\n", ObservedCellNames[i]);
         //fprintf (stderr, "observed data %d \n", *ObservedData[0]);
        
        msa = pll_phylip_load(fileNamePhylip, PLL_FALSE);
       
        int chainNumber=0;
        
        mcmcOptions.numChains=2;
        
        Chain * chains = (Chain *) malloc(mcmcOptions.numChains* ( sizeof(Chain) +  4 * sizeof(TreeNode) + 4 * sizeof(TreeNode) + programOptions.numClones *sizeof(Population)));
        if (!chains)
        {
            fprintf (stderr, "Could not allocate chains");
            exit (1);
        }
       
        InitializeChains(&chains, &programOptions, &mcmcOptions, sampleSizes, &(programOptions.seed), ObservedCellNames, msa);
        
        Chain currentChain;
        int currentIteration;
        int sampleEvery = mcmcOptions.thinning;
        for(chainNumber=0; chainNumber< mcmcOptions.numChains;chainNumber++)
        {
            for (currentIteration = 0; currentIteration < mcmcOptions.Niterations; currentIteration++)
            {
            runChain(&(chains[chainNumber]),    &mcmcOptions,  &(programOptions.seed),  &filePaths, &files, &programOptions,
                     varTimeGMRCA, ObservedCellNames, msa, sampleSizes);
            
            if (currentIteration % sampleEvery == 0 )
            {
              
                PrintTrees(currentIteration, &(chains[chainNumber].root), files.fpTrees, programOptions.mutationRate, programOptions.doUseObservedCellNames);
               PrintTrees2(currentIteration, &(chains[chainNumber].root), files.fpTrees2, programOptions.mutationRate,  ObservedCellNames, programOptions.doUseObservedCellNames);
      
              PrintTimes(currentIteration, files.fpTimes, programOptions.mutationRate, chains[chainNumber].nodes, programOptions.thereisOutgroup);
               PrintTimes2(currentIteration, files.fpTimes2, programOptions.mutationRate, chains[chainNumber].nodes, programOptions.thereisOutgroup);
                
            }
                
            }
            
        }
        pll_msa_destroy(msa);
        
    }//end inference
    
    free (allSites);
    allSites=NULL;
    free (SNVsites);
    SNVsites=NULL;
    free (variantSites);
    variantSites=NULL;
    free (DefaultModelSites);
    DefaultModelSites=NULL;
    free (AltModelSites);
    AltModelSites=NULL;
    //free(cellNames);
    //cellNames=NULL;
    //free(numberDifferences);
    free(SFS);
    SFS=NULL;
    free(populations);
    populations=NULL;
    
  
  //  PrintInfo(noisy, thereisOutgroup, outgroupSelection, outgroupBranchLength_Root1Root2, outgroupBranchLength_RootSample,  numClones, &TotalNumSequences, &TotalN, &TotalBirthRate, mutationRate, &TotalDeathRate, CloneNameBegin, CloneTimeOriginInput, CloneTimeOriginInputSTDPOP, CloneDelta, &expectedPopSize, ClonePopSizeBegin, CloneSampleSizeBegin, CloneBirthRateBegin, CloneDeathRateBegin);
    
   // CheckNumberOfClonesOrigin(numClones,  CloneTimeOriginInputSTDPOP, CloneTimeOriginInput );

  
    /* allocate cell structure */
    //        AllocateCellStructure( cell, numCells, numSites);
    
    
    /* free memory */

    meanNumCA  = cumNumCA /  programOptions.numDataSets;
    meanNumMIG  = cumNumMIG /  programOptions.numDataSets;
    meanNumMU  = cumNumMU /  (double) programOptions.numDataSets;
    fprintf(stderr, "Mean total number of mutations is  %lf \n", meanNumMU );
    if (noisy > 0)
        PrintRunSettings (programOptions, seed,  programOptions.TotalNumSequences,
                          CloneNameBegin,
                          CloneSampleSizeBegin,
                          ClonePopSizeBegin,
                          ClonePopSizeMeffectBegin,
                          CloneBirthRateBegin,
                          CloneDeathRateBegin,
                          CloneTimeOriginInput,
                          CloneGrowthRateBegin,
                          CloneTimeOriginInputSTDPOP,
                          CloneDelta,
                          varTimeGMRCA,
                          programOptions.Nscaling,
                          programOptions.numClones,
                          programOptions.numDataSets,
                          userSeed,
                          programOptions.mutationRate,
                          programOptions.noisy,
                          programOptions.outgroupSelection,
                          programOptions.thereisOutgroup,
                          programOptions.outgroupBranchLength_Root1Root2,
                          programOptions.outgroupBranchLength_RootSample,
                          &expectedPopSize,
                          cumNumCA,
                          meanNumCA,
                          cumNumMIG,
                          meanNumMIG,
                          &numEventsTot);
    //PrintRunSettings (originalSeed);    /* this function writes the value of all variables */
    free (varEvent);
    varEvent=NULL;
    free (varTimeGMRCA);
    varTimeGMRCA=NULL;
    //free(cellNames);
    //cellNames=NULL;
    for (i=0; i< numCells; i++)
    {
       // free(ObservedData[i]) ;
       // ObservedData[i]=NULL;
    }
    if (programOptions.doSimulateData ==NO)
    {
        free(ObservedData);
        ObservedData=NULL;
        
    }
    //free(cellNames);
    //cellNames=NULL;
    //if (doNGS == YES)
   // {
//        for (i=0; i<numCells; i++)
//        {
//            for (j=0; j<numSites; j++)
//            {
//                for (k=0; k<4; k++)
//                {
                  //  free(cell[i].site[j].genLike[k]);
                    //free(cell[i].site[j].scaledGenLike[k]);
                   // free(cell[i].site[j].genLikeDoublet[k]);
                   // free(cell[i].site[j].scaledGenLikeDoublet[k]);
//                }
                //free(cell[i].site[j].readCount);
                //free(cell[i].site[j].readCountDoublet);
                //free(cell[i].site[j].genLike);
                //free(cell[i].site[j].genLikeDoublet);
                //free(cell[i].site[j].scaledGenLike);
                //free(cell[i].site[j].scaledGenLikeDoublet);
//            }
//        }
//        free (cell);
   // }
//    for( i = 0 ; i < numClones; i++)
//          free(*(populations + i));
//     free(populations);
    free(CloneGrowthRateBegin);
    CloneGrowthRateBegin=NULL;
    free(ClonePopSizeMeffectBegin);
    ClonePopSizeMeffectBegin=NULL;
    free(CloneDelta);
    CloneDelta=NULL;
     free(CloneNameBegin);
    CloneNameBegin=NULL;
    free(CloneSampleSizeBegin);
    CloneSampleSizeBegin=NULL;
    free(ClonePopSizeBegin);
    ClonePopSizeBegin=NULL;
    free(CloneBirthRateBegin);
    CloneBirthRateBegin=NULL;
    free(CloneDeathRateBegin);
    CloneDeathRateBegin=NULL;
    free(CloneTimeOriginInput);
    CloneTimeOriginInput=NULL;
    free(proportionsVector);
    proportionsVector=NULL;
 
    
    
    if(programOptions.doPrintSeparateReplicates == NO){
      fprintf(stderr, "\n\nOutput files are in folder \"Results\":");
      if (programOptions.doPrintTrees == YES  )
      {
        fprintf(stderr, "\n Trees printed to files \"%s\"", filePaths.treeFile);
        fclose(files.fpTrees);
        fclose(files.fpTrees2);
      }
     if (programOptions.doPrintTimes == YES)
    {
        fprintf(stderr, "\n Times printed to files  \"%s\"", filePaths.timesFile);
        fclose(files.fpTimes);
        fclose(files.fpTimes2);
    }
    }
    fprintf(stderr, "\n\n*** Simulations finished ***");
    /* ejecution time */
    secs = (double)(clock() - start) / CLOCKS_PER_SEC;
    
    fprintf(stderr, "\n\n_________________________________________________________________");
    fprintf(stderr, "\nTime processing: %G seconds\n", secs);
    fprintf(stderr, "\nIf you need help type '-?' in the command line of the program\n");
    return 0;
    
}
/***************************** comparePopulationsByTimeOrigin*******************************/
int comparePopulationsByTimeOrigin(const void *s1, const void *s2)
{
    struct Population *p1 = *(struct Population **)s1;
    struct Population *p2 = *(struct Population **)s2;
    if (  p1->timeOriginInput  > p2 ->timeOriginInput)
        return 1;
    else if (p1->timeOriginInput < p2->timeOriginInput)
        return -1;
    else
        return 0;
}
/***************************** compare******************************/
int compare (const void * a, const void * b)
{
    double *p1 = (double *)a;
    double *p2 = (double *)b;
    
    if (*p1 > *p2) return 1;
    else if (*p2 > *p1 )return -1;
    else return 0;
}
/***************************** compare1******************************/
// int compare1 (const void * a, const void * b)
//{
//      double *p1 = (double *)a;
//      double *p2 = (double *)b;
//      fprintf (stderr, "compare %d with %d \n", *p1, *p2);
//  if (*p1 > *p2) return 1;
//  else if (*p2 > *p1) return -1;
//  else return 0;
//}

/***************************** comparePopulationsBySTDPopTime*******************************/
static int comparePopulationsBySTDPopTime(const void *s1, const void *s2)
{
    struct Population *p1 = *(struct Population **)s1;
    struct Population *p2 = *(struct Population **)s2;
    if (  p1->timeMigrationSTDCurrentPop  > p2 ->timeMigrationSTDCurrentPop)
        return 1;
    else if (p1->timeMigrationSTDCurrentPop < p2->timeMigrationSTDCurrentPop)
        return -1;
    else
        return 0;
}


/***************************** ListClonesAccordingTimeToOrigin2*******************************/
/* ListClonesAccordingTimeToOrigin2*/
void ListClonesAccordingTimeToOrigin2(Population **populations, int numClones) {
    
    
    qsort(populations, numClones, sizeof(Population* ), comparePopulationsByTimeOrigin);
}
/***************************** ListClonesAccordingTimeToOrigin*******************************/
/* ListClonesAccordingTimeToOrigin*/
void ListClonesAccordingTimeToOrigin(int numClones, int noisy,  int      **CloneNameBeginOrderByTimesOnlyControl, int **CloneNameBeginOrderByTimes, double *CloneTimeOriginInput) {
    int control1;
    int j, m, z;
    *CloneNameBeginOrderByTimes = (int *) calloc(numClones + 1, (long) sizeof(int));
    *CloneNameBeginOrderByTimesOnlyControl = (int *) calloc(numClones + 1, (long) sizeof(int));
    control1 = 0;
    
    for (j = 1; j <= numClones; j++)
    {
        *(*CloneNameBeginOrderByTimes + j) = 0;
        *(*CloneNameBeginOrderByTimesOnlyControl + j) = 0;
        
    }
    
    
    if (noisy > 2)
        fprintf(stderr, "\nClone number according to time to origin (t) (from more recent to more ancestral): ");
    
    
    for (m = 1; m <= numClones; m++) // to get a ranking by times to origin
    {
        for (z = 1; z <= numClones; z++)
        {
            if (m != z)
            {
                if (CloneTimeOriginInput[m] < CloneTimeOriginInput[z])
                    control1++; // number of clones with higher origin times with respect to clone "m"
            }
        }
        printf("\ncontrol1 = %d\n", control1);
        *(*CloneNameBeginOrderByTimesOnlyControl + m) = control1;
        control1 = 0;
    }
    
    for (j = 1; j <= numClones; j++)
        fprintf(stderr, "CloneNameBeginOrderByTimesOnlyControl[%d] = %d\n", j, *(*CloneNameBeginOrderByTimesOnlyControl + j) );
    
    control1 = 0;
    for (j = 1; j <= numClones; j++)
    {
        if (*(*CloneNameBeginOrderByTimesOnlyControl + j) == (numClones - 1)) // youngest clone
            *(*CloneNameBeginOrderByTimes + 1) = j;
        
        if (*(*CloneNameBeginOrderByTimesOnlyControl + j) == 0) // oldest clone
            *(*CloneNameBeginOrderByTimes + numClones) = j;
        
        if (*(*CloneNameBeginOrderByTimesOnlyControl + j) != (numClones - 1) && *(*CloneNameBeginOrderByTimesOnlyControl + j) != 0) // intermediate levels
        {
            for (z = 1; z <= numClones; z++)
            {
                if (*(*CloneNameBeginOrderByTimesOnlyControl + j) == z)
                    *(*CloneNameBeginOrderByTimes + numClones - z) = j;
            }
        }
    }
    
    if (noisy > 2)
    {
        for (j = 1; j <= numClones; j++)
        {
            fprintf(stderr, "%d ",   *(*CloneNameBeginOrderByTimes + j));
        }
        fprintf(stderr, "\n ");
    }
    
    
}
/***************************** CreateFolder*******************************/
/* CreateFolder*/
void CreateFolder(FilePaths filePaths, ProgramOptions programOptions,Files *files//,
//                  FILE  **fpTrees, FILE **fpTimes, FILE **fpTrees2, FILE **fpTimes2,
//                  FILE            **fpSNVgenotypes,   FILE **fpSNVhaplotypes, FILE **fpTrueHaplotypes,
//                  FILE **fpFullGenotypes,FILE **fpFullHaplotypes, FILE **fpVCF, FILE **fpCATG
                  ) {
    char File[MAX_NAME];

    /* Create a Results folder and prepare output files */
    mkdir("Results", S_IRWXU); /* Create "Results" folder (with type S_IRWXU (read, write and execute)) */
    //mkdir("Results",0);
#ifdef MAC
    strcpy (filePaths.dir, ":Results:"); /* Copy the string in char variable dir = Results (char), is different mac vs windows */
#else
    strcpy (filePaths.dir, "Results/");
#endif
    
    if (programOptions.doPrintTrees == YES)    /* if treeFile " " */
    {
        sprintf(File, "%s%s", filePaths.dir, filePaths.treeFile);
       // if ((*fpTrees = fopen(File, "w")) == NULL)
        if (openFile(&files->fpTrees, File) == -1)
        {
            fprintf(stderr, "Can't open %s.\n", File);
            exit(-1);
        }
        
        sprintf(File, "%s%s_2.tre", filePaths.dir, filePaths.treeFile);
        //if ((*fpTrees2 = fopen(File, "w")) == NULL)
       if (openFile(&files->fpTrees2, File) == -1)
        {
            fprintf(stderr, "Can't open %s.\n", File);
            exit(-1);
        }
    }
    if (programOptions.doPrintTimes == YES)    /* if timesFile " " */
    {
        sprintf(File, "%s%s", filePaths.dir, filePaths.timesFile);
       // if ((*fpTimes = fopen(File, "w")) == NULL)
        if (openFile(&files->fpTimes, File) == -1)
        {
            fprintf(stderr, "Can't open %s.\n", File);
            exit(-1);
        }
       
        sprintf(File, "%s%s_2.txt", filePaths.dir, filePaths.timesFile);
      //  if ((*fpTimes2 = fopen(File, "w")) == NULL)
        if (openFile(&files->fpTimes2, File) == -1)
        {
            fprintf(stderr, "Can't open %s.\n", File);
            exit(-1);
        }
    }
    if (programOptions.doSimulateData == YES)
    {
        /* contains SNV genotypes for every cell */
        if (programOptions.doPrintSNVgenotypes == YES)
        {
            sprintf(File,"%s/%s", filePaths.dir, filePaths.SNVgenotypesFile);
           // if ((*fpSNVgenotypes = fopen(File, "w")) == NULL)
            if (openFile(&files->fpTimes2, File) == -1)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        
        /* contains haplotypes for variable sites for every cell */
        if (programOptions.doPrintSNVhaplotypes == YES)
        {
            sprintf(File,"%s/%s", filePaths.dir, filePaths.SNVhaplotypesFile);
            //if ((*fpSNVhaplotypes = fopen(File, "w")) == NULL)
             if (openFile(&files->fpSNVhaplotypes, File) == -1)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        
        /* contains error-free haplotypes for variable sites for every cell */
        if (programOptions.doPrintTrueHaplotypes == YES)
        {
            sprintf(File,"%s/%s", filePaths.dir, filePaths.trueHaplotypesFile);
            //if ((*fpTrueHaplotypes = fopen(File, "w")) == NULL)
            if (openFile(&files->fpTrueHaplotypes, File) == -1)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        
        /* contains all genotypes (variable or invariable) for every cell */
        if (programOptions.doPrintFullGenotypes == YES)
        {
            sprintf(File,"%s/%s", filePaths.dir, filePaths.fullGenotypesFile);
            //if ((*fpFullGenotypes = fopen(File, "w")) == NULL)
             if (openFile(&files->fpFullGenotypes, File) == -1)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        
        /* contains haplotypes for all sites for every cell */
        if (programOptions.doPrintFullHaplotypes == YES)
        {
            sprintf(File,"%s/%s", filePaths.dir, filePaths.fullHaplotypesFile);
           // if ((*fpFullHaplotypes = fopen(File, "w")) == NULL)
               if (openFile(&files->fpFullHaplotypes, File) == -1)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        
        /* contains reads counts and genotype likelihoods for every SNV and cell */
        if (programOptions.doNGS == YES)
        {
            sprintf(File,"%s/%s", filePaths.dir, filePaths.VCFfile);
            //if ((*fpVCF = fopen(File, "w")) == NULL)
            if (openFile(&files->fpVCF, File) == -1)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        
        /* contains reads counts for every SNV and cell */
        if (programOptions.doPrintCATG == YES)
        {
            sprintf(File,"%s/%s", filePaths.dir, filePaths.CATGfile);
            //if ((*fpCATG = fopen(File, "w")) == NULL)
             if (openFile(&files->fpCATG, File) == -1)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        
    }
    
    if (programOptions.doPrintSNVgenotypes == YES)
    {

        
//        //fprintf (*fpSNVgenotypes, "%s - ",PROGRAM_NAME);
//        //PrintDate (fpSNVgenotypes);
//        fprintf (&files->fpSNVgenotypes, "SNV genotypes\n");
//        //PrintCommandLine (fpSNVgenotypes, argc, argv);
//        fprintf (&files->fpSNVgenotypes,"\n%d\n", programOptions.numDataSets);
    }
    
    if (programOptions.doPrintSNVhaplotypes == YES)
    {
 
        
//       // fprintf (*fpSNVhaplotypes, "%s - ",PROGRAM_NAME);
//        //PrintDate (fpSNVhaplotypes);
//        fprintf (&files->fpSNVhaplotypes, "SNV haplotypes\n");
//        //PrintCommandLine (fpSNVhaplotypes, argc, argv);
//        fprintf (&files->fpSNVhaplotypes,"\n%d\n", programOptions.numDataSets);
    }
    
    if (programOptions.doPrintTrueHaplotypes == YES)
    {
   
        
//        fprintf (*fpTrueHaplotypes, "%s - ",PROGRAM_NAME);
//        //PrintDate (fpTrueHaplotypes);
//        fprintf (*fpTrueHaplotypes, "True haplotypes\n");
//        //PrintCommandLine (fpTrueHaplotypes, argc, argv);
//        fprintf (*fpTrueHaplotypes,"\n%d\n", programOptions.numDataSets);
    }
    
    if (programOptions.doPrintFullGenotypes == YES)
    {
 
//       // fprintf (*fpFullGenotypes, "%s - ",PROGRAM_NAME);
//        //PrintDate (fpFullGenotypes);
//        fprintf (&files->fpFullGenotypes, "Full genotypes\n");
//        //PrintCommandLine (fpFullGenotypes, argc, argv);
//        fprintf (&files->fpFullGenotypes,"\n%d\n", programOptions.numDataSets);
    }
    
    if (programOptions.doPrintFullHaplotypes == YES)
    {

        
//        fprintf (&files->fpFullHaplotypes, "%s - ",PROGRAM_NAME);
//        //PrintDate (fpFullHaplotypes);
//        fprintf (&files->fpFullHaplotypes, "Full haplotypes\n");
//        //PrintCommandLine (fpFullHaplotypes, argc, argv);
//        fprintf (&files->fpFullHaplotypes,"\n%d\n", programOptions.numDataSets);
    }
    
    if (programOptions.doNGS == YES)
    {
  
        
//        fprintf (&files->fpVCF, "%s - ",PROGRAM_NAME);
//        //PrintDate (fpVCF);
//        fprintf (&files->fpVCF, "Read counts and genotype likelihoods\n");
//        //PrintCommandLine (fpVCF, argc, argv);
//        fprintf (&files->fpVCF,"\n%d\n", programOptions.numDataSets);
    }
    
    if (programOptions.doPrintCATG == YES)
    {

        
//        fprintf (&files->fpCATG, "%s - ",PROGRAM_NAME);
//        //PrintDate (fpCATG);
//        fprintf (&files->fpCATG, "Read counts\n");
//        //PrintCommandLine (fpCATG, argc, argv);
//        fprintf (&files->fpCATG,"\n%d\n", programOptions.numDataSets);
    }
    
    
    
    
}
/***************************** openFile*******************************/
/* openFile*/
int openFile(FILE **file, char path[MAX_NAME] )
{
    // if ((*fpTrees = fopen(File, "w")) == NULL)
    if ((*file = fopen(path, "w")) == NULL)
    {
        fprintf(stderr, "Can't open %s.\n", path);
        return -1;
    }
    return 0;
}
/***************************** PrintInfo*******************************/
/* PrintInfo*/
void PrintInfo(int noisy, int thereisOutgroup, int outgroupSelection, double outgroupBranchLength_Root1Root2, double outgroupBranchLength_RootSample,  int numClones, int *TotalNumSequences, int *TotalN, double  *TotalBirthRate, double mutationRate, double *TotalDeathRate, int *CloneNameBegin, double *CloneTimeOriginInput, double *CloneTimeOriginInputSTDPOP, double *CloneDelta, double *expectedPopSize, int *ClonePopSizeBegin, int *CloneSampleSizeBegin, double *CloneBirthRateBegin, double *CloneDeathRateBegin  ) {
    int j;
    // print some info
    if (noisy > 1)
        printf("\n\tClone number\t|\tTime to origin\t|\tTime Origin Coal (model units)\t|\tDelta (gr * eps)\n");
    for (j = 0; j < numClones; j++)
    {
        if (noisy > 1)
        {
            printf("\t%d\t\t", CloneNameBegin[j]);
            printf("\t\t\t(-)%lf\t\t", CloneTimeOriginInput[j]);
           // printf("%lf\t\t", CloneTimeOriginInputSTDPOP[j]);
            //printf("\t\t\t\t\t%lf\t\t\n", CloneDelta[j]);
        }
    }
    
    /* expectation */ // The expected number of cells in a birth death process is exp( (birth-death)*(time of origin) )
    if (noisy > 1)
        printf("\nExpected population size per clone");
    for (j = 0; j < numClones; j++)
    {
        *expectedPopSize = 0.0;
      //  *expectedPopSize = exp (CloneTimeOriginInputSTDPOP[j] * CloneDelta[j]);
        //printf("\nClone number %d: Expected initial population size = %lf (exp(%lf)), input initial population size = %d", j, expectedPopSize, CloneTimeOriginInputSTDPOP[j] * CloneDelta[j], ClonePopSizeBegin[j]);
        if (noisy > 1)
            printf("\nClone %d\t-\tExpected initial population size = %lf\t-\tInput initial population size = %d", j, *expectedPopSize, ClonePopSizeBegin[j]);
    }
    if (noisy > 3)
    {
        printf("\n\nTotal number of sequences (all clones) = %d\n", *TotalNumSequences);
        printf("Total population size (all clones) = %d\n", *TotalN);
        printf("Total birth rate (all clones) = %lf\n", *TotalBirthRate);
        printf("Total death rate (all clones) = %lf\n\n", *TotalDeathRate);
    }
    /* Print other parameters */
    if (noisy > 1)
    {
        printf("Mutation rate = %lf\n", mutationRate);
        if (thereisOutgroup == NO)
            printf("Outgroup is not simulated\n");
        else
        {
            if (outgroupSelection == 2)
                printf("Outgroup branch length (Root1 - Root2) = %lf\n\n", outgroupBranchLength_Root1Root2);
            
            printf("Outgroup branch length (Root - sample) = %lf\n", outgroupBranchLength_RootSample);
        }
    }
}

/***************************** InitListClonesTimes*******************************/
/* InitListClonesTimes*/
void InitListClonesTimes(Population **populations, int numClones,  int *doEstimateTimesOriginClones,
                         double *CloneTimeOriginInput
                         ) {
    int z, j;
    z = 0;
    Population *p;
    for (j = 0; j <= (numClones - 1); j++)
    {
        p = *(populations + j);
        p->timeOriginInput = CloneTimeOriginInput[j ];
        p->timeOriginSTD = CloneTimeOriginInput[j ] / p->effectPopSize;
        p->timeMigrationSTDCurrentPop = p->timeOriginSTD;
        p->delta = (double)p->growthRate * p->effectPopSize;
        if ( p->timeOriginSTD == 0) {
            z++;
            p->doEstimateTimeOrigin = 1;
        }
    }
    // active estimation of times
    if (z == numClones)
        *doEstimateTimesOriginClones = YES;
    for (j = 0; j < numClones; j++)
    {
        p = *(populations + j);
        if ( p->delta <= 0)
        {
            fprintf (stderr, "PARAMETER ERROR: The growth rate cannot be lower than the death rate(Delta parameter negative, Delta=(%10.9lf)) for population %d\n\n", p->delta, j);
            PrintUsage();
        }
        if (p->timeOriginInput == 0 && *doEstimateTimesOriginClones == NO)
        {
            fprintf (stderr, "PARAMETER ERROR: Bad time to origin for clone %d (should not be 0; excepting estimation of times where all clones must have a time of 0) (%lf)\n\n", j, p->timeOriginInput);
            PrintUsage();
        }
    }
}
/***************************** InitListClones*******************************/
/* InitListClones*/
void InitListClones(Population **populations, int numClones, int noisy, int *CloneNameBegin, int *CloneSampleSizeBegin, double *CloneBirthRateBegin,  double *CloneDeathRateBegin,  int *ClonePopSizeBegin, int TotalNumSequences  ) {
    int z;
    //struct Population* pops = malloc(numClones * (sizeof(struct Population)+TotalNumSequences * sizeof( int) + numClones * sizeof(double) ));
    struct Population* pops = malloc(numClones * (sizeof(struct Population) ));
    for (z = 0; z <= (numClones - 1); z++)
    {
        //CloneGrowthRateBegin[z ] = 0.0;
        // ClonePopSizeMeffectBegin[z ] = 0.0;
       // CloneGrowthRateBegin[z ] = CloneBirthRateBegin[z ] - CloneDeathRateBegin[z ];
       // ClonePopSizeMeffectBegin[z] = ClonePopSizeBegin[z ] / CloneBirthRateBegin[z ];
        pops[z].index = CloneNameBegin[z ];
        pops[z].order=0;
        pops[z].birthRate = CloneBirthRateBegin[z ];
        pops[z].deathRate = CloneDeathRateBegin[z ];
        pops[z].growthRate = CloneBirthRateBegin[z ] - CloneDeathRateBegin[z ];
        pops[z].sampleSize = CloneSampleSizeBegin[z];
        pops[z].popSize = ClonePopSizeBegin[z ];
        pops[z].effectPopSize =  ClonePopSizeBegin[z ] / CloneBirthRateBegin[z ];
        pops[z].delta = pops[z].growthRate * pops[z].effectPopSize ;
        pops[z].numActiveGametes = CloneSampleSizeBegin[z];
        pops[z].isAlive = 1;
        pops[z].nodeIdAncesterMRCA = 0;
        pops[z].numCompletedCoalescences = 0;
        pops[z].nextAvailableIdInmigrant = 0;
        pops[z].numIncomingMigrations = 0;
         pops[z].numPossibleMigrations = 0;
        pops[z].doEstimateTimeOrigin = NO;
        *(populations + z) = &pops[z];
        if (noisy > 1)
            printf("\t%d\t\t", CloneNameBegin[z ]);
        if ((z + 1) != CloneNameBegin[z ])
        {
            fprintf (stderr, "PARAMETER ERROR: Check order of clones. Clone (%d) in order is different to (%d). (d)\n\n", z, CloneNameBegin[z]);
            PrintUsage();
        }
        if (noisy > 1)
        {
            printf("\t\t\t%d\t\t\t", CloneSampleSizeBegin[z ]);
            printf("\t%d\t\t\t", ClonePopSizeBegin[z ]);
            printf("\t%lf\t\t",  pops[z].effectPopSize);
            printf("\t\t\t%lf\t", CloneBirthRateBegin[z ]);
            printf("\t%lf\t", CloneDeathRateBegin[z ]);
            printf("\t%lf\t\t",pops[z].growthRate);
            //printf("\t%d\t\t", z);
            printf("\n");
        }
    }
}

/***************************** CheckNumberOfClonesOrigin*******************************/
/* CheckNumberofClonesOrigin*/
void CheckNumberOfClonesOrigin(int numClones, double  *CloneTimeOriginInputSTDPOP, double *CloneTimeOriginInput ) {
    int z, j;
    for (z = 0; z < numClones; z++)
    {
        for (j = 0; j < numClones; j++)
        {
            if (j != z && CloneTimeOriginInputSTDPOP[z] == CloneTimeOriginInputSTDPOP[j])
            {
                fprintf (stderr, "PARAMETER ERROR: Clone %d and clone %d present same time of origin (%lf). This program assumes different time origins for different clones. Modify population size, birth rate or death rate for any of these clones.\n", j, z, CloneTimeOriginInputSTDPOP[z]);
                PrintUsage();
            }
            if (j != z && CloneTimeOriginInput[z] == CloneTimeOriginInput[j])
            {
                fprintf (stderr, "PARAMETER ERROR: Clone %d and clone %d present same time of origin (%lf). This program assumes different time origins for different clones. Modify time of origin for any of these clones.\n", j, z, CloneTimeOriginInput[z]);
                PrintUsage();
            }
        }
    }
}


/***************************** GenerateTimesFromPriorsOriginal*******************************/
/* GenerateTimesFromPriorsOriginal */
void GenerateTimesFromPriorsOriginal(int noisy, int numClones,double *proportionsVector, Population **populations,  long int *seed) {
    double    *Uvector;
    int i, j, z, m,  l;
    double    TotalProbability, LocalProbability, aboveTerm, belowTerm, ranHere;
    int  AttemptsAcceptation = 0;
    TotalProbability = 0.0;
    LocalProbability = 0.0;
    aboveTerm = 0.0;
    belowTerm = 0.0;
    ranHere = 0.0;
    if (noisy > 1)
        printf("Estimation of times of origin of clones ..\n");
    Uvector =  malloc((numClones )* (long) sizeof(double));
    if (!Uvector)
    {
        fprintf (stderr, "Could not allocate Xvector (%lu bytes)\n", (numClones ) * (long) sizeof(double));
        exit (-1);
    }
    int doAcceptTimes = NO;
    int doEstimateTimesOriginClones=YES;
    Population *popJ, *popI, *popL;
    double rand;
    if (doEstimateTimesOriginClones == YES)
    {
        while (doAcceptTimes == NO)
        {
            AttemptsAcceptation++;
            // Calculate t and T
            if (noisy > 2)
                printf("\nProposed times of origin of clones ..");
            for (j = 0; j < numClones; j++)
            {
                popJ=*(populations + j);
                Uvector[j] = 0;
                Uvector[j] = RandomUniform (seed);
                rand = RandomUniform(seed);
                popJ->timeOriginSTD= (1/ popJ->delta)*log(1-(popJ->delta / log(rand)));
                popJ->timeMigrationSTDCurrentPop= popJ->timeOriginSTD ;
               // popJ->timeOriginInput=proportionsVector[j]*popJ->timeOriginSTD;
               popJ->timeOriginInput=popJ->effectPopSize*popJ->timeOriginSTD;
            }
            
            if (numClones == 1) {
                doAcceptTimes = YES;
                break;
            }
            ListClonesAccordingTimeToOrigin2(populations, numClones);
            // Calculate probabilities P
            LocalProbability = 0.0;
            TotalProbability = 0.0;
           // printf("Number of clones = %d\n", numClones);
            for (i = 0; i < numClones - 1; i++)
            {
                popI=*(populations + i);
                printf("\ni  = %d ", i);
                LocalProbability = 0.0;
                aboveTerm = 0.0;
                belowTerm = 0.0;
                m = popI->order;
                // printf("\n\nCalculating P for clone = %d ", m);
                for (l = i + 1; l < numClones; l++)
                {
                    popL=*(populations + l);
                    j = popL->order ;
                    aboveTerm = aboveTerm + (popL->popSize * CalculateH(popI->timeOriginSTD * popI->effectPopSize / popL->effectPopSize, popL->timeOriginSTD, popL->delta));
                    belowTerm = belowTerm + popL->popSize;
                    // printf("\nClones %d and %d ", m, j);
                }
                LocalProbability = aboveTerm / belowTerm;
                //printf("\naboveTerm = %lf    belowTerm = %lf", aboveTerm, belowTerm);
                if (i == 1)
                    TotalProbability = 1.0 * LocalProbability;
                else
                    TotalProbability = TotalProbability * LocalProbability;
                // printf("\nLocalProbability = %lf    TotalProbability = %lf", LocalProbability, TotalProbability);
            }
           // printf("\nTotalProbability = %lf \n", TotalProbability);
            // Random number
            ranHere = RandomUniform (seed);
            // printf("\n\nranHere = %lf \n", ranHere);
            // free memory
            //free (CloneNameBeginOrderByTimes);
            // free (CloneNameBeginOrderByTimesOnlyControl);
            // check accept or reject attempt
            if (ranHere <= TotalProbability)
                doAcceptTimes = YES;
            if (noisy > 2)
                printf("\nProbability = %lf (random number = %lf) [#attempt = %d]\n", TotalProbability, ranHere, AttemptsAcceptation);
        }
    }
    if (noisy > 2)
    {
        printf("\nTimes accepted .. \n");
    printf("\n\nDONE! ranHere = %lf / TotalProbability = %lf [total attempts = %d]\n", ranHere, TotalProbability, AttemptsAcceptation);
    }
}
/***************************** GenerateTimesFromPriors*******************************/
/* GenerateTimesFromPriors*/
void GenerateTimesFromPriors(int noisy, int numClones,Population **populations,  long int seed){
    int i, j,  l;
    Population *popI, *popL, *popIPlus1;
    double W, V,P, aboveTerm, belowTerm;
    double sum;
    if (noisy > 1){
        fprintf (stderr, "Estimation of times of origin of clones ..\n");
        fflush (stdout);
    }
    double    *Uvector = (double *) calloc((numClones), (long) sizeof(double));
    if (!Uvector)
    {
        fprintf (stderr, "Could not allocate Xvector (%lu bytes)\n", (numClones + 1) * (long) sizeof(double));
        exit (-1);
    }
    for (j = 0; j < numClones; j++)
    {
        popI=*(populations + j);
        Uvector[j]= (-1/ popI->delta)*log((-1/2)*(popI->delta)+ sqrt(1+0.25*(popI->delta)*(popI->delta)));
    }
    Population *lastPop;
    int iteracion=0;
    do{
        iteracion= iteracion+1;
        //fprintf (stderr, "\n\n iteration = %d \n", iteracion);
        //fflush (stdout);
        // printf("\n\n iteracion = %d \n", iteracion);
        lastPop=*(populations + numClones-1);
        V= RandomUniform (&seed);
        lastPop->timeOriginSTD = (1 / lastPop->delta) * log(1.0 - (lastPop->delta / log(V)));
        lastPop->timeMigrationSTDCurrentPop =lastPop->timeOriginSTD ;
        lastPop->timeOriginInput = lastPop->timeOriginSTD * lastPop->effectPopSize;
         lastPop->timeMigrationSTDCurrentPop =lastPop->timeOriginSTD ;
        for (i = numClones-2; i >= 0; i--)
        {
            do{
                popIPlus1=*(populations + i +1);
                W=RandomUniform (&seed) * popIPlus1->timeOriginInput;
                V=RandomUniform (&seed);
                sum=0.0;
                aboveTerm=0.0;
                belowTerm =0.0;
                for (l = i + 1; l < numClones; l++)
                {
                    popL=*(populations + l);
                    aboveTerm = aboveTerm + (popL->popSize * CalculateH(W / popL->effectPopSize, popL->timeOriginSTD, popL->delta));
                    belowTerm = belowTerm + popL->popSize;
                }
                P = aboveTerm / belowTerm;
            }
            while(V > P);
            //end do-while inside for
            popI=*(populations + i );
            popI->timeOriginSTD= W / popI->effectPopSize;
            popI->timeOriginInput =W;
        }//end for
        // calculate P
        aboveTerm=1;
        belowTerm =1;
        for (l = 0; l < numClones -1; l++)
        {
            popL=*(populations + l);
            aboveTerm = aboveTerm * DensityTime(popL->delta, popL->timeOriginSTD);
            belowTerm = belowTerm * DensityTime(popL->delta, Uvector[l]);
        }
        P = aboveTerm / belowTerm;
        V=RandomUniform (&seed);
    }
    while(V>P);
    fprintf (stderr, "\n\n total iterations = %d \n", iteracion);
    fflush (stdout);
    free(Uvector);
    Uvector=NULL;
}

/***************************** GenerateTimesFromPriors2*******************************/
/* GenerateTimesFromPriors2*/
void GenerateTimesFromPriors2(int noisy, int numClones,Population **populations,  long int seed){
    int i, j, l;
    Population *popI, *popL, *popIPlus1;
    double W, V,P, aboveTerm, belowTerm, temp;
    double sum;
    if (noisy > 1){
        fprintf (stderr, "Estimation of times of origin of clones ..\n");
        fflush (stdout);
    }
    double    *Uvector = (double *) calloc((numClones), (long) sizeof(double));
    if (!Uvector)
    {
        fprintf (stderr, "Could not allocate Xvector (%lu bytes)\n", (numClones ) * (long) sizeof(double));
        exit (-1);
    }
    
    for (j = 0; j < numClones; j++)
    {
        popI=*(populations + j);
        Uvector[j]= (-1/ popI->delta)*log((-1/2)*popI->delta+ sqrt(1+0.25*popI->delta*popI->delta));
    }
    Population *lastPop;
    int iteracion=0;
   
    
        lastPop=*(populations + numClones-1);
        V= RandomUniform (&seed);
        lastPop->timeOriginSTD = (1 / lastPop->delta) * log(1.0 - (lastPop->delta / log(V)));
        lastPop->timeMigrationSTDCurrentPop =lastPop->timeOriginSTD ;
        lastPop->timeOriginInput = lastPop->timeOriginSTD * lastPop->effectPopSize;
         lastPop->timeMigrationSTDCurrentPop =lastPop->timeOriginSTD ;
        for (i = numClones-2; i >= 0; i--)
        {
            do{
                iteracion= iteracion+1;
                fprintf (stderr, "\n\n iteration = %d \n", iteracion);
                //fflush (stdout);
                // printf("\n\n iteracion = %d \n", iteracion);
                popIPlus1=*(populations + i +1);
                W=RandomUniform (&seed) * popIPlus1->timeOriginInput;
                V=RandomUniform (&seed);
                sum=0.0;
                aboveTerm=0.0;
                belowTerm =0.0;
                popI=*(populations + i );
                popI->timeOriginSTD= W / popI->effectPopSize;
                for (l = i + 1; l < numClones; l++)
                {
                    popL=*(populations + l);
                    aboveTerm = aboveTerm + (popL->popSize * CalculateH(W / popL->effectPopSize, popL->timeOriginSTD, popL->delta));
                    belowTerm = belowTerm + popL->popSize;
                    
                }
                 P = (aboveTerm / belowTerm) ;
                temp = DensityTime(popI->delta, popI->timeOriginSTD)/ DensityTime(popI->delta, Uvector[i]);
                P = (aboveTerm / belowTerm) * temp ;
            }
            while(V >P);
            
            popI->timeOriginInput =W;
        }//end for
    fprintf (stderr, "\n\n iterations = %d \n", iteracion);
    fflush (stdout);
    free(Uvector);
    Uvector=NULL;
}

/***************************** GenerateTimesFromPriors3*******************************/
/* GenerateTimesFromPriors3*/
void GenerateTimesFromPriors3(int noisy, int numClones,Population **populations,  long int seed){
    int i, j,  l;
    Population *popI, *popL, *popIPlus1;
    double W, V,P, aboveTerm, belowTerm, temp;
    double sum;
    if (noisy > 1){
        fprintf (stderr, "Estimation of times of origin of clones ..\n");
        fflush (stdout);
    }
    double    *Uvector = (double *) calloc((numClones), (long) sizeof(double));
    if (!Uvector)
    {
        fprintf (stderr, "Could not allocate Xvector (%lu bytes)\n", (numClones ) * (long) sizeof(double));
        exit (-1);
    }
    
    for (j = 0; j < numClones; j++)
    {
        popI=*(populations + j);
        Uvector[j]= (-1/ popI->delta)*log((-1/2)*popI->delta+ sqrt(1+0.25*popI->delta*popI->delta));
    }
    Population *lastPop;
    int iteracion=0;
    lastPop=*(populations + numClones-1);
    V= RandomUniform (&seed);
    lastPop->timeOriginSTD = (1 / lastPop->delta) * log(1.0 - (lastPop->delta / log(V)));
    lastPop->timeMigrationSTDCurrentPop =lastPop->timeOriginSTD ;
    lastPop->timeOriginInput = lastPop->timeOriginSTD * lastPop->effectPopSize;
    lastPop->timeMigrationSTDCurrentPop =lastPop->timeOriginSTD ;
    for (i = numClones-2; i >= 0; i--)
    {
        do{
            iteracion= iteracion+1;
            fprintf (stderr, "\n\n iteration = %d \n", iteracion);
            //fflush (stdout);
            // printf("\n\n iteracion = %d \n", iteracion);
            popIPlus1=*(populations + i +1);
            W=RandomUniform (&seed) * popIPlus1->timeOriginInput;
            V=RandomUniform (&seed);
            sum=0.0;
            aboveTerm=0.0;
            belowTerm =0.0;
            popI=*(populations + i );
            popI->timeOriginSTD= W / popI->effectPopSize;
            for (l = i + 1; l < numClones; l++)
            {
                popL=*(populations + l);
                aboveTerm = aboveTerm + (popL->popSize * CalculateH(W / popL->effectPopSize, popL->timeOriginSTD, popL->delta));
                belowTerm = belowTerm + popL->popSize;
            }
            P = (aboveTerm / belowTerm) ;
            temp = DensityTime(popI->delta, popI->timeOriginSTD)/ DensityTime(popI->delta, Uvector[i]);
            P = (aboveTerm / belowTerm)  ;//here we dont multiply by temp like a previous version
        }
        while(V >P);
        popI->timeOriginInput =W;
    }//end for
    fprintf (stderr, "\n\n iterations = %d \n", iteracion);
    fflush (stdout);
    free(Uvector);
    Uvector=NULL;
}
/***************************** GenerateTimesFromPriors3*******************************/
/* GenerateTimesFromPriors3*/
void GenerateTimesFromPriors4(int noisy, int numClones,double *proportionsVector, Population **populations,  long int seed){
    int j, i, k, l;
    Population* popJ, *popI, *popl, *fatherPop;
    double rand, probIcomeFroml;
    for (j = 0; j < numClones; j++)
    {
        popJ=*(populations + j);
        rand = RandomUniform(&seed);
        popJ->timeOriginSTD= (1/ popJ->delta)*log(1-(popJ->delta / log(rand)));
         popJ->timeMigrationSTDCurrentPop= popJ->timeOriginSTD ;
        popJ->timeOriginInput=proportionsVector[j]*popJ->timeOriginSTD;
    }
     ListClonesAccordingTimeToOrigin2(populations, numClones);
       InitListPossibleMigrations(populations, numClones);
//    for (i = 0; i < numClones -1; i++)
//    {   popI =*(populations + i);
//        fatherPop=  ChooseFatherPopulation(populations, numClones, popI, &seed, noisy);
//           UpdateListMigrants(populations, numClones, popI, fatherPop);
//
//    }
}
/***************************** GeneratePopSizesFromPriors*******************************/
/* GeneratePopSizesFromPriors*/
void GeneratePopSizesFromPriors(int noisy, int numClones,Population **populations,  long int *seed, double gammaParam, int totalPopSize){
    int i, j;
    Population *popI, *popJ, *popL, *popIPlus1;
    double W, V,P, aboveTerm, belowTerm, temp;

    double rand;
    if (noisy > 1){
        fprintf (stderr, "Estimation of pop sizes of clones ..\n");
        fflush (stdout);
    }
    double    *Uvector = (double *) calloc((numClones), (long) sizeof(double));
    if (!Uvector)
    {
        fprintf (stderr, "Could not allocate Xvector (%lu bytes)\n", (numClones ) * (long) sizeof(double));
        exit (-1);
    }
    double    *cumSum = (double *) calloc((numClones +1), (long) sizeof(double));
    if (!cumSum)
    {
        fprintf (stderr, "Could not allocate Xvector (%lu bytes)\n", (numClones +1 ) * (long) sizeof(double));
        exit (-1);
    }
    RandomDirichlet(gammaParam,  numClones, &Uvector, seed );
    cumSum[0]=0.0;
    for (j = 1; j <= numClones; j++)
    {     cumSum[j]=0;
          cumSum[j]=cumSum[j-1]+Uvector[j-1];
          popJ = *(populations + j - 1);
          popJ->effectPopSize=0;
          popJ->popSize=0;
    }
    for (i = 0; i < totalPopSize; i++)
    {
        rand=RandomUniform(seed);
        for (j = 1; j <= numClones;j ++)
        {
            popJ = *(populations + j - 1);
            if (rand <= cumSum[j] && rand > cumSum[j - 1])
            {
                popJ->effectPopSize=popJ->effectPopSize +1;
                break;
            }
       }
    }
    for (j = 0; j < numClones; j++)
    {
        popJ = *(populations + j);
         popJ->growthRate =popJ->delta * popJ->birthRate /popJ->effectPopSize;
        popJ->popSize=popJ->effectPopSize * popJ->growthRate;
       
    }
    free(Uvector);
    Uvector=NULL;
    free(cumSum);
    cumSum=NULL;
}
/***************************** GenerateEffectPopSizesFromPriors2*******************************/
/* GenerateEffectPopSizesFromPriors2 after generating the times of origin*/
void GenerateEffectPopSizesFromPriors2(Chain *chain, int noisy, int numClones,Population **populations,  long int *seed, double gammaParam, int totalPopSize,  int doGenerateProportionsVector){
    int i, j, k, l;
    Population *popI, *popJ, *popL, *popIPlus1;
    double W, V,P, aboveTerm, belowTerm, temp;
    double sum;
    double rand;
  
    if (noisy > 1){
        fprintf (stderr, "Estimation of pop sizes of clones ..\n");
        fflush (stdout);
    }
    double *outputVector;
    if (doGenerateProportionsVector == YES)
    {
        for (j = 0; j < numClones; j++)
        {
            chain->oldproportionsVector[j]=chain->proportionsVector[j] ;
        }

        RandomDirichlet(gammaParam,  numClones, &(chain->proportionsVector), seed );
        //chain->proportionsVector = outputVector;
        
    }
    //oldest population
//    popJ = *(populations + numClones-1);
//    popJ->effectPopSize= oldestPopSize;
//    popJ->growthRate =popJ->delta * popJ->birthRate /popJ->effectPopSize;
//    popJ->popSize=popJ->effectPopSize * popJ->growthRate;
    if (chain->proportionsVector == NULL)
        fprintf (stderr, "ERROR: the proportions vector is null..\n");
    for (j = 0; j < numClones; j++)
    {
        popJ = *(populations + j);
        //temp=*(outputVector + i);
        popJ->oldeffectPopSize = popJ->effectPopSize;
        popJ->effectPopSize = chain->proportionsVector[j] * totalPopSize;
        //popJ->effectPopSize = temp * totalPopSize;
        //  popJ->growthRate =popJ->delta * popJ->birthRate /popJ->effectPopSize;
        //  popJ->popSize=popJ->effectPopSize * popJ->growthRate;
    }
}
/***************************** DensityTime*******************************/
/* DensityTime*/
double DensityTime(double delta, double u){
    double term1=delta * exp(-1*delta*u);
    double term2=1-exp(-1*delta*u);
    return delta * term1 * exp(-1*term1/term2) /(term2 * term2);
}
/***************************** ValidateParameters*******************************/
/* Validate parameters*/

void ValidateParameters(ProgramOptions *programOptions,
                        int *CloneNameBegin , int *CloneSampleSizeBegin, int *ClonePopSizeBegin) {
    if (programOptions->noisy > 1)
    {
        printf("\n>> Settings ..\n");
        
        printf("Number of replicates = %d\n", programOptions->numDataSets);
        if (programOptions->Nscaling == 1)
            printf("Haploid data (%d)\n", programOptions->Nscaling);
        else
            printf("Diploid data (%d)\n", programOptions->Nscaling);
        printf("Number of clones = %d\n", programOptions->numClones);
    }
    for (int j = 0; j < programOptions->numClones; j++)
    {
        // Checking
      //  if (numClones <  CloneNameBegin[j] )
        if (programOptions->numClones <  CloneNameBegin[j] )
        {
            fprintf (stderr, "PARAMETER ERROR: Clon (%d) is higher than the number of clones (%d). (d)\n\n",  CloneNameBegin[j], programOptions->numClones);
            PrintUsage();
        }
        if ( CloneSampleSizeBegin[j] > ClonePopSizeBegin[j] )
        {
            fprintf (stderr, "PARAMETER ERROR: Clone (%d) cannot have sample size (%d) higher than population size (%d). (d)\n\n", j, CloneSampleSizeBegin[j] , ClonePopSizeBegin[j]);
            PrintUsage();
        }
    }
}
/***************** bbinClones *****************/
/* binary search in the probabilities with clones */
int bbinClones (double dat, double *v, int n)
{
    int init, end, middle;
    
    if (dat >= 0 && dat <= v[1])
        return (1); /* first population */
    
    init = 1;
    end = n;
    
    while (init <= end)
    {
        middle = (init + end) / 2;
        
        if (dat > v[middle - 1] && dat <= v[middle])
            return (middle);
        else if (dat > v[middle])
            init = middle + 1;
        else
            end = middle - 1;
    }
    
    fprintf (stderr, "\n Warning in bbinClones function");
    exit (-1);
    return -1;
}


/***************************** ReadParametersFromFile *******************************/
/* Reads parameter values from the parameter file */

void ReadParametersFromFile(ProgramOptions *programOptions, FilePaths *filePaths,
                            int **CloneNameBegin,
                            int **CloneSampleSizeBegin,
                            int **ClonePopSizeBegin,
                            double **CloneBirthRateBegin,
                            double **CloneDeathRateBegin,
                            double **CloneTimeOriginInput,
                            double Mij[4][4],
                             double freq[4]
                            )
{
    int   j, z;
    char  ch;
    float   argument;
    double    sumPi, sum;
    double argumentDouble;
    double argumentDouble1;
    int argumentInt;
    long int argumentLongInt;
    int *pInt;
    double  *pDouble;
    long int *pLongInt;
    
    /* Used: N X C R D M O T K Y # */

    if (feof(stdin))
    {
        fprintf(stderr, "PARAMETER ERROR: Unable to read parameters from stdin\n");
        exit(0);
    }
    
    ch = fgetc(stdin);
    while (isspace(ch))
        ch = fgetc(stdin);
    while (ch == '[')
    {
        ReadUntil(stdin, ']', "closing bracket");
        ch = fgetc(stdin);
        while (isspace(ch))
            ch = fgetc(stdin);
    }
    
    while (!feof(stdin))
    {
        argument = 0;
        ch = toupper(ch);
        switch (ch)
        {
                
            case 'N':
                if (fscanf(stdin, "%d", &argumentInt) != 1 || argumentInt < 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad number of replicates (%d)\n\n", (int)*pInt);
                    PrintUsage();
                }
                programOptions->numDataSets =argumentInt;
                break;
            case '#':
                if (fscanf(stdin, "%lu bytes", &argumentLongInt) != 1 || argumentLongInt < 0)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad seed (#) (%d)\n\n", (int) pLongInt);
                    PrintUsage();
                }
                programOptions->userSeed =argumentLongInt;
                break;
                
                
            case 'X':
                if (fscanf(stdin, "%f", &argument) != 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad haplid/diploid chosen (x) (%d)\n\n", (int) argument);
                    PrintUsage();
                }
                //*Nscaling = (int) argument;
                programOptions->Nscaling =(int) argument;
                if (programOptions->Nscaling < 1 || programOptions->Nscaling > 2)
                {
                    fprintf (stderr, "PARAMETER ERROR: Haploid/diplod option (x) (1-2) (%d)\n\n", programOptions->Nscaling);
                    PrintUsage();
                }
                break;
            case 'C':
                if (fscanf(stdin, "%f", &argument) != 1 || argument <= 0)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad number of clones (must be 1 or higher) (%d)\n\n", (int) argument);
                    PrintUsage();
                }
                //*numClones = (int) argument;
                programOptions->numClones=(int) argument;
                
                *CloneNameBegin =  malloc( programOptions->numClones * (long) sizeof(int));
                *CloneSampleSizeBegin = malloc( programOptions->numClones* (long) sizeof(int));
                *ClonePopSizeBegin = calloc( programOptions->numClones , (long) sizeof(int));
                *CloneBirthRateBegin =  calloc( programOptions->numClones, (long) sizeof(double));
                *CloneDeathRateBegin =  calloc( programOptions->numClones, (long) sizeof(double));
                *CloneTimeOriginInput =  calloc( programOptions->numClones, (long) sizeof(double));
                
                if (*CloneNameBegin == NULL || *CloneSampleSizeBegin == NULL || *ClonePopSizeBegin == NULL || *CloneBirthRateBegin == NULL || *CloneDeathRateBegin == NULL || *CloneTimeOriginInput == NULL)
                {
                    fprintf (stderr, "PARAMETER ERROR: Could not allocate variables for clones\n");
                    exit (1);
                }
                
                for (j = 0; j <  programOptions->numClones; j++)
                {
                    
                    for (z = 1; z <= 6; z++)
                    {
                        if (z == 1)
                        {
                            fscanf(stdin, "%f", &argument);
                            *(*CloneNameBegin + j) = (int) argument;
                            
                            
                            if ( *(*CloneNameBegin + j)  <= 0 ||  *(*CloneNameBegin + j)  >  programOptions->numClones)
                            {
                                fprintf (stderr, "PARAMETER ERROR: Bad number for clone %d (should be higher than 0 and lower than the number of clones %d) (%d)\n\n", j,  programOptions->numClones, *(*CloneNameBegin + j) );
                                PrintUsage();
                            }
                        }
                        if (z == 2)
                        {
                            fscanf(stdin, "%f", &argument);
                            *(*CloneSampleSizeBegin + j) = (int) argument;
                            programOptions->numCells= programOptions->numCells + (int) argument;
                            if (*(*CloneSampleSizeBegin + j)  < 0)
                            {
                                fprintf (stderr, "PARAMETER ERROR: Bad sample size for clone %d (should not be negative) (%d)\n\n", j, *(*CloneSampleSizeBegin + j) );
                                PrintUsage();
                            }
                        }
                        if (z == 3)
                        {
                            fscanf(stdin, "%f", &argument);
                            *(*ClonePopSizeBegin + j)  = (int) argument;
                            if (*(*ClonePopSizeBegin + j) < 0)
                            {
                                fprintf (stderr, "PARAMETER ERROR: Bad population size for clone %d (should be higher than 0) (%d)\n\n", j, *(*ClonePopSizeBegin + j) );
                                PrintUsage();
                            }
                        }
                        if (z == 4)
                        {
                            fscanf(stdin, "%f", &argument);
                            *(*CloneBirthRateBegin + j) = (double) argument;
                        }
                        if (z == 5)
                        {
                            fscanf(stdin, "%f", &argument);
                            *(*CloneDeathRateBegin + j) = (double) argument;
                        }
                        if (z == 6)
                        {
                            fscanf(stdin, "%f", &argument);
                            *(*CloneTimeOriginInput + j) = (double) argument;
                            if (*(*CloneTimeOriginInput + j)  < 0)
                            {
                                fprintf (stderr, "PARAMETER ERROR: Bad time to origin for clone %d (should not be negative) (%lf)\n\n", j,  *(*CloneTimeOriginInput + j) );
                                PrintUsage();
                            }
                        }
                        
                    }
    
                }
                break;
                
                
            case 'U':
                
                if (fscanf(stdin, "%lf", &argumentDouble) != 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad mutation rate (%f) \n\n", argumentDouble);
                    PrintUsage();
                }
                programOptions->mutationRate=  (double) argumentDouble;
                break;
                
            case 'B':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 0  || argument > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad alphabet (%d)\n\n", (int) argument);
                    PrintUsage();
                }
                //*alphabet = (int) argument;
                programOptions->alphabet =(int) argument;
                break;
                
            case 'Y':
                if (fscanf(stdin, "%d",  &argumentInt) != 1 || argumentInt < 0)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad noisy value (%d)\n\n", (int) argumentInt);
                    PrintUsage();
                }
                programOptions->noisy =(int) argumentInt;
                break;
                
            case 'D':
                if (fscanf(stdin, "%lf", &argumentDouble) !=1 ||argumentDouble < 0 ||argumentDouble > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad allelic dropout rate (%f)\n\n", argumentDouble);
                    PrintUsage();
                }
                programOptions->ADOrate =(double) argumentDouble;
                break;
            case 'O':
                if (fscanf(stdin, "%d", &argumentInt) < 0 )
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad selection for outgroup %d (0: No outgroup, 1 outgroup with one branch, 2 outgroup with two branches\n\n", argumentInt);
                    PrintUsage();
                }
                programOptions->outgroupSelection = argumentInt;
                if ( programOptions->outgroupSelection != 0 &&  programOptions->outgroupSelection != 1 &&  programOptions->outgroupSelection != 2)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad selection for outgroup %d (0: No outgroup, 1 outgroup with one branch, 2 outgroup with two branches\n\n",  programOptions->outgroupSelection);
                    PrintUsage();
                }
                
                if ( programOptions->outgroupSelection == 0)
                {
                    programOptions->thereisOutgroup = NO;
                    programOptions->outgroupBranchLength_RootSample = 0.0;
                    programOptions->outgroupBranchLength_Root1Root2 = 0.0;
                }
                else if (programOptions->outgroupSelection == 1)
                {
                    //*thereisOutgroup = YES;
                    programOptions->thereisOutgroup=YES;
                    programOptions->outgroupBranchLength_Root1Root2 = 0.0;
                    
                    if (fscanf(stdin, "%lf", &argumentDouble) < 0)
                    {
                        fprintf(stderr, "PARAMETER ERROR: Bad outgroup branch length (Root-Sample) value (%f)\n\n", argumentDouble);
                        PrintUsage();
                    }
                    programOptions->outgroupBranchLength_RootSample=argumentDouble;
                }
                else if (programOptions->outgroupSelection == 2)
                {
                    //*thereisOutgroup = YES;
                    programOptions->thereisOutgroup=YES;
                    if (fscanf(stdin, "%lf", &argumentDouble) < 0)
                    {
                        fprintf(stderr, "PARAMETER ERROR: Bad outgroup branch length (Root1-Root2) value (%f)\n\n", argumentDouble);
                        PrintUsage();
                    }
                    programOptions->outgroupBranchLength_Root1Root2=argumentDouble;
                    if (fscanf(stdin, "%lf", &argumentDouble) < 0)
                    {
                        fprintf(stderr, "PARAMETER ERROR: Bad outgroup branch length (Root2-Sample) value (%f)\n\n", argumentDouble);
                        PrintUsage();
                    }
                    programOptions->outgroupBranchLength_RootSample=argumentDouble;
                }
                else
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad selection for outgroup %d (0: No outgroup, 1 outgroup with one branch, 2 outgroup with two branches\n\n", argumentInt);
                    PrintUsage();
                }
                //outgroupBranchLength_RootSample = 0;
                break;
            case 'I':
                if (fscanf(stdin, "%lf", &argumentDouble) !=1 || argumentDouble < 0 || argumentDouble > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad allelic imbalance (%f)\n\n", argumentDouble);
                    PrintUsage();
                }
                programOptions->allelicImbalance=argumentDouble;
                break;
            case 'R':
                if (programOptions->doHKY == YES)
                {
                    fprintf (stderr, "PARAMETER ERROR: Cannot specify a mutation matrix (GTR model) and a ti/tv (HKY model) at the same time\n\n");
                    PrintUsage();
                }
                if (fscanf(stdin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf",
                           &Mij[0][0], &Mij[0][1], &Mij[0][2], &Mij[0][3],
                           &Mij[1][0], &Mij[1][1], &Mij[1][2], &Mij[1][3],
                           &Mij[2][0], &Mij[2][1], &Mij[2][2], &Mij[2][3],
                           &Mij[3][0], &Mij[3][1], &Mij[3][2], &Mij[3][3])!=16)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad general rate matrix (-rx x x x x x x x x x x) (AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT)\n\n");
                    PrintUsage();
                }
                
                if (Mij[0][0] != 0  || Mij[1][1] != 0 || Mij[2][2] != 0 || Mij[3][3] != 0)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad general rate matrix: diagonals should be 0 \n\n");
                    PrintUsage();
                }
                //*thereIsMij = YES;
                programOptions->thereIsMij=YES;
                if (CheckMatrixSymmetry (Mij) == YES)
                {
                   
                    //programOptions->doJC = NO;
                    programOptions->doJC = YES;
                    programOptions->doHKY = NO;
                    //programOptions->doGTR = YES;
                    programOptions->doGTR = NO;
                    programOptions->doGTnR = NO;
                    
                }
                else
                {
                 
                    programOptions->doJC = NO;
                    programOptions->doHKY = NO;
                   //  programOptions->doGTR = NO;
                    programOptions->doGTR =  YES;
                   // programOptions->doGTnR = YES;
                    programOptions->doGTnR = NO;
                }
                break;
            case 'P':
                if (fscanf(stdin, "%lf", &argumentDouble) !=1 || argumentDouble < 0 || argumentDouble > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad proportion of alternative model sites (%f)\n\n", argumentDouble);
                    PrintUsage();
                }
                programOptions->propAltModelSites=argumentDouble;
                if (programOptions->propAltModelSites < 0 || programOptions->propAltModelSites > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad proportion of alternative model sites (%f). It has to be between 0 and 1\n\n", programOptions->propAltModelSites);
                    PrintUsage();
                }
                if (programOptions->propAltModelSites > 0 && programOptions->altModel != ISMhap && programOptions->doSimulateFixedNumMutations == YES)
                {
                    fprintf (stderr, "PARAMETER ERROR: You cannot specify a proportion of non-ISM  sites bigger than zero if the number of mutations is fixed\n\n");
                    PrintUsage();
                }
                if (programOptions->alphabet == DNA && programOptions->propAltModelSites > 0)
                {
                    if (programOptions->altModel == Mk)
                    {
                        fprintf (stderr, "PARAMETER ERROR: The DNA alphabet and the alt model (%d) specified are incompatible", (int) programOptions->altModel);
                        PrintUsage();
                    }
                }
                else if (programOptions->alphabet == BINARY && programOptions->propAltModelSites > 0)
                {
                    if (programOptions->altModel == finiteDNA)
                    {
                        fprintf (stderr, "PARAMETER ERROR: The binary alphabet and the alt model (%d) specified are incompatible", (int) programOptions->altModel);
                        PrintUsage();
                    }
                }
                break;
            case 'M':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 0 || argument > 2)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad alternative mutation model (%d)\n\n", (int) argument);
                    PrintUsage();
                }
                programOptions->altModel = (int) argument;
                if (programOptions->alphabet == DNA && programOptions->propAltModelSites > 0)
                {
                    if (programOptions->altModel == Mk)
                    {
                        fprintf (stderr, "PARAMETER ERROR: The DNA alphabet and the model (%d) specified are incompatible", (int) argument);
                        PrintUsage();
                    }
                }
                else if (programOptions->alphabet == BINARY && programOptions->propAltModelSites > 0)
                {
                    if (programOptions->altModel == finiteDNA)
                    {
                        fprintf (stderr, "PARAMETER ERROR: The binary alphabet and the model (%d) specified are incompatible", (int) argument);
                        PrintUsage();
                    }
                }
                break;
            case 'T':
                ch = fgetc(stdin);
                if (isspace(ch))
                {
                    strcpy(filePaths->treeFile, "trees");
                }
                else
                {
                    j = 0;
                    do
                    {
                        filePaths->treeFile[j] = ch;
                        j++;
                        ch = fgetc(stdin);
                    }
                    while (!isspace(ch));
                    filePaths->treeFile[j] = '\0';
                }
                //*doPrintTrees = YES;
                programOptions->doPrintTrees = YES;
                break;
            case 'A':
                if (fscanf(stdin, "%lf %lf %d", &argumentDouble, &argumentDouble1, &argumentInt) != 3)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad mean/var/model amplification error (%f ; %f ; model=%d)\n\n",argumentDouble, argumentDouble1, argumentInt);
                    PrintUsage();
                }
                 programOptions->meanAmplificationError= argumentDouble;
                 programOptions->varAmplificationError= argumentDouble1;
                 programOptions->simulateOnlyTwoTemplates=argumentInt;
                
                if ( programOptions->meanAmplificationError < 0 ||  programOptions->meanAmplificationError > 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad mean amplification error (%f)\n\n",  programOptions->meanAmplificationError);
                    PrintUsage();
                }
                if ( programOptions->varAmplificationError < 0 || ( programOptions->meanAmplificationError > 0 &&  programOptions->varAmplificationError >= ( programOptions->meanAmplificationError * (1.0 -  programOptions->meanAmplificationError))))
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad variance amplification error (%f); it has to be < mean*(1-mean)\n\n",  programOptions->meanAmplificationError);
                    PrintUsage();
                }
                if ( programOptions->simulateOnlyTwoTemplates != 0 &&  programOptions->simulateOnlyTwoTemplates != 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad simulateOnlyTwoTemplates error (%d); it has to be 0 (assume 4 templates) or 1 (assume 2 templates)",  programOptions->simulateOnlyTwoTemplates);
                    PrintUsage();
                }
                break;
            case 'E':
                if (fscanf(stdin, "%lf",  &argumentDouble) !=1 ||  argumentDouble < 0 ||  argumentDouble > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad sequencing error (%f)\n\n",  argumentDouble);
                    PrintUsage();
                }
                 programOptions->sequencingError= argumentDouble;
                break;
            case 'S':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 0 || argument > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad alternative  do simulated data or not (%d)\n\n", (int) argument);
                    PrintUsage();
                }
                programOptions->doSimulateData= (int)argument;
                break;
            case 'K':
                ch = fgetc(stdin);
                if (isspace(ch))
                {
                    strcpy(filePaths->timesFile, "times");
                }
                else
                {
                    j = 0;
                    do
                    {
                        filePaths->timesFile[j] = ch;
                        j++;
                        ch = fgetc(stdin);
                    }
                    while (!isspace(ch));
                    filePaths->timesFile[j] = '\0';
                }
                //*doPrintTimes = YES;
                programOptions->doPrintTimes = YES;
                break;
            case 'J':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad number of mutations (%d)\n\n", (int) argument);
                    PrintUsage();
                }
                programOptions->numFixedMutations = (int) argument;
                //programOptions->doSimulateFixedNumMutations = YES;
                programOptions->doSimulateFixedNumMutations = YES;
                if (programOptions->propAltModelSites > 0 && programOptions->altModel != ISMhap)
                {
                    fprintf (stderr, "PARAMETER ERROR: You cannot specify a fixed number of mutations if there is any non-ISM  site. Set the proportion of non-ISM diploid sites to zero\n\n");
                    PrintUsage();
                }
                break;
            case 'F':
                if (fscanf(stdin, "%lf %lf %lf %lf", &freq[0], &freq[1], &freq[2], &freq[3])!=4)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad Base Frequencies\n\n");
                    PrintUsage();
                }
                else if (freq[0] == freq[1] == freq[2] == freq[3])
                    programOptions->equalBaseFreq = YES;
                else
                    programOptions->equalBaseFreq = NO;
                sumPi = freq[0] + freq[1] + freq[2] + freq[3];
                if (sumPi != 1.0)
                {
                    freq[0]/=sumPi;
                    freq[1]/=sumPi;
                    freq[2]/=sumPi;
                    freq[3]/=sumPi;
                }
                break;
            case '?':
                PrintUsage();
                break;
                /*case 'H':
                 PrintUsage();
                 break;*/
            case 'V':
                if (fscanf(stdin, "%lf", &argumentDouble)!=1 || argumentDouble <= 0)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad coverage dispersion (%f)\n\n", argumentDouble);
                    PrintUsage();
                }
                programOptions->rateVarCoverage = YES;
                break;
            case 'Z':
                if (fscanf(stdin, "%lf", &argumentDouble)!=1 || argumentDouble < 0)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad germline SNP rate (%f)\n\n", argumentDouble);
                    PrintUsage();
                }
                programOptions->SNPrate = YES;
                break;
                
            case 'H':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad sequencing coverage (%d)\n\n", (int) argument);
                    PrintUsage();
                }
                programOptions->coverage = (int) argument;
                if (programOptions->coverage > 0){
                    programOptions->doSimulateReadCounts=YES;
                   // *doSimulateReadCounts = YES;
            
                }
                if (programOptions->genotypingError > 0 && programOptions->doSimulateReadCounts == YES)
                {
                    fprintf (stderr, "PARAMETER ERROR: Cannot specify a coverage larger than 0, which implies read count generation, and a genotyping error at the same time\n\n");
                    PrintUsage();
                }
                break;
            default :
                fprintf(stderr, "PARAMETER ERROR: Incorrect parameter: %c\n\n", ch);
                PrintUsage();
                break;
        }
        ch = fgetc(stdin);
        while (isspace(ch) && !feof(stdin))
            ch = fgetc(stdin);
        while (ch == '[')
        {
            ReadUntil(stdin, ']', "closing bracket");
            ch = fgetc(stdin);
            while (isspace(ch))
                ch = fgetc(stdin);
        }
    }
    
}



/***************************** PrintUsage *******************************/
/* Prints a short description of program usage */
static void PrintUsage()
{
    fprintf (stderr, "\n\nUsage: %s%s [-n# -x# -c# (# # # # # #) -u# -o# -ttrees.tre -ktimes.txt -## -y# -? -h]", PROGRAM_NAME, VERSION_NUMBER);
    fprintf (stderr, "\n-n: number of replicates (e.g. -n1000)");
    fprintf (stderr, "\n-x: haploid/diploid (1 for haploid, 2 for diploid) (e.g. -x2)");
    fprintf (stderr, "\n-c: number of clones; for each clone: ID number, sample size, population size, birth rate, death rate, time to origin (e.g. -c5\n\t1 4 70000 0.3 0.1 56\n\t2 5 60000 0.4 0.3 110\n\t3 5 95000 0.4 0.3 117\n\t4 4 15000 0.5 0.4 95\n\t5 3 44000 0.3 0.1 53)");
    fprintf (stderr, "\n-u: mutation rate (e.g. -u9.1e-6)");
    fprintf (stderr, "\n-o: branch length to the outgroup (root-outgroup) (e.g. -o0.0325) (default is no outgroup)");
    
    /* outputs */
    fprintf (stderr, "\n-t: tree file name (e.g. -ttrees.tre)");
    fprintf (stderr, "\n-k: times file name (e.g. -ktimes.txt)");
    
    /* other settings */
    fprintf (stderr, "\n-#: seed (e.g. -#37864287)");
    fprintf (stderr, "\n-y: noisy (amount of information printed on the screen) (e.g. -y2)");
    fprintf (stderr, "\n-? -h: Print help\n");
    
    exit(-1);
}




/***************************** ReadUntil *******************************/
/* Reading in between [] from input files */
void ReadUntil(FILE *fv, char stopChar, char *what)
{
    char ch;
    
    ch = fgetc(fv);
    while (!feof(fv) && ch != stopChar)
        ch = fgetc(fv);
    
    if (feof(fv) || ch != stopChar)
    {
        fprintf(stderr, "%s missing", what);
        exit(0);
    }
}


/***************************** PrintTitle *******************************/
/* Prints program header */

static void PrintTitle()
{
    fprintf(stderr, "%s - Simulating somatic clones evolutionary histories under diverse evolutionary scenarios", PROGRAM_NAME);
    fprintf(stderr, "\nVersion %s", VERSION_NUMBER);
    fprintf(stderr, "\n___________________________________________________________________________________________________________\n\n");
}



/***************************** PrintDate *******************************/
/* Prints date and time */

static void PrintDate ()
{
    char *date;       /* define date */
    time_t now;
    
    now = time(NULL);
    date = ctime(&now);
    fprintf(stderr, "%s\n", date);
}

/***************************** RandomUniform **********************************/
/* It returns a random uniform variate in range 0..1. It is described in
 Park, S. K. and K. W. Miller. 1988. Random number generators: good
 ones are hard to find. Communications of the ACM, 31(10):1192-1201.
 */

double RandomUniform (long int *seed)
{
    long int  lo, hi, test;
    
    hi = (*seed) / 127773;
    lo = (*seed) % 127773;
    test = 16807 * lo - 2836 * hi;
    if (test > 0)
        *seed = test;
    else
        *seed = test + 2147483647;
    return (double)(*seed) / (double)2147483647;
}

/************************************************************/
/********************* ProbabilityCloneiFromClonej2 ********************/
/* Obtain the probability that clone i is originated from clone j
 */
double ProbabilityCloneiFromClonej2 (Population *PopI, Population* PopJ, Population **populations, int numClones)
{
    double  ProbabilityIJ, AboveTerm, BelowTerm;
    int     l, j;
    double  h, a, b, c, d, e, t, cum;
    Population *p;
    ProbabilityIJ = 0.0;
    AboveTerm = 0.0;
    BelowTerm = 0.0;
    l = 0;
    h = 0.0;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    d = 0.0;
    e = 0.0;
    t = 0.0;
    cum = 0.0;
    // calculate h
    // t = CloneTimeOriginInputSTDPOP[i] * (ClonePopSizeMeffectBegin[i] / ClonePopSizeMeffectBegin[j]);
    t = (PopI->timeOriginSTD ) * (PopI->effectPopSize) / ( PopJ->effectPopSize);
    h = CalculateH(t, PopJ->timeOriginSTD, PopJ->delta);
    AboveTerm = ( PopJ->popSize) * h;
    //fprintf (stderr, "AboveTerm = %lf\n", AboveTerm);
    j=0;
    for (l = PopI->order + 1; l < numClones; l++)
    {    p = *(populations + l);
        //fprintf (stderr, "\ni = %d, j = %d, l = %d\n", i, j, l);
        t = (PopI->timeOriginSTD ) * (PopI->effectPopSize) / ( p->effectPopSize);
        h = CalculateH(t, p->timeOriginSTD, p->delta);
        
        cum = cum + ( ( p->popSize) * h);
    }
    
    BelowTerm = cum;
    //fprintf (stderr, "BelowTerm = %lf\n", BelowTerm);
    
    ProbabilityIJ = AboveTerm / BelowTerm;
    //fprintf (stderr, "ProbabilityIJ = %lf\n", ProbabilityIJ);
    
    return ProbabilityIJ;
}
/************************************************************/
/********************* CalculateH ********************/
/* Calculate H for ProbabilityCloneiFromClonej function.
 */

double CalculateH (double t, double TOrigin, double delta)
{
    double  H, AboveTerm, BelowTerm, firstTerm, secondTerm;
    double  a, b;
    
    AboveTerm = 0.0;
    BelowTerm = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    a = 0.0;
    b = 0.0;
    H = 0.0;
    
    //printf ("\nInput(H), t=%lf T=%lf delta=%lf ", t, T, delta);
    
    a = 1.0 - exp(-1.0 * delta * (TOrigin - t));
    firstTerm = a * a;
    secondTerm = exp(-1.0 * delta * t);
    AboveTerm = firstTerm * secondTerm;
    //printf ("\nAboveTerm(H) = %lf (%lf %lf %lf) / delta = %lf, T = %lf, t = %lf", AboveTerm, a, firstTerm, secondTerm, delta, T, t);
    
    b = 1.0 - exp(-1.0 * delta * TOrigin);
    BelowTerm = b * b;
    //printf ("\nBelowTerm(H) = %lf", BelowTerm);
    
    H = AboveTerm / BelowTerm;
    //printf ("\nH = %lf", H);
    
    return H;
}


/************************************************************/
/********************* RandomExponential ********************/
/* Generates a random number from a Poisson distibution with
 mean lambda.
 */

double RandomExponential (double lambda, long int *seed)
{
    double  exponentialNumber, U;
    
    do
        U = RandomUniform (seed);
    while (U == 0);
    
    exponentialNumber = -log (U) / lambda;
    
    return exponentialNumber;
}

/************************************************************/
/********************* FmodelTstandard ********************/
/* Conversion from model time t>=0 to standard time.
 */

double FmodelTstandard (double t, double TOrigin, double delta)
{
    double  ModelTimeF, firstTerm, secondTerm, thirdTerm;
    double  a, b, c;
    
    ModelTimeF = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    thirdTerm = 0.0;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    
    /*fprintf (stderr, "t = %lf\n", t);
     fprintf (stderr, "T = %lf\n", T);
     fprintf (stderr, "delta = %lf\n", delta);*/
    
    a = 1.0 - exp(-1.0 * delta * TOrigin);
    firstTerm = a * a;
    //fprintf (stderr, "firstTerm = %lf\n", firstTerm);
    
    secondTerm = exp(delta * TOrigin);
    //fprintf (stderr, "secondTerm = %lf\n", secondTerm);
    
    b = 1.0 / (1.0 - exp(-1.0 * delta * (TOrigin - t)));
    c = 1.0 / (1.0 - exp(-1.0 * delta * TOrigin));
    thirdTerm = b - c;
    //fprintf (stderr, "thirdTerm = %lf\n", thirdTerm);
    
    
    ModelTimeF = firstTerm * secondTerm * thirdTerm;
    //fprintf (stderr, "ModelTimeF = %lf\n", ModelTimeF);
    
    // New formula from Carsten!
    ModelTimeF = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    thirdTerm = 0.0;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    
    a = exp(delta * t) - 1.0;
    b = 1.0 - exp(-1.0 * delta * TOrigin);
    c = 1.0 - exp(-1.0 * delta * (TOrigin - t));
    
    ModelTimeF = a * b / (delta * c);
    //fprintf (stderr, "ModelTimeF = %lf\n", ModelTimeF);
    
    
    return ModelTimeF;
}



/************************************************************/
/********************* GmodelTstandard ********************/
/* Conversion from model time t>=0 to standard time.
 */

double GstandardTmodel (double V, double TOrigin, double delta)
{
    double  StandardTimeG, firstTerm, secondTerm, thirdTerm;
    double  a, b, c, d, e;
    
    StandardTimeG = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    thirdTerm = 0.0;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    d = 0.0;
    e = 0.0;
    
    //fprintf (stderr, "\nV = %lf, T = %lf, delta = %lf\n", V, T, delta);
    
    
    firstTerm = TOrigin;
    //fprintf (stderr, "\nfirstTerm = %lf\n", firstTerm);
    
    secondTerm = 1 / delta;
    // fprintf (stderr, "secondTerm = %lf\n", secondTerm);
    
    
    //a = 1.0 - exp(-1.0 * delta * T);
    a =  exp(-1.0 * delta * TOrigin);
    
    //fprintf (stderr, "a = %lf\n", a);
    //b = (a * a) * exp(delta * T);
    b = (1 - a) * (1 - a) * (1.0 / a);
    //fprintf (stderr, "b= %lf\n", b);
    
    //c = 1.0 - exp (-1.0 * delta * T);
    c = 1 - a;
    // fprintf (stderr, "c = %lf\n", c);
    //d = c * exp(delta * T);
    d = (1 - a) * (1.0 / a);
    // fprintf (stderr, "d = %lf\n", d);
    e = V + d;
    //fprintf (stderr, "e = %lf\n", e);
    thirdTerm = log(1 - b / e);
    //fprintf (stderr, "valueOfLog = %lf\n", thirdTerm);
    thirdTerm = log(1 - ((1 - a) * (1 - a) * (1.0 / a)) / (V * delta + (1 - a) * (1.0 / a)));
    //fprintf (stderr, "ArgumentOfLog = %lf\n", 1 - b/e);
    //fprintf (stderr, "valueOfLog = %lf\n", thirdTerm);
    
    
    thirdTerm = log(1 + delta * V - a) - log(1 + (delta * V - 1) * a);
    
    //StandardTimeG = firstTerm + (secondTerm * thirdTerm);
    StandardTimeG = secondTerm * thirdTerm;
    //fprintf (stderr, "StandardTimeG = %lf\n", StandardTimeG);
    
    if ( (1 + delta * V - a) <= 0 ||   (1 + (delta * V - 1)*a ) <= 0 ) // do approximation if required
    {
        fprintf (stderr, "\nApplying approximation of math formula to avoid log(0)\n");
        StandardTimeG = 0.0;
        firstTerm = 0.0;
        secondTerm = 0.0;
        thirdTerm = 0.0;
        a = 0.0;
        b = 0.0;
        c = 0.0;
        d = 0.0;
        e = 0.0;
        
        a = 1 / delta;
        b = log(1 + delta * V);
        firstTerm = a * b;
        //fprintf (stderr, "\nfirstTerm = %lf\n", firstTerm);
        
        
        d = (V * V * delta * exp(-1.0 * delta * TOrigin)) / (1 + V);
        secondTerm =  d;
        //fprintf (stderr, "secondTerm = %lf\n", secondTerm);
        
        StandardTimeG = firstTerm - secondTerm;
        //fprintf (stderr, "StandardTimeG = %lf\n", StandardTimeG);
    }
    
    return StandardTimeG;
}


/***************************** PrintRunSettings *******************************/
/* Prints a summary of run settings */

static void PrintRunSettings (ProgramOptions programOptions, long int seed,  int TotalNumSequences,
                              int *CloneNameBegin,
                              int *CloneSampleSizeBegin,
                              int *ClonePopSizeBegin,
                              double *ClonePopSizeMeffectBegin,
                              double *CloneBirthRateBegin,
                              double *CloneDeathRateBegin,
                              double *CloneTimeOriginInput,
                              double *CloneGrowthRateBegin,
                              double *CloneTimeOriginInputSTDPOP,
                              double *CloneDelta,
                              double  *varTimeGMRCA,
                              int    Nscaling,
                              int    numClones,
                              int  numDataSets,
                              long int userSeed,
                              double  mutationRate,
                              int noisy,
                              int outgroupSelection,
                              int thereisOutgroup,
                              double outgroupBranchLength_Root1Root2,
                              double outgroupBranchLength_RootSample,
                              double          *expectedPopSize,
                              double      cumNumCA,
                              double meanNumCA,
                              double cumNumMIG,
                              double meanNumMIG,
                              double *numEventsTot
                              )
{
    int j, z;
    double meanTotal;
    
    fprintf (stderr, "\n\nRun settings\n------------");
    fprintf (stderr, "\nSeed                                 =   %-3ld", seed);
    fprintf (stderr, "\nNumber replicate data sets           =   %-3d", programOptions.numDataSets);
    fprintf (stderr, "\nNumber of sequences                  =   %-3d", programOptions.TotalNumSequences);
    if (Nscaling == 1)
        fprintf(stderr, "\nHaploid data (%d)", Nscaling);
    else
        fprintf(stderr, "\nDiploid data (%d)", Nscaling);
    fprintf(stderr, "\nNumber of clones = %d\n", numClones);
    
    fprintf(stderr, "Demographic parameters per clone\n");
    fprintf(stderr, "\tClone number\t|\tSample size\t|\tPopulation size\t|\tEfective population size\t|\tBirth rate\t|\tDeath rate\t|\tGrowth rate (br - dr)\n");
    for (z = 1; z <= numClones; z++)
    {
        fprintf(stderr, "\t%d\t\t", CloneNameBegin[z]);
        fprintf(stderr, "\t\t\t%d\t\t\t", CloneSampleSizeBegin[z]);
        fprintf(stderr, "\t%d\t\t\t", ClonePopSizeBegin[z]);
        fprintf(stderr, "\t%lf\t\t", ClonePopSizeMeffectBegin[z]);
        fprintf(stderr, "\t\t\t%lf\t", CloneBirthRateBegin[z]);
        fprintf(stderr, "\t%lf\t", CloneDeathRateBegin[z]);
        fprintf(stderr, "\t%lf\t\t", CloneGrowthRateBegin[z]);
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\tClone number\t|\tTime to origin\t|\tTime Origin Coal (model units)\t|\tDelta (gr * ps)\n");
    for (j = 1; j <= numClones; j++)
    {
        fprintf(stderr, "\t%d\t\t", CloneNameBegin[j]);
        fprintf(stderr, "\t\t\t(-)%lf\t\t", CloneTimeOriginInput[j]);
       // fprintf(stderr, "%lf\t\t", CloneTimeOriginInputSTDPOP[j]);
        //fprintf(stderr, "\t\t\t\t\t%lf\t\t\n", CloneDelta[j]);
    }
    
    /* expectation */ // The expected number of cells in a birth death process is exp( (birth-death)*(time of origin) )
    fprintf(stderr, "Expected population size per clone");
    for (j = 1; j <= numClones; j++)
    {
        //*expectedPopSize = exp (CloneTimeOriginInputSTDPOP[j] * CloneDelta[j]);
        fprintf(stderr, "\nClone %d\t-\tExpected initial population size = %lf\t-\tInput initial population size = %d", j, *expectedPopSize, ClonePopSizeBegin[j]);
    }
    
    
    fprintf(stderr, "\n\nMutation rate                        =   %lf\n", mutationRate);
    if (thereisOutgroup == NO)
        fprintf (stderr, "Simulate outgroup                    =   NO");
    else
    {
        fprintf (stderr, "Simulate outgroup                    =   YES");
        fprintf (stderr, "\n    branch length (root - sample)    =   %lf", outgroupBranchLength_RootSample);
        if (outgroupSelection == 2)
            fprintf (stderr, "\n    branch length (Root1 - Root2)    =   %lf", outgroupBranchLength_Root1Root2);
    }
    fprintf(stderr, "\nSeed                                 =   %ld", userSeed);
    fprintf(stderr, "\nNoisy                                =   %d", noisy);
    
    fprintf (stderr, "\n\nMean number of coalescence events    =   %3.2f", meanNumCA);
    fprintf (stderr, "\nMean number of migration events      =   %3.2f", meanNumMIG);
    *numEventsTot = cumNumCA + cumNumMIG;
    fprintf (stderr, "\nTotal number of events (CA+MIG)      =   %lf", *numEventsTot);
    meanTotal = *numEventsTot / programOptions.numDataSets;
    fprintf (stderr, "\nMean total number of events (CA+MIG) =   %lf", meanTotal);
    meanTotal = 0.0;
    for (j = 0; j < programOptions.numDataSets; j++)
        meanTotal = meanTotal + varTimeGMRCA[j];
    meanTotal = meanTotal / programOptions.numDataSets;
    fprintf (stderr, "\nnumDataSets  =   %d", numDataSets);
    fprintf (stderr, "\nMean time to the MRCA (TMRCA)        =   %lf", meanTotal);
    
}

/************************* SimulatePopulation ************************/
/* simulates the evolution of population until its MRCA and can receive inmigrants  ' */ /* */

void SimulatePopulation( Population *popI, Population** populations,
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
                        )
{
    int  c, d, i, j, w, k, m, cumIndivid, isCoalescence, whichInd,
    firstInd, secondInd, newInd,  foundSuperflousNode,
    isMigration, whichClone, currentNumberAliveClones;
    double     eventTime;
    TreeNode  *p, *q, *r;
    double    ran;
    double    *cumPopulPart;
    *eventNum = 0;
    int numSimClones = 0;
    double ThisRateCA = 0.0;
    double ThisTimeCA_W = 0.0;
    double ThisTimeCA_V1 = 0.0;
    double ThisTimeCA_V2 = 0.0;
    int numParcialActiveGametes;
    int choosePairIndividuals;

    numParcialActiveGametes = popI->numActiveGametes;
//    double CloneDelta = popI->delta;
//    double CloneTimeOriginInputSTDPOP = popI->timeOriginSTD;
    //int numLeftMigrations=popI->numIncomingMigrations;
    int numMigrations = (popI->numIncomingMigrations); //taking into account also the    time of origin as a migration
    
//    int ThisCloneNumber = popI->index;
    double timeNextMigration;
    int indexNextMigration = 0;
    Population *incommingPop;
//    Population *p2;
    //fprintf (stderr, "\n\n> numMigrations= %d \n", numMigrations);
    if (programOptions->noisy > 1)
        fprintf (stderr, "\n\n>> Simulating evolutionary history of clone %d (number active gametes %d, original time to origin %lf)\n", popI->index, popI->numActiveGametes, popI->timeOriginInput);
    if (programOptions->noisy > 1)
        fprintf (stderr, "\n\n> Simulating evolutionary history of clone  or order  %d ..\n", popI->order);
        //fprintf (stderr, "\n\n> Simulating evolutionary history of clone %d ..\n", popI->index);
    *currentTime=0;
    while (indexNextMigration < numMigrations) {
        timeNextMigration = (double)(popI->migrationTimes)[indexNextMigration];
        //fprintf (stderr, "\n\n> numParcialActiveGametes= %d \n", numParcialActiveGametes);
        if ( popI->numActiveGametes >= 2) {
            ThisRateCA = (double)  popI->numActiveGametes * ((double)  popI->numActiveGametes - 1) / 2.0;
            ThisTimeCA_W = RandomExponential (ThisRateCA, seed) ;
            ThisTimeCA_V1 = FmodelTstandard (*currentTime, popI->timeOriginSTD, popI->delta);
            ThisTimeCA_V1 = ThisTimeCA_V1 + ThisTimeCA_W;
            // from standard time to model time, GstandardTmodel(V, T, delta)
            ThisTimeCA_V2 = GstandardTmodel(ThisTimeCA_V1, popI->timeOriginSTD, popI->delta);
        }
        else
        {
            ThisTimeCA_V2 = timeNextMigration + 1.0; // it reached a "provisional" MRCA
        }
        if ( ThisTimeCA_V2 < timeNextMigration)
        {
            //choose randomly two lineages to coalesce
            isCoalescence = YES;
            isMigration = NO;
            *numCA = *numCA + 1;
            *eventNum= *eventNum +1;
            whichClone = popI->index;
            *currentTime = ThisTimeCA_V2; // update current time in model time
            eventTime = *currentTime;
            
            if (programOptions->noisy > 1)
            {
                fprintf (stderr, "\n\n*** Event %3d *** currentTime (model time) = %lf, currentTime (standard time) = %lf\n", *eventNum, ThisTimeCA_V2, ThisTimeCA_V1 );
                fprintf (stderr, "\n\n*** Event %3d *** currentTime (input units) = %lf\n", *eventNum, ThisTimeCA_V2);
            }
            if (programOptions->noisy == 4)
                fprintf (stderr, "* Coalescence *\n");
            MakeCoalescenceEvent(populations, popI, nodes, numClones, seed, programOptions->noisy, numActiveGametes, nextAvailable,
                                 labelNodes, currentTime,  &(programOptions->numNodes));
        }
        else
        {
            if (indexNextMigration < numMigrations - 1) //indexNextMigration corresponds to one of the true migrations
            {
                isCoalescence = NO;
                isMigration = YES;
                *numMIG = *numMIG + 1;
                *eventNum= *eventNum +1;
                *currentTime = timeNextMigration;
                eventTime = *currentTime;
                //update migration times in model time
                if (programOptions->noisy > 1)
                {
                    fprintf (stderr, "\n\n*** Event %3d *** *currentTime (model units) = %lf\n", *eventNum, ThisTimeCA_V2);
                }
                if (programOptions->noisy == 4)
                {fprintf (stderr, "* Migration *\n");}
                newInd = *nextAvailable;
                
                r = *nodes + newInd;   /* new ancester */
                r->index = *nextAvailable;
                r->label = *labelNodes;
                *labelNodes=*labelNodes+1;
                
                r->indexCurrentClone = popI->index;
                r->indexCurrentClone = popI->index;
               r->orderCurrentClone = popI->order;
               r->class = 4;
               
                
                //    p = *nodes + MatrixMigrationIDnodeMRCA[ThisCloneNumber][doAmigration]; // root of younger clone
                incommingPop = *((popI->immigrantsPopOrderedModelTime) + indexNextMigration );
                
                if (!incommingPop){
                    fprintf (stderr, "\nError. The incoming population to  poulation %d is empty \n", popI->index);
                    exit (-1);
                }
                //r->indexOldClone = incommingPop->index;
                p = *nodes + (incommingPop->nodeIdAncesterMRCA); // root of younger clone
                indexNextMigration = indexNextMigration + 1;
                p->indexCurrentClone = popI->index;
                p->indexOldClone = incommingPop->index;
                p->orderCurrentClone = popI->order;
                // link the nodes
                r->left = p;
               
                r->right = NULL;
                //choosePairIndividuals = NO;
                
                //ChooseRandomIndividual(&firstInd, numClones, popI,  &secondInd, seed, choosePairIndividuals);
                //q=*nodes + firstInd;
                //r->right = q;//choose another random living individual of the population
                
                p->anc1 = r;
                //q->anc1 = r;
                
                //connectNodes(p, NULL, r);
                //p->time = *currentTime;
                // p->timePUnits = *currentTime * (popI->effectPopSize);
                
                r->time = *currentTime;// this is equal to the time of the migration
                r->timePUnits = *currentTime * (popI->effectPopSize);
                *nextAvailable=*nextAvailable+1; /* 1 node more is available */
                
                k = p->indexCurrentClone;
                incommingPop->numActiveGametes = incommingPop->numActiveGametes - 1; /* now the other clone has 1 less node */
                // remove node from old clone in list of active gametes and add the new node of the current clone
                //popI->idsActiveGametes[popI->numActiveGametes]=r->index;//adding the superfluos node
                popI->idsActiveGametes[popI->numActiveGametes]=p->index;//adding the superfluos node
                popI->numActiveGametes = popI->numActiveGametes + 1; /* now this clone has 1 more node */
//                if (noisy > 1)
//                    fprintf (stderr, "Migration, creating node %d (clone %d) derived from node %d (clone %d)", newInd, popI->index, incommingPop->nodeIdAncesterMRCA, k);
                //fprintf (stderr, "Migration, creating node %d (clone %d) derived from node %d (clone %d)", newInd, ThisCloneNumber, MatrixMigrationIDnodeMRCA[ThisCloneNumber][doAmigration], k);
                if (programOptions->noisy > 1)
                    fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", p->timePUnits);
                   // fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", r->timePUnits);
                /* memory for number of nodes */
                if (*nextAvailable >= *numNodes)  /* if there aren't enough nodes it go into and it addition more */
                {
                    /* ReallocNodes(&numNodes, activeGametes); */
                    if (programOptions->noisy == 4)
                        fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
                    *numNodes += INCREMENT_NODES;
                    /* realloc */
                    *nodes = (TreeNode *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode));
                    if (!(*nodes))
                    {
                        fprintf (stderr, "Could not reallocate nodes (%lu bytes)\n", *numNodes  * (long) sizeof(TreeNode));
                        exit (-1);
                    }
//                    activeGametes = (int *) realloc (activeGametes, *numNodes * (long) sizeof(int));
//                    if (!activeGametes)
//                    {
//                        fprintf (stderr, "Could not reallocate activeGametes (%lu bytes)\n", *numNodes * (long) sizeof(int));
//                        exit (-1);
//                    }
                }
                
            }
            else {
                //origin reached
                *currentTime = timeNextMigration;
                indexNextMigration = indexNextMigration + 1;
                isCoalescence = NO;
                isMigration = NO;
                eventTime = *currentTime;
                if (programOptions->noisy > 1)
                    fprintf (stderr, "\n\n*** Clone origin ***\n");
                if (programOptions->noisy == 4)
                    fprintf (stderr, "Clone origin %d at time (model units) = %lf\n", popI->index, *currentTime);
                if (popI->order < numClones - 1) // do not do it for the last clone
                {
                    newInd = *nextAvailable;
                    r = *nodes + newInd;    /* new ancestor */
                    r->index = *nextAvailable;
                    r->label = *labelNodes;
                    *labelNodes=*labelNodes+1;
                   r->indexOldClone =r->indexCurrentClone = popI->index;
                    r->orderCurrentClone = popI->order;
                    r->effectPopSize=popI->effectPopSize;
                    popI->nodeIdAncesterMRCA=newInd;
                    
                    r->class = 4;
                 
                    firstInd = *nextAvailable - 1;
                    //p = nodes + activeGametes[firstInd]; // descendant node (previously generated node, nextAvailable - 1)
                    p = *nodes + firstInd; // descendant node (previously generated node, nextAvailable - 1)
                    // link the nodes
                    r->left = p;
                    r->right = NULL;
                    p->anc1 = r;
                    r->time = *currentTime;
                    r->timePUnits = *currentTime * popI->effectPopSize;
                    popI->MRCA = p;
                    
                    //connectNodes(p, NULL, r);
                    //fprintf (stderr, "\n r->index = %d, r->time = %lf\n", r->index, r->time);
                    /* readjust active nodes */
                    *nextAvailable=*nextAvailable+1; /* 1 node more is available */
                    popI->idsActiveGametes[0] = newInd;//always will be in the 0  position because there is only one left
                    
                    //popI->idsGametes[popI->numGametes] = newInd; r is a superflous node and it will be removed so no need to add it
                    //popI->numGametes = popI->numGametes +1;
                    
                    if (programOptions->noisy > 1)
                        fprintf (stderr, "Creating origin node, it creates node %d derived from node %d", newInd, firstInd);
                    if (programOptions->noisy > 1)
                        fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", r->timePUnits);
                    /* memory for number of nodes */
                    if (*nextAvailable >= *numNodes)  /* if there aren't enough nodes it go into and it addition more */
                    {
                        /* ReallocNodes(&numNodes); */
                        if (programOptions->noisy == 4)
                            fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
                        *numNodes += INCREMENT_NODES;
                        /* realloc */
                        *nodes = (TreeNode *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode));
                        if (!(*nodes))
                        {
                            fprintf (stderr, "Could not reallocate nodes (%lu bytes)\n", *numNodes  * (long) sizeof(TreeNode));
                            exit (-1);
                        }

                    }
                }
                else  {//origin of oldest pop reached
                    popI->nodeIdAncesterMRCA=*nextAvailable-1;//for the last population, nodeIdAncesterMRCA is the MRCA instead of ancester of MRCA
                    r = *nodes + *nextAvailable-1;//popI->idsActiveGametes[0]
                    r->indexOldClone = r->indexCurrentClone = popI->index;
                    r->orderCurrentClone = popI->order;
                    popI->MRCA= r;
                    
                
                }
            }
        }
        if (programOptions->noisy > 3)
        {
            fprintf (stderr, "\nActive nodes (%d):",  popI->numActiveGametes);
            for (i = 0; i < popI->numActiveGametes; i++)
                fprintf (stderr, " %d", popI->idsActiveGametes[i]);
            fprintf (stderr, "\t|\tNext node available = %d", *nextAvailable);
        }
    }
    if (programOptions->noisy > 1)
        fprintf (stderr, "\n\nEvolutionary history of clone %d is completed \n", popI->index);
}
/************************* MakeCoalescenceTree2 ************************/
/* Builds a genealogy under the structured  ' */ /* this function go by events CA, MIG, CONV */

void MakeCoalescenceTree2 (long int *seed, Population **populations,
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
                           TreeNode** nodes,
                            TreeNode** treeTips,
                           TreeNode    **treeRootInit,
                           //TreeNode    **treeRootInit,
                           char* ObservedCellNames[],
                           int *sampleSizes
                           ) {
    
    int      c, d, i, j, w, k, m, cumIndivid, *activeGametes = NULL, isCoalescence, whichInd,
    firstInd, secondInd, newInd, eventNum, numActiveGametes, foundSuperflousNode,
    isMigration, whichClone, currentNumberAliveClones;
    int     labelNodes;
    double    currentTime, eventTime;
    TreeNode  *p, *q, *r;
    double    ran;
    double    *cumPopulPart;
    int     *numParcialActiveGametes, *CumSamNodes;
    int         *fromClone;
    int         ThisCloneNumber, ThisOriginCloneNumber;
    double      minValue;
    double      ThisRateCA;
    double      ThisTimeCA_W;
    double      ThisTimeCA_V1;
    double      ThisTimeCA_V2;

    int         doAmigration;
    int         ThisCloneNumberMigrations, ThisM;
    int nextAvailable;

   
    /* defaults */
    isCoalescence = NO;
    isMigration = NO;
    newInd = whichClone = labelNodes = cumIndivid = 0;
    numParcialActiveGametes = NULL;
    fromClone = NULL;
    eventTime = 0.0;
    c = m = 0;
    whichInd = 0;
    minValue = 0.0;
    ThisCloneNumber = 0;
    ThisOriginCloneNumber = 0;
    ThisRateCA = 0.0;
    ThisTimeCA_W = 0.0;
    ThisTimeCA_V1 = 0.0;
    ThisTimeCA_V2 = 0.0;
    doAmigration = -1;
    ThisCloneNumberMigrations = -1;
    ThisM = -1;
    
    
    //*numNodes = 2 * TotalNumSequences * numClones+ 1; // (2 * TotalNumSequences) + numClones(superfluos) - 1, but let's allocate some more..
    currentNumberAliveClones = numClones;
    
    
  //InitListPossibleMigrations(populations, numClones);
    resetMigrationsList( populations,  numClones);
    
    //allocate memory for the treenodes
    *nodes = (TreeNode *) malloc ((programOptions->numNodes + 1)* (sizeof(TreeNode)+ 5* programOptions->numSites * sizeof(int) + 2*MAX_NAME * sizeof(char)+ 3 * sizeof(pll_unode_t) + 3*sizeof(pll_tree_edge_t) )); /* nodes */
    if (!(*nodes))
    {
        fprintf (stderr, "Could not allocate nodes (%lu bytes)\n", (programOptions->numNodes+ 1)  * (long) sizeof(TreeNode));
        exit (1);
    }

    for (i=0; i< programOptions->TotalNumSequences; i++){
        treeTips[i]=NULL;
    }

//    treeRootInit = (TreeNode **) calloc(1, sizeof(TreeNode *)); /* nodes pointers */
//    if (!treeRootInit)
//    {
//        fprintf (stderr, "Could not allocate treeRootInit (%lu bytes)\n", 1  * (long) sizeof(TreeNode));
//        exit (1);
//    }
    /* set everything to null */
    for (i = 0; i < *numNodes; i++)
    {
        p = (*nodes + i);
        p->left = NULL;
        p->right = NULL;
        p->anc1 = NULL;
        p->outgroup = NULL;
        p->nodeRight=NULL;
        p->nodeRight=NULL;
        p->nodeBack=NULL;
        p->time = 0;
        p->timePUnits = 0;
        p->length = 0;
        p->lengthModelUnits = 0;
        p->index = 0;
        p->label = 0;
        p->isOutgroup = NO;
        p->class = 0;
        p->indexOldClone = 0;
        p->indexCurrentClone = 0;
        p->orderCurrentClone = 0;
        p->effectPopSize=0;
        p->isLeaf = NO;
       //  p->cellName =(char *) malloc(MAX_NAME* sizeof(char));
//        if (!(p->cellName))
//        {
//            fprintf (stderr, "Could not allocate p->cellName (%lu bytes)\n", MAX_NAME * (long) sizeof(char));
//            exit (1);
//        }
        if (programOptions->doSimulateData ==YES)
        {
            
            p->maternalSequence= (int*) malloc (programOptions->numSites * sizeof(int));
            if (!(p->maternalSequence))
            {
                fprintf (stderr, "Could not allocate p->maternalSequence (%lu bytes)\n", programOptions->numSites* sizeof(int));
                exit (1);
            }
            
            p->paternalSequence= (int*) malloc (programOptions->numSites* sizeof(int));
            if (!(p->paternalSequence))
            {
                fprintf (stderr, "Could not allocate p->paternalSequence (%lu bytes)\n", programOptions->numSites* sizeof(int));
                exit (1);
            }
            p->numbersMutationsUnderSubtreePerSite=(int*) calloc (programOptions->numSites, sizeof(int));
            if (!(p->numbersMutationsUnderSubtreePerSite))
            {
                fprintf (stderr, "Could not allocate p->NumbersMutationsUnderSubtreePerSite (%lu bytes)\n", programOptions->numSites* sizeof(int));
                exit (1);
            }
            p->numbersMaternalMutationsPerSite=(int*) calloc (programOptions->numSites, sizeof(int));
            if (!(p->numbersMaternalMutationsPerSite))
            {
                fprintf (stderr, "Could not allocate p->NumbersMaternalMutationsPerSite (%lu bytes)\n", programOptions->numSites* sizeof(int));
                exit (1);
            }
            p->numbersPaternalMutationsPerSite=(int*) calloc (programOptions->numSites , sizeof(int));
            if (!(p->numbersMaternalMutationsPerSite))
            {
                fprintf (stderr, "Could not allocate p->NumbersMaternalMutationsPerSite (%lu bytes)\n", programOptions->numSites* sizeof(int));
                exit (1);
            }
            
        }
       
        p->nodeLeft=(pll_unode_t*) calloc (1 , sizeof(pll_unode_t));
        if (!(p->nodeLeft))
        {
            fprintf (stderr, "Could not allocate p->nodeLeft (%lu bytes)\n", 1* sizeof(pll_unode_t));
            exit (1);
        }
        p->nodeRight=(pll_unode_t*) calloc (1 , sizeof(pll_unode_t));
        if (!(p->nodeRight))
        {
            fprintf (stderr, "Could not allocate p->nodeRight (%lu bytes)\n",  sizeof(pll_unode_t));
            exit (1);
        }
        p->nodeBack=(pll_unode_t*) calloc (1 , sizeof(pll_unode_t));
        if (!(p->nodeBack))
        {
            fprintf (stderr, "Could not allocate p->nodeBack (%lu bytes)\n",  sizeof(pll_unode_t));
            exit (1);
        }
        
        p->edgeBack=(pll_tree_edge_t*) calloc (1 , sizeof(pll_tree_edge_t));
        if (!(p->edgeBack))
        {
            fprintf (stderr, "Could not allocate p->edgeBack (%lu bytes)\n",  sizeof(pll_tree_edge_t));
            exit (1);
        }
        p->edgeLeft=(pll_tree_edge_t*) calloc (1 , sizeof(pll_tree_edge_t));
        if (!(p->edgeLeft))
        {
            fprintf (stderr, "Could not allocate p->edgeLeft (%lu bytes)\n",  sizeof(pll_tree_edge_t));
            exit (1);
        }
        p->edgeRight=(pll_tree_edge_t*) calloc (1 , sizeof(pll_tree_edge_t));
        if (!(p->edgeRight))
        {
            fprintf (stderr, "Could not allocate p->edgeRight (%lu bytes)\n",  sizeof(pll_tree_edge_t));
            exit (1);
        }
    }
    AssignCurrentSequencesToPopulation(populations, nodes, programOptions, numClones, *numNodes, programOptions->noisy, programOptions->TotalNumSequences, &numActiveGametes,  &nextAvailable,
                                       &labelNodes, ObservedCellNames, programOptions->doUseObservedCellNames, sampleSizes);
    Population *currentPop;
    Population *fatherPop;
    i=0;
    currentTime=0.0;
    while (i < numClones) {
        currentPop = *(populations + i);
        SimulatePopulation(currentPop, populations,programOptions, seed,
                           &(programOptions->numNodes),
                           numClones,
                           cumNumCA,
                           meanNumCA,
                           cumNumMIG,
                           meanNumMIG,
                           numMIG,
                           numCA,
                           numEventsTot,
                           nodes,
                           &nextAvailable ,
                           &numActiveGametes,
                           &labelNodes,
                           &currentTime,
                           &eventNum);
        if (i< numClones-1)   //if it is not the last one
        {
            //choose the father population from which the population i came
          
          fatherPop= ChooseFatherPopulation(populations, numClones, currentPop, seed,  programOptions->noisy);
            currentPop->FatherPop = fatherPop;
            //update list of migrant times
           UpdateListMigrants(populations, numClones, currentPop, fatherPop);
        }
        i = i + 1;
    }
//    free (CumSamNodes);
    //   free (activeGametes);
    BuildTree(populations,currentPop,
              seed,
            programOptions,
              nodes,
              treeTips,
              treeRootInit,
              //TreeNode    **treeRootInit,
              &nextAvailable,
              &newInd,
              &currentTime,
              &labelNodes
              );

    if (programOptions->noisy > 1)
        fprintf (stderr, "\n\n>> Relabeling nodes on tree... \n\n");
//    if (thereisOutgroup == YES)
//        intLabel = TotalNumSequences + 2;
//    else
//        intLabel = TotalNumSequences;
    
    // RelabelNodes(treeRootInit[0], treeRootInit, &intLabel );
}
/**************** RelabelNodes **************/
/*  After getting rid of superfluos node, we
 need to relabel those so they are consecutive
 Use the indexes as labels when there
 is recombination */
void RelabelNodes(TreeNode *p, TreeNode **treeRootInit, int *intLabel)
{
    if (p != NULL)
    {
        RelabelNodes (p->left, treeRootInit, intLabel);
        RelabelNodes (p->right, treeRootInit, intLabel);
        /*RelabelNodes (p->outgroup);*/
        if (p->left == NULL && p->right == NULL) /* is tip */
        {
           // p->label = intLabel++;
          //  p->label = (*intLabel);
           // *intLabel=*intLabel+1;
            p->label = p->index ;
        }
        else                  /* all ancester */
        {
            //p->label = intLabel++;
            p->label = (*intLabel);
            *intLabel=*intLabel+1;
        }
    }
}


/***************** Index ***************/
/* Returns index for a given node */
int Index (TreeNode *p)
{
    //return (p == NULL) ? -1 : p->index+1; /* If the node haven't got bond => index = -1, else index = index+1 */
    return (p == NULL) ? -1 : p->index; /* If the node haven't got bond => index = -1, else index = index */
}


/***************** Lab ***************/
/* Returns label for a given node */
int Label (TreeNode *p)
{
    return (p->anc1 == NULL && p->left == NULL && p->right == NULL) ? -1 : p->label + 1; /* If the node haven't got ancester and descendants => label = -1, else label = label+1 */
}


/**************** PrintTrees ***************/
/*  Print unrooted trees to treefile in Newick format */
void PrintTrees(int replicate, TreeNode **treeRootInit,   FILE  *fpTrees, double mutationRate, int doUseObservedCellNames)
{
    /* there isn´t recombination */
    /*fprintf(fpTrees,"Tree.%05d = ", replicate+1);*/
//    fprintf(fpTrees, "(");
    WriteTree (treeRootInit[0], mutationRate, fpTrees, doUseObservedCellNames);
//    fprintf(fpTrees, ");\n");
    fprintf (fpTrees,");\n");
}

/**************** PrintTrees2 ***************/
/*  Print unrooted trees to treefile in Newick format */
void PrintTrees2(int replicate, TreeNode **treeRootInit,   FILE   *fpTrees2 , double mutationRate,char * ObservedCellNames[],int doUseObservedCellNames)
{
    int indexCurrentCell=0;
  
    /* there isn´t recombination */
    /*fprintf(fpTrees2,"Tree.%05d = ", replicate+1);*/
//   fprintf(fpTrees2, "(");
    WriteTree2 (treeRootInit[0], mutationRate, fpTrees2, ObservedCellNames, &indexCurrentCell, doUseObservedCellNames);
//     fprintf(fpTrees2, ");\n");
//    long len= strlen(newickString);
//    char *res = malloc(len  + strlen(");\n"));
//    if (res){
//        memcpy(res, newickString, len);
//        memcpy(res + len, ");\n", strlen(");\n")+1);
//    }
    fprintf(fpTrees2, ");\n");
   
    //fprintf (fpTrees2,"\n");
}


/******************* WriteTree ****************/
/* Writes a given (unrooted) tree from PrintTrees */
void WriteTree (TreeNode *p, double mutationRate, FILE    *fpTrees, int doUseObservedCellNames)
{
    char buffer[1024];
  
    if (p != NULL)
    {
        if(p->isOutgroup == YES)            /* Outgroup*/
        {
            /*            fprintf (fpTrees, ",outgroup:%8.6f)",p->length*mutationRate);*/
            //fprintf (fpTrees, ",outgroup:%8.6f",p->length*mutationRate);
               strcpy( p->cellName,"healthycell");
            strcpy( p->observedCellName,"healthycell");
            //p->cellName[MAX_NAME]=0;
//                fprintf (fpTrees, ",outgroup:%10.9lf",p->length);
            //fprintf (fpTrees, ",outgroup:%10.9lf",(p->anc1->time- p->time)*mutationRate);
            fprintf (fpTrees, "healthycell:%10.9lf",(p->anc1->time- p->time)*mutationRate);
        }
        else if (p->left == NULL && p->right == NULL)        /* tip of the tree */
        {
            //fprintf (stderr, "\n\n>> p->index = %d, p->class = %d \n\n", p->index, p->class);
            //fprintf (fpTrees, "samp%05d_C%dR%d:%8.6f", p->index,p->indexOldClone,p->indexOldRegion,(p->anc1->time-p->time)*mutationRate);
           //   snprintf(buffer, sizeof(buffer), "tip_i%05d_C%d_%d", p->index,p->indexOldClone,p->indexCurrentClone);
            snprintf(buffer, sizeof(buffer), "tip_i%05d_C%d_%d", p->index,p->indexOldClone,p->indexCurrentClone);
             strcpy( p->cellName,buffer);
              //strncpy( p->cellName,buffer, sizeof(p->cellName)-1);
              // p->cellName[MAX_NAME]=0;
//            fprintf (fpTrees, "tip_i%05d_C%d_%d:%10.9lf", p->index,p->indexOldClone,p->indexCurrentClone,(p->anc1->time-p->time)*mutationRate);
            //fprintf (fpTrees, "tip_i%05d_C%d_%d:%10.9lf", p->index,p->indexOldClone,p->indexCurrentClone,(p->anc1->time- p->time)*mutationRate);
            fprintf (fpTrees, "tip_i%05d_C%d_%d:%10.9lf", p->index,p->indexOldClone,p->indexCurrentClone,(p->anc1->time- p->time)*mutationRate);
        }
        else                                /* all ancester */
        {
            fprintf (fpTrees, "(");
            WriteTree (p->left, mutationRate, fpTrees, doUseObservedCellNames);
            if (p->right != NULL) // Miguel added this condition to consider an outgroup as this right node that is NULL (see add outgroup)
            {
                fprintf (fpTrees, ",");
                WriteTree (p->right, mutationRate, fpTrees, doUseObservedCellNames);
            }
            if (p->anc1 !=NULL)
            {
                //fprintf (fpTrees, "):%8.6f",(p->anc1->time-p->time)*mutationRate);
                //snprintf(buffer, sizeof(buffer), "int_i%05d_C%d_%d",  p->index,p->indexOldClone,p->indexCurrentClone);
                snprintf(buffer, sizeof(buffer), "int_i%05d_C%d_%d",  p->index,p->indexOldClone,p->indexCurrentClone);
                strcpy( p->cellName,buffer);
                //strncpy( p->cellName,buffer, sizeof(p->cellName)-1);
                //p->cellName[MAX_NAME]=0;
//                 fprintf (fpTrees, "):%10.9lf", (p->anc1->time-p->time)*mutationRate);
                 fprintf (fpTrees, "):%10.9lf", (p->anc1->time- p->time)*mutationRate);
                
//                fprintf (fpTrees, ")int_i%05d_C%d_%d:%10.9lf",p->index, p->indexOldClone, p->indexCurrentClone, (p->anc1->time-p->time)*mutationRate);
            }
             if (p->anc1 ==NULL)  {
                
                //snprintf(buffer, sizeof(buffer), "root_i%05d_C%d_%d",  p->index,p->indexOldClone,p->indexCurrentClone);
                 snprintf(buffer, sizeof(buffer), "root_i%05d_C%d_%d",  p->index,p->indexOldClone,p->indexCurrentClone);
                strcpy( p->cellName,buffer);
                 //strncpy( p->cellName,buffer, sizeof(p->cellName)-1);
                 //p->cellName[MAX_NAME]=0;
//                    fprintf (fpTrees, ")"  );
//                  fprintf (fpTrees, "):0.00"  );
//               fprintf (fpTrees, ")root_i%05d_C%d_%d:0.00", p->index,p->indexOldClone,p->indexCurrentClone );
                
            }
            WriteTree (p->outgroup, mutationRate, fpTrees, doUseObservedCellNames);
        }
    }
}
/******************* WriteTree2 ****************/
/* Writes a given (unrooted) tree from PrintTrees */
void WriteTree2 ( TreeNode *p, double mutationRate, FILE    *fpTrees2, char *cellNames[], int *indexCurrentCell, int doUseObservedCellNames)
{
    char *currentNewick = NULL;
    //asprintf(&currentNewick, *newickString);
    char *temp = NULL;
    if (p != NULL)
    {
        if (p->isOutgroup == YES)     /* Outgroup */
        {
            /*      fprintf (fpTrees2, ",outgroup:%8.6f)",p->length*mutationRate);*/
            //fprintf (fpTrees2, ",outgroup:%8.6f",p->length*mutationRate);
            //fprintf (fpTrees2, ",outgroup:%10.9lf", p->length * mutationRate);
            //fprintf (fpTrees2, "healthycell:%10.9lf", p->length * mutationRate);
            //fprintf (fpTrees2, "healthycell:%10.9lf", p->length * mutationRate);
           fprintf (fpTrees2, "healthycell:%10.9lf", (p->anc1->timePUnits - p->timePUnits) * mutationRate);
        }
        else if (p->left == NULL && p->right == NULL)   /* tip of the tree */
        {
            //fprintf (fpTrees2, "samp%05d_C%dR%d:%8.6f", p->index,p->indexOldClone,p->indexOldRegion,(p->anc1->time-p->time)*mutationRate);
           // fprintf (fpTrees2, "tip_i%05d_C%d_%d:%10.9lf", p->index, p->indexOldClone,p->indexCurrentClone, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
            //fprintf (fpTrees2, "tip_i%05d_C%d_%d:%10.9lf", p->index, p->indexOldClone,p->indexCurrentClone, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
            if (doUseObservedCellNames == YES)
            {
                if (strcmp(cellNames[*indexCurrentCell],"healthycell")==0)
                   *indexCurrentCell =*indexCurrentCell+1;
            }
             if (doUseObservedCellNames == YES)
                 fprintf (fpTrees2, "%s:%10.9lf", p->observedCellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
             else
                 fprintf (fpTrees2, "%s:%10.9lf", p->cellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
           // fprintf (fpTrees2, "%s:%10.9lf", cellNames[*indexCurrentCell], (p->anc1->timePUnits - p->timePUnits)*mutationRate);
            *indexCurrentCell =*indexCurrentCell+1;
        }
        else                /* all ancester */
        {
            fprintf (fpTrees2, "(");
            WriteTree2 (p->left, mutationRate, fpTrees2,  cellNames, indexCurrentCell, doUseObservedCellNames);
            if (p->right != NULL) // Miguel added this condition to consider an outgroup as this right node that is NULL (see add outgroup)
            {
                fprintf (fpTrees2, ",");
                WriteTree2 (p->right, mutationRate, fpTrees2,  cellNames, indexCurrentCell, doUseObservedCellNames);
            }
            if (p->anc1 != NULL)
            {
//                //fprintf (fpTrees2, "):%8.6f",(p->anc1->time-p->time)*mutationRate);
//                fprintf (fpTrees2, ")int_i%05d_C%d:%10.9lf", p->index, p->indexCoalClone, (p->anc1->timePUnits - p->timePUnits)*1);
                   fprintf (fpTrees2, "):%10.9lf",  (p->anc1->timePUnits - p->timePUnits)*mutationRate);
            }
            if (p->anc1 ==NULL)  {
//                  fprintf (fpTrees2, ")root_i%05d_C%d_%d:0.00", p->index,p->indexOldClone,p->indexCurrentClone );
//                 fprintf (fpTrees2, "):0.00" );
            }
            WriteTree2 (p->outgroup, mutationRate, fpTrees2,  cellNames, indexCurrentCell, doUseObservedCellNames);
        }
    }
}
/******************* toNewickString2 ****************/
/*  Build a string with the Newick representation using the left, right, anc1 pointers  */
char * toNewickString2 ( TreeNode *p, double mutationRate,     int doUseObservedCellNames)
{
     char buffer[1024];
    //char *newickString = malloc( size);
    char *newickString =NULL;
    char *left=NULL;
    char *right=NULL;
    char *outgroup =NULL;
    char *temp =NULL;
    if (p != NULL)
    {
        if (p->isOutgroup == YES)     /* Outgroup */
       {
           strcpy( p->cellName,"healthycell");
          
        //strcat(newickString, "healthycell:%10.9lf", (p->anc1->timePUnits - p->timePUnits) * mutationRate);
           if (asprintf(&newickString,  "healthycell:%10.9lf",  (p->anc1->timePUnits - p->timePUnits) * mutationRate)<0)
               return NULL;
        // snprintf(newickString,  size,  "healthycell:%10.9lf",  (p->anc1->timePUnits - p->timePUnits) * mutationRate);
        return newickString;
        }
    else if (p->left == NULL && p->right == NULL)   /* tip of the tree */
        {
            snprintf(buffer, sizeof(buffer), "tip_i%05d_C%d_%d", p->index,p->indexOldClone,p->indexCurrentClone);
            strcpy( p->cellName,buffer);
//            if (p->isOutgroup == YES)     /* Outgroup */
//            {
//                //strcat(newickString, "healthycell:%10.9lf", (p->anc1->timePUnits - p->timePUnits) * mutationRate);
//                snprintf(newickString,  size,  "healthycell:%10.9lf",  (p->anc1->timePUnits - p->timePUnits) * mutationRate);
//                return newickString;
//            }
            //else{
            if  (doUseObservedCellNames == YES)
            {
                if (asprintf(&newickString,   "%s:%10.9lf",  p->observedCellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate)<0)
                  return NULL;
                //snprintf(newickString,  size,  "%s:%10.9lf",  p->observedCellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
                 return newickString;
            }
            else{
                if (asprintf(&newickString,   "%s:%10.9lf",  p->cellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate)<0)
                    return NULL;
                //snprintf(newickString,  size,  "%s:%10.9lf",  p->observedCellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
               
                 return newickString;
            }
          //  }
        }
        else
        {
            // fprintf (fpTrees2, "(");
            if ( p->left != NULL  )
            {
                left = toNewickString2 (p->left, mutationRate,   doUseObservedCellNames);
                //right = toNewickString2 (p->right, mutationRate,   cellNames);
                //snprintf(newickString, size, "(%s,%s):%10.9lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate );
//                free(left);
//                left = NULL;
//                free(right);
//                right=NULL;
//                return newickString;
            }
            if ( p->right != NULL  )
            {
                right = toNewickString2 (p->right, mutationRate,   doUseObservedCellNames);
//                snprintf(newickString, size, "(%s,%s):%10.9lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate );
//                free(left);
//                left = NULL;
//                free(right);
//                right=NULL;
//                return newickString;
            }
            outgroup =toNewickString2 (p->outgroup, mutationRate,   doUseObservedCellNames);
            if(left!= NULL && right!= NULL && p->anc1 != NULL)
            {
                snprintf(buffer, sizeof(buffer), "int_i%05d_C%d_%d",  p->index,p->indexOldClone,p->indexCurrentClone);
                strcpy( p->cellName,buffer);
                
                if (asprintf(&newickString, "(%s,%s):%10.9lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate )<0)
                    return NULL;
                //snprintf(newickString, size, "(%s,%s):%10.9lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate );
                free(left);
                left = NULL;
                free(right);
                right=NULL;
                free(outgroup);
                outgroup=NULL;
            }
            else if (left != NULL &&  right!= NULL  && p->anc1 == NULL)
            {
                snprintf(buffer, sizeof(buffer), "root_i%05d_C%d_%d",  p->index,p->indexOldClone,p->indexCurrentClone);
                strcpy( p->cellName,buffer);
               // left = toNewickString2 (p->left, mutationRate,   cellNames);
                //right = toNewickString2 (p->right, mutationRate,   cellNames);
                if (asprintf(&newickString,  "(%s,%s);", left, right)<0)
                    return NULL;
                //snprintf(newickString, size, "(%s,%s);", left, outgroup);
                free(left);
                left = NULL;
                free(right);
                right=NULL;
                free(outgroup);
                outgroup=NULL;
                // return newickString;
            }
            return newickString;
        }
    }
     return newickString;
}
/******************* toNewickString4 ****************/
/*  Build a string with the Newick representation using the left, right, anc1 pointers  */
char * toNewickString4 ( TreeNode *p, double mutationRate,     int doUseObservedCellNames)
{
    char buffer[1024];
    //char *newickString = malloc( size);
    char *newickString =NULL;
    char *left=NULL;
    char *right=NULL;
    char *outgroup =NULL;
    char *temp =NULL;
    if (p != NULL)
    {
        if (p->isOutgroup == YES)     /* Outgroup */
        {
            strcpy( p->cellName,"healthycell");
            
            //strcat(newickString, "healthycell:%10.9lf", (p->anc1->timePUnits - p->timePUnits) * mutationRate);
            if (asprintf(&newickString,  "healthycell:%10.9lf",  (p->anc1->time - p->time) * p->effectPopSize * mutationRate)<0)
                return NULL;
            // snprintf(newickString,  size,  "healthycell:%10.9lf",  (p->anc1->timePUnits - p->timePUnits) * mutationRate);
            return newickString;
        }
        else if (p->left == NULL && p->right == NULL)   /* tip of the tree */
        {
            snprintf(buffer, sizeof(buffer), "tip_i%05d_C%d_%d", p->index,p->indexOldClone,p->indexCurrentClone);
            strcpy( p->cellName,buffer);
            //            if (p->isOutgroup == YES)     /* Outgroup */
            //            {
            //                //strcat(newickString, "healthycell:%10.9lf", (p->anc1->timePUnits - p->timePUnits) * mutationRate);
            //                snprintf(newickString,  size,  "healthycell:%10.9lf",  (p->anc1->timePUnits - p->timePUnits) * mutationRate);
            //                return newickString;
            //            }
            //else{
            if  (doUseObservedCellNames == YES)
            {
                if (asprintf(&newickString,   "%s:%10.9lf",  p->observedCellName, (p->anc1->time - p->time)* p->effectPopSize* mutationRate)<0)
                    return NULL;
                //snprintf(newickString,  size,  "%s:%10.9lf",  p->observedCellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
                return newickString;
            }
            else{
                if (asprintf(&newickString,   "%s:%10.9lf",  p->cellName, (p->anc1->time - p->time)*  p->effectPopSize* mutationRate)<0)
                    return NULL;
                //snprintf(newickString,  size,  "%s:%10.9lf",  p->observedCellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
                
                return newickString;
            }
            //  }
        }
        else
        {
            // fprintf (fpTrees2, "(");
            if ( p->left != NULL  )
            {
                left = toNewickString2 (p->left, mutationRate,   doUseObservedCellNames);
                //right = toNewickString2 (p->right, mutationRate,   cellNames);
                //snprintf(newickString, size, "(%s,%s):%10.9lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate );
                //                free(left);
                //                left = NULL;
                //                free(right);
                //                right=NULL;
                //                return newickString;
            }
            if ( p->right != NULL  )
            {
                right = toNewickString2 (p->right, mutationRate,   doUseObservedCellNames);
                //                snprintf(newickString, size, "(%s,%s):%10.9lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate );
                //                free(left);
                //                left = NULL;
                //                free(right);
                //                right=NULL;
                //                return newickString;
            }
            outgroup =toNewickString2 (p->outgroup, mutationRate,   doUseObservedCellNames);
            if(left!= NULL && right!= NULL && p->anc1 != NULL)
            {
                snprintf(buffer, sizeof(buffer), "int_i%05d_C%d_%d",  p->index,p->indexOldClone,p->indexCurrentClone);
                strcpy( p->cellName,buffer);
                
                if (asprintf(&newickString, "(%s,%s):%10.9lf", left, right,  (p->anc1->time- p->time)*  p->effectPopSize * mutationRate )<0)
                    return NULL;
                //snprintf(newickString, size, "(%s,%s):%10.9lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate );
                free(left);
                left = NULL;
                free(right);
                right=NULL;
                free(outgroup);
                outgroup=NULL;
            }
            else if (left != NULL &&  right!= NULL  && p->anc1 == NULL)
            {
                snprintf(buffer, sizeof(buffer), "root_i%05d_C%d_%d",  p->index,p->indexOldClone,p->indexCurrentClone);
                strcpy( p->cellName,buffer);
                // left = toNewickString2 (p->left, mutationRate,   cellNames);
                //right = toNewickString2 (p->right, mutationRate,   cellNames);
                if (asprintf(&newickString,  "(%s,%s);", left, right)<0)
                    return NULL;
                //snprintf(newickString, size, "(%s,%s);", left, outgroup);
                free(left);
                left = NULL;
                free(right);
                right=NULL;
                free(outgroup);
                outgroup=NULL;
                // return newickString;
            }
            return newickString;
        }
    }
    return newickString;
}
/******************* toNewickString3 ****************/
/*  Build a string with the Newick representation using  the nodelets  */
char * toNewickString3 ( TreeNode *p, pll_unode_t * pllunodeBack, double mutationRate,   int doUseObservedCellNames)
{
    //char *newickString = malloc( size);
    char *newickString =NULL;
    char *left=NULL;
    char *right=NULL;
    char *outgroup =NULL;
    char *temp =NULL;
     char buffer[1024];
    if (p != NULL && pllunodeBack != NULL )
    {
        if (p->isOutgroup == YES)     /* Outgroup */
        {
            strcpy( p->cellName,"healthycell");
            //strcat(newickString, "healthycell:%10.9lf", (p->anc1->timePUnits - p->timePUnits) * mutationRate);
            if (asprintf(&newickString,  "healthycell:%10.9lf",  (pllunodeBack->length) * mutationRate)<0)
                return NULL;
            // snprintf(newickString,  size,  "healthycell:%10.9lf",  (p->anc1->timePUnits - p->timePUnits) * mutationRate);
            return newickString;
        }
        //else if (p->nodeLeft == NULL && p->nodeRight == NULL)   /* tip of the tree */
        else if (pllunodeBack->next == NULL )   /* tip of the tree */
        {
            snprintf(buffer, sizeof(buffer), "tip_i%05d_C%d_%d", p->index,p->indexOldClone,p->indexCurrentClone);
            strcpy( p->cellName,buffer);
            //            if (p->isOutgroup == YES)     /* Outgroup */
            //            {
            //                //strcat(newickString, "healthycell:%10.9lf", (p->anc1->timePUnits - p->timePUnits) * mutationRate);
            //                snprintf(newickString,  size,  "healthycell:%10.9lf",  (p->anc1->timePUnits - p->timePUnits) * mutationRate);
            //                return newickString;
            //            }
            //else{
            if  (doUseObservedCellNames == YES)
            {
             if (asprintf(&newickString,   "%s:%10.9lf",  p->observedCellName,
                         (pllunodeBack->length) * mutationRate)<0)
                return NULL;
            //snprintf(newickString,  size,  "%s:%10.9lf",  p->observedCellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
             return newickString;
            //  }
            }
            else{
                if (asprintf(&newickString,   "%s:%10.9lf",  p->cellName,
                             (pllunodeBack->length) * mutationRate)<0)
                    return NULL;
                //snprintf(newickString,  size,  "%s:%10.9lf",  p->observedCellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
                return newickString;
                //  }
            }
        }
        else
        {
            // fprintf (fpTrees2, "(");
           // if ( p->left != NULL  )
            if (pllunodeBack->next != NULL)//left
            {
                left = toNewickString3 (p->left, pllunodeBack->next->back, mutationRate,   doUseObservedCellNames);
                //right = toNewickString2 (p->right, mutationRate,   cellNames);
                //snprintf(newickString, size, "(%s,%s):%10.9lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate );
                //                free(left);
                //                left = NULL;
                //                free(right);
                //                right=NULL;
                //                return newickString;
            }
            // if ( p->right != NULL  )
            if ( pllunodeBack->next->next != NULL  )//right
            {
               // right = toNewickString3 (p->right, p->nodeRight->back->back, mutationRate,   cellNames);
                 right = toNewickString3 (p->right, pllunodeBack->next->next->back, mutationRate,   doUseObservedCellNames);
                //                snprintf(newickString, size, "(%s,%s):%10.9lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate );
                //                free(left);
                //                left = NULL;
                //                free(right);
                //                right=NULL;
                //                return newickString;
            }
           // outgroup =toNewickString3 (p->outgroup, pllunodeBack->back ,mutationRate,   cellNames);
            if(left!= NULL && right!= NULL && pllunodeBack->back != NULL)
            {
                snprintf(buffer, sizeof(buffer), "int_i%05d_C%d_%d",  p->index,p->indexOldClone,p->indexCurrentClone);
                strcpy( p->cellName,buffer);
                
                if (asprintf(&newickString, "(%s,%s):%10.9lf", left, right,  (pllunodeBack->length) * mutationRate )<0)
                    return NULL;
                //snprintf(newickString, size, "(%s,%s):%10.9lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate );
                free(left);
                left = NULL;
                free(right);
                right=NULL;
                free(outgroup);
                outgroup=NULL;
                
            }
            else if (left != NULL &&  right!= NULL  && pllunodeBack->back == NULL)
            {
                snprintf(buffer, sizeof(buffer), "root_i%05d_C%d_%d",  p->index,p->indexOldClone,p->indexCurrentClone);
                strcpy( p->cellName,buffer);
                // left = toNewickString2 (p->left, mutationRate,   cellNames);
                //right = toNewickString2 (p->right, mutationRate,   cellNames);
                if (asprintf(&newickString,  "(%s,%s);", left, right)<0)
                    return NULL;
                //snprintf(newickString, size, "(%s,%s);", left, outgroup);
                free(left);
                left = NULL;
                free(right);
                right=NULL;
                free(outgroup);
                outgroup=NULL;
                // return newickString;
            }
            return newickString;
        }
    }
    return newickString;
}
/********************* safe_strcat **********************/
/* safe_strcat:  safe concatenation of strings with strncat*/
char *safe_strcat(char *dest, size_t size, const char *src)
{
    char *p = memchr(dest, '\0', size);
    if (p != NULL) {
        strncat(p, src, size - (p - dest) - 1);
    }
    return dest;
}
/********************* PrintTimes **********************/
/* Prints to timesfile a detailed description of
 the tree: nodes, times, branch lengths */

void PrintTimes(int replicate, FILE   *fpTimes, double mutationRate, TreeNode *nodes,  int thereisOutgroup)
{
    /* there isn't recombination */
    fprintf (fpTimes, "\n\nDataset %d", replicate + 1);
    fprintf (fpTimes, "\n              ------------ Nodes -------------");
    fprintf (fpTimes, "\n    class    | label  index  (left right anc) |         time     time length    branch length");
    fprintf (fpTimes, "\n----------------------------------------------------------------------------------------------\n");
    ListTimes (0, mutationRate, nodes, fpTimes, thereisOutgroup);
}


/********************* PrintTimes2 **********************/
/* Prints to timesfile a detailed description of
 the tree: nodes, times, branch lengths */

void PrintTimes2(int replicate, FILE  *fpTimes2, double mutationRate,  TreeNode *nodes,  int thereisOutgroup)
{
    /* there isn't recombination */
    fprintf (fpTimes2, "\n\nDataset %d", replicate + 1);
    fprintf (fpTimes2, "\n              ------------ Nodes -------------");
    fprintf (fpTimes2, "\n    class    | label  index  (left right anc) |         time     time length    branch length");
    fprintf (fpTimes2, "\n----------------------------------------------------------------------------------------------\n");
    ListTimes2 (0, mutationRate, nodes, fpTimes2, thereisOutgroup);
}



/********************** ListTimes ************************/
/* Writes a given tree description from ListTimes   */

void ListTimes (int j, double mutationRate, TreeNode *nodes, FILE *fpTimes, int thereisOutgroup)
{
    /* It does not list superfluous nodes */
    TreeNode  *p;
    int     i = 0;
    
    do
    {
        p = nodes + i;
        
        if (p->isOutgroup == YES)     /* Outgroup */
            fprintf (fpTimes, "%13s   %4d   %4d  (%4d %4d %4d) |   %10.4lf      %10.4lf       %10.9lf\n",
                     "outgroup", Label(p), p->index, Index(p->left), Index(p->right), Index(p->anc1), p->timePUnits, p->anc1->timePUnits - p->timePUnits, (p->anc1->timePUnits - p->timePUnits) * mutationRate);
        
        else if (p->anc1 != NULL && p->left != NULL && p->right != NULL)        /* No MRCA, no tip (internal ancester) */
            fprintf (fpTimes, "%5s_C%dR%d(f)   %4d   %4d  (%4d %4d %4d) |   %10.4lf      %10.4lf       %10.9lf\n",
                     "int", p->indexOldClone, p->indexCurrentClone, Label(p), p->index, Index(p->left), Index(p->right), Index(p->anc1), p->timePUnits, p->anc1->timePUnits - p->timePUnits, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
        
        else if (p->anc1 != NULL && p->left == NULL && p->right == NULL)        /* tip */
            fprintf (fpTimes, "%8s_C%dR%d   %4d   %4d  (%4d %4d %4d) |   %10.4lf      %10.4lf       %10.9lf\n",
                     "tip", p->indexOldClone, p->indexCurrentClone, Label(p), p->index, Index(p->left), Index(p->right), Index(p->anc1), p->timePUnits, p->anc1->timePUnits - p->timePUnits, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
        
        else if (p->class == 5 || (p->anc1 == NULL && p->left != NULL && p->right != NULL))       /* root, MRCA */
            fprintf (fpTimes, "%8s_C%dR%d   %4d   %4d  (%4d %4d %4d) |   %10.4lf      %10.4lf       %10.9lf\n",
                     "root", p->indexOldClone, p->indexCurrentClone, Label(p), p->index, Index(p->left), Index(p->right), Index(p->anc1), p->timePUnits, 0.0, 0.0);
        
        else
            fprintf (fpTimes, "");
        
        i++;
        
        if (i > 2000)
            exit(-1);
        
    } while    ((thereisOutgroup == NO  && p->anc1  != NULL)    /* no MRCA */
                ||  (thereisOutgroup == NO  && p->left == NULL)   /* tip */
                ||  (thereisOutgroup == YES && p->isOutgroup == NO));
}



/********************** ListTimes2 ************************/
/* Writes a given tree description from ListTimes   */

void ListTimes2 (int j,  double mutationRate, TreeNode *nodes,  FILE *fpTimes2, int thereisOutgroup)
{
    /* It does not list superfluous nodes */
    TreeNode  *p;
    int     i = 0;
    
    do
    {
        p = nodes + i;
        if (p->isOutgroup == YES)     /* Outgroup */
            fprintf (fpTimes2, "%13s   %4d   %4d  (%4d %4d %4d) |   %10.9lf      %10.9lf       %10.9lf\n",
                     "outgroup", Label(p), p->index, Index(p->left), Index(p->right), Index(p->anc1), p->time, p->anc1->time - p->time, (p->anc1->time - p->time) * mutationRate);
        else if (p->anc1 != NULL && p->left != NULL && p->right != NULL)        /* No MRCA, no tip (internal ancester) */
            fprintf (fpTimes2, "%5s_C%dR%d(f)   %4d   %4d  (%4d %4d %4d) |   %10.9lf      %10.9lf       %10.9lf\n",
                     "int", p->indexOldClone, p->indexCurrentClone, Label(p), p->index, Index(p->left), Index(p->right), Index(p->anc1), p->time, p->anc1->time - p->time, (p->anc1->time - p->time)*mutationRate);
        else if (p->anc1 != NULL && p->left == NULL && p->right == NULL)        /* tip */
            fprintf (fpTimes2, "%8s_C%dR%d   %4d   %4d  (%4d %4d %4d) |   %10.9lf      %10.9lf       %10.9lf\n",
                     "tip", p->indexOldClone, p->indexCurrentClone, Label(p), p->index, Index(p->left), Index(p->right), Index(p->anc1), p->time, p->anc1->time - p->time, (p->anc1->time - p->time)*mutationRate);
        else if (p->class == 5 || (p->anc1 == NULL && p->left != NULL && p->right != NULL))       /* root, MRCA */
            fprintf (fpTimes2, "%8s_C%dR%d   %4d   %4d  (%4d %4d %4d) |   %10.9lf      %10.9lf       %10.9lf\n",
                     "root", p->indexOldClone, p->indexCurrentClone, Label(p), p->index, Index(p->left), Index(p->right), Index(p->anc1), p->time, 0.0, 0.0);
        else
            fprintf (fpTimes2, "");
        i++;
        
        if (i > 2000)
            exit(-1);
        
    } while    ((thereisOutgroup == NO  && p->anc1  != NULL)    /* no MRCA */
                ||  (thereisOutgroup == NO  && p->left == NULL)   /* tip */
                ||  (thereisOutgroup == YES && p->isOutgroup == NO));
}

/********************** InitListPossibleMigrations ************************/
/* Initialize the list of possible migrations times in all populations */
void  InitListPossibleMigrations(Population **populations, int numClones)
{
    qsort(populations, numClones, sizeof(Population*), comparePopulationsByTimeOrigin);
    Population *p;
    int i, j;
    struct Population** pops ;
    double* migrationTimes ;
    double d;
    for (i = 0; i < numClones; i++) {
        p = *(populations + i);
        p->order = i;
        p->numPossibleMigrations = i + 1;
        p->numIncomingMigrations = 1; //the time of origin counts as one migration
       // if (!(p->migrationTimes)){
        p->migrationTimes =  malloc( (p->order + 1) * sizeof( double));
        if (!(p->migrationTimes)){
            
            fprintf (stderr, "Could not allocate p->migrationTimes (%lu bytes)\n", (p->order + 1) * (long) sizeof(double));
            exit (1);
            
        }
       // }

        d = (double)(p->timeOriginSTD);
        p->migrationTimes[0] = d;

        
        if(p->order >0 )//&& !(p->immigrantsPopOrderedModelTime) )
            p->immigrantsPopOrderedModelTime = (Population**) malloc(p->order * sizeof(struct Population*)); //besides the other possible immigrants we need to add this populations itself
//        //p->immigrantsPopOrderedModelTime[0] = p;
      
        for (j = 1; j < (p->order +1); j++)
        {
            //p->immigrantsPopOrderedModelTime[j] = (Population *) malloc(sizeof(Population));
            // if (!(p->immigrantsPopOrderedModelTime[j-1]))
            p->immigrantsPopOrderedModelTime[j-1]=NULL;
//            p->immigrantsPopOrderedModelTime[j-1] =  malloc(sizeof(Population));
//            if (!( p->immigrantsPopOrderedModelTime[j-1])){
//
//                fprintf (stderr, "Could not allocate  p->immigrantsPopOrderedModelTime[j-1] (%lu bytes)\n", (p->order + 1) * (long) sizeof(double));
//                exit (1);
//
//            }
            
            d = 2 * ((double)(p->timeOriginSTD)); //  a value greater than time of origin  standarized by  the population
            p->migrationTimes[j] = d;
            //pops[j]->timeMigrationSTDCurrentPop= DBL_MAX;//initialize the rest of migration times to some large values
        }
    }
}
/********************** resetMigrationsList ************************/
/* reset list of migrations*/
void resetMigrationsList(Population **populations, int numClones){
    Population *p;
    int i, j;
    struct Population** pops ;
    double* migrationTimes ;
    double d;
    for (i = 0; i < numClones; i++) {
        p = *(populations + i);
        p->numIncomingMigrations = 1; //the time of origin counts as one migration
        d = (double)(p->timeOriginSTD);
        p->migrationTimes[0] = d;
        for (j = 1; j < (p->order +1); j++)
        {
            //p->immigrantsPopOrderedModelTime[j] = (Population *) malloc(sizeof(Population));
            // if (!(p->immigrantsPopOrderedModelTime[j-1]))
            p->immigrantsPopOrderedModelTime[j-1]=NULL;
            //            p->immigrantsPopOrderedModelTime[j-1] =  malloc(sizeof(Population));
            //            if (!( p->immigrantsPopOrderedModelTime[j-1])){
            //
            //                fprintf (stderr, "Could not allocate  p->immigrantsPopOrderedModelTime[j-1] (%lu bytes)\n", (p->order + 1) * (long) sizeof(double));
            //                exit (1);
            //
            //            }
            d = 2 * ((double)(p->timeOriginSTD)); //  a value greater than time of origin  standarized by  the population
            p->migrationTimes[j] = d;
            //pops[j]->timeMigrationSTDCurrentPop= DBL_MAX;//initialize the rest of migration times to some large values
        }
    }
}
/********************** UpdateListMigrants************************/
/* update the migrations times of the target population(that receives population PopOrigin)  */
void  UpdateListMigrants(Population **populations, int numClones, Population *PopChild, Population *PopFather  ) {
    if (PopChild->index ==  PopFather->index) {
        fprintf (stderr, "\nError. The target population %d for  migration must be different than the population of origin %d \n", PopFather->index, PopChild->index);
        exit (-1);
    }
    if (PopFather->order <= PopChild->order ) {
        fprintf (stderr, "\nError. The target population %d for  migration must be older than the population of origin %d \n", PopFather->index, PopChild->index);
        exit (-1);
    }
    Population *pOrigin, *pTarget, *p;
    double *ptr;
    int i, j;
    int lengthMigrationsArray = (int)(PopFather->order) + 1;
    int updatedNumIncomingMigrations = PopFather->numIncomingMigrations;
   // printf ( "\n lengthMigrationsArray= %d \n", lengthMigrationsArray );
    
    Population **pops  = PopFather->immigrantsPopOrderedModelTime;
    ptr = PopFather->migrationTimes;
    for (j = 0; j < lengthMigrationsArray; j++) {
        if (ptr[j] >  (PopFather->timeOriginSTD) ) {
            break; //first null pointer
        }
    }
    double updatedMigrationTime = (PopChild->timeOriginSTD) * (PopChild->effectPopSize) / (PopFather->effectPopSize);
    PopFather-> migrationTimes[j] = updatedMigrationTime;
    if(updatedNumIncomingMigrations + 1 <= PopFather->numPossibleMigrations){
        updatedNumIncomingMigrations = updatedNumIncomingMigrations + 1;
        PopFather->numIncomingMigrations = updatedNumIncomingMigrations;
        PopFather->immigrantsPopOrderedModelTime[j-1] = PopChild;
    }
    //  fprintf (stderr ,"\n updatedNumIncomingMigrations %d \n",PopFather->numIncomingMigrations);
    //PopFather->immigrantsPopOrderedModelTime[j-1] = PopChild;
    //order immigrant population by time of origin
      if(PopFather->numIncomingMigrations > 1 )
          qsort(PopFather->migrationTimes, PopFather->numIncomingMigrations, sizeof(double), compare);
     if(PopFather->numIncomingMigrations -1 > 1 )
         qsort(PopFather->immigrantsPopOrderedModelTime, PopFather->numIncomingMigrations -1,  sizeof(Population*), comparePopulationsByTimeOrigin);
}
/********************** ChooseFatherPopulation************************/
/* choose probabilistically  the father population of a  population  */
Population* ChooseFatherPopulation(Population **populations, int numClones, Population  *PopChild,  long int *seed, int noisy) {
    
    Population *pOrigin, *pTarget, *p;
    double pij, ran;
    double *ptr;
    int i, j, k;
    double cumProb[numClones  - (int)(PopChild->order)];
    cumProb[0] = 0.0;
    for (j = PopChild->order + 1; j < numClones; j++)
    {
        cumProb[j - PopChild->order] = 0.0;
        p = *(populations + j);
        pij = ProbabilityCloneiFromClonej2(PopChild, p, populations, numClones);
        cumProb[j - PopChild->order] = cumProb[j - 1 - PopChild->order] + pij;
   
    }
    // now selecting the ancestral clone
    ran = RandomUniform(seed);
    //fprintf (stderr, "\n ran = %lf ", ran);
    int w = -1;
    for (k = 1; k < numClones ; k++)
    {
        if (ran <= cumProb[k])
        {
            w = k;
            break;
        }
    }
    Population *result =  *(populations + PopChild->order + w);
    if (noisy > 3)
        fprintf (stderr, "\nClone %d derived from clone %d\n", PopChild->index, result->index); // clone ThisCloneNumber (i) is originated from clone ThisOriginCloneNumber (j)
    /* Update list of migation times considering that clone i comes from clone j */
    if (noisy > 1)
        fprintf (stderr, "\n*** Updating list of migration times (considering that clone %d comes from clone %d) ..\n", PopChild->index, result->index);
    return (result); //the father population has  order  (PopChild->order) + w
}
/********************** AssignCurrentSequencesToPopulation************************/
/* assign current sequences to  population  */
void AssignCurrentSequencesToPopulation(Population **populations, TreeNode **nodes, ProgramOptions* programOptions,
                int numClones, int numNodes, int noisy,  int TotalNumSequences, int *numActiveGametes, int* nextAvailable,
                                        int *labelNodes, char* ObservedCellNames[], int doUseObservedCellNames, int *sampleSizes)
{
    Population *pop;
    TreeNode *p;
    int i, j;
    double currentSampleSize;
    *numActiveGametes = 0;
    int indexFirstObservedCellName;

    int*  CumSamNodes = (int *) malloc ((numClones + 1)* (long) sizeof(int)); /* cumulative of samples in clones */
    if (!CumSamNodes)
    {
        fprintf (stderr, "Could not allocate CumSamNodes (%lu bytes)\n", (numNodes + 1) * (long) sizeof(int));
        exit (1);
    }
    CumSamNodes[0] = 0;
    
    for (i = 1; i <= numClones; i++)
    {
        pop = *(populations + i-1);
        CumSamNodes[i] = 0;
        currentSampleSize = pop->sampleSize;
        //pop->idsActiveGametes= (int*) malloc( (numNodes) * sizeof( int));
        pop->numCompletedCoalescences=0;
        pop->nodeIdAncesterMRCA =0;
         //pop->idsActiveGametes= (int*) malloc( (TotalNumSequences) * sizeof( int));
        pop->idsActiveGametes=NULL;
        pop->idsActiveGametes= (int*) malloc( (pop->sampleSize + pop->numPossibleMigrations) * sizeof( int));
        
        if (!( pop->idsActiveGametes))
        {
            fprintf (stderr, "Could not allocate  pop->idsActiveGametes (%lu bytes)\n", (pop->sampleSize + pop->numPossibleMigrations ) * (long) sizeof(int));
            exit (1);
        }
        pop->idsGametes=NULL;
        pop->idsGametes= (int*) malloc( (2* pop->sampleSize + pop->numPossibleMigrations) * sizeof( int));
        if (!( pop->idsGametes))
        {
            fprintf (stderr, "Could not allocate  pop->idsGametes (%lu bytes)\n", (2* pop->sampleSize + pop->numPossibleMigrations) * (long) sizeof(int));
            exit (1);
        }
        pop ->numActiveGametes=0;
        pop ->numGametes=0;
        CumSamNodes[i] = CumSamNodes[i - 1] + currentSampleSize;
    }
    if (noisy > 1)
        fprintf (stderr, "\n Initial relation nodes-clones:");
    int  cumIndivid = 0;
    int currentPopIndex;
    *numActiveGametes = 0;
   
    for (j = 0; j < TotalNumSequences; j++)
    {
        cumIndivid++;
        p = *nodes + j;
        //activeGametes[*numActiveGametes] = j;
        p->index = j;
        p->label = j;
    
        *labelNodes = j;
        p->class = 1;
        // fprintf (stderr,"\nIn cumIndivid=%d j=%d", cumIndivid, j);
        for (i = 1; i <= numClones; i++)
        {
            pop = *(populations + i - 1);
            currentPopIndex = pop->index;
            indexFirstObservedCellName= pop->indexFirstObservedCellName;
            // Identify to which clone belongs this node (sample)
            if (cumIndivid <= CumSamNodes[i] && cumIndivid > CumSamNodes[i - 1])
            {
                //fprintf (stderr,"\ncumIndivid=%d <= CumSamNodes[i]=%d, in clone %d\n", cumIndivid, CumSamNodes[i], pop->index);
                pop->idsActiveGametes[pop->numActiveGametes]=j;
                
                
                pop->idsGametes[pop->numGametes]=j;
                pop->numGametes=pop->numGametes+1;
                
                p->indexOldClone = currentPopIndex;
                p->indexCurrentClone = currentPopIndex;
                p->effectPopSize= pop->effectPopSize;
                p->orderCurrentClone = pop->order;
                
               
                
                if (programOptions->doUseObservedCellNames)
                   strcpy( p->observedCellName,ObservedCellNames[indexFirstObservedCellName + pop->numActiveGametes ]);
                pop->numActiveGametes=pop->numActiveGametes+1;
                
                break;
            }
        }
        
//        if(doUseObservedCellNames == YES)
//            strcpy( p->observedCellName,ObservedCellNames[j]);
        
        
        if (noisy > 1)
            fprintf (stderr,"\n > The node %d(%d) belongs to clone %d", p->index, cumIndivid, p->indexOldClone);
        *numActiveGametes = *numActiveGametes + 1;
    }
    //AssignObservedCellNamestoTips(nodes, populations, sampleSizes,ObservedCellNames,  programOptions);
   // AssignObservedCellNamestoTips2(nodes, populations, sampleSizes,ObservedCellNames,  programOptions);
    free(CumSamNodes);
    CumSamNodes=NULL;
    *nextAvailable = *numActiveGametes;
    *labelNodes = *labelNodes + 1;
}
/********************** ChooseRandomIndividual************************/
/* ChooseRandomIndividual  */
void ChooseRandomIndividual(int *firstInd,   int numClones, Population *popI,  int *secondInd, long *seed, int choosePairIndividuals)
{
    double random;
    int k, w;
    double *cumPopulPart = (double *) malloc((popI->numActiveGametes + 1)* (long) sizeof(double));
    if (!cumPopulPart)
    {
        fprintf (stderr, "Could not allocate cumPopulPart (%lu bytes)\n", (popI->numActiveGametes + 1) * (long) sizeof(double));
        exit (-1);
    }
    cumPopulPart[0] = 0;
    for (k = 1; k <= popI->numActiveGametes; k++)
        cumPopulPart[k] = 0;
    for (k = 1; k <= popI->numActiveGametes; k++)
        cumPopulPart[k] = cumPopulPart[k - 1] + 1.0 / (popI->numActiveGametes);

    w = 0;
//    for (k = 0; k < numClones; k++)
//    {
//        pop = *(populations + k);
//        w = w + (pop->numActiveGametes);
//        //fprintf (stderr, "\nClone %d with %d gametes. w=%d ", k, pop->numActiveGametes, w);
//    }
//    if (w != numActiveGametes)
//    {
//        fprintf (stderr, "\nError. The sum of partial active gametes is different to the total number of gametes, w %d != numActiveGametes %d. In Coalescence.", w, numActiveGametes);
//        exit (-1);
//    }
    random = RandomUniform(seed);
    //fprintf (stderr, "\nran = %lf ", ran);
    *firstInd = bbinClones(random, cumPopulPart, popI->numActiveGametes)-1;
    w = 0;
    
    if (*firstInd >= popI->numActiveGametes || *firstInd < 0 ) /* checking */
    {
        fprintf (stderr, "\n\nERROR: firstInd out of range!\n");
        exit (-1);
    }
    
    if (choosePairIndividuals== YES && popI->numActiveGametes > 1) {
        
        do//choose randomly another individual to coalesce
        {
            random = RandomUniform(seed);
            *secondInd = bbinClones(random, cumPopulPart, popI->numActiveGametes)-1;
            
        } while (*firstInd == *secondInd  );
    }
    free (cumPopulPart);
    cumPopulPart=NULL;
}

/********************** MakeCoalescenceEvent************************/
/*  choose 2  active individuals  to make coalescent  */
void MakeCoalescenceEvent(Population **populations, Population *popI, TreeNode **nodes, int numClones, long int* seed, int noisy,   int *numActiveGametes, int* nextAvailable,
                           int*labelNodes, double *currentTime, int *numNodes)
{
    int k, w;
    double rand, ran;
    TreeNode  *p, *q, *r;
    int firstInd, i, j, secondInd=0, newInd=0;
    int choosePairIndividuals = YES;
    
    ChooseRandomIndividual(&firstInd, numClones, popI,  &secondInd, seed, choosePairIndividuals);
    
    newInd = *nextAvailable;
    if (noisy > 1)
        fprintf (stderr, "Coalescence involving %d and %d to create node %d (in clone %d)", popI->idsActiveGametes[firstInd], popI->idsActiveGametes[secondInd], newInd, popI->index);
    /*  set pointers between nodes */
    p = *nodes + popI->idsActiveGametes[firstInd];
    q = *nodes + popI->idsActiveGametes[secondInd];
    r = *nodes + newInd;    /* new ancester */
    r->index = *nextAvailable;
    r->label = *labelNodes;
    *labelNodes=*labelNodes+1;
    r->indexOldClone = r->indexCurrentClone = popI->index;//here the clone number is updated
    // r->indexCurrentClone = p->indexCurrentClone;
    // r->orderCurrentClone = p->orderCurrentClone;
    r->orderCurrentClone = popI->order;
    r->effectPopSize=popI->effectPopSize;
    r->class = 4;
    // link the nodes
    r->left = p;
    r->right = q;
    p->anc1 = r;
    q->anc1 = r;
    r->time = *currentTime;
    r->timePUnits = *currentTime * (popI->effectPopSize);

    //fprintf (stderr, "\n r->index = %d, r->time = %lf, ClonePopSizeMeffectBegin[ThisCloneNumber] = %lf, ThisCloneNumber = %d\n", r->index, r->time, ClonePopSizeMeffectBegin[ThisCloneNumber], ThisCloneNumber);
    if (noisy > 1)
        fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", r->timePUnits);
    /* readjust active nodes */
  
    popI->idsActiveGametes[firstInd] = newInd;
    popI->idsActiveGametes[secondInd] = popI->idsActiveGametes[popI->numActiveGametes - 1];;
    *numActiveGametes = *numActiveGametes - 1; /* less 1 active node */
    
    //update list ids nodes
    popI->idsGametes[popI->numGametes] = newInd;
    popI->numGametes = popI->numGametes +1;
    
    *nextAvailable=*nextAvailable+1; /* 1 node more is available */
    
    popI->CoalescentEventTimes[ popI->numCompletedCoalescences]=  r->time;
    popI->numActiveGametes = popI->numActiveGametes - 1; /* now this clone
                                                          has 1 less node */
   
    popI->numCompletedCoalescences= popI->numCompletedCoalescences+1;
    /* memory for number of nodes */
    if (*nextAvailable >= *numNodes)  /* if there aren't enough nodes it go into and it addition more */
    {
        /* ReallocNodes(&numNodes, activeGametes); */
        if (noisy == 4)
            fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
        *numNodes += INCREMENT_NODES;
        /* realloc */
        *nodes = (TreeNode *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode) );
        if (!(*nodes))
        {
            fprintf (stderr, "Could not reallocate nodes (%lu bytes)\n", *numNodes  * (long) sizeof(TreeNode));
            exit (-1);
        }
        popI->idsActiveGametes = (int *) realloc (popI->idsActiveGametes, *numNodes * (long) sizeof(int));
        if (!(popI->idsActiveGametes))
        {
            fprintf (stderr, "Could not reallocate idsActiveGametes for the current population(%lu bytes)\n", *numNodes * (long) sizeof(int));
            exit (-1);
        }
    }
}
/********************** BuildTree************************/
/*  build tree */
void BuildTree(Population **populations,Population *CurrentPop,
               long int *seed,
               ProgramOptions *programOptions,
               TreeNode    **nodes,
               TreeNode   **treeTips,
               TreeNode    **treeRootInit,
               //TreeNode    **treeRootInit,
               int *nextAvailable,
               int *newInd,
               double *currentTime,
               int *labelNodes
               )
{
    int i, j, k;
    int indexCurrentTip;
    int  foundSuperflousNode;
    TreeNode *p, *q, *r;
    /********** BUILDING TREES ***********/
    if (programOptions->noisy > 1)
    {
        fprintf (stderr, "\n>> Building trees ..");
    }
    
    j=0;
    indexCurrentTip=0;
    /* get rid of superflous nodes */
    foundSuperflousNode = YES;
    while (foundSuperflousNode == YES)
    {
        foundSuperflousNode = NO;
        for (i = 0; i < *nextAvailable; i++) // available all the nodes
        {
            p = *nodes + i;

            //fprintf (stderr, "\n\np->index = %d", p->index);
            
            if (p->left == NULL && p->right == NULL && p->anc1 == NULL)
            {
                // nothing to do with this node because it is not connected to anything
                //fprintf (stderr, "\n * nothing to do with this node because it is not connected to anything");
            }
            else if (p->left == NULL && p->right == NULL && p->anc1 != NULL)
            {
               // (*treeTips)[indexCurrentTip]=p;
                if(indexCurrentTip <  programOptions->TotalNumSequences){
                    treeTips[indexCurrentTip]=p;
                    indexCurrentTip++;
                }
               connectNodes(NULL, NULL, p);
                // do not do anything with this node because it is a tip
                //fprintf (stderr, "\n * do not do anything with this node because it is a tip");
            }
            else if (p->left != NULL && p->right == NULL && p->anc1 != NULL)
            {
                // this is a superflous node and can be removed(this superfluos nodes are the MRCA nodes of the demes
                foundSuperflousNode = YES;
                q = p->left;
                r = p->anc1;
                if (p->anc1->left == p)  // p->anc up, p->left down, total: up and down for left
                {
                    r->left = q;
                    q->anc1 = r;
                    p->left = NULL;
                    p->anc1 = NULL;
                    
                    connectNodes(q, r->right, r);
                }
                else
                {
                    r->right = q;
                    q->anc1 = r;
                    p->left = NULL;
                    p->anc1 = NULL;
                    connectNodes(r->left, q, r);
                }
               
                //fprintf (stderr, "\n - this is a superflous node and can be removed (1)");
            }
            else if (p->left == NULL && p->right != NULL && p->anc1 != NULL)
            {
                // this is a superflous node and can be removed
                foundSuperflousNode = YES;
                q = p->right;
                r = p->anc1;
                if (p->anc1->left == p)
                {
                    r->left = q;
                    q->anc1 = r;
                    p->right = NULL;
                    p->anc1 = NULL;
                    connectNodes(q, r->right, r);
                }
                else
                {
                    r->right = q;
                    q->anc1 = r;
                    p->right = NULL;
                    p->anc1 = NULL;
                    connectNodes(r->left, q, r);
                }
                
                //fprintf (stderr, "\n - this is a superflous node and can be removed (2)");
            }
            else if (p->left != NULL && p->right != NULL && p->anc1 != NULL)
            {
                connectNodes(p->left, p->right, p);
                // this is an internal node formed by a coalescence event, do not touch
                //fprintf (stderr, "\n * this is an internal node formed by a coalescence event, do not touch");
            }
            else if (p->left != NULL && p->right != NULL && p->anc1 == NULL)
            {
                connectNodes(p->left, p->right, p);
                // this is the last (coalescence event) in the tree, MRCA
                //fprintf (stderr, "\n * this is the last (coalescence event) in the tree, MRCA");
            }
            else if (p->left != NULL && p->right == NULL && p->anc1 == NULL)
            {
                // Seems to be the last coalescent event among sequences with non-ancestral material
                // it is not superfluous, we just remove it
                p->left->anc1 = NULL;
                //fprintf (stderr, "\n - this is a superflous node and can be removed (3)");
            }
            else if (p->left == NULL && p->right != NULL && p->anc1 == NULL)
            {
                // not clear what this node could be doing, but we will remove it anyway
                fprintf (stderr, "strange\n");
                p->left = NULL;
                p->right->anc1 = NULL;
                //fprintf (stderr, "\n - this is a superflous node and can be removed (4)");
            }
            else
            {
                fprintf (stderr, "You should not be here, I think\n");
                fprintf (stderr, "%d %d-- %d %d %d\n", Index(p), j, Index(p->left), Index(p->right), Index(p->anc1));
            }
            if (p->anc1 != NULL)
            {//update length field
                
                //   p->length = p->anc1->time- p->time;
                p->length = (p->anc1->timePUnits- p->timePUnits);
                //*mutationRate;
                p->lengthModelUnits = (p->anc1->time- p->time);
                //*mutationRate;
                setLength(p);
            }
        }
        //fprintf (stderr, "\n");
    }//while
    
    /* about the MRCA */
    *newInd=*nextAvailable-1;
    p = *nodes + *newInd; /* because the last one event is the last coalescence */
    //fprintf (stderr, "\n\n\n>> newInd = %d\n", newInd);
    
    if (programOptions->thereisOutgroup == NO)
    {
        p = *nodes + *newInd;
        p->class = 5;
        *treeRootInit = p;
        //treeRootInit[0] = p;
        p->anc1 = NULL;
    }
    if (programOptions->thereisOutgroup == YES && programOptions->outgroupSelection > 0)  /*** Root and outgroup ***/
    {
        p = *nodes + *newInd; // MRCA
        p->class = 4;
        
        if (programOptions->noisy > 1)
            fprintf (stderr, "\n\n>> Attaching outgroup .. ");
        
        //fprintf (stderr, "\n>> ThisCloneNumber = %d, ListMigrationTimesInitial[ThisCloneNumber] = %lf, ClonePopSizeMeffectBegin[ThisCloneNumber] = %lf \n", ThisCloneNumber, ListMigrationTimesInitial[ThisCloneNumber], ClonePopSizeMeffectBegin[ThisCloneNumber]);
        
        if (programOptions->outgroupSelection == 1)  /*** Root 2 times and outgroup ***/
            *currentTime = CurrentPop->timeOriginSTD; // origin of the clone; // currentTime + (outgroupBranchLength_Root1Root2 / mutationRate); // set time of the new root (from which the MRCA and outgroup nodes are derived)
        else if (programOptions->outgroupSelection == 2) { /*** Root 2 times and outgroup ***/
            *currentTime = CurrentPop->timeOriginSTD + (programOptions->outgroupBranchLength_Root1Root2 / CurrentPop->effectPopSize) ; // origin of the clone + time given by the user
            
        }
        else
        {
            fprintf (stderr, "\n\nError simulationg the outgroup. Check input settings\n");
            PrintUsage();
        }
       TreeNode*       healthyRoot = *nodes + *nextAvailable;
       healthyRoot->index = *nextAvailable;
        healthyRoot->label = *labelNodes;
        healthyRoot->effectPopSize= p->effectPopSize;
        *labelNodes=*labelNodes+1;
        healthyRoot->left = p;//coalTreeMRCA;
//        coalTreeMRCA->anc = healthyRoot;
        p->anc1 = healthyRoot;
       
        healthyRoot->timePUnits = p->timePUnits * healthyRoot->effectPopSize;
        healthyRoot->class = 5;
//        coalTreeMRCA->length = transformingBranchLength/mutationRate;
        p->length = 0;
        
//        coalTreeMRCA->branchLength = transformingBranchLength;
        p->lengthModelUnits = 0;
        
//        healthyRoot->time = currentTime +  transformingBranchLength/mutationRate;
        healthyRoot->time = *currentTime  ;
        
        int transformingBranchLength=1.001;
       // healthyRoot->time = p->time * transformingBranchLength ;
        healthyRoot->timePUnits = *currentTime * healthyRoot->effectPopSize;
        p->length = (p->anc1->timePUnits- p->timePUnits);
        //*mutationRate;
        p->lengthModelUnits = (p->anc1->time- p->time);
        //*mutationRate;
        
         healthyRoot->length = 0;
//        healthyRoot->length = 0;
    
//        if (noisy > 2)
//            fprintf (stderr, "DONE");
//
       (*nextAvailable)++;
//
//        /* connect the healthy ancestral cell with the tip healthy cell*/
//        if (noisy > 2)
//            fprintf (stderr, "\n>> Adding healthy tip ... ");
      TreeNode* healthyTip = *nodes + *nextAvailable;
        healthyTip->left = NULL;
        healthyTip->right = NULL;
         healthyTip->effectPopSize= healthyRoot->effectPopSize;
        
        connectNodes(NULL, NULL, healthyTip);
        
       healthyTip->anc1 = healthyRoot;
       healthyRoot->right = healthyTip;
        healthyTip->time = 0;
         healthyTip->timePUnits = 0;
        double  healthyTipBranchLengthRatio =1;
       // healthyTip->time = (healthyRoot->time - healthyTip->time) * healthyTipBranchLengthRatio;
        //healthyTip->timePUnits = (healthyRoot->timePUnits - healthyTip->timePUnits) * healthyTipBranchLengthRatio
   
       // healthyTip->timePUnits = (healthyRoot->time - healthyTip->time) * healthyTipBranchLengthRatio * healthyTip->effectPopSize;
        healthyTip->length = (healthyTip->anc1->timePUnits- healthyTip->timePUnits);
        //*mutationRate;
        healthyTip->lengthModelUnits = (healthyTip->anc1->time- healthyTip->time);
        //*mutationRate;
        //healthyTip->length = healthyTip->length * mutationRate;
        healthyTip->isOutgroup= YES;
        
        connectNodes(p, healthyTip, healthyRoot);
         setLength(p);
        setLength(healthyTip);
//      healthyTip->time = 0;
////        healthyTip->length = healthyTipBranchLength/mutationRate;
//         healthyTip->lengthModelUnits =p->lengthModelUnits;
//          healthyTip->length =p->length;
       
       // treeRootInit[0]=healthyRoot;// this is the older version
      //  treeRootInit[0]=healthyRoot;
        *treeRootInit=healthyRoot;
        
//        healthyTip->branchLength = healthyTipBranchLength;
//        healthyTip->isHealthyTip = YES;
//        healthyTip->index = nextAvailableNode;
        /* Root node: r; outgroup node: q */
//        r = *nodes + *nextAvailable; // root
//        r->index = *nextAvailable;
//        r->label = *labelNodes;
////        r->left=p;
//        *labelNodes=*labelNodes+1;
//        r->indexOldClone = r->indexCurrentClone = p->indexCurrentClone;
//        //r->indexCurrentClone = p->indexCurrentClone = 0;
//        r->indexCoalClone = r->indexCurrentClone = p->indexCurrentClone;
//        r->effectPopSize=CurrentPop->effectPopSize;
//        nextAvailable++; /* 1 node more is available */
//        r->class = 5;
//        r->left = p;
//        r->anc1 = NULL;
//        r->length = 0;
//          r->lengthModelUnits = 0;
//        r->time = *currentTime;
//        r->timePUnits = *currentTime * CurrentPop->effectPopSize;
//        treeRootInit[0] = r;
//        p->anc1 = r;
        
        /* q is the outgroup node */
//        q = *nodes + *nextAvailable;
//        r->right = NULL; // q??, let's implement it as outgroup, not as right, but be carefull with superfluos nodes!!
//        q->left = NULL;
//        q->right = NULL;
//        q->anc1 = treeRootInit[0];
//        q->time = 0;
//        q->timePUnits = 0;
//       q->length = outgroupBranchLength_RootSample / mutationRate;
//        q->length = treeRootInit[0]->timePUnits  / mutationRate;
//        q->isOutgroup = YES;
//        q->index = *nextAvailable + 1;
        /* p->label = p->index; */
//        q->label = TotalNumSequences + 1;
//        treeRootInit[0]->outgroup = q;
        //fprintf (stderr, "\n>> r->time = %lf, r->timePUnits = %lf\n", r->time, r->timePUnits);
    }
    
//    *counterTimeInit = *counterTimeInit + *currentTime;

    //fprintf ( stderr, "\n numCA = %d\n", *numCA);
    //fprintf ( stderr, "\n numMIG = %d\n", *numMIG);
    int intLabel;
    if (programOptions->noisy > 1)
        fprintf (stderr, "\n\n>> Relabeling nodes on tree... \n\n");
    if (programOptions->thereisOutgroup == YES)
        intLabel = programOptions->TotalNumSequences + 1;
    else
        intLabel = programOptions->TotalNumSequences;
   // intLabel = 0;
   // RelabelNodes(treeRootInit[0], treeRootInit, &intLabel );
    RelabelNodes(*treeRootInit, &treeRootInit, &intLabel );
    
}

/********************************** AddGermlineVariation ***********************************/
/* Introduce germline SNPs in the healthy root genome for DNA or 0/1 data  */
void AddGermlineVariation (TreeNode *treeRoot,long int *seed, int numSites, double SNPrate, SiteStr* allSites, int alphabet, int*** data,  int HEALTHY_ROOT, double cumMij[4][4])
{
    int        j, site, genome, current_state;
    double    ran;
    int cell=-1;
    for (site=0; site<numSites; site++)
    {
        if(RandomUniform(seed) < SNPrate) //there is a SNP at this site
        {
            allSites[site].isSNP = YES;
          if (treeRoot != NULL)
          {
                if (treeRoot->isOutgroup  == YES || treeRoot->anc1 == NULL)
                {
                     cell = treeRoot->label;
                }
            
            if(RandomUniform(seed) < 0.5)  // in maternal/paternal genome
                genome = MATERNAL;
            else
                genome = PATERNAL;
            
            if (alphabet == BINARY){
                
                //data[genome][HEALTHY_ROOT][site] = 1;
                data[genome][cell][site] = 1;
                
            }
            
            else // DNA mutation according to Mij matrix
            {
                //current_state = data[genome][HEALTHY_ROOT][site];
                current_state = data[genome][cell][site];
                ran = RandomUniform(seed) * cumMij[current_state][3];
                for (j=0; j<4; j++)
                {
                    if (ran <= cumMij[current_state][j])
                    {
                        //data[genome][HEALTHY_ROOT][site] = j;
                        data[genome][cell][site] = j;
                        break;
                    }
                }
            }
          }// if (treeRoot != NULL)
            
        }
    }
}
/********************* RandomPoisson ********************/
/* Generates a random number from a Poisson distibution with
 mean lambda.
 */

int RandomPoisson (double lambda, long int *seed)
{
    int        poissonNumber;
    double    sum;
    
    sum = 0;
    poissonNumber = -1;
    
    while (sum <= 1.0)
    {
        sum += RandomExponential (lambda, seed);
        poissonNumber++;
    }
    
    return poissonNumber;
}


/*************************** SimulateISM **********************************/
/*  Simulates mutations under an infinite diploid/haploid sites model (ISM).
 
 Diploid here means that if site 33 mutates in the female genome,
 site 33 in the male genome will not mutate, and viceversa.
 
 Haploid here means that site 33 in the maternal genome and site 33 in the paternal
 genome will be considered different sites, so both can mutate. The rest is as
 in the diploid, standard ISM
 
 The number of mutations will be distributed as:
 
 (1) a Poisson with parameter the sum of
 branch lengths over all sites  Times are already scaled in
 2N. We will force it to get just one mutation per site
 
 (2) A used-defined number of mutations (segregating sites in the ISM diploid)
 
 Approach 1, which conditions on theta and takes in account
 the branch lengths, is preferred than Hudson's approach 2, which selects per each
 site a random branch and puts a mutation there, and therefore it is
 conditioned on the number of segregating sites. With the coalescent
 you get large or short trees, which imply different number of mutations
 but if you do Hudson's approach you always put the same number of mutations.
 */

void SimulateISM (TreeNode *treeRoot, int genome, int doISMhaploid, long int *seed,  int *DefaultModelSites, int numDefaultModelSites, int* AltModelSites, int numAltModelSites, double totalTreeLength , int *numISMmutations, int numFixedMutations, int numSNVmaternal, int doSimulateFixedNumMutations,  int alphabet,  int ***data, SiteStr* allSites, int  *numMU, double cumMij[4][4], double mutationRate)
{
    int   i, trials, numMutations, mutationsSoFar;
    double  totalBranchSum;
    int   *modelSites, numModelSites;
    double cumBranchLength =0;
    double uniform =0;
    int mutationAdded;
    double ran=0;
//    if (doISMhaploid == NO)
//    {
        modelSites = DefaultModelSites;
        numModelSites = numDefaultModelSites;
//    }
//    else
//    {
//        modelSites = AltModelSites;
//        numModelSites = numAltModelSites;
//    }
    
   
    
    //totalBranchSum = totalTreeLength * numModelSites;//this only works if totalTreeLength <1(time is measured in 2*N generations
//    totalBranchSum = totalTreeLength * mutationRate*numModelSites;
    totalBranchSum = totalTreeLength *numModelSites;
//    if (doSimulateFixedNumMutations == YES) /* conditioned on a specific number of mutations */
//    {
////        if (genome == MATERNAL)
////            numMutations = numSNVmaternal =  RandomBinomial (0.5, numFixedMutations, seed);
////        else
////            numMutations = numFixedMutations - numSNVmaternal;
//    }
   // else
   // {
    int ploidyFactor;
    
    if (doISMhaploid==NO){
        ploidyFactor=1;
    }
    else{
         ploidyFactor=2;
    }
    
    do{
         numMutations = RandomPoisson (totalBranchSum, seed);
    }while( *numISMmutations + numMutations > (ploidyFactor *numModelSites));
  
       // fprintf(stderr, "num mutations  %d\n",numMutations );
    
  //  }/* the number of mutations will be distributed as a Poisson with parameter totalBranchSum */
    
    /* if the number of ISM mutations is bigger than the number of ISM sites quit the program and warn the user about violation of ISM */
//    if (doISMhaploid == NO && *numISMmutations + numMutations > ( numModelSites))
//    {
//        fprintf (stderr, "\n\nERROR: The diploid infinite sites model (ISM) has been violated. There will be");
//        fprintf (stderr, "\nmore mutations (%d existing + %d proposed]) than available sites under this model (%d).", *numISMmutations, numMutations, numDefaultModelSites);
//        fprintf (stderr, "\nTry using a smaller mutation rate or effective population size");
//        fprintf (stderr, "\n(i.e. smaller theta) or check the proportion of JC sites.");
//        fprintf (stderr, "\nIf using an user tree, try with smaller branch lengths (including the transforming and healthy branches).\n");
//        exit (-1);
//    }
//    else if (doISMhaploid == YES && *numISMmutations + numMutations > (2*numModelSites))
//    {
//        fprintf (stderr, "\n\nERROR: The haploid infinite sites model (ISM) has been violated. There will be");
//        fprintf (stderr, "\nmore mutations (%d existing + %d proposed]) than available sites under this model (%d).", *numISMmutations, numMutations, 2*numModelSites);
//        fprintf (stderr, "\nTry using a smaller mutation rate or effective population size");
//        fprintf (stderr, "\n(i.e. smaller theta) or check the proportion of JC sites.");
//        fprintf (stderr, "\nIf using an user tree, try with smaller branch lengths (including the transforming and healthy branches).\n");
//        exit (-1);
//    }
//
    *numISMmutations = *numISMmutations + numMutations;
    trials = 0;
    mutationsSoFar = 0;
    mutationAdded=NO;
    while (mutationsSoFar < numMutations)
    {
        i = RandomUniformTo(numModelSites, seed);
        
        /* choose a site at random */
//        while (((doISMhaploid == NO)  &&(allSites[i].numMutations != 0))
//               ||
//               ((doISMhaploid == YES) && (genome == MATERNAL) && (allSites[i].numMutationsMaternal != 0)) ||
//              ((doISMhaploid == YES) && (genome == PATERNAL) && (allSites[i].numMutationsPaternal != 0))
//               )
        while (((doISMhaploid == NO)  && (allSites[i].numMutations != 0))  ||
               ((doISMhaploid == YES) && (genome == MATERNAL) && (allSites[i].numMutationsMaternal != 0)) ||
               ((doISMhaploid == YES) && (genome == PATERNAL) && (allSites[i].numMutationsPaternal != 0)))
        {
//            fprintf (stderr, "\n\n site %d again!", i);
            i = RandomUniformTo(numModelSites, seed);
            if (trials++ > 1000*numModelSites)
            {
                fprintf (stderr, "\n\n ERROR: after %d trials cannot find an unmuted site",1000*numModelSites);
                fprintf (stderr, "\nmutations = %d   mutations so far = %d\n\n",numMutations, mutationsSoFar);
                for (i=0; i<numModelSites; i++)
                    fprintf (stderr, "\nsite %d has %d mutations",i,allSites[i].numMutations);
                exit(-1);
            }
        }
        if (alphabet == DNA)
            SimulateISMDNAforSite (treeRoot, genome, i, doISMhaploid, seed, totalTreeLength, data, allSites, numMU,cumMij, mutationRate,&uniform, &cumBranchLength,  &ran);
        else{
            mutationAdded=NO;
            SimulateISMforSite (treeRoot, genome, i, doISMhaploid, seed, totalTreeLength, data, allSites, numMU,cumMij,mutationRate ,    &cumBranchLength,  &uniform,  &mutationAdded);
        }
        mutationsSoFar++;
#ifdef MYDEBUG
        fprintf (stderr, "\nmutations = %d   mutations so far = %d\n",numMutations, mutationsSoFar);
        fprintf (stderr, "\n position = %d ",i);
        if (allSites[i].numMutations > 1)
        {
            fprintf (stderr, "\n\n ERROR: %d mutations in site %d",allSites[i].numMutations, i);
            exit(-1);
        }
        for (i=0; i<numDefaultModelSites; i++)
        {
          if (allSites[i].numMutations > 0)
            fprintf (stderr, "%2d[%d,%d,%d ] ",i, allSites[i].numMutations, allSites[i].numMutationsMaternal, allSites[i].numMutationsPaternal);
        }
#endif
        
    }
 }

///******************** ReadParametersFromCommandLine *************************/
//static void ReadParametersFromCommandLine (int argc,char **argv)
//    {
//    int   i, z, j;
//    char  flagb;
//    float argument;
//
//    z = 0;
//
//    for (i=1; i<argc; i++)
//        {
//        argv[i]++;
//        flagb=*argv[i];
//        argv[i]++;
//        argument = -9999;
//
//        /* Used: N X C R D M O T K Y # */
//
//        switch (toupper(flagb))
//            {
//
//            case 'N':
//                argument = atof(argv[i]);
//                numDataSets = (int) argument;
//                if (numDataSets <1)
//                    {
//                    fprintf (stderr, "PARAMETER ERROR: Bad number of replicates (%d)\n\n", numDataSets);
//                    PrintUsage();
//                    }
//                break;
//
//
//            case '#':
//                argument = atof(argv[i]);
//                userSeed = (int) argument;
//                if (userSeed < 0)
//                    {
//                    fprintf (stderr, "PARAMETER ERROR: Bad user seed for random number (%lu bytes)\n\n", userSeed);
//                    PrintUsage();
//                    }
//                break;
//
//
//            case 'X':
//                argument = atof(argv[i]);
//                Nscaling = (int) argument;
//                if (Nscaling < 1 || Nscaling > 2)
//                    {
//                    fprintf (stderr, "PARAMETER ERROR: Bad haplid/diploid chosen (x) (%d)\n\n", Nscaling);
//                    PrintUsage();
//                    }
//                break;
//
//
//            case 'Y':
//                argument = atof(argv[i]);
//                noisy = (int) argument;
//                if (noisy < 0)
//                    {
//                    fprintf (stderr, "PARAMETER ERROR: Bad noisy value (%d)\n\n", noisy);
//                    PrintUsage();
//                    }
//                break;
//
//
//            //[outgroup (optional), branch length (root-sample) branch length (root-mrca)] o0.04 0.06
//            case 'O':
//                outgroupSelection = atof(argv[i]);
//                if (outgroupSelection < 0)
//                    {
//                    fprintf (stderr, "PARAMETER ERROR: Bad selection for outgroup %d (0: No outgroup, 1 outgroup with one branch, 2 outgroup with two branches\n\n", outgroupSelection);
//                    PrintUsage();
//                    }
//                if (outgroupSelection != 0 && outgroupSelection != 1 && outgroupSelection != 2)
//                    {
//                    fprintf (stderr, "PARAMETER ERROR: Bad selection for outgroup %d (0: No outgroup, 1 outgroup with one branch, 2 outgroup with two branches\n\n", outgroupSelection);
//                    PrintUsage();
//                    }
//
//
//                if (outgroupSelection == 0)
//                    {
//                    thereisOutgroup = NO;
//                    outgroupBranchLength_RootSample = 0.0;
//                    outgroupBranchLength_Root1Root2 = 0.0;
//                    }
//                else if (outgroupSelection == 1)
//                    {
//                    thereisOutgroup = YES;
//                    outgroupBranchLength_Root1Root2 = 0.0;
//
//                    outgroupBranchLength_RootSample = atof(argv[++i]);
//                    if (outgroupBranchLength_RootSample < 0)
//                        {
//                        fprintf(stderr, "PARAMETER ERROR: Bad outgroup branch length (Root-Sample) value (%f)\n\n", outgroupBranchLength_RootSample);
//                        PrintUsage();
//                        }
//                    }
//                else if (outgroupSelection == 2)
//                    {
//                    thereisOutgroup = YES;
//
//                    outgroupBranchLength_Root1Root2 = atof(argv[++i]);
//                    if (outgroupBranchLength_RootSample < 0)
//                        {
//                        fprintf(stderr, "PARAMETER ERROR: Bad outgroup branch length (Root1-Root2) value (%f)\n\n", outgroupBranchLength_Root1Root2);
//                        PrintUsage();
//                        }
//                    outgroupBranchLength_RootSample = atof(argv[++i]);
//                    if (outgroupBranchLength_RootSample < 0)
//                        {
//                        fprintf(stderr, "PARAMETER ERROR: Bad outgroup branch length (Root2-Sample) value (%f)\n\n", outgroupBranchLength_RootSample);
//                        PrintUsage();
//                        }
//                    }
//                else
//                    {
//                    fprintf (stderr, "PARAMETER ERROR: Bad selection for outgroup %d (0: No outgroup, 1 outgroup with one branch, 2 outgroup with two branches\n\n", outgroupSelection);
//                    PrintUsage();
//                    }
//                //outgroupBranchLength_RootSample = 0;
//            break;
//
//
//            case 'C':
//                argument = atof(argv[i]);
//                numClones = (int) argument;
//                if (numClones < 1)
//                    {
//                    fprintf (stderr, "PARAMETER ERROR: Bad number of clones (c) (%d)\n\n", numClones);
//                    PrintUsage();
//                    }
//
//
//                CloneNameBegin = (int *) calloc(numClones+1,(long) sizeof(int));
//                CloneSampleSizeBegin = (int *) calloc(numClones+1,(long) sizeof(int));
//                ClonePopSizeBegin = (int *) calloc(numClones+1,(long) sizeof(int));
//                CloneBirthRateBegin = (double *) calloc(numClones+1,(long) sizeof(double));
//                CloneDeathRateBegin = (double *) calloc(numClones+1,(long) sizeof(double));
//                CloneTimeOriginInput = (double *) calloc(numClones+1,(long) sizeof(double));
//                if (CloneNameBegin == NULL || CloneSampleSizeBegin == NULL || ClonePopSizeBegin == NULL || CloneBirthRateBegin == NULL || CloneDeathRateBegin == NULL || CloneTimeOriginInput == NULL)
//                    {
//                    fprintf (stderr, "PARAMETER ERROR: Could not allocate variables for clones\n");
//                    exit (1);
//                    }
//
//                for (j=1; j<=numClones; j++)
//                    {
//                    for (z=1; z<=6; z++)
//                        {
//                        if (z == 1)
//                            {
//                            argument = atof(argv[++i]);
//                            CloneNameBegin[j] = (int) argument;
//                            if (CloneNameBegin[j] <= 0 || CloneNameBegin[j] > numClones)
//                                {
//                                fprintf (stderr, "PARAMETER ERROR: Bad number for clone %d (should be higher than 0 and lower than the number of clones %d) (%d)\n\n", j, numClones, CloneNameBegin[j]);
//                                PrintUsage();
//                                }
//                            }
//                        if (z == 2)
//                            {
//                            argument = atof(argv[++i]);
//                            CloneSampleSizeBegin[j] = (int) argument;
//                            if (CloneSampleSizeBegin[j] < 0)
//                                {
//                                fprintf (stderr, "PARAMETER ERROR: Bad sample size for clone %d (should not be negative) (%d)\n\n", j, CloneSampleSizeBegin[j]);
//                                PrintUsage();
//                                }
//                            }
//                        if (z == 3)
//                            {
//                            argument = atof(argv[++i]);
//                            ClonePopSizeBegin[j] = (int) argument;
//                            if (ClonePopSizeBegin[j] < 0)
//                                {
//                                fprintf (stderr, "PARAMETER ERROR: Bad population size for clone %d (should be higher than 0) (%d)\n\n", j, ClonePopSizeBegin[j]);
//                                PrintUsage();
//                                }
//                            }
//                        if (z == 4)
//                            {
//                            argument = atof(argv[++i]);
//                            CloneBirthRateBegin[j] = (double) argument;
//                            }
//                        if (z == 5)
//                            {
//                            argument = atof(argv[++i]);
//                            CloneDeathRateBegin[j] = (double) argument;
//                            }
//                        if (z == 6)
//                            {
//                            argument = atof(argv[++i]);
//                            CloneTimeOriginInput[j] = (double) argument;
//                            if (CloneTimeOriginInput[j] < 0)
//                                {
//                                fprintf (stderr, "PARAMETER ERROR: Bad time to origin for clone %d (should not be negative) (%lf)\n\n", j, CloneTimeOriginInput[j]);
//                                PrintUsage();
//                                }
//                            }
//                        }
//                    }
//                break;
//
//
//            case 'U':
//                mutationRate = atof(argv[i]);
//                if (mutationRate < 0)
//                    {
//                    fprintf (stderr, "PARAMETER ERROR: Bad mutation rate (%f)\n\n", mutationRate);
//                    PrintUsage();
//                    }
//                break;
//
//
//            case 'J':
//                strcpy(treeFile, "trees");
//                doPrintTrees = YES;
//                break;
//
//
//            case 'K':
//                strcpy(timesFile, "times");
//                doPrintTimes = YES;
//                break;
//
//
//            case '?':
//                PrintUsage();
//                break;
//
//
//            default :
//                fprintf(stderr, "PARAMETER ERROR: Incorrect parameter: %c\n\n", flagb);
//                PrintUsage();
//                break;
//            }
//        }
//    }



/********************* RandomBinomial ********************/
/* Generates a random number from a Binomial distibution with sucess probabilty p
 and n trials using the direct method (sum of Bernoulli variables).
 */

int RandomBinomial (double prob, int numTrials, long int *seed)
{
    int i, sum;
    sum = 0;
    
    for(i=0; i<numTrials; i++)
    {
        if(RandomUniform(seed) < prob)
            sum++;
    }
    return sum;
}

/********************** RandomUniformTo ****************************/
/* it returns random uniform in range 0...max-1          */

static int RandomUniformTo (int max, long int *seed)
{
    double    rd;
    rd = RandomUniform (seed);
    return (floor(rd*max));
}
/********************************** SimulateISMDNAforSite ***********************************/
/*    Simulates a ACGT mutation under an infinite sites model (ISM) for a given site. The branch
 where this mutation is placed is chosen according to its length.
 The reference (healthy) allele for each site will be determined by the nucleotide frequencies
 A mutation matrix will define the probability of changing to other nucleotides given
 the healthy alelle chose
 */
void SimulateISMDNAforSite (TreeNode *p, int genome, int site, int doISMhaploid, long int *seed, double totalTreeLength, int ***data, SiteStr* allSites, int  *numMU, double cumMij[4][4],double mutationRate, double *uniform, double *cumBranchLength, double* ran )
{
//    static double    cumBranchLength, uniform, ran;
    int             j, cell, anccell, ancstate;
    
    if (p != NULL)
    {
        cell = p->label;
        
        if ( p->anc1 == NULL)
        {
            *cumBranchLength = 0;
            *uniform = RandomUniform(seed) * totalTreeLength;
        }
        else
        {
            anccell = p->anc1->label;
             if(genome == MATERNAL )
                 ancstate = p->anc1->maternalSequence[site];
            else
                ancstate = p->anc1->paternalSequence[site];
//            cumBranchLength += p->length;// ->branchLength;
            *cumBranchLength =*cumBranchLength+ p->length;
          
            if ((*cumBranchLength < *uniform) || /* => there will be no change */
                ((doISMhaploid == NO)  && (allSites[site].numMutations > 0))  ||
                ((doISMhaploid == YES) && (genome == MATERNAL) && (allSites[site].numMutationsMaternal > 0)) ||
                ((doISMhaploid == YES) && (genome == PATERNAL) && (allSites[site].numMutationsPaternal > 0)))
            {
                if(genome == MATERNAL )
                    p->maternalSequence[site]=p->anc1->maternalSequence[site];
                else
                    p->paternalSequence[site]=p->anc1->paternalSequence[site];
                
                p->numbersMaternalMutationsPerSite[site]=p->anc1->numbersMaternalMutationsPerSite[site];
                p->numbersPaternalMutationsPerSite[site]=p->anc1->numbersPaternalMutationsPerSite[site];
                p->numbersMutationsUnderSubtreePerSite[site]=p->anc1->numbersMutationsUnderSubtreePerSite[site];
            }
            else /* => there will be change */
            {
                *ran = RandomUniform(seed) * cumMij[ancstate][3];
                for (j=0; j<4; j++)
                {
                    if (*ran <= cumMij[ancstate][j])
                    {
                        //data[genome][cell][site] = j;
                        if(genome == MATERNAL )
                            p->maternalSequence[site]=j;
                        else
                            p->paternalSequence[site]=j;
                        break;
                    }
                }
                if (genome == MATERNAL){
                    allSites[site].numMutationsMaternal++;
                    p->numbersMaternalMutationsPerSite[site]=p->anc1->numbersMaternalMutationsPerSite[site]+1;
                    p->numbersPaternalMutationsPerSite[site]=p->anc1->numbersPaternalMutationsPerSite[site];
                    p->numbersMutationsUnderSubtreePerSite[site]=p->anc1->numbersMutationsUnderSubtreePerSite[site]+1;
                }
                else if (genome == PATERNAL){
                    allSites[site].numMutationsPaternal++;
                    p->numbersMaternalMutationsPerSite[site]=p->anc1->numbersMaternalMutationsPerSite[site];
                    p->numbersPaternalMutationsPerSite[site]=p->anc1->numbersPaternalMutationsPerSite[site]+1;
                    p->numbersMutationsUnderSubtreePerSite[site]=p->anc1->numbersMutationsUnderSubtreePerSite[site]+1;
                }
                allSites[site].numMutations++;
                (*numMU)++;
            }
        }
        SimulateISMDNAforSite (p->left, genome, site, doISMhaploid, seed,  totalTreeLength, data, allSites, numMU,cumMij,  mutationRate, uniform, cumBranchLength,  ran);
        SimulateISMDNAforSite (p->right, genome, site, doISMhaploid, seed, totalTreeLength, data, allSites, numMU,cumMij,  mutationRate, uniform, cumBranchLength,  ran);
    }
}
/********************************** SimulateISMForSite ***********************************/
/*    Simulates a 0/1 mutation under an infinite sites model (ISM) for a given site. The branch
 where this mutation is placed is chosen according to its length.
 0 is the reference (healthy) allele
 */
void SimulateISMforSite (TreeNode *p, int genome, int site, int doISMhaploid, long int *seed, double totalTreeLength, int ***data, SiteStr* allSites, int  *numMU, double cumMij[4][4], double mutationRate, double*    cumBranchLength, double* uniform, int* mutationAdded)
{
//    static double    cumBranchLength, uniform;
//    static double    cumBranchLength, uniform;
    int             cell, anccell;
//    if (*mutationAdded==YES)
//        return;
    if (p != NULL)
    {
        cell = p->label;
        
//        if (p->isOutgroup  == YES || p->anc1 == NULL)
         if ( p->anc1 == NULL)
        {
            *cumBranchLength = 0;
//            double rUniform=RandomUniform(seed) * totalTreeLength;
            *uniform = RandomUniform(seed) * totalTreeLength;
            
        }
        else
        {
            anccell = p->anc1->label;
//            *cumBranchLength =*cumBranchLength+ p->length;// ->branchLength;
             // *cumBranchLength =*cumBranchLength+ p->lengthModelUnits;// ->branchLength;
            *cumBranchLength =*cumBranchLength+ p->length;
//            if ((*cumBranchLength < *uniform) || /* => there will be no change */
//                ((doISMhaploid == NO)  && (allSites[site].numMutations > 0))
//               ||
//               ((doISMhaploid == YES) && (genome == MATERNAL) && (allSites[site].numMutationsMaternal > 0)) ||
//              ((doISMhaploid == YES) && (genome == PATERNAL) && (allSites[site].numMutationsPaternal > 0))
//                )
            if ((*cumBranchLength < *uniform) ||// (*mutationAdded==YES)/* => there will be no change */
                ((doISMhaploid == NO)  && (allSites[site].numMutations > 0))  ||
                ((doISMhaploid == YES) && (genome == MATERNAL) && (allSites[site].numMutationsMaternal > 0)) ||
                ((doISMhaploid == YES) && (genome == PATERNAL) && (allSites[site].numMutationsPaternal > 0)))
            {
                if(genome == MATERNAL )
                    p->maternalSequence[site]=p->anc1->maternalSequence[site];
                else
                     p->paternalSequence[site]=p->anc1->paternalSequence[site];
          
                p->numbersMaternalMutationsPerSite[site]=p->anc1->numbersMaternalMutationsPerSite[site];
                 p->numbersPaternalMutationsPerSite[site]=p->anc1->numbersPaternalMutationsPerSite[site];
                p->numbersMutationsUnderSubtreePerSite[site]=p->anc1->numbersMutationsUnderSubtreePerSite[site];
            }
            else /* => there will be change */
            {
                //if (data[genome][anccell][site] == 0)  /* checking all this might be excessive */
                if(genome == MATERNAL && p->anc1->maternalSequence[site]==0 )
                {
                    p->maternalSequence[site]=1;
                    p->numbersMaternalMutationsPerSite[site]=p->anc1->numbersMaternalMutationsPerSite[site]+1;
                     p->numbersPaternalMutationsPerSite[site]=p->anc1->numbersPaternalMutationsPerSite[site];
                       p->numbersMutationsUnderSubtreePerSite[site]=p->anc1->numbersMutationsUnderSubtreePerSite[site]+1;
                    *numMU=*numMU+1;
                    *mutationAdded=YES;
                }
                 else if(genome == PATERNAL && p->anc1->paternalSequence[site]==0 )
                {
                    p->paternalSequence[site]=1;
                    p->numbersPaternalMutationsPerSite[site]=p->anc1->numbersPaternalMutationsPerSite[site]+1;
                    p->numbersMaternalMutationsPerSite[site]=p->anc1->numbersMaternalMutationsPerSite[site];
                    p->numbersMutationsUnderSubtreePerSite[site]=p->anc1->numbersMutationsUnderSubtreePerSite[site]+1;
                    *numMU=*numMU+1;
                    *mutationAdded=YES;
                }
//                if (data[genome][anccell][site] == 0)
//                {
//                    data[genome][cell][site] = 1;
//                  if (genome == MATERNAL)
//                    { allSites[site].numMutationsMaternal++;
//                       }
//                    else // (genome == PATERNAL)
//                 { allSites[site].numMutationsPaternal++;
//                        }
//                    allSites[site].numMutations++;
//                    *numMU=*numMU+1;
//                    *mutationAdded=YES;
//                    return;
//
//
//
//                }
                else if (genome == MATERNAL && p->anc1->maternalSequence[site]==1)
                {
                    
//                    data[genome][cell][site] = 0;
                    fprintf(stderr,"\n\nERROR: anccell=%d, doISMhaploid:%d , site %d in genome %d of cell %d cannot mutate twice under the ISM model",anccell, doISMhaploid ,site, genome, anccell);
                    fprintf (stderr, "\n\nnumber of mutations %d  maternal: %d and paternal:  %d ", p->numbersMutationsUnderSubtreePerSite[site], p->numbersMaternalMutationsPerSite[site], p->numbersPaternalMutationsPerSite[site]);
                 exit(-1);
                }
                else if (genome == PATERNAL && p->anc1->paternalSequence[site]==1)
                {
                    
                    //                    data[genome][cell][site] = 0;
                    fprintf(stderr,"\n\nERROR: anccell=%d, doISMhaploid:%d , site %d in genome %d of cell %d cannot mutate twice under the ISM model",anccell, doISMhaploid ,site, genome, anccell);
                    fprintf (stderr, "\n\nnumber of mutations %d  maternal: %d and paternal:  %d ", p->numbersMutationsUnderSubtreePerSite[site], p->numbersMaternalMutationsPerSite[site], p->numbersPaternalMutationsPerSite[site]);
                    exit(-1);
                }
                else
                {
                    fprintf (stderr, "\n\nERROR: site %d in genome %d of cell %d has an unknow state %d, %d under the ISM model", site, genome, anccell, p->anc1->maternalSequence[site], p->anc1->paternalSequence[site]);
                    exit(-1);
                }
                if (genome == MATERNAL)
                    allSites[site].numMutationsMaternal++;
                else if (genome == PATERNAL)
                    allSites[site].numMutationsPaternal++;
                allSites[site].numMutations++;
                
               
            }
        }
//        double goLeftFirst=RandomUniform(seed);
//        if(goLeftFirst < 0.5)
//        {
            SimulateISMforSite (p->left, genome, site, doISMhaploid, seed, totalTreeLength, data, allSites, numMU, cumMij, mutationRate, cumBranchLength, uniform, mutationAdded);
            SimulateISMforSite (p->right, genome, site, doISMhaploid, seed,totalTreeLength, data, allSites, numMU,cumMij, mutationRate, cumBranchLength,  uniform, mutationAdded);
            
//        }
//        else{
//            SimulateISMforSite (p->right, genome, site, doISMhaploid, seed, totalTreeLength, data, allSites, numMU, cumMij, mutationRate, cumBranchLength, uniform, mutationAdded);
//            SimulateISMforSite (p->left, genome, site, doISMhaploid, seed,totalTreeLength, data, allSites, numMU,cumMij, mutationRate, cumBranchLength,  uniform, mutationAdded);
//
//
//        }
    }
}

/* Returns the sum of the branch lengths for a given tree */
double SumBranches (TreeNode *p, double mutationRate)
{
    static double sum;
    
    if (p != NULL)
    {
        if (p->anc1 == NULL)
            sum = 0;
        else{
            //sum += (p->anc1->time- p->time)* mutationRate;//p->lengthModelUnits;
            sum += p->length;
//            sum += p->length;
        }
//            sum += p->lengthModelUnits;//length;
        SumBranches (p->left,  mutationRate);
        SumBranches (p->right,   mutationRate);
    }
    
    return sum;
}

/********************************** EvolveSitesOnTree ***********************************/
/* Evolves all sites (maternal and paternal genomes) on the given tree
 We assume that a site will be ISM or Mk in both maternal and paternal genome
 */
void EvolveSitesOnTree (TreeNode *treeRoot, int genome, long int *seed, int rateVarAmongSites, int numSites, SiteStr* allSites, int doGeneticSignatures, int alphaSites, int  propAltModelSites , int numDefaultModelSites, int numAltModelSites, int* DefaultModelSites, int* AltModelSites,  double totalTreeLength , int *numISMmutations, int numFixedMutations, int numSNVmaternal, int doSimulateFixedNumMutations,  int alphabet,  int ***data,  int  *numMU, double cumMij[4][4], int altModel, double altModelMutationRate, int doUserTree,int doJC, int doHKY, int doGTR, int doGTnR, double freqR, double freqY, double freqAG, double freqCT, double titv, double freq[4], double Mij[4][4],  double Root[], double Cijk[])
{
    int i = 0;

//    numDefaultModelSites = numSites;
//    numAltModelSites = 0;
//    for (i=0; i<numSites; i++)
//        DefaultModelSites[i] = i;
    
//    SimulateISM (treeRoot, genome, NO, seed,
//                 DefaultModelSites,numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet,   data, allSites, numMU, cumMij,altModelMutationRate
//                 );
    
    if (rateVarAmongSites == YES)
       for (i=0; i<numSites; i++)
           allSites[i].rateMultiplier = RandomGamma (alphaSites, seed) / alphaSites;
//
//    if (doGeneticSignatures == YES)
//    {
//   //     SimulateSignatureISM(treeRoot, genome, seed); /* we assume diploid ISM for genetic signatures */
//    }
//    else
//    {
//      if (propAltModelSites == 0)  /* only default model (ISM diploid) sites */
//        {
//            numDefaultModelSites = numSites;
//            numAltModelSites = 0;
//            for (i=0; i<numSites; i++)
//                DefaultModelSites[i] = i;
//            SimulateISM (treeRoot, genome, NO, seed,
//                          DefaultModelSites,numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet,   data, allSites, numMU, cumMij,altModelMutationRate
//                         );
//        }
//        else if (propAltModelSites == 1.0)  /* only alternative model sites */
//        {
//            numDefaultModelSites = 0;
//            numAltModelSites = numSites;
//            for (i=0; i<numSites; i++)
//                AltModelSites[i] = i;
//            if (altModel == ISMhap)
//                fprintf(stderr, "pop");
////                SimulateISM (treeRoot, genome, YES, seed,  DefaultModelSites,numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet,   data, allSites, numMU, cumMij,  altModelMutationRate);
//            else if (altModel == Mk2)
//                SimulateMk2 (treeRoot, genome, seed, AltModelSites,   numAltModelSites, doUserTree, rateVarAmongSites,  altModelMutationRate, allSites, data, numMU);
//            else if (altModel == finiteDNA)
//                SimulateFiniteDNA (treeRoot, genome, seed, doJC, doHKY,  doGTR, doGTnR,  freqR, freqY,  freqAG,  freqCT,  titv,  freq,  Mij, numAltModelSites);
//            else
//            {
//                fprintf (stderr, "\n\nERROR: Sorry, the specified model is unknown for me");
//                exit(0);
//            }
//        }
//        else /* both ISM and non-ISM sites */
//        {
//            numDefaultModelSites = 0;
//            numAltModelSites = 0;
//
//            for (i=0; i<numSites; i++)
//            {
//                if (RandomUniform (seed) < propAltModelSites)
//                    AltModelSites[numAltModelSites++] = i;
//                else
//                    DefaultModelSites[numDefaultModelSites++] = i;
//            }
    if (propAltModelSites == 0)/* only default model (ISM diploid) sites */
    {
        numDefaultModelSites = numSites;
        numAltModelSites = 0;
        for (i=0; i<numSites; i++)
            DefaultModelSites[i] = i;
        SimulateISM (treeRoot, genome, NO, seed,
                     DefaultModelSites,numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet,   data, allSites, numMU, cumMij,altModelMutationRate
                     );
        
        
    }
    else if (propAltModelSites == 1)
    {
        numDefaultModelSites = 0;
        numAltModelSites = numSites;
        for (i=0; i<numSites; i++)
            AltModelSites[i] = i;
        
        if (altModel == ISMhap)
        { fprintf(stderr, "only non ISM sites");
            SimulateISM (treeRoot, genome, YES, seed,DefaultModelSites,numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet,   data, allSites, numMU, cumMij,    altModelMutationRate);
        }
        else if (altModel == Mk2)
        { SimulateMk2 (treeRoot, genome, seed, AltModelSites,   numAltModelSites, doUserTree, rateVarAmongSites,  altModelMutationRate, allSites, data, numMU);
        }
        else if (altModel == finiteDNA)
        {
            SimulateFiniteDNA (treeRoot, genome, seed,doJC,  doHKY, doGTR,  doGTnR,  freqR, freqY,  freqAG,  freqCT, titv,  freq, Mij, numAltModelSites, AltModelSites, allSites,   rateVarAmongSites, altModelMutationRate,  numMU,   Root,  Cijk);
            
        }
        else
        {
            fprintf (stderr, "\n\nERROR: Sorry, the specified model is unknown for me");
            exit(-1);
        }
        
    }
    else {/* both ISM and non-ISM sites */
        numDefaultModelSites = 0;
        numAltModelSites = 0;
        
        for (i=0; i<numSites; i++)
        {
            if (RandomUniform (seed) < propAltModelSites)
                AltModelSites[numAltModelSites++] = i;
            else
                DefaultModelSites[numDefaultModelSites++] = i;
        }
//            /* Evolve ISM sites */
           if (numDefaultModelSites > 0)
               {
                   SimulateISM (treeRoot, genome, NO, seed,
                                DefaultModelSites,numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet,   data, allSites, numMU, cumMij,altModelMutationRate
                                );
                   
               }
//            /* Evolve non-SIM sites */
           if (numAltModelSites > 0)
            {
                if (altModel == ISMhap)
                { fprintf(stderr, "pip");
                    SimulateISM (treeRoot, genome, YES, seed,DefaultModelSites,numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet,   data, allSites, numMU, cumMij,    altModelMutationRate);
                }
                else if (altModel == Mk2)
                { SimulateMk2 (treeRoot, genome, seed, AltModelSites,   numAltModelSites, doUserTree, rateVarAmongSites,  altModelMutationRate, allSites, data, numMU);
                }
                else if (altModel == finiteDNA)
                {
                    SimulateFiniteDNA (treeRoot, genome, seed,doJC,  doHKY, doGTR,  doGTnR,  freqR, freqY,  freqAG,  freqCT, titv,  freq, Mij, numAltModelSites,AltModelSites,  allSites,    rateVarAmongSites,  altModelMutationRate, numMU,   Root,  Cijk);
                  }
                else
                {
                    fprintf (stderr, "\n\nERROR: Sorry, the specified model is unknown for me");
                    exit(-1);
                }
           }
        }
//    }
}
/**************************** RandomGamma *************************/
/*    Generates a gamma number using routines in Ziheng's
 Yang tools.h in PAML
 
 Random standard gamma (Mean=Var=s,  with shape par=s, scale par=1)
 r^(s-1)*exp(-r)
 
 J. Dagpunar (1988) Principles of random variate generation,
 Clarendon Press, Oxford
 
 Calling rndgamma1() if s<1 or rndgamma2() if s>1 or exponential if s=1
 */

double    RandomGamma (double shape, long int *seed)
{
    double gammaNumber = 0;
    
    if (shape <= 0)
        fprintf (stderr, "ERROR: problems with gamma variable generation, shape < 0");
    else if (shape < 1)
        gammaNumber = RandomGamma1 (shape, seed);
    else if (shape > 1)
        gammaNumber = RandomGamma2 (shape, seed);
    else
        gammaNumber = -log (RandomUniform(seed));
    return (gammaNumber);
}
/*************** RandomDirichlet ***************/
//first generates random samples from a gamma and then divide each value by the sum
void  RandomDirichlet (double s, int vectorSize, double **outputVector, long int *seed)
{   int i;
    double sum=0.0;
    double current;
   // *outputVector = malloc(vectorSize * sizeof(double));
   // if (*outputVector == NULL)
   //     return;
    for (i=0; i < vectorSize; i++){
        current = RandomGamma(s, seed);
        (*outputVector)[i] = current;
        //*(*outputVector + i)=current;
        sum=sum+current;
    }
    for (i=0; i < vectorSize; i++){
        (*outputVector)[i] =  (*outputVector)[i] / sum;
       // *(*outputVector + i)= *(*outputVector + i)/sum;
    }
}
/*************** RandomDirichlet2 ***************/
//first generates random samples from a gamma with different parameters and then divide each value by the sum
void  RandomDirichlet2 (double *vectorGamma, int vectorSize, double *outputVector, long int *seed)
{   int i;
    double sum=0.0;
    double current;
    for (i=0; i < vectorSize; i++){
        current = RandomGamma(vectorGamma[i], seed);
        *(outputVector + i)=current;
        sum=sum+current;
    }
    for (i=0; i < vectorSize; i++){
        *(outputVector + i)= *(outputVector + i)/sum;
    }
}
/*************** RandomGamma1 ***************/
double RandomGamma1 (double s, long int *seed)
{
    /* Random standard gamma for s<1
     switching method
     */
    double            r, x=0.0, small=1e-37, w;
    static double   a, p, uf, ss=10.0, d;
    
    if (s!=ss)
    {
        a  = 1.0-s;
        p  = a/(a+s*exp(-a));
        uf = p*pow(small/a,s);
        d  = a*log(a);
        ss = s;
    }
    for (;;)
    {
        r = RandomUniform(seed);
        if (r > p)
        {
            x = a-log((1.0-r)/(1.0-p));
            w=a*log(x)-d;  /* this was with comma in line above before 270917*/
        }
        else if (r>uf)
        {
            x = a*pow(r/p,1/s);
            w=x; /* this was with comma in line above before 270917*/
        }
        else
            return (0.0);
        r = RandomUniform(seed);
        if (1.0-r <= w && r > 0.0)
            if (r*(w+1.0) >= 1.0 || -log(r) <= w)
                continue;
        break;
    }
    return (x);
}
/*************** RandomGamma2 ***************/
double RandomGamma2 (double s, long int *seed)
{
    /* Random standard gamma for s>1
     Best's (1978) t distribution method
     */
    double            r ,d, f, g, x;
    static double    b, h, ss=0;
    
    if (s!=ss)
    {
        b  = s-1.0;
        h  = sqrt(3.0*s-0.75);
        ss = s;
    }
    for (;;)
    {
        r = RandomUniform(seed);
        g = r-r*r;
        f = (r-0.5)*h/sqrt(g);
        x = b+f;
        if (x <= 0.0)
            continue;
        r = RandomUniform(seed);
        d = 64*r*r*g*g*g;
        if (d*x < x-2.0*f*f || log(d) < 2*(b*log(x/b)-f))
            break;
    }
    return (x);
}

/************************************* SimulateMk2 **********************************************/
/* Simulates the nucleotide substitution process under the Mk2 model (see Lewis 2001),
 also called Cavender-Farris-Neyman CFN  model or Jukes-Cantor (1969) model for two alleles */

void SimulateMk2 (TreeNode *p, int genome, long int *seed, int* AltModelSites, int  numAltModelSites,int doUserTree, int rateVarAmongSites, double altModelMutationRate, SiteStr* allSites, int ***data, int* numMU)
{
    int     i;
    
    for (i=0; i<numAltModelSites; i++)
        SimulateMk2forSite (p, genome, AltModelSites[i], seed,  doUserTree,  rateVarAmongSites,  altModelMutationRate, allSites, data, numMU);
    }
/************************************* SimulateMk2ForSite ***************************************/
/* Simulates the nucleotide substitution process for a given site under Mk2 model (see Lewis 2001)
 with equal rates. 0 is the reference (healthy) allele */

void SimulateMk2forSite (TreeNode *p, int genome, int site, long int *seed, int doUserTree, int rateVarAmongSites, double altModelMutationRate, SiteStr* allSites, int ***data, int* numMU )
{
    double    probOfChange, uniform, branchLength;
    int     cell, anccell;
    
    if (p != NULL)
    {
        if (p->isOutgroup == NO)
        {
            cell = p->label;
            anccell = p->anc1->label;
            
            if (doUserTree == YES){
                branchLength = p->lengthModelUnits;//>branchLength;
//                 branchLength = p->length;//>branchLength;
                
            }
            else
            {
                if (rateVarAmongSites == YES)
                    branchLength = altModelMutationRate * p->length * allSites[site].rateMultiplier;
                else
                    branchLength = altModelMutationRate * p->length;
            }
            
            probOfChange = 0.5 - 0.5 * exp (-2.0 * branchLength);
            
            uniform = RandomUniform(seed);
            if (uniform >= probOfChange) /* => no change */
                data[genome][cell][site] = data[genome][anccell][site];
            else /* => there will be change */
            {
                if (data[genome][anccell][site] == 0)
                    data[genome][cell][site] = 1;
                else
                    data[genome][cell][site] = 0;
                
                if (genome == MATERNAL)
                    allSites[site].numMutationsMaternal++;
                else if (genome == PATERNAL)
                    allSites[site].numMutationsPaternal++;
                allSites[site].numMutations++;
                numMU=numMU+1;
            }
        }
        SimulateMk2forSite (p->left,  genome, site, seed,  doUserTree,  rateVarAmongSites,  altModelMutationRate, allSites, data, numMU);
        SimulateMk2forSite (p->right, genome, site, seed,  doUserTree,  rateVarAmongSites,  altModelMutationRate, allSites, data, numMU);
    }
}

/************************************* SimulateFiniteDNA **********************************************/
/* Simulates the nucleotide substitution process under a 4-state Markov model including JC, HKY, GTR and GTRnr */
/* Note that beta is set such that mean substitution rate will be 1.
 E.g., for J-C model, beta=4/3, where 12(1/4)(1/4)(4/3) = 1.      */

void SimulateFiniteDNA (TreeNode *p, int genome, long int *seed, int doJC, int doHKY, int doGTR, int doGTnR, double freqR, double freqY, double freqAG, double freqCT, double titv, double freq[4], double Mij[4][4], int numAltModelSites, int *AltModelSites,SiteStr* allSites,  int rateVarAmongSites, double altModelMutationRate, int *numMU, double Root[], double Cijk[])
{
    int     i, j;
    double beta, kappa;
    //double Qij[16];
    //double mr;
    
    if (doJC == YES)
    {
        beta = 4./3;
    }
    else if (doHKY == YES)
    {
        freqR = freq[A] + freq[G];
        freqY = freq[C] + freq[T];
        freqAG = freq[A] * freq[G];
        freqCT = freq[C] * freq[T];
        kappa = (titv*freqR*freqY)/(freqAG+freqCT);
        beta = 0.5 / (freqR*freqY + kappa*(freqAG+freqCT));
    }
    else if (doGTR == YES || doGTnR == YES)
    {
        for (i=0; i<4; i++)
            for (j=0; j<4; j++)
                Qij[i*4+j] = Mij[i][j] / Mij[2][3] * freq[j];
        mr=0;
        for (i=0; i<4; i++)
        {
            Qij[i*4+i]=0;
            Qij[i*4+i]=-(Qij[i*4]+Qij[i*4+1]+Qij[i*4+2]+Qij[i*4+3]);
            mr-=freq[i]*Qij[i*4+i];
        }
        EigenREV(Root, Cijk);
    }
    
    for (i=0; i<numAltModelSites; i++)
        SimulateFiniteDNAforSite (p,  genome, AltModelSites[i], allSites,  seed,  rateVarAmongSites,  altModelMutationRate, numMU,  doJC,  doHKY,  doGTR,  doGTnR,    beta,    kappa,   freqR,   freqY,   freq,  Root,  Cijk);
        
}
/************************************* SimulateFiniteDNAforSite **********************************************/
/* Simulates JC, HKY, GTR or GTRnr for a given site */
void SimulateFiniteDNAforSite (TreeNode *p, int genome, int site,SiteStr* allSites,  long int *seed, int rateVarAmongSites, double altModelMutationRate, int *numMU, int doJC, int doHKY, int doGTR, int doGTnR, double beta,  double kappa, double freqR, double freqY, double freq[4], double Root[], double Cijk[])
{
    double    branchLength, Pij[4][4];
    int     cell, anccell, ancstate, newstate;
    
    if (p != NULL)
    {
       // if (p->isHealthyRoot == NO)
        if (p->isOutgroup== NO)
        {
            cell = p->label;
            if (p->anc1 !=NULL)
            {
                
                anccell = p->anc1->label;
                 if (genome == MATERNAL)
                     ancstate = p->anc1->maternalSequence[site];
                else
                    ancstate = p->anc1->paternalSequence[site];
           // }
           // if (doUserTree == YES)
                //branchLength = p->branchLength;
          
                
                if (rateVarAmongSites == YES){
                   // branchLength = altModelMutationRate * p->length * allSites[site].rateMultiplier;
                     branchLength =  p->length * allSites[site].rateMultiplier;
                }
                else{
                    branchLength =  p->length;
                    // branchLength = altModelMutationRate * p->length;
                    
                }
                
            
            FillSubstitutionMatrix (Pij, branchLength,
                                    doJC, doHKY,  doGTR,  doGTnR, beta,  kappa,  freqR,  freqY,  freq,  Root,  Cijk);
            if (genome ==MATERNAL)
                newstate =p->maternalSequence[site]=ChooseUniformState (Pij[ancstate], seed);
            else// paternal;
                 newstate =p->paternalSequence[site]=ChooseUniformState (Pij[ancstate], seed);
            //newstate = data[genome][cell][site] = ChooseUniformState (Pij[ancstate], seed);
            
            if (newstate != ancstate)
            {
                if (genome == MATERNAL)
                {
                 allSites[site].numMutationsMaternal++;
                    
                }
                else if (genome == PATERNAL){
                    
                    allSites[site].numMutationsPaternal++;
                }
                
                allSites[site].numMutations++;
                (*numMU)++;
            }
          }
        }
        SimulateFiniteDNAforSite (p->left,  genome, site,allSites, seed,  rateVarAmongSites, altModelMutationRate, numMU,  doJC,  doHKY, doGTR,  doGTnR,  beta,    kappa,   freqR,   freqY,   freq,  Root,  Cijk);
        SimulateFiniteDNAforSite (p->right, genome, site,allSites, seed, rateVarAmongSites, altModelMutationRate, numMU,  doJC,  doHKY,  doGTR,  doGTnR,  beta,    kappa,   freqR,   freqY,   freq ,Root,  Cijk);
    }
}


/********************* FillSubstitutionMatrix **********************/
/* Sets the apropriate model of nucleotide substitution   */
void FillSubstitutionMatrix (double ch_prob[4][4], double branchLength, int doJC, int doHKY, int doGTR, int doGTnR, double beta, double kappa, double freqR, double freqY, double freq[4], double Root[], double Cijk[])
{
    int i, j;
    
    if (branchLength<1e-6)
    {
        for (i=0; i<4; i++)
        {
            for (j=0; j<4; j++)
            {
                if (i == j)
                    ch_prob[i][j] = 1.0;
                else
                    ch_prob[i][j] = 0.0;
            }
        }
    }
    else if (doJC == YES)
        JCmodel (ch_prob, branchLength,  beta);
    else if (doHKY == YES)
        HKYmodel (ch_prob, branchLength, kappa,  freqR,  freqY,  beta, freq);
    else if (doGTR == YES)
        GTRmodel (ch_prob, branchLength,  Root,  Cijk);
    else if (doGTnR == YES)
        GTRmodel (ch_prob, branchLength, Root,  Cijk);
}



/*********************************** JC **************************************/
/*    JC performs Jukes-Cantor 69 correction */
/* Note that beta was set such that mean substitution rate will be 1.
 for the JC model, beta=4/3, where 12(1/4)(1/4)(4/3) = 1.      */
void JCmodel (double Pij[4][4], double branchLength, double beta )
{
    int i, j;
    
    for (i=0; i<4; i++)
    {
        for (j=0; j<4; j++)
        {
            if (i == j)
                Pij[i][j] = 0.25 + 0.75*exp(beta*-branchLength);
            else
                Pij[i][j] = 0.25 - 0.25*exp(beta*-branchLength);
        }
    }
}

/*********************************** HKY **************************************/
/*    HKY performs Hasegawa-Kishino-Yano 85 correction */

void HKYmodel (double Pij[4][4], double branchLength, double kappa, double freqR, double freqY, double beta, double freq[4])
{
    int            i, j;
    double        AA1, t, PIj;
    
    t = branchLength;
    
    for (i=0; i<4; i++)
    {
        for (j=0; j<4; j++)
        {
            if (j == A || j == G)    /* purine */
                PIj = freqR;
            else
                PIj = freqY; /* pyrimidine */
            
            AA1 = 1 + PIj*(kappa-1.0);
            
            if (i==j)
                Pij[i][j] = freq[j] + freq[j]*(1/PIj - 1)*exp(-beta*t) + ((PIj-freq[j])/PIj)*exp(-beta*t*AA);
            else if ((i==A && j==G) || (i==C && j==T) || (i==G && j==A) || (i==T && j==C)) /* transition */
                Pij[i][j] = freq[j] + freq[j]*(1/PIj - 1)*exp(-beta*t) - (freq[j]/PIj)*exp(-beta*t*AA);
            else /* transversion */
                Pij[i][j] = freq[j]*(1-exp(-beta*t));
        }
    }
}


/*************** GTR **********************/
void GTRmodel (double Pij[4][4], double branchLength, double Root[], double Cijk[])
{
    int     i, j, k;
    double    t, expt[4];
    
    t = branchLength;
    
    /* P(t)ij = SUM Cijk * exp{Root*t} */
    for (k=1; k<4; k++)
        expt[k]=exp(t*Root[k]);
    for (i=0; i<4; i++)
        for (j=0; j<4; j++)
        {
            Pij[i][j]=Cijk[i*4*4+j*4+0];
            for (k=1; k<4; k++)
                Pij[i][j]+=Cijk[i*4*4+j*4+k]*expt[k];
        }
}



/********************* AllocateCellStructure  ************************/
/*
 This cell structure keeps information basically about read count functions. It was
 build initially for making doublets.
 */
void AllocateCellStructure(CellStr* cell, int numCells, int numSites)
{
    int i, j, k;
    
    cell = (CellStr*) malloc((numCells+1) * sizeof(CellStr));
    if (!cell)
    {
        fprintf (stderr, "Could not allocate cell (%ld)\n", (numCells+1) * sizeof(CellStr));
        exit (-1);
    }
    
    for (i=0; i<numCells+1; i++)
    {
        cell[i].site  = (CellSiteStr*) malloc (numSites* sizeof(CellSiteStr));
        if (!cell[i].site)
        {
            fprintf (stderr, "Could not allocate the cell[%d].site structure\n", i);
            exit (-1);
        }
        for (j=0; j<numSites; j++)
        {
            cell[i].site[j].readCount = (int*) calloc (4, sizeof(int));
            if (!cell[i].site[j].readCount)
            {
                fprintf (stderr, "Could not allocate the cell[%d].site[%d].readCount structure\n", i, j);
                exit (-1);
            }
            
            cell[i].site[j].readCountDoublet = (int*) calloc (4, sizeof(int));
            if (!cell[i].site[j].readCountDoublet)
            {
                fprintf (stderr, "Could not allocate the cell[%d].site[%d].readCountDoublet structure\n", i, j);
                exit (-1);
            }
            
            cell[i].site[j].genLike = (double**) calloc (4, sizeof(double*));
            if (!cell[i].site[j].genLike)
            {
                fprintf (stderr, "Could not allocate the cell[%d].site[%d].genLike structure\n", i, j);
                exit (-1);
            }
            for (k=0; k<4; k++)
            {
                 cell[i].site[j].genLike[k] = (double*) calloc (4, sizeof(double));
                if (! cell[i].site[j].genLike[k] )
                {
                    fprintf (stderr, "Could not allocate the cell[%d].site[%d].genLike[%d] structure\n", i,j,k);
                    exit (-1);
                }
            }
            
            cell[i].site[j].scaledGenLike = (double**) calloc (4, sizeof(double*));
            if (!cell[i].site[j].scaledGenLike)
            {
                fprintf (stderr, "Could not allocate the cell[%d].site[%d].scaledGenLike structure\n", i, j);
                exit (-1);
            }
            for (k=0; k<4; k++)
            {
                cell[i].site[j].scaledGenLike[k] = (double*) calloc (4, sizeof(double));
                if (!cell[i].site[j].scaledGenLike[k] )
                {
                    fprintf (stderr, "Could not allocate the cell[%d].site[%d].scaledGenLike[%d] structure\n", i,j,k);
                    exit (-1);
                }
            }
            
            cell[i].site[j].genLikeDoublet = (double**) calloc (4, sizeof(double*));
            if (!cell[i].site[j].genLikeDoublet)
            {
                fprintf (stderr, "Could not allocate the cell[%d].site[%d].genLikeDoublet structure\n", i, j);
                exit (-1);
            }
            for (k=0; k<4; k++)
            {
                cell[i].site[j].genLikeDoublet[k] = (double*) calloc (4, sizeof(double));
                if (!cell[i].site[j].genLikeDoublet[k] )
                {
                    fprintf (stderr, "Could not allocate the cell[%d].site[%d].genLikeDoublet[%d] structure\n", i,j,k);
                    exit (-1);
                }
            }
            
            cell[i].site[j].scaledGenLikeDoublet = (double**) calloc (4, sizeof(double*));
            if (!cell[i].site[j].scaledGenLikeDoublet)
            {
                fprintf (stderr, "Could not allocate the cell[%d].site[%d].scaledGenLikeDoublet structure\n", i, j);
                exit (-1);
            }
            for (k=0; k<4; k++)
            {
                cell[i].site[j].scaledGenLikeDoublet[k] = (double*) calloc (4, sizeof(double));
                if (!cell[i].site[j].scaledGenLikeDoublet[k] )
                {
                    fprintf (stderr, "Could not allocate the cell[%d].site[%d].scaledGenLikeDoublet[%d] structure\n", i,j,k);
                    exit (-1);
                }
            }
        }
    }
}

/************************ CheckMatrixSymmetry **************************/
/* Checks whether a given matrix is symmetric */

int CheckMatrixSymmetry(double matrix[4][4])
{
    int i,j;
    
    for(i=0; i<4; i++)
        for(j=0; j<4; j++)
            if(matrix[i][j] != matrix[j][i])
                return NO;
    return YES;
}
/********************************** InitializeGenomes ***********************************/
/* Initialize all genomes with the reference states  */
void InitializeGenomes (TreeNode *p, long int *seed,  int alphabet, int doUserGenome, int numSites, SiteStr* allSites, int doGeneticSignatures, double cumfreq[4],double *triNucFreq, char **cellNames)
{
    int     i, cell, anccell, site;
    double  ran;
    
    if (p != NULL)
    {
        cell = p->label;
       // if (p->cellName !=NULL)
       //    strcpy(cellNames[cell],p->cellName);
//        p->maternalSequence= (int*) calloc (numSites, sizeof(int));
//        p->paternalSequence= (int*) calloc (numSites, sizeof(int));
//        p->NumbersMutationsUnderSubtreePerSite=(int*) calloc (numSites, sizeof(int));
//        p->NumbersMaternalMutationsPerSite=(int*) calloc (numSites, sizeof(int));
//        p->NumbersPaternalMutationsPerSite=(int*) calloc (numSites, sizeof(int));
       
        if (alphabet == DNA)
        {
//            if (p->isOutgroup == YES)
             if (p->anc1 == NULL)
            {
                 if (doGeneticSignatures == NO) /* initialize genome with mononucleotide frequencies */
                {
                    for (site=0; site<numSites; site++)
                    {
                        ran = RandomUniform(seed);
                        for (i=0; i<4; i++)
                        {
                            if (ran <= cumfreq[i])
                            {
                                p->maternalSequence[site]=p->paternalSequence[site]=i;
                        
                            allSites[site].referenceAllele =p->maternalSequence[site];//   // then allSites[site].referenceAllele hosts the reference genome
                                break;
                            }
                        }
                    }
                }
                else /* initialize genome with trinucleotide frequencies */
                {
                    SimulateTriNucFreqGenome (cell, seed, p, alphabet,  doUserGenome,  numSites, allSites,  doGeneticSignatures,  cumfreq,triNucFreq);
                }
            }
            else
            {
                anccell = p->anc1->label;
                for (site=0; site<numSites; site++)
                {
                    p->maternalSequence[site]=p->anc1->maternalSequence[site];
                    p->paternalSequence[site]=p->anc1->paternalSequence[site];
                }
            }
        }
        else{
            for (site=0; site<numSites; site++){
                p->maternalSequence[site]=0;
                p->paternalSequence[site]=0;
                p->numbersMutationsUnderSubtreePerSite[site]=0;
                p->numbersMaternalMutationsPerSite[site]=0;
                p->numbersPaternalMutationsPerSite[site]=0;
               
            }
        }
        InitializeGenomes (p->left, seed,   alphabet,  doUserGenome,  numSites, allSites,  doGeneticSignatures,  cumfreq,triNucFreq, cellNames);
        InitializeGenomes (p->right, seed,  alphabet,  doUserGenome,  numSites, allSites,  doGeneticSignatures,  cumfreq,triNucFreq, cellNames);
    }
}

/********************* WhichNucChar ************************/
/* Returns integer representation for character nucleotudes */

int WhichNucChar (char nucleotide)
{
    if (nucleotide == 'A')
        return (A);
    else if (nucleotide == 'C')
        return (C);
    else if (nucleotide == 'G')
        return (G);
    else if (nucleotide == 'T')
        return (T);
    else if (nucleotide == '?')
        return (ADO);
    else if (nucleotide == '-')
        return (DELETION);
    else if (nucleotide == 'N')
        return (N);
    else if (nucleotide == 'R')
        return (R);
    else
    {
        fprintf (stderr, "\nERROR in WhichNucChar: nucleotide = %c\n",  nucleotide);
        exit(-1);
    }
}
/********************* WhichGenotypeChar ************************/
/* Returns integer representation for character nucleotudes */

int WhichGenotypeChar (char nucleotide)
{
    if (nucleotide == 'A')
        return (AA);
    else if (nucleotide == 'C')
        return (CC);
    else if (nucleotide == 'G')
        return (GG);
    else if (nucleotide == 'T')
        return (TT);
    else if (nucleotide == '?')
        return (__);
    else if (nucleotide == '-')
        return (__);
    else if (nucleotide == 'N')
        return (N);
    else if (nucleotide == 'R')
        return (AG);
    else if (nucleotide == 'M')
        return (AC);
    else if (nucleotide == 'W')
        return (AT);
    else if (nucleotide == 'S')
        return (CG);
    else if (nucleotide == 'Y')
        return (CT);
    else if (nucleotide == 'K')
        return (GT);
    else if (nucleotide == 'a')
        return (A_);
    else if (nucleotide == 'c')
        return (C_);
    else if (nucleotide == 'g')
        return (G_);
    else if (nucleotide == 't')
        return (T_);
    else
    {
        fprintf (stderr, "\nERROR in WhichGenotypeChar: nucleotide = %c\n",  nucleotide);
        exit(-1);
    }
}
/********************************** SimulateTriNucFreqGenome ***********************************/
/*
 Simulate a homozygous diploid genome with the trinucleotide frequencies of the human genome
 
 We simulate first a random trinucleotide.
 Then we go site by site taking into account the last two letters  of the
 previous trinucleotide.  if the first trinucleotide is, for example, CAT
 then we simulate trinucleotide starting in AT_  (i.e. ATA, ATC, ATG or ATT),
 and keep going  taking into account the last two letters of the new trinucleotide,
 and so on.
 */

void SimulateTriNucFreqGenome (int cell, long int *seed, TreeNode *p, int alphabet, int doUserGenome, int numSites, SiteStr* allSites, int doGeneticSignatures, double cumfreq[4], double *triNucFreq )
{
    int         chosenTriNucleotide, nextNucleotide;
    int            k, n1, n2, n3, rest, site;
    double         *prob4, sum;
    
    /* memory allocations */
    prob4 = (double *) calloc (4, sizeof(double));
    if (!prob4)
    {
        fprintf (stderr, "Could not allocate the prob4 vector\n");
        exit (-1);
    }
    
    /* choose first trinucleotide */
    chosenTriNucleotide = ChooseUniformState(triNucFreq, seed);
    
    /* find bases of the selected trinucleotide */
    n1 = chosenTriNucleotide/16;
    rest = chosenTriNucleotide%16;
    n2 = rest/4;
    n3 = rest%4;
    //fprintf (stderr, "\n%2d %c%c%c ", chosenTriNucleotide, WhichNuc(n1), WhichNuc(n2), WhichNuc(n3));
    
    site = 0;
    p->maternalSequence[site]=p->paternalSequence[site]=n1;
//    data[MATERNAL][cell][site] = data[PATERNAL][cell][site] = allSites[site].referenceAllele = n1;
    site++;
     p->maternalSequence[site]=p->paternalSequence[site] = allSites[site].referenceAllele = n2;
    site++;
     p->maternalSequence[site]=p->paternalSequence[site]= allSites[site].referenceAllele = n3;
    
    /* fill the rest of the genome */
    /* choose next nucleotide given the last two bases of the previous trinucleotide  */
    for (site=3; site<numSites; site++)
    {
        /* normalize frequencies given the last two bases */
        sum = 0;
        for (k=0; k<4; k++)
            sum += triNucFreq[trinuc(n2,n3,k)];
        for (k=0; k<4; k++)
            prob4[k] = triNucFreq[trinuc(n2,n3,k)] / sum;
        
        nextNucleotide = ChooseUniformState(prob4, seed);
         p->maternalSequence[site]=p->paternalSequence[site]= allSites[site].referenceAllele = nextNucleotide;
        
        /* move downstream one position */
        n1 = n2;
        n2 = n3;
        n3 = nextNucleotide;
    }
    
    free(prob4);
    prob4=NULL;
}
/************************************* ChooseUniformState **********************************************/
/* Chooses uniformy a random state according to a vector of state probabilities */

int ChooseUniformState (double *prob, long int *seed)
{
    int            chosenState;
    double        ran, cumProb;
    
    chosenState = 0;
    cumProb = prob[chosenState];
    ran = RandomUniform(seed);
    
    while (ran > cumProb)
        cumProb += prob[++chosenState];
    
    return chosenState;
}
/********************* PrepareGlobalFiles **********************/
/* Open global files to output results */
void PrepareGlobalFiles(int argc, char **argv, int doPrintTree, FILE *fpTrees, char   resultsDir[MAX_NAME] ,  char        treeFile[MAX_NAME],  char        timesFile[MAX_NAME], char        SNVgenotypesFile[MAX_NAME],char SNVhaplotypesFile[MAX_NAME],char trueHaplotypesFile[MAX_NAME], char fullGenotypesFile[MAX_NAME], char fullHaplotypesFile[MAX_NAME],char VCFfile[MAX_NAME],
                        char CATGfile[MAX_NAME] ,FILE            *fpSNVgenotypes, int doPrintTimes, int doSimulateData, FILE  *fpTimes, FILE *fpSNVhaplotypes, FILE *fpTrueHaplotypes,
    FILE *fpFullGenotypes,FILE *fpFullHaplotypes, FILE *fpVCF, FILE *fpCATG, int doPrintSNVgenotypes, int  doPrintSNVhaplotypes, int doPrintTrueHaplotypes, int doPrintFullGenotypes, int doPrintFullHaplotypes, int doNGS,  int doPrintCATG, int numDataSets)
{
    char File[MAX_NAME];
    
    /* contains the simulated tree in Newick format */
//    if (doPrintTree == YES)
//    {
//        sprintf(File,"%s/%s", resultsDir, treeFile);
//        if ((fpTrees = fopen(File, "w")) == NULL)
//        {
//            fprintf (stderr, "Can't open \"%s\"\n", File);
//            exit(-1);
//        }
//    }
    
    /* contains a list of times and branch lenghts for all nodes in the simulated trees */
//    if (doPrintTimes == YES)
//    {
//        sprintf(File,"%s/%s", resultsDir, timesFile);
//        if ((fpTimes = fopen(File, "w")) == NULL)
//        {
//            fprintf (stderr, "Can't open \"%s\"\n", File);
//            exit(-1);
//        }
//    }
    
    if (doSimulateData == YES)
    {
        /* contains SNV genotypes for every cell */
        if (doPrintSNVgenotypes == YES)
        {
            sprintf(File,"%s/%s", resultsDir, SNVgenotypesFile);
            if ((fpSNVgenotypes = fopen(File, "w")) == NULL)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        
        /* contains haplotypes for variable sites for every cell */
        if (doPrintSNVhaplotypes == YES)
        {
            sprintf(File,"%s/%s", resultsDir, SNVhaplotypesFile);
            if ((fpSNVhaplotypes = fopen(File, "w")) == NULL)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        
        /* contains error-free haplotypes for variable sites for every cell */
        if (doPrintTrueHaplotypes == YES)
        {
            sprintf(File,"%s/%s", resultsDir, trueHaplotypesFile);
            if ((fpTrueHaplotypes = fopen(File, "w")) == NULL)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        
        /* contains all genotypes (variable or invariable) for every cell */
        if (doPrintFullGenotypes == YES)
        {
            sprintf(File,"%s/%s", resultsDir, fullGenotypesFile);
            if ((fpFullGenotypes = fopen(File, "w")) == NULL)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        
        /* contains haplotypes for all sites for every cell */
        if (doPrintFullHaplotypes == YES)
        {
            sprintf(File,"%s/%s", resultsDir, fullHaplotypesFile);
            if ((fpFullHaplotypes = fopen(File, "w")) == NULL)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        
        /* contains reads counts and genotype likelihoods for every SNV and cell */
        if (doNGS == YES)
        {
            sprintf(File,"%s/%s", resultsDir, VCFfile);
            if ((fpVCF = fopen(File, "w")) == NULL)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        
        /* contains reads counts for every SNV and cell */
        if (doPrintCATG == YES)
        {
            sprintf(File,"%s/%s", resultsDir, CATGfile);
            if ((fpCATG = fopen(File, "w")) == NULL)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        
    }
    
    if (doPrintSNVgenotypes == YES)
    {
        //fprintf (fpSNVgenotypes, "%s - ",PROGRAM_NAME);
        //PrintDate (fpSNVgenotypes);
       // fprintf (fpSNVgenotypes, "SNV genotypes\n");
        //PrintCommandLine (fpSNVgenotypes, argc, argv);
        fprintf (fpSNVgenotypes,"\n%d\n", numDataSets);
    }
    
    if (doPrintSNVhaplotypes == YES)
    {
        fprintf (fpSNVhaplotypes, "%s - ",PROGRAM_NAME);
        //PrintDate (fpSNVhaplotypes);
        fprintf (fpSNVhaplotypes, "SNV haplotypes\n");
        //PrintCommandLine (fpSNVhaplotypes, argc, argv);
        fprintf (fpSNVhaplotypes,"\n%d\n", numDataSets);
    }
    
    if (doPrintTrueHaplotypes == YES)
    {
//        fprintf (fpTrueHaplotypes, "%s - ",PROGRAM_NAME);
//        //PrintDate (fpTrueHaplotypes);
//        fprintf (fpTrueHaplotypes, "True haplotypes\n");
//        //PrintCommandLine (fpTrueHaplotypes, argc, argv);
//        fprintf (fpTrueHaplotypes,"\n%d\n", numDataSets);
    }
    
    if (doPrintFullGenotypes == YES)
    {
        fprintf (fpFullGenotypes, "%s - ",PROGRAM_NAME);
        //PrintDate (fpFullGenotypes);
        fprintf (fpFullGenotypes, "Full genotypes\n");
        //PrintCommandLine (fpFullGenotypes, argc, argv);
        fprintf (fpFullGenotypes,"\n%d\n", numDataSets);
    }
    
    if (doPrintFullHaplotypes == YES)
    {
       // fprintf (fpFullHaplotypes, "%s - ",PROGRAM_NAME);
        //PrintDate (fpFullHaplotypes);
        fprintf (fpFullHaplotypes, "Full haplotypes\n");
        //PrintCommandLine (fpFullHaplotypes, argc, argv);
        fprintf (fpFullHaplotypes,"\n%d\n", numDataSets);
    }
    
    if (doNGS == YES)
    {
       // fprintf (fpVCF, "%s - ",PROGRAM_NAME);
        //PrintDate (fpVCF);
        fprintf (fpVCF, "Read counts and genotype likelihoods\n");
        //PrintCommandLine (fpVCF, argc, argv);
        fprintf (fpVCF,"\n%d\n", numDataSets);
    }
    
    if (doPrintCATG == YES)
    {
        //fprintf (fpCATG, "%s - ",PROGRAM_NAME);
        //PrintDate (fpCATG);
        fprintf (fpCATG, "Read counts\n");
        //PrintCommandLine (fpCATG, argc, argv);
        fprintf (fpCATG,"\n%d\n", numDataSets);
    }
}
/***************************** PrintFullHaplotypes *******************************/
/* Prints observed/ML haplotypes for all sites (variable + invariable) to a file */

static void PrintFullHaplotypes (FILE *fp, int ***data, int doPrintIUPAChaplotypes, int doPrintAncestors, int doNGS, int numCells, int numSites, int alphabet, int doUserTree, char **cellNames, CellStr            *cell, int        HEALTHY_ROOT, int TUMOR_ROOT)
{
    int         i, j;
    
//    if (doPrintIUPAChaplotypes == NO)
//    {
//        if (doPrintAncestors == YES && doNGS == NO)
//            fprintf (fp,"%d %d\n",2*(numCells+3), numSites);
//        else
//            fprintf (fp,"%d %d\n",2*(numCells+1), numSites);
//    }
//    else
//    {
        if (doPrintAncestors == YES && doNGS == NO)
            fprintf (fp,"%d %d\n", numCells+3, numSites);
        else
            fprintf (fp,"%d %d\n", numCells+1, numSites);
//            fprintf (fp,"%d %d\n", numCells+1, numSites);
//    }
    
    if (alphabet == DNA)
    {
        if (doPrintIUPAChaplotypes == YES)
        {
            for (i=0; i<numCells+1; i++)
            {
                /* print IUPAC haplotype */
                if (i == numCells)
                    fprintf (fp,"healthycell  ");
                else
                {
                    if (doUserTree == NO)
                        fprintf (fp,"tumcell%04d  ", i+1);
                    else
                        fprintf (fp,"%-12s ", cellNames[i]);
                }
                for (j=0; j<numSites; j++)
                {
                    if (doNGS == NO)
                        fprintf (fp, "%c", WhichIUPAC(data[MATERNAL][i][j],data[PATERNAL][i][j]));
                    else if (cell[i].hasDoublet == NO)
                        fprintf (fp, "%c", WhichIUPAC(cell[i].site[j].MLmatAllele, cell[i].site[j].MLpatAllele));
                    else
                        fprintf (fp, "%c", WhichIUPAC(cell[i].site[j].MLmatAlleleDoublet, cell[i].site[j].MLpatAlleleDoublet));
                }
                fprintf (fp,"\n");
            }
            
            if (doPrintAncestors == YES && doNGS == NO)
            {
                fprintf (fp,"hearoot%04d  ", 0);
                for (j=0; j<numSites; j++)
                    fprintf (fp, "%c", WhichIUPAC(data[MATERNAL][HEALTHY_ROOT][j],data[PATERNAL][HEALTHY_ROOT][j]));
                
                fprintf (fp,"\ntumroot%04d  ", 0);
                for (j=0; j<numSites; j++)
                    fprintf (fp, "%c", WhichIUPAC(data[MATERNAL][TUMOR_ROOT][j],data[PATERNAL][TUMOR_ROOT][j]));
                fprintf (fp,"\n");
            }
        }
        else // print maternal and paternal DNA haplotypes
        {
            for (i=0; i<numCells+1; i++)
            {
                /* print maternal haplotype */
                if (i == numCells)
                    fprintf (fp,"healthycellm ");
                else
                {
                    if (doUserTree == NO)
                        fprintf (fp,"tumcell%04dm ", i+1);
                    else
                        fprintf (fp,"m%-12s", cellNames[i]);
                }
                
                for (j=0; j<numSites; j++)
                {
                    if (doNGS == NO)
                        fprintf (fp, "%c", WhichNuc(data[MATERNAL][i][j]));
                    else if (cell[i].hasDoublet == NO)
                        fprintf (fp, "%c", WhichNuc(cell[i].site[j].MLmatAllele));
                    else
                        fprintf (fp, "%c", WhichNuc(cell[i].site[j].MLmatAlleleDoublet));
                }
                fprintf (fp,"\n");
                
                /* print paternal haplotype */
                if (i == numCells)
                    fprintf (fp,"healthycellp ");
                else
                {
                    if (doUserTree == NO)
                        fprintf (fp,"tumcell%04dp ", i+1);
                    else
                        fprintf (fp,"p%-12s", cellNames[i]);
                }
                
                for (j=0; j<numSites; j++)
                {
                    if (doNGS == NO)
                        fprintf (fp, "%c", WhichNuc(data[PATERNAL][i][j]));
                    else if (cell[i].hasDoublet == NO)
                        fprintf (fp, "%c", WhichNuc(cell[i].site[j].MLpatAllele));
                    else
                        fprintf (fp, "%c", WhichNuc(cell[i].site[j].MLpatAlleleDoublet));
                }
                fprintf (fp,"\n");
            }
            
            if (doPrintAncestors == YES && doNGS == NO)
            {
                fprintf (fp,"hearoot%04dm ", 0);
                for (j=0; j<numSites; j++)
                    fprintf (fp, "%c", WhichNuc(data[MATERNAL][HEALTHY_ROOT][j]));
                fprintf (fp,"\nhearoot%04dp ", 0);
                for (j=0; j<numSites; j++)
                    fprintf (fp, "%c", WhichNuc(data[PATERNAL][HEALTHY_ROOT][j]));
                
                fprintf (fp,"\ntumroot%04dm ", 0);
                for (j=0; j<numSites; j++)
                    fprintf (fp, "%c", WhichNuc(data[MATERNAL][TUMOR_ROOT][j]));
                fprintf (fp,"\ntumroot%04dp ", 0);
                for (j=0; j<numSites; j++)
                    fprintf (fp, "%c", WhichNuc(data[PATERNAL][TUMOR_ROOT][j]));
                fprintf (fp,"\n");
            }
        }
    }
    else  //print binary haplotypes
    {
        if (doPrintIUPAChaplotypes == YES) // print binary consensus haplotypes
        {
            for (i=0; i<numCells+1; i++)
            {
                if (i == numCells)
                    fprintf (fp,"healthycell  ");
                else
                {
                    if (doUserTree == NO)
                        fprintf (fp,"tumcell%04d ", i+1);
                    else
                        fprintf (fp,"%-12s ", cellNames[i]);
                }
                for (j=0; j<numSites; j++)
                    fprintf (fp, "%c", WhichConsensusBinary(data[MATERNAL][i][j],data[PATERNAL][i][j]));
                fprintf (fp,"\n");
            }
            
            if (doPrintAncestors == YES)
            {
                fprintf (fp,"hearoot%04d  ", 0);
                for (j=0; j<numSites; j++)
                    fprintf (fp, "%c", WhichConsensusBinary(data[MATERNAL][HEALTHY_ROOT][j],data[PATERNAL][HEALTHY_ROOT][j]));
                
                fprintf (fp,"\ntumroot%04d  ", 0);
                for (j=0; j<numSites; j++)
                    fprintf (fp, "%c", WhichConsensusBinary(data[MATERNAL][TUMOR_ROOT][j],data[PATERNAL][TUMOR_ROOT][j]));
                fprintf (fp,"\n");
            }
        }
        else // print maternal and paternal binary haplotypes
        {
            for (i=0; i<numCells+1; i++)
            {
                /* print maternal haplotype */
                if (i == numCells)
                    fprintf (fp,"healthycellm ");
                else
                {
                    if (doUserTree == NO)
                        fprintf (fp,"tumcell%04dm ", i+1);
                    else
                        fprintf (fp,"m%-12s", cellNames[i]);
                }
                for (j=0; j<numSites; j++)
                    fprintf (fp, "%c", WhichMut(data[MATERNAL][i][j]));
                fprintf (fp,"\n");
                
                /* print paternal haplotype */
                if (i == numCells)
                    fprintf (fp,"healthycellp ");
                else
                {
                    if (doUserTree == NO)
                        fprintf (fp,"tumcell%04dp ", i+1);
                    else
                        fprintf (fp,"p%-12s", cellNames[i]);
                }
                for (j=0; j<numSites; j++)
                    fprintf (fp, "%c", WhichMut(data[PATERNAL][i][j]));
                fprintf (fp,"\n");
            }
            
            if (doPrintAncestors == YES)
            {
                fprintf (fp,"hearoot%04dm ", 0);
                for (j=0; j<numSites; j++)
                    fprintf (fp, "%c", WhichMut(data[MATERNAL][HEALTHY_ROOT][j]));
                fprintf (fp,"\nhearoot%04dp ", 0);
                for (j=0; j<numSites; j++)
                    fprintf (fp, "%c", WhichMut(data[PATERNAL][HEALTHY_ROOT][j]));
                
                fprintf (fp,"\ntumroot%04dm ", 0);
                for (j=0; j<numSites; j++)
                    fprintf (fp, "%c", WhichMut(data[MATERNAL][TUMOR_ROOT][j]));
                fprintf (fp,"\ntumroot%04dp ", 0);
                for (j=0; j<numSites; j++)
                    fprintf (fp, "%c", WhichMut(data[PATERNAL][TUMOR_ROOT][j]));
                fprintf (fp,"\n");
            }
        }
    }
}

/***************************** PrintTrueFullHaplotypes *******************************/
/* Prints observed/ML haplotypes for all sites (variable + invariable) to a file */

 void PrintTrueFullHaplotypes (FILE *fp, TreeNode* nodes,TreeNode* treeRoot, int numNodes, int doPrintIUPAChaplotypes, int doPrintAncestors, int numSites, int numCells, int alphabet, int doUserTree , int ***data,   int doNGS,   char **cellNames, CellStr            *cell, int        HEALTHY_ROOT, int TUMOR_ROOT , char *cellnames[], int doUseObservedCellName)
{
    int         i, j;
    char *temp;
    TreeNode *p;
//    if (doPrintIUPAChaplotypes == NO)
//    {
//        if (doPrintAncestors == YES)
//            fprintf (fp,"%d %d\n",2*(numCells+3), numSites);
//        else
//           fprintf (fp,"%d %d\n",numCells+1, numSites);
//    }
//    else
//    {
//        if (doPrintAncestors == YES)
//            fprintf (fp,"%d %d\n", numCells+3, numSites);
//        else
//            fprintf (fp,"%d %d\n", numCells+1, numSites);
//    }
    TreeNode * healthyTip= getHealthyTip(treeRoot);
    if (alphabet == DNA)
    {
        
        if (doPrintIUPAChaplotypes == YES)
        {
              fprintf (fp,"%d %d\n",numCells +1, numSites);
            for (i=0; i<numCells; i++){
                p = (nodes + i);
                /* print IUPAC haplotype */
                if (p !=NULL){
                    
                    if (p->left==NULL && p->right ==NULL){
                        if (doUseObservedCellName == YES)
                            temp=p->observedCellName;
                        else
                           temp=p->cellName;
                        fprintf (fp,"%-12s ", temp);
                        for (j=0; j<numSites; j++)
                            fprintf (fp, "%c", WhichIUPAC(p->maternalSequence[j],p->paternalSequence[j]));
                        fprintf (fp,"\n");
                        
                    }
                }
            }
            
            fprintf (fp,"%-12s ", healthyTip->observedCellName);
           for (j=0; j<numSites; j++)
              fprintf (fp, "%c", WhichIUPAC(healthyTip->maternalSequence[j],healthyTip->paternalSequence[j]));
           fprintf (fp,"\n");
        }
        else // print maternal and paternal DNA haplotypes
        {
            fprintf (fp,"%d %d\n",2*(numCells+1), numSites);
            for (i=0; i<numCells; i++){
                p = (nodes + i);
                if (p !=NULL){
                    
                    if (p->left==NULL && p->right ==NULL){
                        if (doUseObservedCellName == YES)
                            temp=p->observedCellName;
                        else
                            temp=p->cellName;
                        fprintf (fp,"m%-12s ", temp);
                        for (j=0; j<numSites; j++)
                            fprintf (fp, "%c", WhichNuc(p->maternalSequence[j]));
                        fprintf (fp,"\n");
                        fprintf (fp,"p%-12s ", temp);
                        for (j=0; j<numSites; j++)
                            fprintf (fp, "%c", WhichNuc(p->paternalSequence[j]));
                        fprintf (fp,"\n");
                    }
                }
            }
            if (doUseObservedCellName == YES)
                fprintf (fp,"m%-12s ", healthyTip->observedCellName);
            else
                fprintf (fp,"m%-12s ", healthyTip->cellName);
           
            for (j=0; j<numSites; j++)
                fprintf (fp, "%c", WhichNuc(healthyTip->maternalSequence[j]));
            fprintf (fp,"\n");
            if (doUseObservedCellName == YES)
                fprintf (fp,"p%-12s ", healthyTip->observedCellName);
            else
                fprintf (fp,"p%-12s ", healthyTip->cellName);
    
            for (j=0; j<numSites; j++)
                fprintf (fp, "%c", WhichNuc(healthyTip->paternalSequence[j]));
            fprintf (fp,"\n");
//            for (i=0; i<numCells+1; i++)
//            {
//                /* print maternal haplotype */
//                if (i == numCells)
//                    fprintf (fp,"healthycellm ");
//                else
//                {
////                    if (doUserTree == NO)
////                        fprintf (fp,"tumcell%04dm ", i+1);
////                    else
//                        fprintf (fp,"m%-12s", cellNames[i]);
//                }
//
//                for (j=0; j<numSites; j++)
//                {
//                    fprintf (fp, "%c", WhichNuc(data[MATERNAL][i][j]));
//                }
//                fprintf (fp,"\n");
//
//                /* print paternal haplotype */
//                if (i == numCells)
//                    fprintf (fp,"healthycellp ");
//                else
//                {
////                    if (doUserTree == NO)
////                        fprintf (fp,"tumcell%04dp ", i+1);
////                    else
//                        fprintf (fp,"p%-12s", cellNames[i]);
//                }
//
//                for (j=0; j<numSites; j++)
//                {
//                    fprintf (fp, "%c", WhichNuc(data[PATERNAL][i][j]));
//                }
//                fprintf (fp,"\n");
//            }
            
//            if (doPrintAncestors == YES)
//            {
//                fprintf (fp,"hearoot%04dm ", 0);
//                for (j=0; j<numSites; j++)
//                    fprintf (fp, "%c", WhichNuc(data[MATERNAL][HEALTHY_ROOT][j]));
//                fprintf (fp,"\nhearoot%04dp ", 0);
//                for (j=0; j<numSites; j++)
//                    fprintf (fp, "%c", WhichNuc(data[PATERNAL][HEALTHY_ROOT][j]));
//
//                fprintf (fp,"\ntumroot%04dm ", 0);
//                for (j=0; j<numSites; j++)
//                    fprintf (fp, "%c", WhichNuc(data[MATERNAL][TUMOR_ROOT][j]));
//                fprintf (fp,"\ntumroot%04dp ", 0);
//                for (j=0; j<numSites; j++)
//                    fprintf (fp, "%c", WhichNuc(data[PATERNAL][TUMOR_ROOT][j]));
//                fprintf (fp,"\n");
//            }
        }
    }
    else  //print binary haplotypes
    {
        if (doPrintIUPAChaplotypes == YES) // print binary consensus haplotypes
        {
            for (i = 0; i < numCells; i++)
            {
                p = (nodes + i);
                
                if (p->left==NULL && p->right ==NULL){
                    if (doUseObservedCellName == YES)
                        temp=p->observedCellName;
                    else
                        temp=p->cellName;
                    
                    fprintf (fp,"%-12s", temp);
                    for (j=0; j<numSites; j++)
                        fprintf (fp, "%c", WhichConsensusBinary(p->maternalSequence[j],p->paternalSequence[j]));
                    fprintf (fp,"\n");
                }
            }
//            {
////                if (i == numCells)
////                    fprintf (fp,"healthycell  ");
////                else
////                {
////                if (doUserTree == NO){
////                    fprintf (fp,"tumcell%04d ", i+1);
////                    fprintf (fp,"tumcell%04d ", i+1);
////                   }
//
////                    else
//                        temp=cellNames[i];
//
//                        fprintf (fp,"%-12s", cellNames[i]);
////                }
//                for (j=0; j<numSites; j++)
//                    fprintf (fp, "%c", WhichConsensusBinary(data[MATERNAL][i][j],data[PATERNAL][i][j]));
//                fprintf (fp,"\n");
//            }
            
//            if (doPrintAncestors == YES)
//            {
//                fprintf (fp,"hearoot%04d  ", 0);
//                for (j=0; j<numSites; j++)
//                    fprintf (fp, "%c", WhichConsensusBinary(data[MATERNAL][HEALTHY_ROOT][j],data[PATERNAL][HEALTHY_ROOT][j]));
//
//                fprintf (fp,"\ntumroot%04d  ", 0);
//                for (j=0; j<numSites; j++)
//                    fprintf (fp, "%c", WhichConsensusBinary(data[MATERNAL][TUMOR_ROOT][j],data[PATERNAL][TUMOR_ROOT][j]));
//                fprintf (fp,"\n");
//            }
        }
        else // print maternal and paternal binary haplotypes
        {
            int i=0;
            int numAddedTips=0;
            fprintf (fp,"%d %d\n",(numCells+1), numSites);
            for (i=0; i<numCells; i++){
                p = (nodes + i);
                if (p !=NULL){
                    
                    if (p->left==NULL && p->right ==NULL){
                        if (doUseObservedCellName == YES)
                            temp=p->observedCellName;
                        else
                            temp=p->cellName;
                        numAddedTips++;
                        fprintf (fp,"%-12s ", temp);
                        for (j=0; j<numSites; j++)
                            fprintf (fp, "%c", WhichMut(p->maternalSequence[j]+p->paternalSequence[j]));
                        fprintf (fp,"\n");
                    }
                }
            }
            
          //this next part is for printing the root/healthy cell
            if (doUseObservedCellName == YES)
                fprintf (fp,"%-12s ", healthyTip->observedCellName);
            else
                 fprintf (fp,"%-12s ", healthyTip->cellName);

            for (j=0; j<numSites; j++)
                 fprintf (fp, "%c", WhichMut(healthyTip->maternalSequence[j]+healthyTip->paternalSequence[j]));
            fprintf (fp,"\n");

            
            
////            for (i=0; i<numCells+1; i++)
//            for (i=0; i<numCells; i++)
//            {
//                /* print maternal haplotype */
////                if (i == numCells)
////                    fprintf (fp,"healthycellm ");
////                else
////                {
////                    if (doUserTree == NO)
////                        fprintf (fp,"tumcell%04dm ", i+1);
////                    else
//                        temp=cellNames[i];
//
//                        fprintf (fp,"%-12s ", cellNames[i]);
////                }
//                for (j=0; j<numSites; j++)
//                    fprintf (fp, "%c", WhichMut(data[MATERNAL][i][j]+data[PATERNAL][i][j]));
//                fprintf (fp,"\n");
//
                /* print paternal haplotype */
//                if (i == numCells)
//                    fprintf (fp,"healthycellp ");
//                else
//                {
//                    if (doUserTree == NO)
//                        fprintf (fp,"tumcell%04dp ", i+1);
//                    else
//                        fprintf (fp,"p%-12s ", cellNames[i]);
//                }
//                for (j=0; j<numSites; j++)
//                    fprintf (fp, "%c", WhichMut(data[PATERNAL][i][j]));
//                fprintf (fp,"\n");
//            }
//            now tumor  root
//            fprintf (fp,"%-12s ", cellNames[HEALTHY_ROOT]);
//            for (j=0; j<numSites; j++)
//               fprintf (fp, "%c", WhichMut(data[MATERNAL][HEALTHY_ROOT][j]+data[PATERNAL][HEALTHY_ROOT][j]));
//            fprintf (fp,"\n");
//            fprintf (fp,"p%-12s ", cellNames[HEALTHY_ROOT]);
//            for (j=0; j<numSites; j++)
//                fprintf (fp, "%c", WhichMut(data[PATERNAL][HEALTHY_ROOT][j]));
//            fprintf (fp,"\n");
            
//            if (doPrintAncestors == YES)
//            {
//                fprintf (fp,"hearoot%04dm ", 0);
//                for (j=0; j<numSites; j++)
//                    fprintf (fp, "%c", WhichMut(data[MATERNAL][HEALTHY_ROOT][j]));
//                fprintf (fp,"\nhearoot%04dp ", 0);
//                for (j=0; j<numSites; j++)
//                    fprintf (fp, "%c", WhichMut(data[PATERNAL][HEALTHY_ROOT][j]));
//
//                fprintf (fp,"\ntumroot%04dm ", 0);
//                for (j=0; j<numSites; j++)
//                    fprintf (fp, "%c", WhichMut(data[MATERNAL][TUMOR_ROOT][j]));
//                fprintf (fp,"\ntumroot%04dp ", 0);
//                for (j=0; j<numSites; j++)
//                    fprintf (fp, "%c", WhichMut(data[PATERNAL][TUMOR_ROOT][j]));
//                fprintf (fp,"\n");
//            }
        }
    }
}
/********************* WhichIUPAC ************************/
/* Returns the IUPAC representation of the genotype */
/*
 UPAC nucleotide code    Base
 A    Adenine
 C    Cytosine
 G    Guanine
 T (or U)    Thymine (or Uracil)
 R    A or G
 Y    C or T
 S    G or C
 W    A or T
 K    G or T
 M    A or C
 B    C or G or T
 D    A or G or T
 H    A or C or T
 V    A or C or G
 N    unknown state
 . or -    gap
 
 This is what we do:
 
 A/A => A
 A/C => M
 A/G => R
 A/T => W
 A/_ => a
 
 C/A => M
 C/C => C
 C/G => S
 C/T => Y
 C/_ => c
 
 G/A => R
 G/C => S
 G/G => G
 G/T => K
 G/_ => g
 
 T/A => W
 T/C => Y
 T/G => K
 T/T => T
 T/_ => t
 
 _/A => a
 _/C => c
 _/G => g
 _/T => t
 _/_ => -
 
 */

char WhichIUPAC (int allele1, int allele2)
{
    if (allele1 == 0)
    {
        if (allele2 == 0)        //AA
            return ('A');
        else if (allele2 == 1)    //AC
            return ('M');
        else if (allele2 == 2)    //AG
            return ('R');
        else if (allele2 == 3)    //AT
            return ('W');
        else if (allele2 == ADO)    //A?
            return ('a');
        else if (allele2 == DELETION)    //A–
            return ('a');
        else
            return ('N');
    }
    else if (allele1 == 1)
    {
        if (allele2 == 0)        //CA
            return ('M');
        else if (allele2 == 1)    //CC
            return ('C');
        else if (allele2 == 2)    //CG
            return ('S');
        else if (allele2 == 3)    //CT
            return ('Y');
        else if (allele2 == ADO)    //C?
            return ('c');
        else if (allele2 == DELETION)    //C–
            return ('c');
        else
            return ('N');
    }
    else if (allele1 == 2)
    {
        if (allele2 == 0)        //GA
            return ('R');
        else if (allele2 == 1)    //GC
            return ('S');
        else if (allele2 == 2)    //GG
            return ('G');
        else if (allele2 == 3)    //GT
            return ('K');
        else if (allele2 == ADO)    //G?
            return ('g');
        else if (allele2 == DELETION)    //G–
            return ('g');
        else
            return ('N');
    }
    else if (allele1 == 3)
    {
        if (allele2 == 0)        //TA
            return ('W');
        else if (allele2 == 1)    //TC
            return ('Y');
        else if (allele2 == 2)    //TG
            return ('K');
        else if (allele2 == 3)    //TT
            return ('T');
        else if (allele2 == ADO)    //T?
            return ('t');
        else if (allele2 == DELETION)    //T–
            return ('t');
        else
            return ('N');
    }
    else if (allele1 == ADO)
    {
        if (allele2 == 0)        //?A
            return ('a');
        else if (allele2 == 1)    //?C
            return ('c');
        else if (allele2 == 2)    //?G
            return ('g');
        else if (allele2 == 3)    //?T
            return ('t');
        else if (allele2 == ADO)    //??
            return ('-');
        else if (allele2 == DELETION)    //?-
            return ('-');
        else
            return ('N');
    }
    else if (allele1 == DELETION)
    {
        if (allele2 == 0)        //-A
            return ('a');
        else if (allele2 == 1)    //-C
            return ('c');
        else if (allele2 == 2)    //-G
            return ('g');
        else if (allele2 == 3)    //-T
            return ('t');
        else if (allele2 == ADO)    //-?
            return ('-');
        else if (allele2 == DELETION)    //--
            return ('-');
        else
            return ('N');
    }
    else
        return ('N');
}
/********************* WhichGenotypeFromIUPAC ************************/
/* Returns the genotypes from the IUPAC code */
//int  WhichGenotypeFromIUPAC (char  iupac)
//{
//    int result;
// 
//        if (iupac == A)        //AA
//            result= AA;
//        else if (iupac == 'M')    //AC
//             result= AC;
//        else if (iupac == 'R')    //AG
//             result='AG';
//        else if (iupac == 'W')    //AT
//             result='AT';
//        else if (iupac == 'a')    //A?
//             result='A?';
////        else if (iupac == DELETION)    //A–
////            return ('A-');
//        else if (iupac == 'N')
//            result='--';
//        else if (iupac == 'M')        //CA
//            result='CA';
//        else if (iupac == 'C')    //CC
//            result='CC';
//        else if (iupac == 'S')    //CG
//            result='CG';
//        else if (iupac == 'Y')    //CT
//           result='CT';
//        else if (iupac == 'c')    //C?
//            result='C?';
////        else if (iupac == DELETION)    //C–
////            return ('C–');
//       else if (iupac == 'R')        //GA
//            result='GA';
//        else if (iupac == 'S')    //GC
//            result='GC';
//        else if (iupac == 'G')    //GG
//            result='GG';
//        else if (iupac == 'K')    //GT
//            result='GT';
//        else if (iupac == 'g')    //G?
//            result='G?';
//    //        else if (iupac == DELETION)    //C–
//    //            return ('G–');
//        if (iupac == 'W')        //TA
//            result='TA';
//        else if (iupac == 'Y')    //TC
//            result='TC';
//        else if (iupac == 'K')    //TG
//            result='TG';
//        else if (iupac == 'T')    //TT
//            result='TT';
//        else if (iupac == 't')    //T?
//            result='T?';
////        else if (iupac == DELETION)    //T–
////            return ('t');
//        if (iupac == 'a')        //?A
//           result='?A';
//        else if (iupac == 'c')    //?C
//            result='?C';
//        else if (iupac == 'g')    //?G
//            result='?G';
//        else if (iupac == 't')    //?T
//            result='?T';
//        else if (iupac == '-')    //??
//            result='??';
////        else if (allele2 == DELETION)    //?-
////            return ('-');
//        if (iupac == 'a')        //-A
//            result='-A';
//        else if (iupac == 'c')    //-C
//            result='-C';
//        else if (iupac == 'g')    //-G
//            result='-G';
//        else if (iupac == 't')    //-T
//            result='-T';
//        else if (iupac == ADO)    //-?
//            result='-?';
////        else if (iupac == DELETION)    //--
////            return ('-');
//       else
//            result='--';
//    
//    return(result);
//}

/********************* WhichNuc ************************/
/* Returns character representation for nucleotides */

char WhichNuc (int nucleotide)
{
    if (nucleotide == A)
        return ('A');
    else if (nucleotide == C)
        return ('C');
    else if (nucleotide == G)
        return ('G');
    else if (nucleotide == T)
        return ('T');
    else if (nucleotide == ADO)
        return ('?');
    else if (nucleotide == DELETION)
        return ('-');
    else
        return ('N');
}

/********************* WhichConsensusBinary ************************/
/* Returns a consensus representation of the binary genotype */
/*
 0/0 => 0
 0/1 => 1
 1/0 => 1
 1/1 => 2
 
 0/_ => 0
 _/0 => 0
 
 1/_ => 2
 _/1 => 2
 
 _/_ => -
 */

char WhichConsensusBinary (int allele1, int allele2)
{
    if (allele1 == 0)
    {
        if (allele2 == 0)        //00
            return ('0');
        else if (allele2 == 1)    //01
            return ('1');
        else if (allele2 == ADO)    //0?
            return ('0');
        else if (allele2 == DELETION)    //0-
            return ('0');
        else
            return ('N');
    }
    else if (allele1 == 1)
    {
        if (allele2 == 0)        //10
            return ('1');
        else if (allele2 == 1)    //11
            return ('2');
        else if (allele2 == ADO)    //1?
            return ('2');
        else if (allele2 == DELETION)    //0-
            return ('2');
        else
            return ('N');
    }
    else if (allele1 == ADO)
    {
        if (allele2 == 0)        //?0
            return ('0');
        else if (allele2 == 1)    //?1
            return ('2');
        else if (allele2 == ADO)    //??
            return ('-');
        else if (allele2 == DELETION)    //?-
            return ('-');
        else
            return ('N');
    }
    else if (allele1 == DELETION)
    {
        if (allele2 == 0)        //-0
            return ('0');
        else if (allele2 == 1)    //-1
            return ('2');
        else if (allele2 == ADO)    //-?
            return ('-');
        else if (allele2 == DELETION)    //--
            return ('-');
        else
            return ('N');
    }
    else
        return ('N');
}

/********************* WhichMut ************************/
/* Returns character representation for binary data */

char WhichMut (int state)
{
    if (state == 0)
        return ('0');
    else if (state == 1)
        return ('1');
    else if (state == ADO)
        return ('?');
    else if (state == DELETION)
        return ('-');
    else
        return ('N');
}
/********************* CompareGenotypes  ************************/
/*
 Compares two unphased genotypes (i.e, A/T = T/A)
 */
int CompareGenotypes (int a1, int a2, int b1, int b2)
{
    int temp, equal;
    
    temp = 0;
    equal = YES;
    
    /* order the alleles to facilitate the comparisons */
    if (a1 > a2)
    {
        temp = a1;
        a1 = a2;
        a2 = temp;
    }
    
    if (b1 > b2)
    {
        temp = b1;
        b1 = b2;
        b2 = temp;
    }
    
    if (a1 != b1 || a2 != b2)
        equal = NO;
    
    return equal;
}

/************************* CountTrueVariants  ************************/
/* Count number of variants in a given data set, assuming there is no ADO (yet).
 A variant is defined as no-deletion genotype different from the healthy root genotype */

int CountTrueVariants (TreeNode *nodes,  int numSites, int numCells, TreeNode *HEALTHY_ROOT, SiteStr* allSites, int* variantSites, int *SNVsites )
{
    int        cell, site;
    int        nVariants = 0;
    int numberDiff=0;
    int isVariant=NO;
    TreeNode* p;
    for (site=0; site<numSites; site++)
    {   isVariant=NO;
        numberDiff=0;
        for (cell=0; cell<numCells; cell++)
        {    p= (nodes +cell);
            if (p->maternalSequence[site] != DELETION && p->maternalSequence[site] != HEALTHY_ROOT->maternalSequence[site])
            {
               isVariant=YES;
                numberDiff++;
            }
            if (p->paternalSequence[site] != DELETION && p->paternalSequence[site] != HEALTHY_ROOT->paternalSequence[site])
            {
                isVariant=YES;
               numberDiff++;
            }
        }
        if (isVariant){
            allSites[site].isVariant = YES;
            SNVsites[ nVariants] = site;
            nVariants++;
        }
        allSites[site].numberDiffReference=numberDiff;
    }
    return nVariants;
}
/************************* ComputeAvgHeterocigocity  ************************/
/* Identify reference and alternate alleles plus SNVs, in ingroup plus outgroup cell genotypes
 after introducing mutations, cnLOH, deletions, ADO and genotype errors
 A SNV is defined as no-deletion genotype different from the reference genotype */

double ComputeAvgHeterocigocity (int numSites, int numCells, TreeNode **treeTips, TreeNode *treeRoot, SiteStr* allSites, int* SNVsites)
{
    int        numAltAlleles, cell, site;
    int        countA, countC, countG, countT, countADO, countDEL, nSNVs;
    TreeNode  *p;
    nSNVs = 0;
    double sum=0;
    double HSS=0.0;
    for (site=0; site<numSites; site++)
    {
        countA = countC = countG = countT = countADO = countDEL = 0;
        for (cell=0; cell<numCells; cell++)
        {
            p= *(treeTips +cell);
            
            if (p->maternalSequence[site]  == A)
                countA++;
            else if (p->maternalSequence[site]  == C)
                countC++;
            else if (p->maternalSequence[site]  == G)
                countG++;
            else if (p->maternalSequence[site]  == T)
                countT++;
            else if (p->maternalSequence[site]  == ADO)
                countADO++;
            else if (p->maternalSequence[site]  == DELETION)
                countDEL++;
            
            if (p->paternalSequence[site]  == A)
                countA++;
            else if (p->paternalSequence[site]  == C)
                countC++;
            else if (p->paternalSequence[site]  == G)
                countG++;
            else if (p->paternalSequence[site]  == T)
                countT++;
            else if (p->paternalSequence[site]  == ADO)
                countADO++;
            else if (p->paternalSequence[site]  == DELETION)
                countDEL++;
        }
        if (treeRoot->paternalSequence[site]  == A)
        countA++;
        else if (treeRoot->paternalSequence[site]  == C)
        countC++;
        else if (treeRoot->paternalSequence[site]  == G)
        countG++;
        else if (treeRoot->paternalSequence[site]  == T)
        countT++;
        else if (treeRoot->paternalSequence[site]  == ADO)
        countADO++;
        else if (treeRoot->paternalSequence[site]  == DELETION)
        countDEL++;
        
        allSites[site].countA = countA;
        allSites[site].countC = countC;
        allSites[site].countG = countG;
        allSites[site].countT = countT;
        
        sum = sum + (1 -  (allSites[site].countA / (numCells +1)) ^2);
        sum = sum + (1 -  (allSites[site].countC / (numCells +1)) ^2);
        sum = sum + (1 -  (allSites[site].countG / (numCells +1)) ^2);
        sum = sum + (1 -  (allSites[site].countT / (numCells +1)) ^2);
    
        allSites[site].countACGT = countA + countC + countG + countT;
        allSites[site].countDropped = countADO + countDEL;
        
        /* count number of alternate alleles, ignoring ADO or DELETION */
        numAltAlleles = 0;
        if (countA > 0 && allSites[site].referenceAllele != A)
            allSites[site].alternateAlleles[numAltAlleles++] = A;
        if (countC > 0 && allSites[site].referenceAllele != C)
            allSites[site].alternateAlleles[numAltAlleles++] = C;
        if (countG > 0 && allSites[site].referenceAllele != G)
            allSites[site].alternateAlleles[numAltAlleles++] = G;
        if (countT > 0 && allSites[site].referenceAllele != T)
            allSites[site].alternateAlleles[numAltAlleles++] = T;
        allSites[site].numAltAlleles = numAltAlleles;
        
        /* find out whether this site is a SNV */
        if (numAltAlleles > 0)
        {
            allSites[site].isSNV = YES;
            SNVsites[nSNVs++] = site;
        }
    }
    HSS= (numCells +1)*  sum    /(numCells * nSNVs);
    return HSS;
}
/************************ PrintSNVGenotypes ***********************/
/* Prints oberved/ML genotypes at variable sites (SNVs) to a file */
 void PrintSNVGenotypes (FILE *fp, int numClones, Population **populations, TreeNode* nodes,TreeNode* treeRoot, int doPrintAncestors, int doNGS, int numCells, int numSNVs, int* SNVsites, int alphabet, int doUserTree, int ***data, char **cellNames, int TUMOR_ROOT, int HEALTHY_ROOT, CellStr* cell, int* SFS, int *numberDifferences, int TajimaD, double mutationRate)
{
    int        i, j, k;
    TreeNode *p;
    char *temp;
    Population *popI;
//    if (doPrintAncestors == YES && doNGS == NO)
//        fprintf (fp, "%d %d\n", numCells+3, numSNVs);
//    else
//        fprintf (fp, "%d %d\n", numCells+1, numSNVs);
    
    /* site information */
    /* fprintf (fp, "%14s", ""); */
//
//    for (i=0; i<numSNVs; i++){
//
//        fprintf (fp, "%d ", SNVsites[i]+1);
//
//    }
//    fprintf (fp, "\n");
//    for (i=0; i<numSNVs; i++){
//        posSNV =numberDifferences[i];
//        fprintf (fp, "%d ", posSNV);
//        //        fprintf (fp, "%d ", SFS[ SNVsites[i]]);
//
//    }
//    fprintf (fp, "\n");
//     qsort(SFS, numSNVs, sizeof(int), compareIntDescending);
    //for( i = 0 ; i < numClones; i++)
    for( i = 0 ; i < numClones; i++)
    {
        popI=*(populations + i);
         fprintf (fp, "Delta=%lf \ntOriginSTD=%lf \ngrowthRate=%lf \nPopSize=%d \neffectPopSize=%lf \nmutationRate=%10.7lf", popI->delta, popI->timeOriginSTD, popI->growthRate,popI->popSize, popI->effectPopSize, mutationRate);
        fprintf (fp, "\n");
    }
    
   
    //for (i=1; i<numSNVs; i++){
    for (i=1; i< numCells; i++){
//        if (SFS[i]>0)
           fprintf (fp, "%d ", SFS[i]);
//        else
//            break;
        //        fprintf (fp, "%d ", SFS[ SNVsites[i]]);
    }
//    fseek (fp, -1, SEEK_CUR);
   fprintf (fp, "\n");
    fprintf (fp, "Avg Heterocigocity=%d", TajimaD);
    fprintf (fp, "\n");
    k=0;
    if (alphabet == DNA)
    {
////        fprintf (fp,"%d %d\n",2*(numCells+1), numSites);
//            for (i=0; i<numCells; i++){
//                p = (nodes + i);
//                /* print IUPAC haplotype */
//                if (p !=NULL){
//                    if (p->left==NULL && p->right ==NULL){
//                        temp=p->cellName;
//                        fprintf (fp,"%-12s ", temp);
//                        for (j=0; j<numSNVs; j++){
//                             fprintf (fp, " %c%c", WhichNuc(p->maternalSequence[SNVsites[j]]), WhichNuc(p->paternalSequence[SNVsites[j]]));
//                        }
//                        fprintf (fp,"\n");
//                    }
//                }
//            }
//            fprintf (fp,"%-12s ", treeRoot->cellName);
//            for (j=0; j<numSNVs; j++)
//                  fprintf (fp, " %c%c", WhichNuc(treeRoot->maternalSequence[SNVsites[j]]), WhichNuc(treeRoot->paternalSequence[SNVsites[j]]));
//            fprintf (fp,"\n");
        }
//        for (i=0; i<numCells+1; i++)
//        {
//            if (i == numCells)
//                fprintf (fp,"healthycell ");
//            else
//            {
//                if (doUserTree == NO)
//                    fprintf (fp,"tumcell%04d ", i+1);
//                else
//                    fprintf (fp,"%-12s", cellNames[i]);
//            }
//            for (j=0; j<numSNVs; j++)
//            {
//                if (doNGS == NO)
//                    fprintf (fp, " %c%c", WhichNuc(data[MATERNAL][i][SNVsites[j]]), WhichNuc(data[PATERNAL][i][SNVsites[j]]));
//                else if (cell[i].hasDoublet == NO)
//                    fprintf (fp, " %c%c", WhichNuc(cell[i].site[SNVsites[j]].MLmatAllele), WhichNuc(cell[i].site[SNVsites[j]].MLpatAllele));
//                else
//                    fprintf (fp, " %c%c", WhichNuc(cell[i].site[SNVsites[j]].MLmatAlleleDoublet), WhichNuc(cell[SNVsites[j]].site[j].MLpatAlleleDoublet));
//            }
//            fprintf (fp,"\n");
//        }
        
//        if (doPrintAncestors == YES & doNGS == NO)
//        {
//            fprintf (fp,"hearoot%04d ", i+1);
//            for (j=0; j<numSNVs; j++)
//                fprintf (fp, " %c%c", WhichNuc(data[MATERNAL][HEALTHY_ROOT][SNVsites[j]]),WhichNuc(data[PATERNAL][HEALTHY_ROOT][SNVsites[j]]));
//            fprintf (fp,"\ntumroot%04d ", i+1);
//            for (j=0; j<numSNVs; j++)
//                fprintf (fp, " %c%c", WhichNuc(data[MATERNAL][TUMOR_ROOT][SNVsites[j]]),WhichNuc(data[PATERNAL][TUMOR_ROOT][SNVsites[j]]));
//            fprintf (fp,"\n");
//        }
//    }
   else // binary data
    {
        for (i=0; i<numCells; i++){
            p = (nodes + i);
            if (p !=NULL){
                
                if (p->left==NULL && p->right ==NULL){
                    temp=p->cellName;
                    fprintf (fp,"%-12s ", temp);
                    for (j=0; j<numSNVs; j++)
                        fprintf (fp, " %c%c", WhichMut(p->maternalSequence[SNVsites[j]]),WhichMut(p->paternalSequence[SNVsites[j]]));
                    fprintf (fp,"\n");
                }
            }
        }
        fprintf (fp,"%-12s ", treeRoot->cellName);
        for (j=0; j<numSNVs; j++)
            fprintf (fp, " %c%c", WhichMut(treeRoot->maternalSequence[SNVsites[j]]),WhichMut(treeRoot->paternalSequence[SNVsites[j]]));
        fprintf (fp,"\n");
//        for (i=0; i<numCells; i++)
//        {
//            if (i == numCells)
//                fprintf (fp,"healthycell ");
//            else
//            {
//                if (doUserTree == NO)
//                    fprintf (fp,"tumcell%04d ", i+1);
//                else
//                   fprintf (fp,"%-12s", cellNames[i]);//uncomment this to write names of cells
//            }
//            for (j=0; j<numSNVs; j++)
//                 fprintf (fp, " %c%c", WhichMut(data[MATERNAL][i][SNVsites[j]]),WhichMut(data[PATERNAL][i][SNVsites[j]]));

//            fprintf (fp,"\n");
//        }
//        fprintf (fp,"%-12s", cellNames[HEALTHY_ROOT]);
//        for (j=0; j<numSNVs; j++)
//            fprintf (fp, " %c%c", WhichMut(data[MATERNAL][HEALTHY_ROOT][SNVsites[j]]),WhichMut(data[PATERNAL][HEALTHY_ROOT][SNVsites[j]]));
//
//        fprintf (fp,"\n");
        
//        if (doPrintAncestors == YES)
//        {
//            fprintf (fp,"hearoot%04d ", i+1);
//            for (j=0; j<numSNVs; j++)
//                fprintf (fp, " %d%d", data[MATERNAL][HEALTHY_ROOT][SNVsites[j]],data[PATERNAL][HEALTHY_ROOT][SNVsites[j]]);
//            fprintf (fp,"\ntumroot%04d ", i+1);
//            for (j=0; j<numSNVs; j++)
//                fprintf (fp, " %d%d", data[MATERNAL][TUMOR_ROOT][SNVsites[j]],data[PATERNAL][TUMOR_ROOT][SNVsites[j]]);
//            fprintf (fp,"\n");
//        }
    }
}
/************************ computeUnfoldedISMSFS ***********************/
/* compute unfolded version of SFS */
void computeUnfoldedISMSFS(int numSites,SiteStr* allSites,int numSNVs, int* SNVsites, int* SFS, int *numberDifferences){
    int numberDiff,j;
    int count=0;
    int site;
    int posSNV;
    for (site=0; site<numSites; site++){
        SFS[site]=0;
    }
    for (site=0; site<numSNVs; site++){
       // if (allSites[site].isVariant == YES){
             posSNV =SNVsites[site];
             numberDiff=allSites[posSNV].numberDiffReference;
             SFS[numberDiff]=SFS[numberDiff]+1;
             numberDifferences[site]=numberDiff;
        //}
    }
}
/************************ computeTajimaD ***********************/
/* computeTajimaD */
int computeTajimaD(TreeNode **treeTips,  int numSites, int numCells){
    int        cell, cell1, site;
    int numberPairDiff=0;
    int numberDiff=0;
    int isVariant=NO;
    TreeNode* p;
    TreeNode* p1;
    for (site=0; site<numSites; site++)
    {   isVariant=NO;
        numberPairDiff=0;
        for (cell=0; cell< numCells -1 ; cell++)
        {
             p= *(treeTips +cell);
             for (cell1=cell+1; cell1< numCells; cell1++)
              {
                  p1 = *(treeTips +cell1);
                  if (p->maternalSequence[site] != DELETION &&
                      p1->maternalSequence[site] != DELETION && p->maternalSequence[site] != p1->maternalSequence[site])
                  {
                      isVariant=YES;
                      numberPairDiff++;
                  }
                  if (p->paternalSequence[site] != DELETION
                      && p1->paternalSequence[site] != DELETION && p->paternalSequence[site] != p1->paternalSequence[site])
                  {
                      isVariant=YES;
                      numberPairDiff++;
                  }
              }
        }
        numberDiff= numberDiff+numberPairDiff;
    
    }
    int tajimaD= 2 * 2 * numberDiff/ (numCells * (numCells-1));
   // int tajimaD=  2 * numberDiff/ (numCells * (numCells-1));
    return tajimaD;
}
/************************ compareIntDescending ***********************/
/* function to compare two integers in descending order */
int compareIntDescending(const void *a, const void *b)
{
    const int *ia = (const int *)a;
    const int *ib = (const int *)b;
    return *ib  - *ia;

}
/********************* RandomNegativeBinomial ********************/
/*
 *    Random variates from the negative binomial distribution.
 *
 *  NOTES
 *
 *    x = the number of failures before the n-th success
 *
 *  REFERENCE
 *
 *    Devroye, L. (1986).
 *    Non-Uniform Random Variate Generation.
 *    New York:Springer-Verlag.  Pages 488 and 543.
 *
 *  METHOD
 *
 *    Generate lambda as gamma with shape parameter "dispersion" (aka size) and scale
 *    parameter "mean/dispersion".  Return a Poisson deviate with mean lambda.
 
 **** NOTE: Extracted from rnbinom.c R code:
 rpois(rgamma(size, (1 - prob) / prob));
 rpois(rgamma(size, mu / size));
 
 The negative binomial distribution with dispersion = n and prob = p has density
 
 p(x) = Gamma(x+n)/(Gamma(n) x!) p^n (1-p)^x
 
 for x = 0, 1, 2, ..., n > 0 and 0 < p <= 1.
 
 A negative binomial distribution can arise as a mixture of Poisson distributions with mean distributed as a Γ (pgamma) distribution with scale parameter (1 - prob)/prob and shape parameter dispersion. In this model prob = scale/(1+scale), and the mean is dispersion * (1 - prob)/prob. The variance in this parametrization is n (1-p)/p^2.
 
 The alternative parameterization, often used in ecology, and the one used here, is by the mean mu, and the dispersion parameter, where prob = dispersion/(dispersion+mu). The variance is mu + mu^2/dispersion in this parametrization.
 */

int RandomNegativeBinomial (double mean, double dispersion, long int *seed)
{
    int        poissonRand;
    double    gammaRand;
    /* the RandomGamma function here has mean 1, so we need to scale it ourselves */
    gammaRand = mean/dispersion * RandomGamma (dispersion, seed);
    poissonRand = RandomPoisson (gammaRand, seed);
    return poissonRand;
}


/************************RandomLogNormal ***********************/
/* generate random log normal */
double RandomLogNormal( double mean, double dispersion, long int *seed){
    
    return(exp(RandomNormal(mean, dispersion, seed)));
}
/************************RandomLogUniform ***********************/
/* generate random log uniform */
double RandomLogUniform( double from, double to, long int *seed){
    
    return(exp(from + RandomUniform(seed)*(to -from)));
}


/************************ generateGrowthRatefromPrior ***********************/
/* generate random growth rate from a prior */
double generateTimeofOriginfromPrior( double mean, double dispersion, long int *seed){
    // we would use a log normal prior
    return(exp(RandomNormal(mean, dispersion, seed)));
}
/************************ RandomNormal ***********************/
/* generate random from normal distribution */
double RandomNormal(double mean, double stddev,  long int* seed)
    {//Box muller method
        static double n2 = 0.0;
        static int n2_cached = 0;
        if (!n2_cached)
        {
            double x, y, r;
            do
            {
               // x = 2.0*rand()/RAND_MAX - 1;
              //  y = 2.0*rand()/RAND_MAX - 1;
                 x = 2.0*RandomUniform(seed) - 1;
                  y = 2.0*RandomUniform(seed)- 1;
                
                r = x*x + y*y;
            }
            while (r == 0.0 || r > 1.0);
            {
                double d = sqrt(-2.0*log(r)/r);
                double n1 = x*d;
                n2 = y*d;
                double result = n1*stddev + mean;
                n2_cached = 1;
                return result;
            }
        }
        else
        {
            n2_cached = 0;
            return n2*stddev + mean;
        }
    }
/************************ AvgHeterozygosity ***********************/
/* Average heterozygosity of all segregating sites */
double AvgHeterozygosity(int numSites,SiteStr* allSites,int numSNVs, int* SNVsites, int* SFS, TreeNode** treeTips, int numCells, int numSNPs){
    double result = 0;
    int i, j, cell, site;
    TreeNode *p;
    for (site=0; site<numSites; site++)
    {
        for (cell=0; cell<numCells; cell++)
        {    p= *(treeTips +cell);
        }
    }
    return result;
}
/************************ DistanceBetweenSFS ***********************/
/* distance between SFS */
double DistanceBetweenSFS(int* SFS1, int*SFS, int numSNVs , int numSites)
{
    int i;
    int sum=0;
    return sum;
}
/************************ ComputeESS ***********************/
/*  Computes the Effective Sample Size */
double ComputeESS(double *weights, int numberWeights){
    double ESS=0.0;
    int i;
    for (i=0; i<numberWeights; i++)
    {   ESS= ESS + 1.0 / (weights[i] * weights[i]);
    }
    return(ESS);
}
/************************ normalizeVector ***********************/
/*  normalizeVector */
void  normalizeVector(double *vector, int length){
    double sum=0.0;
    int i;
    for (i=0; i<length; i++)
    {   sum= sum +  vector[i] ;
    }
    for (i=0; i<length; i++)
    {    vector[i]= vector[i]/ sum ;
    }
}
/************************ ReadParametersFromFastaFile ***********************/
/*  ReadParametersFromFastaFile */
void ReadParametersFromFastaFile(char *fileName, ProgramOptions *programOptions){
    //read fasta
    FILE *fastaFile;
    kseq_t *seq;
    int l1;
    int current=0;
    char *currentSeq;
    char *currentQual;
    int  seqlength=0;
    char *t;
    int index,i;
    int max_length=0.0;
    int numberSeq;
    //    char script[80];
    //    strcpy(script,"grep -c" ) ;
    //    strcpy(script," \">\" ") ;
    //    strcpy(script,fileName) ;
    //
    //    numberSeq =system(script);
    
    if ((fastaFile = freopen(fileName, "r", stdin)) != NULL)
    {
        seq = kseq_init(fileno(fastaFile));
        numberSeq=0;
        while ((l1 = kseq_read(seq)) >= 0 )
        {
            //printf("name: %s\n", seq->name.s);
            //if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
            //printf("seq: %s\n", seq->seq.s);
            numberSeq=numberSeq+1;
            if ( l1 > max_length)
            {
                max_length=l1;
            }
        }
        //*numSites=max_length;
        programOptions->numSites=max_length;
        if (numberSeq >=1){
           // *numCells =numberSeq;//not counting the healthy cell
            programOptions->numCells=numberSeq;
        }
        else{
         //   *numCells =0;
             programOptions->numCells=numberSeq;
            programOptions->TotalNumSequences=numberSeq;
        }
        kseq_destroy(seq);
        fclose(fastaFile);
    }
    else{
        
        fprintf (stderr, "\nERROR: Can't read parameters file.");
        PrintUsage();
        
    }
}
/************************ ReadFastaFile ***********************/
/*  ReadFastaFile */
void ReadFastaFile(char *fileName, int** ObservedData,  char **ObservedCellNames, ProgramOptions *programOptions){
    FILE *fastaFile;
    kseq_t *seq;
    int l1;
    int current=0;
    char *currentSeq;
    char *currentQual;
    int  seqlength=0;
    char *t;
    int index,i;
    int max_length=0.0;
    int numberSeq;
   
    if ((fastaFile = freopen(fileName, "r", stdin)) != NULL){
        seq = kseq_init(fileno(fastaFile));
        while ((l1 = kseq_read(seq)) >= 0 ) {
            //printf("name: %s\n", seq->name.s);
            //strcpy( cellNames[current] , seq->name.s);
            ObservedCellNames[current] = malloc(MAX_NAME);
            if (ObservedCellNames[current] != NULL)
            {
                strcpy(ObservedCellNames[current], seq->name.s);
                //memcpy(ObservedCellNames[current], seq->name.s, sizeof(seq->name.s));
                
            }
            
            //if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
            //printf("seq: %s\n", seq->seq.s);
            seqlength=0;
            //if(current < *numSites){
            currentSeq=seq->seq.s;
               // for ( t= currentSeq; *t != '\0'; t++) {
            for ( index= 0; index < l1; index++) {
                    t= currentSeq+index;
                    if (programOptions->doUseGenotypes == NO) // use sequences
                        ObservedData[current][index]= WhichNucChar(*t);
                    else // use genotypypes
                        ObservedData[current][index]= WhichGenotypeChar(*t);
                        
                     seqlength++;
             //   }
                //seqlength= sizeof( seq->seq.s) / sizeof(char);
                 //dataFromFile[current++]= WhichNucChar(*temp);
            }
            current++;
            if (seq->qual.l){
                // printf("qual: %s\n", seq->qual.s);
                currentQual=seq->qual.s;
            }
        }
        //printf("return value: %d\n", l1);
        kseq_destroy(seq);
         //  gzclose(fastaFile);
        fclose(fastaFile);
    }
}
/************************ LogConditionalLikelihoodTree ***********************/
/*  LogConditionalLikelihoodTree */
double LogConditionalLikelihoodTree(TreeNode  *tree, TreeNode *nodes, Population **populations, int numClones)
{
    Population* popI;
    Population* popJ;
    Population* fatherPop;
    double product=0;
    int i, j;
    double temp;
     for ( i = 0; i < numClones; i++)
     {
         popI=*(populations + i );
         product = product + log( DensityTime(popI->delta, popI->timeOriginSTD));
     }
    for ( j = 0; j < numClones - 1; j++)
    {
        popJ = *(populations + j );
        product = product + log( popJ->popSize);
        fatherPop = popJ -> FatherPop;
        temp=popJ->timeOriginSTD * popJ->effectPopSize / fatherPop->effectPopSize;
        temp=CalculateH(popJ->timeOriginSTD * popJ->effectPopSize / fatherPop->effectPopSize, fatherPop->timeOriginSTD, fatherPop->delta);
        product = product  + log( temp);
    }
    //for ( i = 0; i < numClones; i++)
  //  {
   //     popI=*(populations + i );
        product = product + LogDensityCoalescentTimesForPopulation(tree, nodes, populations, numClones);
    //}
    return product;
}
/************************ LogDensityCoalescentTimesForPopulation ***********************/
/*  LogDensityCoalescentTimesForPopulation */
double LogDensityCoalescentTimesForPopulation(TreeNode  *tree, TreeNode *nodes,  Population **populations, int numClones)
{
    double result =0;
    int i, k;
    Population *popI;
    int numberLeftCoalescences;
    int numberLeftMigrations;
    int numberAliveCells;
    int currentCoalescentEvent=0;
    int currentMigrationEvent=0;
    double temp;
    for ( i = 0; i < numClones; i++){
        popI = *(populations + i );
        currentCoalescentEvent=0;
        currentMigrationEvent=0;
        numberAliveCells= popI->sampleSize;
       // numberLeftCoalescences =  popI->numCompletedCoalescences; //in the numCompletedCoalescences we are considering also the migrations
        numberLeftCoalescences =  popI->numCompletedCoalescences - (popI->numIncomingMigrations-1);
        numberLeftMigrations = popI->numIncomingMigrations-1;
        //we are not counting time of origin as a migration
        if (numberLeftCoalescences ==0)
            return  result;
        while(numberLeftMigrations > 0)
        {
            while(popI->CoalescentEventTimes[currentCoalescentEvent] < popI->migrationTimes[currentMigrationEvent] && numberAliveCells > 1)
            {
                temp=log(numberAliveCells * (numberAliveCells-1)/2);
                result= result + temp;
                temp = log(1/CalculateH(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta));
                result= result + temp;
                temp =  (numberAliveCells/2)* (numberAliveCells-1)*(FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta)-FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent+1], popI->timeOriginSTD, popI->delta));
                result= result -temp;
                currentCoalescentEvent++;
                numberLeftCoalescences--;
                numberAliveCells--;
            }
            //if (numberLeftMigrations > 0 && numberAliveCells > 1)// if there are migrations
            if (numberLeftMigrations > 0 )// if there are migrations
            {  temp= LogProbNoCoalescentEventBetweenTimes(popI,popI->CoalescentEventTimes[currentCoalescentEvent],popI-> migrationTimes[currentMigrationEvent], numberAliveCells );
                result= result+ temp;
                numberLeftMigrations--;
                currentMigrationEvent++;
            }
        }
        //here there are only coalescents events left(al least one event)
        while(numberLeftCoalescences > 0 && numberAliveCells > 1)
        {   temp = log(numberAliveCells * (numberAliveCells-1)/2);
            result= result + temp;
            temp = log(1/CalculateH(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta));
            result= result + temp;
            temp=( numberAliveCells/2)* (numberAliveCells-1)*(FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta)-FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent+1], popI->timeOriginSTD, popI->delta));
            result= result -  temp;
            currentCoalescentEvent++;
            numberLeftCoalescences--;
            numberAliveCells--;
        }
    }
   
    return result;
}
/************************ LogProbNoCoalescentEventBetweenTimes ***********************/
/*  LogProbNoCoalescentEventBetweenTimes */
double LogProbNoCoalescentEventBetweenTimes(Population *popI,double from, double to, int numberActiveInd)
{   int j=numberActiveInd;
    double result=0.0;
//    result= (j * (j-1)/2) * exp(-1 * j* (j-1)*(FmodelTstandard(to,popI->timeOriginSTD, popI->delta)-FmodelTstandard(from, popI->timeOriginSTD, popI->delta))/2);
    result=  -1 * j* (j-1)*(FmodelTstandard(to,popI->timeOriginSTD, popI->delta)-FmodelTstandard(from, popI->timeOriginSTD, popI->delta))/2;
    return result;
}
/************************ InitPopulationSampleSizes ***********************/
/*  InitPopulationSampleSizes */
static void InitPopulationSampleSizes(Population **populations, int TotalSampleSize, int numClones, double *proportionsVector, long int *seed)
{
    int i,j;
    Population *popI, *popJ;
    double rand;
    double    *cumSum = (double *) calloc((numClones +1), (long) sizeof(double));
    if (!cumSum)
    {
        fprintf (stderr, "Could not allocate Xvector (%lu bytes)\n", (numClones +1 ) * (long) sizeof(double));
        exit (-1);
    }
    cumSum[0]=0.0;
    for (j = 1; j <= numClones; j++)
    {
        cumSum[j]=0;
        cumSum[j]=cumSum[j-1]+proportionsVector[j-1];
        popJ = *(populations + j - 1);
        popJ->sampleSize=0;
    }
    for (i = 0; i < TotalSampleSize; i++)
    {
        rand=RandomUniform(seed);
        for (j = 1; j <= numClones;j ++)
        {
            popJ = *(populations + j - 1);
            if (rand <= cumSum[j] && rand > cumSum[j - 1])
            {
                popJ->sampleSize=popJ->sampleSize +1;
                break;
            }
        }
    }
    free(cumSum);
    cumSum=NULL;
}
/************************ LogConditionalLikelihoodSequences ***********************/
/*  LogConditionalLikelihoodSequences */


static double  LogConditionalLikelihoodSequences(Chain *chain, pll_msa_t * msa, char* NewickString, ProgramOptions *programOptions){

      FILE    *outputShell;
      char script[80];
      char buf[1000];

    //pll_state_t pll_map_gt10_2[256];
    if (NewickString == NULL)
    {
        fprintf (stderr, "\nERROR: The newick representation of the tree cannot be empty\n\n");
         PrintUsage();
        return 0;
    }
    
    unsigned int i;
    unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
    pll_partition_t * partition;
    //pll_utree_t *unrootedTree = pll_utree_parse_newick_unroot(TreFileName);
    pll_utree_t *unrootedTree = pll_utree_parse_newick_string_unroot(NewickString);
 
    /* compute node count information */
    tip_nodes_count = unrootedTree->tip_count;
    inner_nodes_count = unrootedTree->inner_count;
    nodes_count = inner_nodes_count + tip_nodes_count;
    branch_count = unrootedTree->edge_count;

      pllmod_treeinfo_t * treeinfo = pllmod_treeinfo_create(unrootedTree->vroot,
                                                          tip_nodes_count,
                                                          1,
                                                          PLLMOD_COMMON_BRLEN_LINKED);
    pll_unode_t ** tipnodes = unrootedTree->nodes;
    
    /* create a libc hash table of size tip_nodes_count */
    hcreate(tip_nodes_count);
    
    /* populate a libc hash table with tree tip labels */
    unsigned int * data = (unsigned int *)malloc(tip_nodes_count *
                                                sizeof(unsigned int));
    char *label;
    for (i = 0; i < tip_nodes_count; ++i)
    {
        data[i] = tipnodes[i]->clv_index;
        ENTRY entry;
        label= tipnodes[i]->label;
        entry.key = tipnodes[i]->label;
         entry.data = (void *)(data+i);
        hsearch(entry, ENTER);
    }
  
  
    if (msa !=NULL)
       printf("Original sequence (alignment) length : %d\n", msa->length);
    else{
        fprintf (stderr, "\nERROR: The multiple sequence alignment is empty\n\n");
        PrintUsage();
        return 0;
    }
   
   

    pllmod_subst_model_t * model = pllmod_util_model_info_genotype(GT_MODEL);
    
    /* create the PLL partition instance
     
     tip_nodes_count : the number of tip sequences we want to have
     inner_nodes_count : the number of CLV buffers to be allocated for inner nodes
     model->states : the number of states that our data have
     1 : number of different substitution models (or eigen decomposition)
     to use concurrently (i.e. 4 for LG4)
     branch_count: number of probability matrices to be allocated
     RATE_CATS : number of rate categories we will use
     inner_nodes_count : how many scale buffers to use
     PLL_ATTRIB_ARCH_AVX : list of flags for hardware acceleration
     */
    partition = pll_partition_create(tip_nodes_count,
                                     inner_nodes_count,
                                     model->states,
                                     (unsigned int)(msa->length),
                                     1,
                                     branch_count,
                                     RATE_CATS,
                                     inner_nodes_count,
                                     PLL_ATTRIB_ARCH_AVX);
    
    set_partition_tips(chain, partition, msa, programOptions);
    
    treeinfo = pllmod_treeinfo_create(unrootedTree->vroot,
                                      tip_nodes_count,
                                      1,
                                      PLLMOD_COMMON_BRLEN_LINKED);
      int params_to_optimize = 0;
      unsigned int params_indices[RATE_CATS] = {0};
    
    int retval = pllmod_treeinfo_init_partition(treeinfo,
                                                0,
                                                partition,
                                                params_to_optimize,
                                                PLL_GAMMA_RATES_MEAN,
                                                1.0, /* alpha*/
                                                params_indices, /* param_indices */
                                                model->rate_sym /* subst matrix symmetries*/
                                                );
    
     double * empirical_frequencies = pllmod_msa_empirical_frequencies1(partition);
    
    unsigned int * weight = pll_compress_site_patterns(msa->sequence,
                                                       pll_map_gt10,
                                                       tip_nodes_count,
                                                       &(msa->length));
    printf("Number of unique site patterns: %d\n\n", msa->length);
    
    /* initialize the array of base frequencies  AA CC GG TT AC/CA AG/GA AT/TA CG/GC CT/TC GT/TG  */
    double user_freqs[10] = { 0.254504, 0.157238, 0.073667, 0.287452, 0.027336,
        0.058670, 0.041628, 0.024572, 0.064522, 0.010412 };
    
    
    
    //computeGenotypesFreq(user_freqs,  msa);
    
    
    /* substitution rates: for GTR4 model those are 6 "regular" DNA susbt. rates + 1 rate
     * for "unlikely" double substitutions (eg A/A -> C/T) */
    double unique_subst_rates[7] = { 0.001000, 0.101223, 0.001000, 0.001000, 1.000000,
        0.001000, 0.447050 };
    
    /* get full above-diagonal half-matrix */
    double * user_subst_rates = expand_uniq_rates(model->states, unique_subst_rates,
                                                  model->rate_sym);
    
    double rate_cats[RATE_CATS] = {0};
    
    /* compute the discretized category rates from a gamma distribution
     with alpha shape 1 and store them in rate_cats  */
    pll_compute_gamma_cats(1, RATE_CATS, rate_cats, PLL_GAMMA_RATES_MEAN);
    
    /* set frequencies at model with index 0 (we currently have only one model) */
    //pll_set_frequencies(partition, 0, model->freqs ? model->freqs : user_freqs);
    pll_set_frequencies(partition, 0, model->freqs ? model->freqs : empirical_frequencies);
    
    /* set substitution parameters at model with index 0 */
    pll_set_subst_params(partition, 0, model->rates ? model->rates : user_subst_rates);
    free(user_subst_rates);
    
    /* set rate categories */
    pll_set_category_rates(partition, rate_cats);
    
    /* set pattern weights and free the weights array */
    pll_set_pattern_weights(partition, weight);
    free(weight);
    
    //set_partition_tips(chain, partition, msa, programOptions);
    
    // pll_msa_destroy(msa);
    
    /* destroy hash table */
    hdestroy();
    /* we no longer need these two arrays (keys and values of hash table... */
    free(data);
    
    /* update matrix_count probability matrices using the rate matrix with
     index 0. The i-th matrix (i ranges from 0 to matrix_count - 1) is
     generated using branch length branch_lengths[i] and rate matrix
     (substitution rates + frequencies) params_indices[i], and can be refered
     to with index matrix_indices[i] */
    //unsigned int params_indices[RATE_CATS] = {0};
    
    /* we do not want to optimize anything */
    //int params_to_optimize = 0;
    
    /* create treeinfo structure */
//     treeinfo = pllmod_treeinfo_create(unrootedTree->vroot,
//                                                          tip_nodes_count,
//                                                          1,
//                                                          PLLMOD_COMMON_BRLEN_LINKED);
//
//    int retval = pllmod_treeinfo_init_partition(treeinfo,
//                                                0,
//                                                partition,
//                                                params_to_optimize,
//                                                PLL_GAMMA_RATES_MEAN,
//                                                1.0, /* alpha*/
//                                                params_indices, /* param_indices */
//                                                model->rate_sym /* subst matrix symmetries*/
//                                                );
    
    if (!retval)
        fprintf(stderr, "Error initializing partition!");
    /* Compute initial LH of the starting tree */
    double loglh = pllmod_treeinfo_compute_loglh(treeinfo, 0);
    pllmod_treeinfo_destroy(treeinfo);
    double * clv ;
    int s, clv_index, j, k;
    int scaler_index = PLL_SCALE_BUFFER_NONE;
    unsigned int * scaler = (scaler_index == PLL_SCALE_BUFFER_NONE) ?
    NULL : partition->scale_buffer[scaler_index];
    unsigned int states = partition->states;
    unsigned int states_padded = partition->states_padded;
    unsigned int rates = partition->rate_cats;
    double prob;
    unsigned int *site_id = 0;
    int index;
    unsigned int float_precision;
    
//    for (clv_index = 0 ; clv_index < 2*tip_nodes_count -1 ; ++clv_index)
//    {
//        clv = partition->clv[clv_index];
//        if (pll_repeats_enabled(partition) && partition->repeats->pernode_ids[clv_index]) {
//            site_id = partition->repeats->pernode_site_id[clv_index];
//        }
//        if ((clv_index < partition->tips) &&
//            (partition->attributes & PLL_ATTRIB_PATTERN_TIP))
//            return 0;
//        for (s = 0; s < partition->sites; ++s)
//        {
//            i = site_id ? site_id[s] : s;
//            for (j = 0; j < rates; ++j)
//            {
//                for (k = 0; k < states-1; ++k)
//                {
//                    prob = clv[i*rates*states_padded + j*states_padded + k];
//                    if (scaler) unscale(&prob, scaler[i]);
//                }
//                index = i*rates*states_padded + j*states_padded + k;
//               // if (index >=0 && index <= states-1)
//                  prob = clv[index];
//
//                if (scaler) unscale(&prob, scaler[i]);
//                // if (j < rates - 1) printf(",");
//            }
//        }
//    }
//
    /* destroy all structures allocated for the concrete PLL partition instance */
    pll_partition_destroy(partition);
    /* we will no longer need the tree structure */
    //pll_utree_destroy(unrootedTree, NULL);
    destroyTree(unrootedTree, NULL);
    pllmod_util_model_destroy(model);
//    if ((fp = freopen("genotype", "r", stdin)) == NULL)
//    {
//        strcpy(script,"./genotype ") ;
//        strcpy(script, TreFileName) ;
//        strcpy(script, SeqFileName) ;
//        outputShell =  popen( script, "r");
//        if (NULL == (outputShell = popen(script, "r")))
//        {
//
//        }
//
//        while(fgets(buf, sizeof(buf), outputShell) != NULL)
//        {
//
//        }
//        // fscanf(outputShell,....);
//        pclose(outputShell);
//    }
    return loglh;
}
/************************ destroyTree ***********************/
/*  destroyTree */
void  destroyTree(pll_utree_t * tree, void (*cb_destroy)(void *))
{
    
        unsigned int i;
        
        /* deallocate tip nodes */
        for (i = 0; i < tree->tip_count; ++i)
        {
            dealloc_data(tree->nodes[i], cb_destroy);
         
            free(tree->nodes[i]);
        }
        /* deallocate inner nodes */
        for (i = tree->tip_count; i < tree->tip_count + tree->inner_count; ++i)
        {
            pll_unode_t * first = tree->nodes[i];
            assert(first);
            if (first->label)
                free(first->label);
            
            pll_unode_t * node = first;
            do
            {
                pll_unode_t * next = node->next;
                dealloc_data(node, cb_destroy);
                free(node);
                node = next;
            }
            while(node && node != first);
        }
        
        /* deallocate tree structure */
        free(tree->nodes);
        free(tree);
   
}
/************************ unscale ***********************/
/*  unscale  from Pll_module library*/
static void unscale(double * prob, unsigned int times)
{
    unsigned int i;
    
    for (i = 0; i < times; ++i)
        *prob *= PLL_SCALE_THRESHOLD;
}
/************************ dealloc_data ***********************/
/*  dealloc_data  from Pll_module library*/
void dealloc_data(pll_unode_t * node, void (*cb_destroy)(void *))
{
    if (node->data)
    {
        if (cb_destroy)
            cb_destroy(node->data);
    }
}
/************************ compute_asc_bias_correction ***********************/
/*  compute_asc_bias_correction */
static double compute_asc_bias_correction(double logl_base,
                                          unsigned int sum_w,
                                          unsigned int sum_w_inv,
                                          int asc_bias_type)
{
    double logl_correction = 0.0;
    switch (asc_bias_type)
    {
        case PLL_ATTRIB_AB_LEWIS:
            logl_correction = -(sum_w*log(1 - logl_base));
            break;
        case PLL_ATTRIB_AB_STAMATAKIS:
            /* no need to add anything here */
            logl_correction = logl_base;
            break;
        case PLL_ATTRIB_AB_FELSENSTEIN:
            logl_correction = sum_w_inv * log(logl_base);
            break;
        default:
           // pll_errno = PLL_ERROR_AB_INVALIDMETHOD;
            //snprintf(pll_errmsg, 200, "Illegal ascertainment bias algorithm");
            return -INFINITY;
    }
    return logl_correction;
}
/********************* PrepareSeparateFiles **********************/
/* Open individual files to output results */

void PrepareSeparateFiles(int ChainNumber, int paramSetNumber, int replicate,const FilePaths *filePaths, const ProgramOptions *programOptions,Files *files)
{
    char File[MAX_NAME];
    char dir[MAX_NAME];
    /* contains the simulated tree in Newick format
     */
    mkdir("Results", S_IRWXU); /* Create "Results" folder (with type S_IRWXU (read, write and execute)) */
    //mkdir("Results",0);
    #ifdef MAC
    strcpy (dir, ":Results:"); /* Copy the string in char variable dir = Results (char), is different mac vs windows */
    #else
      strcpy (dir, "Results/");
    #endif
    //strcpy (resultsDir, dir);
    if (programOptions->doPrintTrees == YES)
    {
        if (programOptions->doSimulateData ==YES)
        {
            sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->treeDir );
        }
        else
        {
            sprintf(File,"%s/%s_Chain%d", filePaths->resultsDir, filePaths->treeDir, ChainNumber );
            
        }
        mkdir(File,S_IRWXU);
        
        if (programOptions->doSimulateData ==YES)
        {
        sprintf(File,"%s/%s/%s_%04d_%04d.tre", filePaths->resultsDir, filePaths->treeDir, filePaths->treeFile, paramSetNumber+1, replicate+1);
        }
        else
        {
           sprintf(File,"%s/%s_Chain%d/%s_%04d_%04d.tre", filePaths->resultsDir, filePaths->treeDir,ChainNumber , filePaths->treeFile, paramSetNumber+1, replicate+1);
        }
        //if ((*fpTrees = fopen(File, "w")) == NULL)
         if (openFile(&files->fpTrees, File) == -1)
        {
            fprintf (stderr, "Can't open \"%s\"\n", File);
            exit(-1);
        }
        
        sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->treeDir);
        
        if (programOptions->doSimulateData ==YES)
        {
        sprintf(File,"%s/%s/%s_2_%04d_%04d.tre", filePaths->resultsDir, filePaths->treeDir, filePaths->treeFile, paramSetNumber+1, replicate+1);
        }
        else
        {
             sprintf(File,"%s/%s_Chain%d/%s_2_%04d_%04d.tre", filePaths->resultsDir, filePaths->treeDir, ChainNumber, filePaths->treeFile, paramSetNumber+1, replicate+1);
            
        }
         //sprintf(File,"%s/%s/%s_2_%04d.tre", resultsDir, treeDir, treeFile, replicate+1);
        //if ((*fpTrees2 = fopen(File, "w")) == NULL)
        if (openFile(&files->fpTrees2, File) == -1)
        {
            fprintf(stderr, "Can't open %s.\n", File);
            exit(-1);
        }
    }
    
    if (programOptions->doPrintTimes == YES)
    {
        if (programOptions->doSimulateData ==YES)
        {
        sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->timesDir);
        }
        else
        {
            sprintf(File,"%s/%s_Chain%d", filePaths->resultsDir, filePaths->timesDir, ChainNumber );
        }
        mkdir(File,S_IRWXU);
        if (programOptions->doSimulateData ==YES)
        {
        sprintf(File,"%s/%s/%s_%04d_%04d.txt", filePaths->resultsDir, filePaths->timesDir, filePaths->timesFile, paramSetNumber+1, replicate+1);
        }
        else
        {
            sprintf(File,"%s/%s_Chain%d/%s_%04d_%04d.txt", filePaths->resultsDir, filePaths->timesDir, ChainNumber, filePaths->timesFile, paramSetNumber+1, replicate+1);
            
        }
        //if ((*fpTimes = fopen(File, "w")) == NULL)
    if (openFile(&files->fpTimes, File) == -1)
        {
            fprintf (stderr, "Can't open \"%s\"\n", File);
            exit(-1);
        }
        sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->timesDir);
        if (programOptions->doSimulateData ==YES)
        {
        sprintf(File,"%s/%s/%s_2_%04d_%04d.txt", filePaths->resultsDir, filePaths->timesDir, filePaths->timesFile, paramSetNumber+1, replicate+1);
        }
        else
        {
             sprintf(File,"%s/%s_Chain%d/%s_2_%04d_%04d.txt", filePaths->resultsDir, filePaths->timesDir, ChainNumber,filePaths->timesFile, paramSetNumber+1, replicate+1);
            
        }
        // sprintf(File,"%s/%s/%s_2_%04d.txt", resultsDir, timesDir, timesFile, replicate+1);
        
        //if ((*fpTimes2 = fopen(File, "w")) == NULL)
        if (openFile(&files->fpTimes2, File) == -1)
        {
            fprintf(stderr, "Can't open %s.\n", File);
            exit(-1);
        }
    }
    
//    if (doSimulateData == YES)
//    {
//        /* contains SNV genotypes for every cell */
//        if (doPrintSNVgenotypes == YES)
//        {
//            sprintf(File,"%s/%s", resultsDir, SNVgenotypesDir);
//            mkdir(File,S_IRWXU);
//            sprintf(File,"%s/%s/%s_%04d_%04d.txt", resultsDir, SNVgenotypesDir, SNVgenotypesFile, paramSetNumber+1, replicate+1);
//            //sprintf(File,"%s/%s/%s.%04d", resultsDir, SNVgenotypesDir, SNVgenotypesFile, replicate+1);
//            if ((*fpSNVgenotypes = fopen(File, "w")) == NULL)
//            {
//                fprintf (stderr, "Can't open \"%s\"\n", File);
//                exit(-1);
//            }
//        }
//
//        /* contains haplotypes for variable sites for every cell */
//        if (doPrintSNVhaplotypes == YES)
//        {
//            sprintf(File,"%s/%s", resultsDir, SNVhaplotypesDir);
//            mkdir(File,S_IRWXU);
//            sprintf(File,"%s/%s/%s_%04d_%04d.txt", resultsDir, SNVhaplotypesDir, SNVhaplotypesFile, paramSetNumber+1, replicate+1);
//           // sprintf(File,"%s/%s/%s.%04d", resultsDir, SNVhaplotypesDir, SNVhaplotypesFile, replicate+1);
//            if ((*fpSNVhaplotypes = fopen(File, "w")) == NULL)
//            {
//                fprintf (stderr, "Can't open \"%s\"\n", File);
//                exit(-1);
//            }
//        }
//
//        /* contains reference haplotypes (before errors)  for every cell */
//        if (doPrintTrueHaplotypes == YES)
//        {
//            sprintf(File,"%s/%s", resultsDir, trueHaplotypesDir);
//            mkdir(File,S_IRWXU);
//             sprintf(File,"%s/%s/%s_%04d_%04d.txt", resultsDir, trueHaplotypesDir, trueHaplotypesFile, paramSetNumber+1, replicate+1);
//
//
//           // sprintf(File,"%s/%s/%s.%04d", resultsDir, trueHaplotypesDir, trueHaplotypesFile, replicate+1);
//            if ((*fpTrueHaplotypes = fopen(File, "w")) == NULL)
//            {
//                fprintf (stderr, "Can't open \"%s\"\n", File);
//                exit(-1);
//            }
//        }
//
//        /* contains ML haplotypes  for every cell */
//        if (doPrintMLhaplotypes == YES)
//        {
//            sprintf(File,"%s/%s", resultsDir, MLhaplotypesDir);
//            mkdir(File,S_IRWXU);
//            sprintf(File,"%s/%s/%s_%04d_%04d.txt", resultsDir, MLhaplotypesDir, MLhaplotypesFile, paramSetNumber+1, replicate+1);
//           // sprintf(File,"%s/%s/%s.%04d", resultsDir, MLhaplotypesDir, MLhaplotypesFile, replicate+1);
//            if ((*fpMLhaplotypes = fopen(File, "w")) == NULL)
//            {
//                fprintf (stderr, "Can't open \"%s\"\n", File);
//                exit(-1);
//            }
//        }
//
//        /* contains all genotypes (variable or invariable) for every cell */
//        if (doPrintFullGenotypes == YES)
//        {
//            sprintf(File,"%s/%s", resultsDir, fullGenotypesDir);
//            mkdir(File,S_IRWXU);
//
//             sprintf(File,"%s/%s/%s_%04d_%04d.txt", resultsDir, fullGenotypesDir, fullGenotypesFile, paramSetNumber+1, replicate+1);
//            //sprintf(File,"%s/%s/%s.%04d", resultsDir, fullGenotypesDir, fullGenotypesFile, replicate+1);
//            if ((*fpFullGenotypes = fopen(File, "w")) == NULL)
//            {
//                fprintf (stderr, "Can't open \"%s\"\n", File);
//                exit(-1);
//            }
//        }
//
//        /* contains haplotypes for all sites for every cell */
//        if (doPrintFullHaplotypes == YES)
//        {
//            sprintf(File,"%s/%s", resultsDir, fullHaplotypesDir);
//            mkdir(File,S_IRWXU);
//
//              sprintf(File,"%s/%s/%s_%04d_%04d.txt", resultsDir, fullHaplotypesDir, fullHaplotypesFile, paramSetNumber+1, replicate+1);
//            //sprintf(File,"%s/%s/%s.%04d", resultsDir, fullHaplotypesDir, fullHaplotypesFile, replicate+1);
//            if ((*fpFullHaplotypes = fopen(File, "w")) == NULL)
//            {
//                fprintf (stderr, "Can't open \"%s\"\n", File);
//                exit(-1);
//            }
//        }
//
//        /* contains reads counts and log10 normalized genotype likelihoods for every SNV and cell */
//        if (doSimulateReadCounts == YES)
//        {
//            sprintf(File,"%s/%s", resultsDir, VCFdir);
//            mkdir(File,S_IRWXU);
//
//            sprintf(File,"%s/%s/%s_%04d_%04d.txt", resultsDir, VCFdir, VCFfile, paramSetNumber+1, replicate+1);
//            //sprintf(File,"%s/%s/%s.%04d", resultsDir, VCFdir, VCFfile, replicate+1);
//            if ((*fpVCF = fopen(File, "w")) == NULL)
//            {
//                fprintf (stderr, "Can't open \"%s\"\n", File);
//                exit(-1);
//            }
//        }
//
//        /* contains reads counts for every SNV and cell */
//        if (doPrintCATG == YES)
//        {
//            sprintf(File,"%s/%s", resultsDir, CATGdir);
//            mkdir(File,S_IRWXU);
//
//            sprintf(File,"%s/%s/%s_%04d_%04d.txt", resultsDir, CATGdir, CATGfile, paramSetNumber+1, replicate+1);
//           // sprintf(File,"%s/%s/%s.%04d", resultsDir, CATGdir, CATGfile, replicate+1);
//            if ((*fpCATG = fopen(File, "w")) == NULL)
//            {
//                fprintf (stderr, "Can't open \"%s\"\n", File);
//                exit(-1);
//            }
//        }
//    }
}
/********************* PrepareSeparateFilesGenotypes **********************/
/* Open individual genotypes files to output results */
void PrepareSeparateFilesGenotypes(int paramSetNumber, int TreeNum,int MutationAssignNum,
       const FilePaths *filePaths, const ProgramOptions *programOptions,Files *files)
{
    char File[MAX_NAME];
    char dir[MAX_NAME];
    /* contains the simulated tree in Newick format
     */
    mkdir("Results", S_IRWXU); /* Create "Results" folder (with type S_IRWXU (read, write and execute)) */
    //mkdir("Results",0);
    #ifdef MAC
       strcpy (dir, ":Results:"); /* Copy the string in char variable dir = Results (char), is different mac vs windows */
    #else
       strcpy (dir, "Results/");
    #endif
    //strcpy (resultsDir, dir);
    
    if (programOptions->doSimulateData == YES)
    {
        /* contains SNV genotypes for every cell */
        if (programOptions->doPrintSNVgenotypes == YES)
        {
            sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->SNVgenotypesDir);
            mkdir(File,S_IRWXU);
            sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths->resultsDir, filePaths->SNVgenotypesDir, filePaths->SNVgenotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
            //sprintf(File,"%s/%s/%s.%04d", resultsDir, SNVgenotypesDir, SNVgenotypesFile, replicate+1);
            //if ((*fpSNVgenotypes = fopen(File, "w")) == NULL)
            if (openFile(&files->fpSNVgenotypes, File) == -1)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        
        /* contains haplotypes for variable sites for every cell */
        if (programOptions->doPrintSNVhaplotypes == YES)
        {
            sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->SNVhaplotypesDir);
            mkdir(File,S_IRWXU);
            sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths->resultsDir, filePaths->SNVhaplotypesDir, filePaths->SNVhaplotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
            // sprintf(File,"%s/%s/%s.%04d", resultsDir, SNVhaplotypesDir, SNVhaplotypesFile, replicate+1);
            //if ((*fpSNVhaplotypes = fopen(File, "w")) == NULL)
            if (openFile(&files->fpSNVhaplotypes, File) == -1)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        
        /* contains reference haplotypes (before errors)  for every cell */
        if (programOptions->doPrintTrueHaplotypes == YES)
        {
            sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->trueHaplotypesDir);
            mkdir(File,S_IRWXU);
            sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths->resultsDir, filePaths->trueHaplotypesDir, filePaths->trueHaplotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
            
            
            // sprintf(File,"%s/%s/%s.%04d", resultsDir, trueHaplotypesDir, trueHaplotypesFile, replicate+1);
            //if ((*fpTrueHaplotypes = fopen(File, "w")) == NULL)
               if (openFile(&files->fpTrueHaplotypes, File) == -1)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        
        /* contains ML haplotypes  for every cell */
        if (programOptions->doPrintMLhaplotypes == YES)
        {
            sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->MLhaplotypesDir);
            mkdir(File,S_IRWXU);
            sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths->resultsDir, filePaths->MLhaplotypesDir, filePaths->MLhaplotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
            // sprintf(File,"%s/%s/%s.%04d", resultsDir, MLhaplotypesDir, MLhaplotypesFile, replicate+1);
           // if ((*fpMLhaplotypes = fopen(File, "w")) == NULL)
        if (openFile(&files->fpMLhaplotypes, File) == -1)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        
        /* contains all genotypes (variable or invariable) for every cell */
        if (programOptions->doPrintFullGenotypes == YES)
        {
            sprintf(File,"%s/%s", filePaths->resultsDir, filePaths
                    ->fullGenotypesDir);
            mkdir(File,S_IRWXU);
            
            sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths->resultsDir, filePaths->fullGenotypesDir, filePaths->fullGenotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
            //sprintf(File,"%s/%s/%s.%04d", resultsDir, fullGenotypesDir, fullGenotypesFile, replicate+1);
            //if ((*fpFullGenotypes = fopen(File, "w")) == NULL)
            if (openFile(&files->fpFullGenotypes, File) == -1)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        /* contains haplotypes for all sites for every cell */
        if (programOptions->doPrintFullHaplotypes == YES)
        {
            sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->fullHaplotypesDir);
            mkdir(File,S_IRWXU);
            
            sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths->resultsDir, filePaths->fullHaplotypesDir, filePaths->fullHaplotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
            //sprintf(File,"%s/%s/%s.%04d", resultsDir, fullHaplotypesDir, fullHaplotypesFile, replicate+1);
            //if ((*fpFullHaplotypes = fopen(File, "w")) == NULL)
            if (openFile(&files->fpFullHaplotypes, File) == -1)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        /* contains reads counts and log10 normalized genotype likelihoods for every SNV and cell */
        if (programOptions->doSimulateReadCounts == YES)
        {
            sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->VCFdir);
            mkdir(File,S_IRWXU);
            
            sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths->resultsDir, filePaths->VCFdir, filePaths->VCFfile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
            //sprintf(File,"%s/%s/%s.%04d", resultsDir, VCFdir, VCFfile, replicate+1);
            //if ((*fpVCF = fopen(File, "w")) == NULL)
         if (openFile(&files->fpVCF, File) == -1)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
        /* contains reads counts for every SNV and cell */
        if (programOptions->doPrintCATG == YES)
        {
            sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->CATGdir);
            mkdir(File,S_IRWXU);
            
            sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths->resultsDir, filePaths->CATGdir, filePaths->CATGfile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
            // sprintf(File,"%s/%s/%s.%04d", resultsDir, CATGdir, CATGfile, replicate+1);
            //if ((*fpCATG = fopen(File, "w")) == NULL)
            if (openFile(&files->fpCATG, File) == -1)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
    }
}
/********************* load_tree **********************/
/* load_tree(Function that uses code from genotype.c from pll-modules/examples) */
pll_utree_t * load_tree(const char *fname)
{
    /* parse the unrooted binary tree in newick format, and store the number
     of tip nodes in tip_nodes_count */
    //pll_utree_t * tree = pll_utree_parse_newick(fname);
    pll_utree_t * tree = pll_utree_parse_newick_unroot(fname);
    if (!tree)
        fprintf(stderr, "Tree must be an unrooted binary tree");
    
    /* fix all missing branch lengths (i.e. those that did not appear in the
     newick) to 0.000001 */
    pllmod_utree_set_length_recursive(tree, BRLEN_MIN, 1);
    
    return tree;
}
/********************* set_partition_tips **********************/
/* set_partition_tips(Function that uses code from genotype.c from pll-modules/examples) */
void set_partition_tips(Chain *chain, pll_partition_t * partition, pll_msa_t * msa, ProgramOptions *programOptions)
{
    //pll_state_t pll_map_gt10_2[256];
     int states =10;
    int i, currentState;
    int from, to;
  
 
    unsigned int state;
    //double * _freqs;
   
    /* find sequences in hash table and link them with the corresponding taxa */
    for (i = 0; i < msa->count; ++i)
    {
        ENTRY query;
        query.key = msa->label[i];
        ENTRY * found = NULL;
        
        found = hsearch(query,FIND);
        
        if (!found)
            fprintf(stderr,"Sequence with header %s does not appear in the tree", msa->label[i]);
        
        unsigned int tip_clv_index = *((unsigned int *)(found->data));
        
        if (programOptions->doUseGenotypes == NO)
        {
            pll_set_tip_states(partition, tip_clv_index, pll_map_gt10, msa->sequence[i]);
        }
        else
        {
            
//            for ( currentState = 0; currentState < states; currentState++)
//            {
//                from=1;
//                to=1;
                set_tipclv1( partition,
                            tip_clv_index,
                            pll_map_gt10,
                            msa->sequence[i],
                            chain->seqErrorRate,
                            chain->dropoutRate
                            );
//                 compute_state_probs( msa->sequence[i],  &(partition->clv[tip_clv_index]),  states, from, to,
//                                chain->seqErrorRate,
//                                 chain->dropoutRate);
      //      }
            
            //pll_set_tip_clv(partition, tip_clv_index, tipCLV, PLL_FALSE);
            
           
        }
    }
  
}
/********************* expand_uniq_rates **********************/
/* expand_uniq_rates(Function that uses code from genotype.c from pll-modules/examples) */
double * expand_uniq_rates(int states, const double * uniq_rates, const int * rate_sym)
{
    unsigned int i;
    unsigned int num_rates = states * (states-1) / 2;
    double * subst_rates = calloc(num_rates, sizeof(double));
    for (i = 0; i < num_rates; ++i)
        subst_rates[i] = rate_sym ? uniq_rates[rate_sym[i]] : uniq_rates[i];
    return subst_rates;
}
/********************* AssignTreeNodestoObservedCellNames **********************/
/* AssignTreeNodestoObservedCellNames*/
void AssignTreeNodestoObservedCellNames(TreeNode **treeRootInit, char* ObservedCellNames[])
{
    int currentCellName=0;
    AssignTreeNodestoObservedCellNamesRecursive ( treeRootInit[0],  ObservedCellNames, &currentCellName);
    
}
/********************* AssignTreeNodestoObservedCellNamesRecursive **********************/
/* AssignTreeNodestoObservedCellNames*/
void AssignTreeNodestoObservedCellNamesRecursive(TreeNode *p, char* ObservedCellNames[], int *currentCellName)
{
    if (p != NULL)
    {
        if(p->isOutgroup == YES)            /* Outgroup*/
        {
            strcpy( p->observedCellName,"healthycell");
        }
        else if (p->left == NULL && p->right == NULL)        /* tip of the tree */
        {
            strcpy( p->observedCellName,ObservedCellNames[*currentCellName]);
            *currentCellName = *currentCellName +1;
        }
        else
        {
            if(p->left != NULL)
               AssignTreeNodestoObservedCellNamesRecursive (p->left, ObservedCellNames, currentCellName);
            if(p->right != NULL)
               AssignTreeNodestoObservedCellNamesRecursive (p->right, ObservedCellNames, currentCellName);
        }
    }
}
/********************* getHealthyTip **********************/
/* getHealthyTip*/
TreeNode *getHealthyTip(TreeNode *treeRootInit)
{
    if (treeRootInit !=NULL && treeRootInit->right!=NULL)
        return treeRootInit->right;
    else
       return NULL;
}
/********************* SaveCurrentPartialTreeForPopulation **********************/
/* SaveCurrentPartialTreeForPopulation*/
void SaveCurrentPartialTreeForPopulation(Population *popI, FILE * fpFile)
{
    int i=0;
    int idTemp;
    for(int i=0; i< popI->numActiveGametes ; i ++)
    {
        idTemp = popI->idsActiveGametes[i];
    
        
    }
    
    //save the every sub tree into the file fpFile in newick format
}
/********************* runChain **********************/
/* runChain*/
void runChain(Chain *chain,  const  MCMCoptions *opt,  long int *seed, const FilePaths *filePaths, Files *files, const ProgramOptions *programOptions,
              double  *varTimeGMRCA,char* ObservedCellNames[], pll_msa_t * msa, int *sampleSizes
              ){
    Population** populations;
    TreeNode** nodes; TreeNode** treeTips;
    if (chain->numClones <= 0)
    {
        fprintf (stderr, "\nERROR: The number of clones cannot be negative.");
        PrintUsage();
        
    }
    Population *popI;
    double      cumNumCA, meanNumCA, cumNumMIG, meanNumMIG, numEventsTot, countTMRCA;
    int        numCA, numMIG;
    int i, j, k;
    double    *proportionsVector;
    long numCells;
    int currentIteration;
    double randomDelta;
    double meanNumSNVs, meanNumMU, meanNumDEL, meanNumCNLOH;
    double   cumNumSNVs, cumNumMU, cumNumDEL, cumNumCNLOH, cumCountMLgenotypeErrors;
    double cumNumMUSq, cumNumSNVsSq, cumNumDELSq, cumNumCNLOHSq;
    double varNumMU, varNumSNVs, varNumDEL, varNumCNLOH;
    double expNumMU, expVarNumMU;
    double logConditionalLikelihoodTree;
    double totalTreeLength;
    cumNumCA=0;
    meanNumCA=0;
    cumNumMIG=0;
    meanNumMIG=0;
    int numNodes=0;

    /* Computations per iteration in the chain */
   // for (currentIteration = 0; currentIteration < opt->Niterations; currentIteration++)
   // {
        //generate proposals for the variables
        
        if (programOptions->doPrintSeparateReplicates == YES)
            PrepareSeparateFiles(chain->chainNumber, currentIteration, 1,  filePaths, programOptions, files);

        numCA = numMIG = 0;
        numEventsTot=0;
        countTMRCA = 0.0;

        // order of proposals:
        fprintf (stderr, "\n>> Current log conditional Likelihood tree of the chain %d is = %lf  \n", chain->chainNumber,chain->currentlogConditionalLikelihoodTree );

        // 1- Update the topology using Wilson Balding move
        if (programOptions->doUseFixedTree == NO)
          WilsonBaldingMove(chain, seed,  chain->currentlogConditionalLikelihoodSequences, programOptions,  ObservedCellNames,  msa);
       
        //2- Update M_T
        newTotalEffectivePopulationSizeMove(chain, seed,  programOptions, ObservedCellNames,  msa,  opt, sampleSizes);
        
        //3-Update the vector of proportions \theta(this will change M_i)
        
        newProportionsVectorMove(chain, seed,  programOptions, ObservedCellNames, msa, opt, sampleSizes);
        
        //        4- Update the Delta_i
        for( i = 0 ; i < programOptions->numClones; i++)
        {
            
            popI=*(chain->populations + i);
            newScaledGrowthRateMoveforPopulation(chain, popI, seed,  programOptions,ObservedCellNames, msa, opt, sampleSizes);
            
        }
        //        5-Update T_i conditional on Delta_i
    
    
    
         //        6-Update the sub tree_i for sub population i by changing the topology of the sub tree of the time of some internal nodes.
        
 //   }
    
}
/********************* WilsonBaldingMove **********************/
/* WilsonBaldingMove*/
void WilsonBaldingMove(Chain *chain, long int *seed, double currentLogLikelihoodSequences, ProgramOptions *programOptions, char *ObservedCellNames[], pll_msa_t * msa){
    //select 2 random nodes to do the move
    
    if (chain->numClones > 1)
    {
        int i,j;
        double currentProb, proposalProb;
        pll_utree_rb_t * rb = (pll_utree_rb_t *)  malloc( (long)sizeof(pll_utree_rb_t));;
       int min =0;
        int max=chain->numClones -2;// except the oldest population
       int orderPopulationToDisconnect = rand() % (max + 1 - min) + min;// generate random integer between 0 and max -min = max-0=max= numClones -2
        int orderPopulationToReAttach;
        int indexPopulationToReAttach;
       int indexNodeToReAttach =0;
       TreeNode *tempNode;
       Population *popI=*(chain->populations + orderPopulationToDisconnect);
       Population *proposalPop;
        Population *currentFatherPop;
       TreeNode * MRCA=  popI->MRCA;
        
       currentFatherPop= popI->FatherPop;
        
       double timeOrigin = popI->timeOriginSTD;
       TreeNode *nodeReAttach;
       Population *tempPopulation;
       double timeOriginScaledOtherUnits;
        double distanceModelTimeProposalEdge;
       double acceptanceProbability=0;
       int numCandidateNodes=0;
       double sumNumerators=0;
       double sumDenominators=0;
    
       TreeNode *listCandidateNodes[chain->numNodes];
       //finding branches candidates
       for( i = orderPopulationToDisconnect+1 ; i < chain->numClones; i++)
       {
           tempPopulation = *(chain->populations + i);
           for (j = tempPopulation->sampleSize  ; j < tempPopulation->numGametes ; j++)
           {
              // tempNode= *(nodes) + tempPopulation->idsGametes[j];
               tempNode= chain->nodes + tempPopulation->idsGametes[j];
               timeOriginScaledOtherUnits = (timeOrigin * popI->effectPopSize ) /(tempPopulation->effectPopSize);
               if ( tempNode->index != MRCA->anc1->index && tempNode->index != MRCA->index   && timeOriginScaledOtherUnits > tempNode->time )
                   // if tempNode is older in time backwards and different from MRCA
               {
                   if (( tempNode->anc1 != NULL ) && (tempNode->anc1->orderCurrentClone == tempNode->orderCurrentClone)  && (tempNode->anc1->index != MRCA->anc1->index) && ( timeOriginScaledOtherUnits <= tempNode->anc1->time ))
                       //if tempNode has father node and the father node is older than  the rescaled time of origin
                   {
                       listCandidateNodes[numCandidateNodes]= tempNode;
                       numCandidateNodes ++;
                   }
               }
           }
       }
        if (numCandidateNodes ==0)
        {
              fprintf (stderr, "\n: No candidate edges  to Wilson Balding move.\n");
            return;
        }
        
        indexNodeToReAttach = rand() % (numCandidateNodes  - min) + min;
        
        nodeReAttach = listCandidateNodes[indexNodeToReAttach];
        indexPopulationToReAttach = nodeReAttach->indexCurrentClone;
        orderPopulationToReAttach = nodeReAttach->orderCurrentClone;
        
        if (orderPopulationToReAttach <= (chain->numClones- 1))
           proposalPop =  *(chain->populations + orderPopulationToReAttach);
        else{
            fprintf (stderr, "\n: Index of population outside range.\n");
            return;
        }
       // proposalPop = *(populations + orderPopulationToReAttach);
        
        pll_tree_edge_t  * edgeReAttach = nodeReAttach->edgeBack;
        
        pll_tree_rollback_t * rollback_info= (pll_tree_rollback_t *)  malloc( (long)sizeof(pll_tree_rollback_t));
        
        if (currentFatherPop !=NULL && proposalPop !=NULL )
        {
            if (currentFatherPop->order !=proposalPop->order)
            {
                proposalProb=ProbabilityCloneiFromClonej2(popI, proposalPop, chain->populations, chain->numClones);
                currentProb=ProbabilityCloneiFromClonej2(popI, currentFatherPop , chain->populations, chain->numClones);
                sumNumerators= sumNumerators + log(proposalProb);
                sumDenominators= sumDenominators + log(currentProb);
                
                popI->oldFatherPop =currentFatherPop;
                popI->FatherPop = proposalPop;
                
            }
        }
        distanceModelTimeProposalEdge=(nodeReAttach->anc1->time - (timeOrigin * popI->effectPopSize ) /(proposalPop->effectPopSize));
        double time_from_r_to_p = RandomUniform(seed)*(distanceModelTimeProposalEdge *proposalPop->effectPopSize)  ;
        double effect_pop_size =proposalPop->effectPopSize ;
        
        if (distanceModelTimeProposalEdge >0)
          sumNumerators= sumNumerators + log(distanceModelTimeProposalEdge);
        
        double distanceModelTimeCurrentEdge =MRCA->anc1->time- (timeOrigin * popI->effectPopSize ) /(currentFatherPop->effectPopSize);
        if (distanceModelTimeCurrentEdge >0)
            sumDenominators= sumDenominators +  log(distanceModelTimeCurrentEdge);
        
         char *oldNewickString=NULL;
        oldNewickString = toNewickString2 ( chain->root, chain->mutationRate, programOptions->doUseObservedCellNames);
       printf("\n before newick = %s  \n", oldNewickString);
        oldNewickString=NULL;
        //oldNewickString = toNewickString3 ( treeRootInit[0], treeRootInit[0]->nodeBack, mutationRate, ObservedCellNames);
       //  oldNewickString=NULL;
       // printf("\n  before newick = %s  \n", oldNewickString);
        TreeNode * container_of_u;
        TreeNode * container_of_v;
        
        if (! pllmod_utree_spr1(MRCA->nodeBack->back, nodeReAttach->nodeBack->back, rollback_info,rb,  &(chain->nodes), time_from_r_to_p, currentFatherPop->effectPopSize,
                                proposalPop->effectPopSize, MRCA->anc1, nodeReAttach->anc1,
                                MRCA,//mrca
                                nodeReAttach,//the other end point of the edge to reconnect
                                container_of_u,
                                container_of_v
                                ))//first argument is the back nodelet of the ancester of the MRCA and
        {
            fprintf (stderr, "\nERROR: Wilson Balding move cannot be applied.\n");
        
            popI->FatherPop = popI->oldFatherPop;
            
            free(oldNewickString);
            oldNewickString= NULL;
            return;
            //exit (1);
        }
         char *newNewickString=NULL;
     //   newNewickString=NULL;
      //   newNewickString = toNewickString3 ( treeRootInit[0], treeRootInit[0]->nodeBack, mutationRate, ObservedCellNames);
       // printf("\n after newick 1 = %s  \n", newNewickString);
      newNewickString = toNewickString2 ( chain->root, chain->mutationRate, programOptions->doUseObservedCellNames);
        printf("\n after newick 2 = %s  \n", newNewickString);
      
        //compute the new likelihood
        double newlogConditionalLikelihoodSequences;
        if (newNewickString !=NULL)
           newlogConditionalLikelihoodSequences = LogConditionalLikelihoodSequences(chain, msa,  newNewickString, programOptions);
        
        
        free(newNewickString);
        newNewickString=NULL;
//        logl = pllmod_utree_compute_lk(partition,
//                                       tree,
//                                       params_indices,
//                                       1,
//                                       1);
        
        /* compute marginal likelihoods */
//        printf("\nMarginal likelihoods:\n");
//        logl = pll_compute_root_loglikelihood (partition,
//                                               tree->clv_index,
//                                               tree->scaler_index,
//                                               params_indices,
//                                               NULL);
//        printf ("  Log-L Partial at %s: %f\n", tree->label, logl);
//        logl = pll_compute_root_loglikelihood (partition,
//                                               tree->back->clv_index,
//                                               tree->back->scaler_index,
//                                               params_indices,
//                                               NULL);
//        printf ("  Log-L Partial at %s: %f\n", tree->back->label, logl);
//
//        /* compute global likelihood */
        
        sumNumerators = sumNumerators + newlogConditionalLikelihoodSequences;
        sumDenominators = sumDenominators + currentLogLikelihoodSequences;
        
      
           acceptanceProbability =sumNumerators- sumDenominators;
        
        double randomNumber= RandomUniform (seed);
        
        double LogAcceptanceRate = (sumNumerators - sumDenominators) >0? (sumNumerators - sumDenominators) :0;
        
        if (log(randomNumber) < LogAcceptanceRate )
         {
            // accept the WB move
               printf("\n Accepted Wilson Balding move\n");
            //int result= pll_utree_tbr(MRCA->nodeBack,edgeReAttach,rollback_info);
             chain->currentlogConditionalLikelihoodSequences =newlogConditionalLikelihoodSequences;
             double  currentLikelihoodTree= chain->currentlogConditionalLikelihoodTree;
             //update the coalescentEventTimes
             double newLogConditionalLikelihoodTree= LogConditionalLikelihoodTree(chain->root, chain->nodes, chain->populations,  chain->numClones);
             chain->currentlogConditionalLikelihoodTree =newLogConditionalLikelihoodTree;
             
         }
        else {
            //rollback the move
            if (rollback_info !=NULL){
                 printf("\n Rejected Wilson Balding move\n");
                
                pllmod_utree_rollback_spr1(rollback_info,rb, &(chain->nodes),
                                           currentFatherPop->effectPopSize,
                                           proposalPop->effectPopSize,
                                           MRCA->anc1,
                                           nodeReAttach->anc1,
                                           MRCA,//mrca
                                           nodeReAttach,
                                            container_of_u,
                                           container_of_v
                                           );
                free(rollback_info);
                free(rb);
                rollback_info = NULL;
                rb=NULL;
            }
            
        }

    }
}
/********************* ScalingMove **********************/
/* ScalingMove*/
double ScalingMove(double var, double beta, long int seed ){// beta is greater than 1
    
    double delta = (1.0/beta) + RandomUniform(&seed)*(beta - (1.0/ beta)) ;
    return(var * delta );
}
/********************* proposalScalingMoveProportionsVector **********************/
/* proposalScalingMoveProportionsVector*/
void proposalScalingMoveProportionsVector(double *proportionsVector, int numClones, double omega, long int *seed )
{
    int i;
    if (omega >= 0)
    {
       double displacement =(-1 * omega) + RandomUniform(seed)*(2* omega) ;
       for(i=0; i< numClones; i++)
       {
          proportionsVector[i]=proportionsVector[i] +displacement;
       }
    }
}
/********************* proposeMutationRate **********************/
/* proposeMutationRate*/
void proposalMutationRate(double mutationRate, int numClones, double omega, long int *seed )
{
    int i;
    if (omega >= 0)
    {
        double displacement =(-1 * omega) + RandomUniform(seed)*(2* omega) ;
        mutationRate=mutationRate +displacement;
    }
}
/********************* RangeNodeTree **********************/
/* RangeNodeTree*/
// gives the time range  of a node in a tree between the ancestral node and the oldest child.
void RangeNodeTree(  TreeNode *node, double *from, double* to)
{
    if (node !=NULL)
   {
      if (node->anc1 == NULL)
       {
        *from =node->timePUnits ;
       }
      else
        *from = node->anc1->timePUnits;
        if (node->left !=NULL && node->right !=NULL)
        {
            *to = (node->left->timePUnits >= node->right->timePUnits)? node->left->timePUnits : node->right->timePUnits ;
        }
   }
}
/********************* proposalSlidingWindow **********************/
/* proposalSlidingWindow*/
void proposalSlidingWindow(  double  *newvalue, double oldvalue,  double windowSize, long int *seed)
{
   *newvalue = oldvalue + RandomUniform(seed)* windowSize ;
    if (*newvalue <0)
        *newvalue = -  *newvalue;
}
/********************* proposalChangeCoalTimeInternalNodePopulation **********************/
/* proposalChangeCoalTimeInternalNodePopulation*/
void proposalChangeCoalTimeInternalNodePopulation(TreeNode **root, TreeNode **nodes, Population *pop, int numNodes, long int seed, double mutationRate)
{
    if (pop !=NULL)
    {
        int i;
        TreeNode *p;
        TreeNode *candidateNode;
        double from,  to;
        double newTime;
      
        TreeNode *MRCA =pop->MRCA;
        int populationSample= pop->sampleSize ;
        int numberInternalNodes = populationSample -1 + pop->numIncomingMigrations -1; // the last -1 is because the index start with 0.
    
        double rand = RandomUniform(&seed);
        
        int indexCandidateNode = round(rand * numberInternalNodes);
    
        int numberVisitedNode= 0;
        if (numberInternalNodes > 1  && MRCA != NULL)
        {
            candidateNode=getInternalNodebyIndexSubTree(MRCA, indexCandidateNode , pop);
            RangeNodeTree( candidateNode, &from, &to);
            proposalSlidingWindow(&newTime, candidateNode->timePUnits, from - to, seed);
            if (newTime >= to && newTime >= from){
                
                candidateNode->timePUnits = newTime;
                candidateNode->time =candidateNode->timePUnits / pop->effectPopSize;
                candidateNode->length = (candidateNode->anc1->timePUnits- candidateNode->timePUnits);
                //* mutationRate;
                candidateNode->lengthModelUnits = (candidateNode->anc1->time- candidateNode->time);
               // * mutationRate;
            }
        }
    }
}
///********************* getInternalNodebyIndexSubTree **********************/
///* getInternalNodebyIndexSubTree*/
//TreeNode* getInternalNodebyIndexSubTree(TreeNode *currentNode,int indexNode, int *numberAlreadyVisited)
//{
//
//
////    if (indexNode == *numberAlreadyVisited && currentNode != NULL && isLeaf(currentNode) ==  NO)
////    {
////        return currentNode;
////    }
////    else
////    {
////        if (currentNode != NULL )
////        {
////            if (currentNode->left != NULL && isLeaf(currentNode->left) == NO)
////           //if (isLeaf(currentNode->left) == NO)
////            {
////                *numberAlreadyVisited = *numberAlreadyVisited + 1;
////                return getInternalNodebyIndexSubTree(currentNode->left, indexNode, numberAlreadyVisited);
////            }
////            else if (currentNode->right != NULL && isLeaf(currentNode->right) == NO)
////           //else if (isLeaf(currentNode->right) == NO)
////            {
////                *numberAlreadyVisited = *numberAlreadyVisited + 1;
////                return getInternalNodebyIndexSubTree(currentNode->right, indexNode, numberAlreadyVisited);
////            }
////
////        }
////        else
////        {
////            return NULL;
////
////        }
////    }
//    return NULL;
//}
/********************* getInternalNodebyIndexSubTree **********************/
/* getInternalNodebyIndexSubTree*/
TreeNode* getInternalNodebyIndexSubTree(TreeNode *currentNode,int indexNode, Population * pop)
{
    int numberInternalNodes = pop->sampleSize -1 + pop->numIncomingMigrations ;
    //int numberNodes = pop->sampleSize  +  numberInternalNodes ;
    TreeNode *current;
    TreeNode *queue_internalNodes[numberInternalNodes];
    TreeNode *queue_allNodes[numberInternalNodes];
    int numberAlreadyVisited=0;

    queue_allNodes[0]=currentNode;
    int nextNodeinQueue=0;
    int numberNodesinQueue =1;
    int isInternal =NO;
    
    while(numberNodesinQueue >0)
    {
        current = queue_allNodes[nextNodeinQueue];
        nextNodeinQueue++;
        isInternal =NO;
        
        if (current->left !=NULL)
        {
            isInternal = YES;
            queue_allNodes[numberNodesinQueue]=current->left;
            numberNodesinQueue ++;
        }
        if (current->right!=NULL)
        {
            isInternal = YES;
            queue_allNodes[numberNodesinQueue]=current->right;
            numberNodesinQueue ++;
        }
        if (isInternal == YES)
        {
            queue_internalNodes[numberAlreadyVisited]=current;
            if (indexNode == numberAlreadyVisited)
                return queue_internalNodes[indexNode];
            numberAlreadyVisited++;
        }
    }
    if (indexNode < numberInternalNodes)
        return queue_internalNodes[indexNode];
    else return NULL;
}
/********************* isLeaf **********************/
/* isLeaf*/
bool isLeaf(TreeNode *node){
    if (node != NULL)
    {
      if (node->left == NULL && node->right == NULL && node->anc1 != NULL )
          return YES;
     else
         return NO;
    }
    return NO;
}
/********************* SetPopulationsBirthRate **********************/
/* SetPopulationsBirthRate*/
void SetPopulationsBirthRate(Population** populations, double lambda, int numClones){
    int i;
    Population *popI;
    for( i = 0 ; i < numClones; i++)
    {
        popI=*(populations + i);
        if (popI != NULL)
           popI->birthRate = lambda;
    }
}
/********************* proposalProportionsVector **********************/
/* proposalProportionsVector*/
void proposalProportionsVector(double *newProportionsvector, double *currentProportionVector, double tuningParameter, int numClones,  long int *seed)
{
    int i=0;
    double vectorForGamma[numClones];
    for( i = 0 ; i < numClones; i++)
    {
        vectorForGamma[i]= tuningParameter * currentProportionVector[i];
    }
     RandomDirichlet2 (vectorForGamma, numClones, newProportionsvector, seed);
}
/********************* AcceptanceRateTheta **********************/
/* AcceptanceRateTheta*/
double AcceptanceRateTheta(double *newProportionsvector, TreeNode  *root, TreeNode *nodes, Population **populations, int numClones, int totalPopSize, double totalDelta){
    Population *popI;
    int i;
    double currentProportionsVector[numClones];
    double vectorOnes[numClones];
    double logConditionalLikelihoodTreeCurrentValues= LogConditionalLikelihoodTree(root, nodes, populations,  numClones);
    for( i = 0 ; i < numClones; i++)
    {
        popI=*(populations + i);
        vectorOnes[i]=1;
        if (popI != NULL){
            currentProportionsVector[i]= popI->effectPopSize;
            popI->oldeffectPopSize = popI->effectPopSize;
            popI->effectPopSize = newProportionsvector[i] * totalPopSize;
            popI->delta = newProportionsvector[i]* totalDelta;
            popI->growthRate =popI->delta  /popI->effectPopSize;
            popI->popSize=popI->effectPopSize * popI->birthRate;
            popI->deathRate= popI->birthRate -popI->growthRate;
        }
        
    }
    double logConditionalLikelihoodTreeNewValues= LogConditionalLikelihoodTree(root, nodes, populations,  numClones);
    double logNumerator=logConditionalLikelihoodTreeCurrentValues + DirichletDensity( newProportionsvector, vectorOnes , numClones)+DirichletDensity(newProportionsvector, currentProportionsVector, numClones);
    double logDenominator=logConditionalLikelihoodTreeNewValues + DirichletDensity( currentProportionsVector, vectorOnes , numClones)+DirichletDensity(currentProportionsVector, newProportionsvector, numClones);
    return(logNumerator - logDenominator);
}
/********************* DirichletDensity **********************/
/* DirichletDensity*/
double DirichletDensity(double *proportionsVector, double *concentrationVector, int sizeVector)
{
    int i;
    double sum=0;
    double logResult=0;
    Population *popI;
    for( i = 0 ; i < sizeVector; i++)
    {
        sum = sum +concentrationVector[i];
        logResult= logResult+(concentrationVector[i]-1)*log(proportionsVector[i]);
        logResult= logResult-lgamma(concentrationVector[i]);
    }
    logResult = logResult+lgamma(concentrationVector[i]);
    return logResult;
}
/********************* connectNodes **********************/
/* connectNodes*/
void connectNodes(TreeNode *left, TreeNode *right, TreeNode *ancester  ){
    if (left!=NULL && right!= NULL && ancester!=NULL )
    {
        connectNodelets(left);
        connectNodelets(right);
        connectNodelets(ancester);
        //connect the child nodes
        left->nodeBack->back =ancester->nodeLeft;
        right->nodeBack->back =ancester->nodeRight;
        
        //connect the ancester node
        ancester->nodeLeft->back =left->nodeBack;
        ancester->nodeRight->back =right->nodeBack;
        
        //connect the edges
        left->edgeBack->edge.utree.parent =ancester->nodeLeft;
        right->edgeBack->edge.utree.parent=ancester->nodeRight;
       
        ancester->edgeLeft->edge.utree.child =left->nodeBack;
        ancester->edgeRight->edge.utree.child=right->nodeBack;

        ancester->isLeaf=NO;
    }
    else if(left==NULL && right== NULL && ancester!=NULL )
    { // the ancester node is a leaf
        connectNodelets(ancester);
        //connect the child nodes
        ancester->nodeLeft =NULL;
        ancester->nodeRight =NULL;
        ancester->isLeaf=YES;
        
        ancester->edgeLeft=NULL;
        ancester->edgeRight=NULL;
    }
    else if(left!=NULL && right== NULL && ancester!=NULL )
    {
        connectNodelets(left);
        connectNodelets(ancester);
         //connect the child nodes
        left->nodeBack->back =ancester->nodeLeft;
        //connect the ancester node
        ancester->nodeLeft->back =left->nodeBack;

        //connect the edges
        left->edgeBack->edge.utree.parent =ancester->nodeLeft;
        
        ancester->edgeLeft->edge.utree.child =left->nodeBack;

        ancester->isLeaf=NO;
    }
    else if(left==NULL && right!= NULL && ancester!=NULL )
    {
        connectNodelets(right);
        connectNodelets(ancester);

        //connect the child nodes
        right->nodeBack->back =ancester->nodeRight;
        //connect the ancester node
        ancester->nodeRight->back =right->nodeBack;
        //connect the edges
        right->edgeBack->edge.utree.parent=ancester->nodeRight;

        ancester->edgeRight->edge.utree.child=right->nodeBack;
        
        ancester->isLeaf=NO;
    }
}
/********************* connectNodelets **********************/
/* connectNodelets*/
void connectNodelets(TreeNode *node )
{
    if (node != NULL)
    {
        if (node->left == NULL && node->right== NULL)
        {
            char * temp;
            node->isLeaf=YES;
            node->nodeLeft= NULL;
            node->nodeRight= NULL;
            node->nodeBack->next = NULL;
            
            node->nodeBack->node_index= node->index;
            if (asprintf(&temp,  "%d_back",  node->label)<0)
                return;
            node->nodeBack->label=temp;
        }
        else
        {
            char * temp1;
            char * temp2;
            char *temp3;
            node->isLeaf=NO;
            node->nodeBack->next=node->nodeLeft;
            node->nodeLeft->next=node->nodeRight;
            node->nodeRight->next =node->nodeBack;
            
            node->nodeLeft->node_index= node->index;
            
            node->nodeRight->node_index= node->index;
            
            node->nodeBack->node_index= node->index;
            if (asprintf(&temp1,  "%d_back",  node->label)<0)
                return;
            node->nodeBack->label=temp1;
            if (asprintf(&temp2,  "%d_left",  node->label)<0)
                return;
            node->nodeLeft->label= temp2;
            if (asprintf(&temp3,  "%d_right",  node->label)<0)
                return;
            node->nodeRight->label= temp3;
      
        }
    }
}
/********************* setLength **********************/
/* setLength*/
void setLength(TreeNode *node )
{
    double lengthEdge;
    if(node->anc1!=NULL)
    {
        //node->nodeBack->length= node->lengthModelUnits;//this takes into account  the mutation rate
        // lengthEdge = node->anc1->timePUnits - node->timePUnits; // this doesnt take into account the mutation rate
        lengthEdge = node->anc1->timePUnits - node->timePUnits;//not in model time, already includes the effect pop size
        node->nodeBack->length= lengthEdge;
        
        if (node->isLeaf==NO)
        {
          //node->nodeLeft->length=  node->timePUnits - node->left->timePUnits;
           //node->nodeRight->length= node->timePUnits - node->right->timePUnits;
            
            node->nodeLeft->length=  node->timePUnits - node->left->timePUnits;
            node->nodeRight->length= node->timePUnits - node->right->timePUnits;
        
        //node->edgeLeft->length  =node->timePUnits - node->left->timePUnits;
        //node->edgeRight->length  =node->timePUnits - node->right->timePUnits;
            node->edgeLeft->length  =node->timePUnits - node->left->timePUnits;
            node->edgeRight->length  =node->timePUnits - node->right->timePUnits;
        }
        
        if(node->edgeBack!=NULL)
        {
            node->edgeBack->length=lengthEdge;
        
        }
       
    }
}
/***************************** pllmod_utree_spr1*******************************/
/* pllmod_utree_spr1*/
PLL_EXPORT int pllmod_utree_spr1(pll_unode_t * p_edge,
                                pll_unode_t * r_edge,
                                pll_tree_rollback_t * rollback_info,
                                 pll_utree_rb_t * rb,
                                 TreeNode **nodes,
                                 double time_from_r_to_p,
                                 double old_father_pop_effect_pop_size,
                                 double new_father_pop_effect_pop_size,
                                 TreeNode *container_of_p,
                                 TreeNode *container_of_r,
                                 TreeNode *container_of_p_prime,//mrca
                                 TreeNode *container_of_r_prime,//the other end point of the edge to reconnect
                                 TreeNode * container_of_u,
                                 TreeNode * container_of_v
                                 )
{
    int retval;
    
    if (pllmod_utree_is_tip(p_edge))
    {
        /* invalid move */
        pllmod_set_error(PLLMOD_TREE_ERROR_SPR_INVALID_NODE,
                         "Attempting to prune a leaf branch");
        return PLL_FAILURE;
    }
    
    /* save rollback information */
    if (rollback_info)
    {
        rollback_info->rearrange_type     = PLLMOD_TREE_REARRANGE_SPR;
        rollback_info->rooted             = 0;
        rollback_info->SPR.prune_edge     = (void *) p_edge;
        rollback_info->SPR.regraft_edge   = (void *) p_edge->next->back;
        rollback_info->SPR.prune_bl       = p_edge->length;
        rollback_info->SPR.prune_left_bl  = p_edge->next->length;
        rollback_info->SPR.prune_right_bl = p_edge->next->next->length;
        rollback_info->SPR.regraft_bl     = r_edge->length;
    }
    
   // pll_utree_rb_t * rb;
    
    retval = pll_utree_spr1(p_edge,
                           r_edge,
                           rb, 0, 0,
                            nodes,
                            time_from_r_to_p,
                            old_father_pop_effect_pop_size,
                            new_father_pop_effect_pop_size,
                            container_of_p,
                            container_of_r,
                            container_of_p_prime,
                            container_of_r_prime,
                            container_of_u,
                             container_of_v
                            );
    
    return retval;
}
/***************************** pll_utree_spr1*******************************/
/* pll_utree_spr1*/
PLL_EXPORT int pll_utree_spr1(pll_unode_t * p,//ancester MRCA
                             pll_unode_t * r,//ancester of the edge to reconnect
                             pll_utree_rb_t * rb,
                             double * branch_lengths,
                             unsigned int * matrix_indices,
                               TreeNode **nodes,
                               double time_from_r_to_p,
                              double old_father_pop_effect_pop_size,
                              double new_father_pop_effect_pop_size,
                              TreeNode *container_of_p,
                              TreeNode *container_of_r,
                              TreeNode *container_of_p_prime,//mrca
                              TreeNode *container_of_r_prime,//the other end point of the edge to reconnect,
                               TreeNode * container_of_u,
                               TreeNode * container_of_v
                              )
{
    /* given nodes p and r, perform an SPR move in the following way,
     i.e. prune subtree C and make it adjacent to subtree D:
     
     A           B          C             D           A          B
     ____        ____       ____          ____        ____       ____
     \  /        \  /       \  /          \  /        \  /       \  /
     \/          \/         \/            \/          \/         \/
     *          *          * p'           *          *          *
     \         |     q   /                \         |         /
     *'*_____.*._____*'* p     --->       *'*_____.*._____*'*
     '*'     *.*     '*'                  '*'     *.*     '*'
     / r       u    q' \                  /                 \
     r' *                   * v              *                   *
     /\                   /\              /\                   /\
     /__\                 /__\            /__\                 /__\
     
     D                    E               C                    E
     
     node p must be part of an inner node (i.e. node with ->next set). The
     procedure prunes the subtree rooted at the opposite end-point of p
     (subtree C in our case) and regrafts it on the edge r'<->r. It is done
     in the following way:
     
     (a) prune the subtree rooted at the opposite end-point of p (p' on figure)
     by breaking the edges q<->u and q'<->v
     
     (b) connect node u with node v
     
     (c) break edge r<->r' by connecting node r with node q, and node r' with
     node q'
     
     Node r must not be part of the subtree to be pruned (C in this case). Note
     that for speed reasons, the function *does not* check this property to save
     a tree traversal. A safer (albeit slower) function that checks this
     property is pll_utree_spr_safe
     */
    
    int k = 0;
    
    if ((!branch_lengths && matrix_indices) ||
        (branch_lengths && !matrix_indices))
    {
        pll_errno = PLL_ERROR_PARAM_INVALID;
        snprintf(pll_errmsg, 200, "Parameters 4,5 must be both NULL or both set");
        return PLL_FAILURE;
    }
    
    /* if p is a tip node then prompt an error */
    if (!p->next)
    {
        pll_errno = PLL_ERROR_SPR_TERMINALBRANCH;
        snprintf(pll_errmsg, 200, "Prune edge must be defined by an inner node");
        return PLL_FAILURE;
    }
    
    /* check whether the move will result in the same tree */
    if (r == p || r == p->back ||
        r == p->next || r == p->next->back ||
        r == p->next->next || r == p->next->next->back)
    {
        pll_errno = PLL_ERROR_SPR_NOCHANGE;
        snprintf(pll_errmsg, 200, "Proposed move yields the same tree");
        return PLL_FAILURE;
    }
    
    /* check if rollback buffer is provided, and fill it up */
    if (rb)
    {
        rb->move_type = PLL_UTREE_MOVE_SPR;
        rb->spr.p = p;
        rb->spr.r = r;
        rb->spr.rb = r->back;
        rb->spr.r_len = r->length;
        rb->spr.pnb = p->next->back;
        rb->spr.pnb_len = p->next->length;
        rb->spr.pnnb = p->next->next->back;
        rb->spr.pnnb_len = p->next->next->length;
    }
    
    /* (b) connect u and v */
    pll_unode_t * u = p->next->back;
    pll_unode_t * v = p->next->next->back;
    utree_link(u,
               v,
               u->length + v->length,
               u->pmatrix_index);
    
    /*connect the TreeNodes */
    //TreeNode * container_of_p =  *nodes + p->node_index;
   container_of_u = container_of_p->anc1;
   // TreeNode * container_of_v =  *nodes + v->node_index;
   // TreeNode * container_of_v ;
    if (container_of_p_prime->anc1->left == container_of_p_prime)//container_of_p_prime(MRCA) is left child
        container_of_v=container_of_p->right;
    else //container_of_p_prime is right child
        container_of_v=container_of_p->left;
    
    // TreeNode * container_of_p_prime =  *nodes + p->back->node_index;
    //disconnect container_of_p
  
    if (container_of_v->anc1->left == container_of_v)//container_of_v is left child
            container_of_p->left=NULL;
    else //container_of_v is right child
            container_of_p->right=NULL;
 
    container_of_v->anc1=NULL;
    
    //connect  container_of_u and container_of_v
    if (container_of_u->left == container_of_p)
    {//container_of_p is left child
        container_of_u->left=NULL;
       container_of_u->left=container_of_v;
    }
    else{ //container_of_p is right child
        container_of_u->right=NULL;
        container_of_u->right=container_of_v;
    }
    
    container_of_v->anc1 =container_of_u;
    container_of_v->length =u->length + v->length;
    container_of_v->lengthModelUnits =(u->length + v->length) / old_father_pop_effect_pop_size;
    
    
    /* if requested, store the new branch length for the corresponding
     pmatrix index */
    if (branch_lengths)
    {
        branch_lengths[k] = u->length;
        matrix_indices[k] = u->pmatrix_index;
    }
    
    /* (a) prune subtree C */
    p->next->back = p->next->next->back = NULL;
    //now for the TreeNodes
    container_of_p->anc1 = NULL;
    
    
    
    /* (c) regraft C at r<->r' */
    double length_from_r_prime_to_p;
    //length = r->length / 2;
    length_from_r_prime_to_p = r->length  -   time_from_r_to_p  ;
    
    /* r' <-> q' */
    utree_link(r->back,//r_prime
               p->next->next,//q_prime
               length_from_r_prime_to_p,
               p->next->next->pmatrix_index);
    //regraft C using TreeNodes
   // TreeNode * container_of_r =  *nodes + r->node_index;
    //TreeNode * container_of_r_prime =  *nodes + r->back->node_index;
    // TreeNode * MRCA =  *nodes + p->back->node_index;
    
    
    if (container_of_r_prime->anc1->left == container_of_r_prime)
    {//container_of_r_prime is left child of container_of_r
        container_of_r->left=NULL;
        container_of_r->left=container_of_p;
    }
    else
    {//container_of_r_prime is right child of container_of_r
        container_of_r->right=NULL;
        container_of_r->right=container_of_p;
    }
    
    container_of_r_prime->anc1=NULL;
    container_of_r_prime->anc1 =container_of_p;
    container_of_p->anc1=NULL;
    container_of_p->anc1 =container_of_r;
    
    container_of_r_prime->length =length_from_r_prime_to_p;
    container_of_r_prime->lengthModelUnits = container_of_r_prime->length / new_father_pop_effect_pop_size;
    
    
    container_of_p->lengthModelUnits =time_from_r_to_p / new_father_pop_effect_pop_size;
    
     container_of_p->length =time_from_r_to_p ;
    
    container_of_p->time = container_of_r->time - (time_from_r_to_p / new_father_pop_effect_pop_size) ;
    
    container_of_p->timePUnits =  container_of_p->time * new_father_pop_effect_pop_size;
    
    //container_of_p->right->time = (r->length -   time_from_r_to_p)/effect_pop_size;
    if  (container_of_p->right ==container_of_p_prime)//the MRCA is on the right
    {
        container_of_p->left =container_of_r_prime;

    }
    else if (container_of_p->left ==container_of_p_prime)//the MRCA is on the left
    {
        container_of_p->right =container_of_r_prime;
    
    }
    
    container_of_p_prime->time = (container_of_p_prime->timePUnits)/new_father_pop_effect_pop_size;
    container_of_p_prime->timePUnits =    container_of_p->right->time * new_father_pop_effect_pop_size;
    
    
    container_of_p_prime->length = container_of_p->timePUnits - container_of_p_prime->timePUnits;
    container_of_p_prime->lengthModelUnits = container_of_p->time - container_of_p_prime->time;
  
    
    /* if requested, store the new branch length for the corresponding
     pmatrix index */
    if (branch_lengths)
    {
        ++k;
        branch_lengths[k] = length_from_r_prime_to_p;
        matrix_indices[k]   = p->next->next->pmatrix_index;
    }
    
    /* r<->q */
    utree_link(r,
               p->next,
               time_from_r_to_p,
               //length,
               r->pmatrix_index);
    /* if requested, store the new branch length for the corresponding
     pmatrix index */
    if (branch_lengths)
    {
        ++k;
        branch_lengths[k] = length_from_r_prime_to_p;
        matrix_indices[k]   = r->pmatrix_index;
    }
    
    return PLL_SUCCESS;
}
/***************************** utree_link*******************************/
/* utree_link*/
static void utree_link(pll_unode_t * a,
                       pll_unode_t * b,
                       double length,
                       unsigned int pmatrix_index)
{
    a->back = b;
    b->back = a;
    a->length = length;
    b->length = length;
    
    a->pmatrix_index = b->pmatrix_index = pmatrix_index;
}
/***************************** Initialize*******************************/
/* Initialize*/
static void Initialize( double (*Eij)[4], double (*Mij)[4], double *freq,  ProgramOptions *programOptions ) {
    programOptions->numDataSets = 10;            /* the number of samples to simulate */
    programOptions->numCells = 8;                /* number of cells in each data set */
    programOptions->ploidy = 2;                 /* we assume diploid genomes */
    programOptions->numSites = 10000;                /* number of sites (markers, loci) to simulate = length of the chromosomes */
    //  N = 1000;                    /* effective population size */
    programOptions->numPeriods = 0;                /* number of distinct demographic periods */
    programOptions->doDemographics = NO;        /* whether to implement demographics */
    programOptions->doExponential = YES;            /* whether to do exponential growth */
    programOptions->growthRate = 0;                /* rate for the exponential population growth */
    programOptions->mutationRate = 1.0e-7;        /* nucleotide mutation rate per site per generation */
    programOptions->rateVarAmongLineages = NO;    /* modify rate variation among branches  (to non-clock) */
    programOptions->alphabet = BINARY;          /* alphabet 0/1 or DNA") */
    programOptions->altModel = 0;                /* by default the alternative model will be ISM haploid */
    programOptions->propAltModelSites = 0;        /* proportion of sites that will mutate according to alternative model */
    programOptions->nonISMRelMutRate = 1.0;        /* relative rate alternative/default model for sites */
    programOptions->equalBaseFreq = YES;        /* DNA base frequencies */
    freq[0] = freq[1] = freq[2] = freq[3] = 0.25;
    programOptions->titv = 0.5;                    /* transition/transversion rate ratio */
    programOptions->thereIsMij = NO;            /* mutation rate matrix*/
    Mij[0][0] = Mij[1][1] = Mij[2][2] = Mij[3][3] = 0;  /* mutation probabilities */
    Mij[0][1] = Mij[0][2] = Mij[0][3] = 1.0/3;
    Mij[1][0] = Mij[1][2] = Mij[1][3] = 1.0/3;
    Mij[2][0] = Mij[2][1] = Mij[2][3] = 1.0/3;
    Mij[3][0] = Mij[3][1] = Mij[3][2] = 1.0/3;
    programOptions->thereIsEij = NO;            /* error rate matrix*/
    Eij[0][0] = Eij[1][1] = Eij[2][2] = Eij[3][3] = 0;  /* sequencing error probabilities */
    Eij[0][1] = Eij[0][2] = Eij[0][3] = 1.0/3;
    Eij[1][0] = Eij[1][2] = Eij[1][3] = 1.0/3;
    Eij[2][0] = Eij[2][1] = Eij[2][3] = 1.0/3;
    Eij[3][0] = Eij[3][1] = Eij[3][2] = 1.0/3;
    /*                                        AA CC GG TT AC AG AT CG CT GT            */
    static const char mut_dist[10][10] = {  {  0, 2, 2, 2, 1, 1, 1, 2, 2, 2 },   /* AA */
        {  2, 0, 2, 2, 1, 2, 2, 1, 1, 2 },   /* CC */
        {  2, 2, 0, 2, 2, 1, 2, 1, 2, 1 },   /* GG */
        {  2, 2, 2, 0, 2, 2, 1, 2, 1, 1 },   /* TT */
        {  1, 1, 2, 2, 0, 1, 1, 1, 1, 2 },   /* AC */
        {  1, 2, 1, 2, 1, 0, 1, 1, 2, 1 },   /* AG */
        {  1, 2, 2, 1, 1, 1, 0, 2, 1, 1 },   /* AT */
        {  2, 1, 1, 2, 1, 1, 2, 0, 1, 1 },   /* CG */
        {  2, 1, 2, 1, 1, 2, 1, 1, 0, 1 },   /* CT */
        {  2, 2, 1, 1, 2, 1, 1, 1, 1, 0 }    /* GT */
    };
    
    programOptions->doJC = YES;
    programOptions->doHKY = NO;
    programOptions->doGTR = NO;
    programOptions->doGTnR = NO;
    programOptions->rateVarAmongSites = NO;         /* rate variation among different sites along the genome */
    programOptions->alphaSites = infinity;          /* alpha shape of the gamma distribution for rate variation among sites */
    programOptions->alphaBranches = infinity;       /* alpha shape of the gamma distribution for rate variation among lineages */
    programOptions->alphaCoverage = infinity;       /* alpha shape of the gamma distribution for coverage */
    programOptions->doSimulateFixedNumMutations = NO;    /* whether to simulate a fixed number of mutations */
    programOptions->doUserTree = NO;                /* whether to assume a user tree instead od making the coalescent */
    programOptions->doUserGenome = NO;                /* whether to use a user genome instead of a simulated one */
    programOptions->doPrintSNVgenotypes = NO;        /* whether to print SNVs */
    programOptions->doPrintSNVhaplotypes = NO;      /* whether to print haplotypes */
    programOptions->doPrintTrueHaplotypes = NO;      /* whether to print haplotypes without errors */
    programOptions->doPrintMLhaplotypes = NO;          /* whether to print ML haplotypes */
    programOptions->doPrintFullHaplotypes = NO;        /* whether to print sequences */
    programOptions->doPrintFullGenotypes = NO;        /* whether to print all genotypes (variable + invariable) */
    programOptions->doPrintTree = NO;                /* whether to print the coalescent tree */
    programOptions->doPrintTimes = NO;                /* whether to print coalescent times */
    programOptions->doPrintAncestors = NO;          /* whether to print data for ancestral cells */
    programOptions->doSimulateReadCounts = NO;      /* do not produce reads by default */
    programOptions->doPrintCATG = NO;                /* whether to print read counts for SNVs in CATG format*/
    programOptions->doSimulateData = NO;            /* whether to simulate any data or do inference from real data */
    programOptions->doPrintSeparateReplicates = YES; /* whether to put every replica in its own file */
    
    programOptions->doPrintIUPAChaplotypes = YES;    /* whether to print IUPAC halotypes */
    programOptions->doGeneticSignatures = NO;        /* whether to use a genetic signature to model trinucleotide mutations */
    programOptions->numUserSignatures = 0;            /* by default we do not use a genetic signature */
    programOptions->healthyTipBranchLength = 0;     /* length of the branch leading to the healthy cell */
    programOptions->transformingBranchLength = 0;     /* length of the transforming branch leading to the healthy ancestral cell */
    programOptions->coverage = 0;                    /* NGS  depth for read counts */
    programOptions->rateVarCoverage = NO;            /* there is coverage dispersion */
    programOptions->ADOrate = 0;                    /* allelic dropout */
    programOptions->sequencingError = 0;            /* NGS error rate */
    programOptions->genotypingError = 0;            /* add errors directly in the genotypes */
    programOptions->SNPrate = 0.0;                    /* germline variation rate for rooth healthy genome */
    programOptions->meanAmplificationError = 0;     /* mean of beta distribution for WGA errors */
    programOptions->varAmplificationError = 0;      /* variance of beta distribution for WGA errors */
    programOptions->simulateOnlyTwoTemplates = NO;    /* whether simualate maximum of two templates after single-cell amplification, or there can be all four */
    programOptions->haploidCoverageReduction = 0.5; /* proportion of reads produced when a single allele is present */
    programOptions->allelicImbalance = 0.5;            /* proportion of maternal/ paternal reads */
    programOptions->doubletRate = 0.0;                /* no doublets by default */
    programOptions->numNodes = 3000;                /* initial number of nodes allocated to build the coalescent trees */
    programOptions->seed = time(NULL);                 /* seed for random numbers */
    programOptions->userSeed = 0;                    /* seed entered by the user */
    programOptions->noisy = 1;                        /* level of information to be printed in the screen (see below) */
    programOptions->doNGS = NO;
}
/***************************** InitChainPopulations*******************************/
/* InitChainPopulations*/
void InitChainPopulations(Population **populations, int numClones, int noisy,  int TotalNumSequences  ) {
    int z;
    //struct Population* pops = malloc(numClones * (sizeof(struct Population)+TotalNumSequences * sizeof( int) + numClones * sizeof(double) ));
    struct Population* pops = malloc(numClones * (sizeof(struct Population) ));
    for (z = 0; z <= (numClones - 1); z++)
    {
        pops[z].index =z;
        pops[z].order=0;
        pops[z].birthRate =0.0;
        pops[z].deathRate = 0;
        pops[z].growthRate = 0;
        pops[z].sampleSize = 0;
        pops[z].popSize =0;
        pops[z].effectPopSize =0;
        pops[z].delta =0 ;
        pops[z].numActiveGametes = 0;
        pops[z].isAlive = 1;
        pops[z].nodeIdAncesterMRCA = 0;
        pops[z].numCompletedCoalescences = 0;
        pops[z].nextAvailableIdInmigrant = 0;
        pops[z].numIncomingMigrations = 0;
        pops[z].numPossibleMigrations = 0;
        pops[z].doEstimateTimeOrigin = NO;
        pops[z].timeOriginInput = 0;
        pops[z].timeOriginSTD =0;
        pops[z].timeMigrationSTDCurrentPop = 0;
        pops[z].FatherPop =NULL;
     
        *(populations + z) = &pops[z];

    }
}
/***************************** InitializeChains*******************************/
/* InitializeChains*/
void InitializeChains(Chain **chains,  const ProgramOptions *programOptions, const MCMCoptions *mcmcOptions,int * sampleSizes, long int *seed, char* ObservedCellNames[], pll_msa_t * msa)
{
    int chainNumber;
    Chain *currentChain=NULL;
    double totalTreeLength;
     int        numCA, numMIG;
    double      cumNumCA, meanNumCA, cumNumMIG, meanNumMIG, numEventsTot;
     char *newickString2;
    for(chainNumber=0; chainNumber< mcmcOptions->numChains;chainNumber++)
    {
       
        chains[chainNumber] = (Chain *)malloc(sizeof(Chain));
        
        chains[chainNumber]->chainNumber = chainNumber;
        chains[chainNumber]->currentNumberIerations=0;
        chains[chainNumber]->numClones = programOptions->numClones;
        chains[chainNumber]->mutationRate =programOptions->mutationRate;
        chains[chainNumber]->gammaParam=1;
        chains[chainNumber]->seqErrorRate= programOptions->seqErrorRate;
        chains[chainNumber]->dropoutRate =programOptions->dropoutRate;
        if (programOptions->numberClonesKnown)
        {
            chains[chainNumber]->numNodes = programOptions->numNodes;
            
            
        }
        
        chains[chainNumber]->proportionsVector = (double *) calloc((programOptions->numClones), (long) sizeof(double));
        if (!(chains[chainNumber]->proportionsVector))
        {
            fprintf (stderr, "Could not allocate proportions vector (%lu bytes)\n", (programOptions->numClones ) * (long) sizeof(double));
            exit (-1);
        }
        chains[chainNumber]->oldproportionsVector = (double *) calloc((programOptions->numClones), (long) sizeof(double));
        if (!(chains[chainNumber]->oldproportionsVector))
        {
            fprintf (stderr, "Could not allocate proportions vector (%lu bytes)\n", (programOptions->numClones ) * (long) sizeof(double));
            exit (-1);
        }
        chains[chainNumber]->populations = (Population**)malloc (sizeof(struct Population*)  * programOptions->numClones);
        if (!(chains[chainNumber]->populations))
        {
            fprintf (stderr, "Could not allocate populations (%lu bytes)\n", (programOptions->numClones)  * (long) sizeof(Population*));
            exit (1);
        }
     
//        currentChain->treeRootInit = (TreeNode **) calloc(1, sizeof(TreeNode *)); /* nodes pointers */
//        if (!(currentChain->treeRootInit))
//        {
//            fprintf (stderr, "Could not allocate treeRootInit (%lu bytes)\n", 1  * (long) sizeof(TreeNode));
//            exit (1);
//        }
//        currentChain->oldtreeRootInit = (TreeNode **) calloc(1, sizeof(TreeNode *)); /* nodes pointers */
//        if (!(currentChain->oldtreeRootInit))
//        {
//            fprintf (stderr, "Could not allocate treeRootInit (%lu bytes)\n", 1  * (long) sizeof(TreeNode));
//            exit (1);
//        }
        //TreeNode  **treeTips;
        
        chains[chainNumber]->treeTips =  malloc (programOptions->TotalNumSequences* sizeof(TreeNode*));
        if (!(chains[chainNumber]->treeTips))
        {
            fprintf (stderr, "Could not allocate the treeTips array\n");
            exit (-1);
        }
        chains[chainNumber]->oldtreeTips =  malloc (programOptions->TotalNumSequences* sizeof(TreeNode*));
        if (!(chains[chainNumber]->oldtreeTips))
        {
            fprintf (stderr, "Could not allocate the treeTips array\n");
            exit (-1);
        }
        InitChainPopulations(chains[chainNumber]->populations, chains[chainNumber]->numClones, programOptions->noisy, programOptions->TotalNumSequences );
        FillChainPopulationsFromPriors(chains[chainNumber]->populations, chains[chainNumber], programOptions,mcmcOptions, sampleSizes, seed );
        
        MakeCoalescenceTree2 (seed, chains[chainNumber]->populations,
                              &(chains[chainNumber]->numNodes),
                              chains[chainNumber]->numClones,
                              programOptions,
                              cumNumCA,
                              meanNumCA,
                              cumNumMIG,
                              meanNumMIG,
                              &numMIG,
                              &numCA,
                              &numEventsTot,
                              &(chains[chainNumber]->nodes),
                              chains[chainNumber]->treeTips,
                              &(chains[chainNumber]->root),
                              ObservedCellNames, sampleSizes
                              ) ;
        //        cumNumCA += numCA;
        //        cumNumMIG += numMIG;
        //        countTMRCA = treeRootInit[0]->timePUnits;
        //
        //        varTimeGMRCA[currentIteration] = countTMRCA;
        
        totalTreeLength = SumBranches(chains[chainNumber]->root, chains[chainNumber]->mutationRate);
        //        cumNumMUperTree=0;
        
        newickString2=NULL;
        newickString2 = toNewickString2 ( chains[chainNumber]->root, chains[chainNumber]->mutationRate,     programOptions->doUseObservedCellNames);
        printf("\n newick = %s  \n", newickString2);
        
        
        chains[chainNumber]->currentlogConditionalLikelihoodTree= LogConditionalLikelihoodTree(chains[chainNumber]->root, chains[chainNumber]->nodes, chains[chainNumber]->populations,  chains[chainNumber]->numClones);
        printf ( "Initial likelihood of the tree of chain %d is:  %lf \n",chainNumber, chains[chainNumber]->currentlogConditionalLikelihoodTree );
        
        chains[chainNumber]->currentlogConditionalLikelihoodSequences= LogConditionalLikelihoodSequences(chains[chainNumber], msa,  newickString2, programOptions);
        
        fprintf (stderr, "Initial likelihood of the sequences of chain %d  is = %lf  \n", chainNumber,chains[chainNumber]->currentlogConditionalLikelihoodSequences );
        
        free(newickString2);
        newickString2=NULL;
        
    }
}
/***************************** FillChainPopulationsFromPriors*******************************/
/* FillChainPopulationsFromPriors*/
void FillChainPopulationsFromPriors(Population **populations,Chain *chain,const ProgramOptions *programOptions, const MCMCoptions *mcmcOptions, int *sampleSizes, long int *seed)
{
    if (chain !=NULL && populations !=NULL)
    {
        int i;
        Population *popI;
        double randomDelta;
        chain->mutationRate = RandomLogUniform(mcmcOptions->MutRatefrom, mcmcOptions->MutRateto, seed);
        
        chain->totalPopSize= RandomLogUniform(mcmcOptions->totalEffectPopSizefrom, mcmcOptions->totalEffectPopSizeto, seed);
        
        double lambda = 1;
        
        SetPopulationsBirthRate(populations,  lambda, chain->numClones);
 
        GenerateEffectPopSizesFromPriors2(chain, programOptions->noisy,  chain->numClones,populations,  seed,  chain->gammaParam, chain->totalPopSize, YES);
        if (programOptions->populationSampleSizesKnown == NO)
        {
            InitPopulationSampleSizes(populations, programOptions->numCells, programOptions->numClones, chain->proportionsVector, seed);
        }
        //else fill the sample  sizes
        else{
            setChainPopulationSampleSizes(chain, sampleSizes, programOptions);
        }
        //generate the population  sizes
        
        InitPopulationsCoalescentEvents( chain->numClones,  populations) ;
        if (programOptions->populationSampleSizesKnown ==YES)
        {
            for( i = 0 ; i < chain->numClones; i++)
            {
                popI=*(populations + i);
                do {
                    randomDelta = RandomLogUniform(mcmcOptions->Deltafrom, mcmcOptions->Deltato, seed);
                    //popI->delta = chain->proportionsVector[i] * randomDelta;
                    popI->delta =  randomDelta;
                    popI->growthRate =popI->delta  / popI->effectPopSize;
                    popI->popSize=popI->effectPopSize * popI->birthRate;
                    popI->deathRate= popI->birthRate - popI->growthRate;
                }
                while(popI->popSize < popI->sampleSize);
            }
        }
        ListClonesAccordingTimeToOrigin2(populations, chain->numClones);
        GenerateTimesFromPriorsOriginal(programOptions->noisy, chain->numClones,chain->proportionsVector, populations, seed);
        ListClonesAccordingTimeToOrigin2(populations, chain->numClones);
        InitListPossibleMigrations(populations, chain->numClones);
    }
}
/***************************** pllmod_utree_rollback_spr1*******************************/
/* pllmod_utree_rollback_spr1 rollback SPR move*/
void pllmod_utree_rollback_spr1(pll_tree_rollback_t *rollback_info,pll_utree_rb_t *rb, TreeNode   **nodes,
                                double old_father_pop_effect_pop_size,
                                double new_father_pop_effect_pop_size,
                                TreeNode *container_of_p,
                                TreeNode *container_of_r,
                                TreeNode *container_of_p_prime,//mrca
                                TreeNode *container_of_r_prime,
                                TreeNode *container_of_u,
                                TreeNode *container_of_v
                                )
{
    int i, j, k;
    //rollback_info->SPR.
    //                rollback_info->rearrange_type     = PLLMOD_TREE_REARRANGE_SPR;
    //                rollback_info->rooted             = 0;
    //                rollback_info->SPR.prune_edge     = (void *) p_edge;
    //                rollback_info->SPR.regraft_edge   = (void *) p_edge->next->back;
    //                rollback_info->SPR.prune_bl       = p_edge->length;
    //                rollback_info->SPR.prune_left_bl  = p_edge->next->length;
    //                rollback_info->SPR.prune_right_bl = p_edge->next->next->length;
    //                rollback_info->SPR.regraft_bl     = r_edge->length;
    pll_unode_t * p_edge;
    pll_unode_t *regraft_edge ;
    pll_unode_t * r_edge;
    
    double time_from_r_to_p;
    
    pll_unode_t* p = rb->spr.p;// p;
    pll_unode_t* r = rb->spr.r;// r;
    pll_unode_t* r_prime= rb->spr.rb ;//= r->back;
    double r_len = rb->spr.r_len; //= r->length;
    pll_unode_t*  u= rb->spr.pnb; //= p->next->back;
    double length_to_u = rb->spr.pnb_len; //= p->next->length;
    pll_unode_t* v = rb->spr.pnnb ;// p->next->next->back;
    double v_len = rb->spr.pnnb_len; //= p->next->next->length;
  
    p_edge =(pll_unode_t*) rollback_info->SPR.prune_edge;
    regraft_edge =  rollback_info->SPR.regraft_edge  ;// = (void *) p_edge->next->back;
    double prune_bl   = rollback_info->SPR.prune_bl  ;//     = p_edge->length;
    double prune_left_bl = rollback_info->SPR.prune_left_bl ; //= p_edge->next->length;
    double prune_right_bl = rollback_info->SPR.prune_right_bl  ; //p_edge->next->next->length;
    double regraft_bl  = rollback_info->SPR.regraft_bl ;  //  = r_edge->length;
    
    /* (c) connect r<->r' */
    utree_link(r,
               r_prime,
               r_len,
               r->pmatrix_index);
    
    double length_from_r_prime_to_p;
    length_from_r_prime_to_p = r->length  -   time_from_r_to_p  ;

    container_of_r_prime->anc1=container_of_r;
    
    if (container_of_r_prime->anc1->left == container_of_r_prime)
    {//container_of_r_prime is left child of container_of_r
        container_of_r->left=NULL;
        container_of_r->left=container_of_r_prime;
    }
    else
    {//container_of_r_prime is right child of container_of_r
        container_of_r->right=NULL;
        container_of_r->right=container_of_r_prime;
    }
//    pll_unode_t * u = p->next->back;
//    pll_unode_t * v = p->next->next->back;
    utree_link(u,
               p->next,
               length_to_u ,
               u->pmatrix_index);
    
    utree_link(v,
               p->next->next,
               v_len ,
               u->pmatrix_index);
    
    /*connect the TreeNodes */

   // TreeNode * container_of_u = container_of_p->anc1;
    container_of_p->anc1 = container_of_u;
  
    //TreeNode * container_of_v ;
    if (container_of_v->anc1->left == container_of_v)//container_of_v is left child
        container_of_u->left=container_of_p;
    else //container_of_v is right child
        container_of_u->right=container_of_p;
    
    // TreeNode * container_of_p_prime =  *nodes + p->back->node_index;
    //disconnect container_of_p
    
   if (container_of_p_prime->anc1->left == container_of_p)//container_of_p_prime is left child
       container_of_p->right=container_of_v;
   else //container_of_p_prime is right child
       container_of_p->left=container_of_v;
}
/***************************** newTotalEffectivePopulationSizeMove*******************************/
/* newTotalEffectivePopulationSizeMove*/
void newTotalEffectivePopulationSizeMove(Chain *chain, long int *seed,  ProgramOptions *programOptions, char *ObservedCellNames[], pll_msa_t * msa, MCMCoptions * mcmcOptions, int *sampleSizes)
{
    int i,j;
    Population *popJ;
    double totalTreeLength;
    int        numCA, numMIG;
    double      cumNumCA, meanNumCA, cumNumMIG, meanNumMIG, numEventsTot;
    double newlogLikelihoodTree;
    Population *popI;
    // save the current values of
    chain->oldtotalPopSize = chain->totalPopSize;
    double newTotalPopulationSize;
     proposalSlidingWindow(  &newTotalPopulationSize,  chain->oldtotalPopSize,  mcmcOptions->slidingWindowSizeTotalEffectPopSize, seed);
    
    chain->totalPopSize = (int)newTotalPopulationSize;
    // RandomLogUniform(mcmcOptions->totalEffectPopSizefrom,mcmcOptions->totalEffectPopSizeto, seed);
    // chain->totalPopSize= RandomLogUniform(7,13, seed);
     GenerateEffectPopSizesFromPriors2(chain, programOptions->noisy,  chain->numClones, chain->populations,  seed,  chain->gammaParam, chain->totalPopSize,   NO);
    
    //setChainPopulationSampleSizes(chain, sampleSizes, programOptions);
    InitPopulationsCoalescentEvents( chain->numClones,  chain->populations) ;
  
    if (programOptions->populationSampleSizesKnown ==YES)
    {
        for( i = 0 ; i < chain->numClones; i++)
        {
            popI=*(chain->populations + i);
            do {
                popI->growthRate =popI->delta  / popI->effectPopSize;
                popI->popSize=popI->effectPopSize * popI->birthRate;
                popI->deathRate= popI->birthRate - popI->growthRate;
            }
            while(popI->popSize < popI->sampleSize);
        }
    }

    ListClonesAccordingTimeToOrigin2(chain->populations, chain->numClones);

//    InitListPossibleMigrations(chain->populations, chain->numClones);
//
//    MakeCoalescenceTree2 (seed, chain->populations,
//                          &(chain->numNodes),
//                          chain->numClones,
//                          programOptions,
//                          cumNumCA,
//                          meanNumCA,
//                          cumNumMIG,
//                          meanNumMIG,
//                          &numMIG,
//                          &numCA,
//                          &numEventsTot,
//                          &(chain->oldnodes),
//                          chain->oldtreeTips,
//                          &(chain->oldroot),
//                          ObservedCellNames,
//                          sampleSizes
//                          ) ;
    
   // double newLogConditionalLikelihoodTree= LogConditionalLikelihoodTree(chain->oldroot, chain->oldnodes, chain->populations,  chain->numClones);
    
     double newLogConditionalLikelihoodTree= LogConditionalLikelihoodTree(chain->root, chain->nodes, chain->populations,  chain->numClones);

    fprintf (stderr, "\n>> New log conditional Likelihood tree after the new total effective population size,  of the chain %d is = %lf  \n", chain->chainNumber,newLogConditionalLikelihoodTree );
    
    char *newickString2;
    newickString2=NULL;
    newickString2 = toNewickString2 ( chain->oldroot, programOptions->mutationRate,     programOptions->doUseObservedCellNames);
    printf("\n newick after move= %s  \n", newickString2);
    newickString2 = toNewickString4 ( chain->oldroot, programOptions->mutationRate,     programOptions->doUseObservedCellNames);
    printf("\n newick after move= %s  \n", newickString2);
    free(newickString2);
    newickString2=NULL;
    
    double priorDensityNewTotalEffectivePopulationSize= LogUniformDensity(newTotalPopulationSize, mcmcOptions->totalEffectPopSizefrom, mcmcOptions->totalEffectPopSizeto);
   double priorDensityCurrentTotalEffectivePopulationSize= LogUniformDensity(chain->oldtotalPopSize, mcmcOptions->totalEffectPopSizefrom, mcmcOptions->totalEffectPopSizeto);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewTotalEffectivePopulationSize;
   
    double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensityCurrentTotalEffectivePopulationSize;
    
    double randomNumber= RandomUniform (seed);
    
    double LogAcceptanceRate = (sumLogNumerators - sumLogDenominators) >0? (sumLogNumerators - sumLogDenominators) :0;
    
    if (log(randomNumber) < LogAcceptanceRate )
    {//accept the move
        printf("\n Accepted new total effective population size move\n");
        //free nodes, treeTips, treeRootInit and make them point to the new ones
       // free (chain->nodes);
       // chain->nodes=NULL;
        //free (chain->root);
       // chain->root=NULL;
        chain->nodes = chain->oldnodes;
        chain->root = chain->oldroot;
    }
    else
    {
        //generate a random number
        double randomNumber= RandomUniform (seed);
       
        //reject the move
        
        
        printf("\n Rejected new total effective population size move\n");
        chain->totalPopSize = chain->oldtotalPopSize;
        for (j = 0; j < chain->numClones; j++)
        {
            popJ = *(chain->populations + j);
            popJ->effectPopSize = popJ->oldeffectPopSize;
        }
        //free oldnodes, oldtreeTips, oldtreeRootInit
        free (chain->oldnodes);
        chain->oldnodes=NULL;
        free (chain->oldroot);
      chain->oldroot=NULL;
    }
}
/***************************** setPopulationSampleSizes*******************************/
/* setPopulationSampleSizes*/
void setChainPopulationSampleSizes(Chain *chain, int *sampleSizes,  ProgramOptions *programOptions)
{
    int i;
    Population *popI;
    int cumSampleSize=0;
    
    if (sampleSizes== NULL)
        fprintf (stderr, "ERROR: the sample sizes  vector is null..\n");
    if (chain->numClones >0 && chain->populations!=NULL)
    {
        for( i = 0 ; i < chain->numClones; i++)
        {
            popI=*(chain->populations + i);
            popI->sampleSize=sampleSizes[i];
            popI->indexFirstObservedCellName =cumSampleSize;
            cumSampleSize +=sampleSizes[i];
        }
    }
}
/***************************** densityLogUniform*******************************/
/* densityLogUniform*/
double LogUniformDensity(double value, double from, double to)
{
    double result;
    if (value >= exp(from) && value<= exp(to) )
    {
       result = 1 / (value * (to-from));
        return log(result);
    }
    else
    {
        result=0;
        return result;
    }
}
/***************************** InitPopulationsCoalescentEvents*******************************/
/* InitPopulationsCoalescentEvents*/
void InitPopulationsCoalescentEvents( int numClones,  Population **populations) {
    int i;
    Population *popI;
    for( i = 0 ; i < numClones; i++)
    {
        popI=*(populations + i);
        popI->CoalescentEventTimes=(double *) calloc((popI->sampleSize -1)+ numClones -1, (long) sizeof(double));
        if (!(popI->CoalescentEventTimes))
        {
            fprintf (stderr, "CoalescentEventTimes (%lu bytes)\n", (popI->sampleSize -1 +numClones -1) * (long) sizeof(double));
            exit (-1);
        }
    }
}
/***************************** compute_state_probs*******************************/
/* compute_state_probs*/
void compute_state_probs(unsigned int state, double  **clvp,  unsigned int _states, int from, int to,
                         double _seq_error_rate,
                         double _dropout_rate)
{
    unsigned int state_id = __builtin_ctz(state);
    //unsigned int state_id = state;
    static const double one_3 = 1. / 3.;
    static const double one_6 = 1. / 6.;
    static const double one_8 = 1. / 8.;
    static const double three_8 = 3. / 8.;
    static const double one_12 = 1. / 12.;
    
    // TODO: move it out of here
    unsigned int undef_state = (unsigned int) (pow(2, _states)) - 1;
    
    double sum_lh = 0.;
    
    for (size_t k = 0; k < _states; ++k)
    {
        if (state == undef_state)
            (*clvp)[k] = 1.;
        else
        {
            if (k == state_id)
            {
                /* 0 letters away */
                if (HOMO(state_id))
                    (*clvp)[k] = 1. - _seq_error_rate + 0.5 * _seq_error_rate * _dropout_rate;
                else
                    (*clvp)[k] =  (1. - _dropout_rate ) * (1. - _seq_error_rate) + one_12 * _seq_error_rate * _dropout_rate;
            }
            else if (mut_dist[state_id][k] == 1)
            {
                /* 1 letter away */
                if (HOMO(k))
                {
                    (*clvp)[k] = one_12 * _seq_error_rate * _dropout_rate +
                    one_3  * (1. - _dropout_rate) * _seq_error_rate;
                }
                else
                {
                    if (HOMO(state_id))
                    {
                        (*clvp)[k] = 0.5 * _dropout_rate + one_6 * _seq_error_rate -
                        three_8 * _seq_error_rate * _dropout_rate;
                    }
                    else
                    {
                        (*clvp)[k] = one_6 * _seq_error_rate -
                        one_8 * _seq_error_rate * _dropout_rate;
                    }
                }
            }
            else
            {
                /* 2 letters away */
                if (HOMO(state_id))
                    (*clvp)[k] = one_12 * _seq_error_rate * _dropout_rate;
                else
                    (*clvp)[k] = 0.;
            }
            sum_lh += (*clvp)[k];
        }
    }
}
/***************************** computeGenotypesFreq*******************************/
/* computeGenotypesFreq*/
void computeGenotypesFreq(double freqs[10], pll_msa_t * msa)
{
    int i, j;
    int sum=0;;
    char genotypeIUPAC;
    int index;
    char * currentSequence;
    
     for (j= 0; j < 10; j++)
         freqs[j]=0;
  
      for (j= 0; j < msa->count; j++)
      {
          currentSequence=msa->sequence[j];
          for (i= 0; i < msa->length; i++)
          {
              genotypeIUPAC =currentSequence[i];
              if (genotypeIUPAC == 'A')
                  index=0;
              else if (genotypeIUPAC == 'C')
                  index=1;
              else if (genotypeIUPAC == 'G')
                   index=2;
              else if (genotypeIUPAC == 'T')
                   index=3;
              else if (genotypeIUPAC == 'R')
                   index=4;
              else if (genotypeIUPAC == 'M')
                   index=5;
              else if (genotypeIUPAC == 'W')
                   index=6;
              else if (genotypeIUPAC == 'S')
                   index=7;
              else if (genotypeIUPAC == 'Y')
                   index=8;
              else if (genotypeIUPAC == 'K')
                   index=9;
              else{
                  fprintf (stderr, "\nERROR in computeGenotypesFreq: nucleotide = %c\n",  genotypeIUPAC);
                  exit(-1);
              }
             freqs[index]= freqs[index]+1;
              sum = sum +1;
          }
      }
    for (j= 0; j < 10; j++)
        freqs[j]=freqs[j] /sum;
    
}
/***************************** AssignObservedCellNamestoNodes*******************************/
/* AssignObservedCellNamestoNodes*/

void AssignObservedCellNamestoTips(TreeNode **nodes, Population ** populations, int *sampleSizes, char* ObservedCellNames[],   ProgramOptions *programOptions){
    
    int i, j,k,  currentSampleSize, cumSumSampleSizes;
    TreeNode *p;
    Population *popI;
    cumSumSampleSizes=0;
    for (i = 0; i < programOptions->numClones; i++)
    {
        popI = *(populations + i);
  
       popI->CellAssignationCompleted =NO;
      
    }
    
    for (j = 0; j < programOptions->numClones; j++)
    {
         currentSampleSize= sampleSizes[j];
        
        for (i = 0; i < programOptions->numClones; i++)
        {
             popI = *(populations + i);
            if (currentSampleSize ==  popI->sampleSize && popI->CellAssignationCompleted==NO)
            {
                 popI->CellAssignationCompleted =YES;
                
                 break;
            }
        }
        
        for (k = cumSumSampleSizes; k < cumSumSampleSizes + currentSampleSize; k++)
        {
            p = *nodes + k;
           
            strcpy( p->observedCellName,ObservedCellNames[k]);
            
        }
        cumSumSampleSizes = cumSumSampleSizes + currentSampleSize;
    }

}
/***************************** AssignObservedCellNamestoNodes2*******************************/
/* AssignObservedCellNamestoNodes2*/

void AssignObservedCellNamestoTips2(TreeNode **nodes, Population ** populations, int *sampleSizes, char* ObservedCellNames[],   ProgramOptions *programOptions){
    
    int i, j,k,  currentSampleSize, cumSumSampleSizes;
    TreeNode *p;
    Population *popI;
    cumSumSampleSizes=0;
    int indexFirstObservedCellName;
    
    for (i = 0; i < programOptions->numClones; i++)
    {
        popI = *(populations + i);
        
        for (j = 0; j < popI->numActiveGametes ; j++)
        {
            k=popI->idsActiveGametes[j];
            p = *nodes + k;
            indexFirstObservedCellName= popI->indexFirstObservedCellName;
            strcpy( p->observedCellName,ObservedCellNames[indexFirstObservedCellName + j ]);
            
        }
    }
    
}
double * pllmod_msa_empirical_frequencies1(pll_partition_t * partition)
{
    unsigned int i, j, k, n;
    unsigned int states         = partition->states;
    unsigned int states_padded  = partition->states_padded;
    unsigned int sites          = partition->sites;
    unsigned int rate_cats      = partition->rate_cats;
    unsigned int tips           = partition->tips;
    const pll_state_t * tipmap  = partition->tipmap;
    const unsigned int * w      = partition->pattern_weights;
    double * frequencies;
    
    if ((frequencies = (double *) calloc ((size_t) states, sizeof(double)))
        == NULL)
    {
        pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                         "Cannot allocate memory for empirical frequencies");
        return NULL;
    }
    
    if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
    {
        if (states == 4)
        {
            for (i = 0; i < tips; ++i)
            {
                const unsigned char *tipchars = partition->tipchars[i];
                for (n = 0; n < sites; ++n)
                {
                    unsigned int state = (unsigned int) tipchars[n];
                    double sum_site = 1.0 * PLL_POPCNT32(state);
                    for (k = 0; k < states; ++k)
                    {
                        if (state & 1)
                            frequencies[k] += w[n] / sum_site;
                        state >>= 1;
                    }
                }
            }
        }
        else
        {
            for (i = 0; i < tips; ++i)
            {
                const unsigned char *tipchars = partition->tipchars[i];
                for (n = 0; n < sites; ++n)
                {
                    pll_state_t state = tipmap[(int) tipchars[n]];
                    double sum_site = 1.0 * PLL_STATE_POPCNT(state);
                    for (k = 0; k < states; ++k)
                    {
                        if (state & 1)
                            frequencies[k] += w[n] / sum_site;
                        state >>= 1;
                    }
                }
            }
        }
    }
    else
    {
        for (i = 0; i < tips; ++i)
        {
            unsigned int *site_to_id = 0;
            if ((partition->attributes & PLL_ATTRIB_SITE_REPEATS)
                && partition->repeats->pernode_ids[i]) {
                site_to_id = partition->repeats->pernode_site_id[i];
            }
            
            for (n = 0; n < sites; ++n)
            {
                j = site_to_id ? site_to_id[n] : n;
                j *= (states_padded * rate_cats);
                double sum_site = 0.0;
                for (k = 0; k < states; ++k)
                    sum_site += partition->clv[i][j + k];
                
                for (k = 0; k < states; ++k)
                    frequencies[k] += w[n] * partition->clv[i][j + k] / sum_site;
            }
        }
    }
    
    /* IMPORTANT: we must use the original number of sites in alignment (before pattern compression),
     * since we previously multiplied our base counts with the respective column weights! */
    unsigned int uncomp_sites = 0;
    for (i = 0; i < sites; ++i)
        uncomp_sites += w[i];
    
    for (k = 0; k < states; ++k)
        frequencies[k] /= uncomp_sites * tips;
    
#ifndef NDEBUG
    double sum_test = 0.0;
    for (k = 0; k < states; ++k)
    {
        sum_test += frequencies[k];
    }
    assert(fabs (sum_test - 1) < 1e-6);
#endif
    
    return frequencies;
}
//int pll_set_tip_states1(pll_partition_t * partition,
//                                  unsigned int tip_index,
//                                  const pll_state_t * map,
//                                  const char * sequence)
//{
//    int rc;
//
//    if (pll_repeats_enabled(partition))
//    {
//        if (PLL_FAILURE == pll_update_repeats_tips(partition, tip_index, map, sequence))
//            return PLL_FAILURE;
//    }
//    if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
//    {
//        /* create (or update) character map for tip-tip precomputations */
//        if (partition->tipchars)
//        {
//            update_charmap(partition,map);
//        }
//        else
//        {
//            if (!create_charmap(partition,map))
//            {
//                dealloc_partition_data(partition);
//                return PLL_FAILURE;
//            }
//        }
//
//        if (partition->states == 4)
//            rc = set_tipchars_4x4(partition, tip_index, map, sequence);
//        else
//            rc = set_tipchars(partition, tip_index, map, sequence);
//    }
//    else
//        rc = set_tipclv1(partition, tip_index, map, sequence);
//
//    return rc;
//}
//
//
//
static int set_tipclv1(pll_partition_t * partition,
                      unsigned int tip_index,
                      const pll_state_t * map,
                      const char * sequence,
                       double _seq_error_rate,
                       double _dropout_rate
                       )
{
    pll_state_t c;
    unsigned int i,j;
    double * tipclv = partition->clv[tip_index];
    unsigned int state_id ;

    pll_repeats_t * repeats = partition->repeats;
    unsigned int use_repeats = pll_repeats_enabled(partition);
    unsigned int ids = use_repeats ?
    repeats->pernode_ids[tip_index] : partition->sites;
    
    static const double one_3 = 1. / 3.;
    static const double one_6 = 1. / 6.;
    static const double one_8 = 1. / 8.;
    static const double three_8 = 3. / 8.;
    static const double one_12 = 1. / 12.;
    
    // TODO: move it out of here
    unsigned int undef_state = (unsigned int) (pow(2, partition->states)) - 1;
    
    double sum_lh = 0.;
    
    /* iterate through sites */
    for (i = 0; i < ids; ++i)
    {
        unsigned int index = use_repeats ?
        repeats->pernode_id_site[tip_index][i] : i;
        if ((c = map[(int)sequence[index]]) == 0)
        {
            pll_errno = PLL_ERROR_TIPDATA_ILLEGALSTATE;
            snprintf(pll_errmsg, 200, "Illegal state code in tip \"%c\"", sequence[index]);
            return PLL_FAILURE;
        }

        /* decompose basecall into the encoded residues and set the appropriate
         positions in the tip vector */
         state_id = __builtin_ctz(c);
        for (j = 0; j < partition->states; ++j)
        {
            if (c == undef_state)
                tipclv[j] = 1.;
            else
            {
                if (j == state_id)
                {
                    /* 0 letters away */
                    if (HOMO(state_id))
                        tipclv[j] = 1. - _seq_error_rate + 0.5 * _seq_error_rate * _dropout_rate;
                    else
                        tipclv[j] =  (1. - _dropout_rate ) * (1. - _seq_error_rate) + one_12 * _seq_error_rate * _dropout_rate;
                }
                else if (mut_dist[state_id][j] == 1)
                {
                    /* 1 letter away */
                    if (HOMO(j))
                    {
                        tipclv[j] = one_12 * _seq_error_rate * _dropout_rate +
                        one_3  * (1. - _dropout_rate) * _seq_error_rate;
                    }
                    else
                    {
                        if (HOMO(state_id))
                        {
                            tipclv[j] = 0.5 * _dropout_rate + one_6 * _seq_error_rate -
                            three_8 * _seq_error_rate * _dropout_rate;
                        }
                        else
                        {
                            tipclv[j]= one_6 * _seq_error_rate -
                            one_8 * _seq_error_rate * _dropout_rate;
                        }
                    }
                }
                else
                {
                    /* 2 letters away */
                    if (HOMO(state_id))
                        tipclv[j] = one_12 * _seq_error_rate * _dropout_rate;
                    else
                        tipclv[j] = 0.;
                }
                sum_lh += tipclv[j];
            }
           // tipclv[j] = c & 1;
          //  c >>= 1;
        }

        /* fill in the entries for the other gamma values */
        tipclv += partition->states_padded;
        for (j = 0; j < partition->rate_cats - 1; ++j)
        {
            memcpy(tipclv, tipclv - partition->states_padded,
                   partition->states * sizeof(double));
            tipclv += partition->states_padded;
        }
    }

    /* if asc_bias is set, we initialize the additional positions */
    if (partition->asc_bias_alloc)
    {
        for (i = 0; i < partition->states; ++i)
        {
            for (j = 0; j < partition->states; ++j)
            {
                tipclv[j] = j==i;
            }

            /* fill in the entries for the other gamma values */
            tipclv += partition->states_padded;
            for (j = 0; j < partition->rate_cats - 1; ++j)
            {
                memcpy(tipclv, tipclv - partition->states_padded,
                       partition->states * sizeof(double));
                tipclv += partition->states_padded;
            }
        }
    }

  return PLL_SUCCESS;
}
/***************************** newProportionsVectorMove*******************************/
/* newProportionsVectorMove*/
void newProportionsVectorMove(Chain *chain, long int *seed,  ProgramOptions *programOptions, char *ObservedCellNames[], pll_msa_t * msa, MCMCoptions *mcmcOptions, int *sampleSizes)
{
    int i,j,k;
    Population *popJ;
    double totalTreeLength;
    int        numCA, numMIG;
    double      cumNumCA, meanNumCA, cumNumMIG, meanNumMIG, numEventsTot;
    double newlogLikelihoodTree;
    Population *popI;
    proposalProportionsVector(chain->oldproportionsVector, chain->proportionsVector, mcmcOptions->tuningParameter, chain->numClones, seed);
    
    GenerateEffectPopSizesFromPriors2(chain, programOptions->noisy,  chain->numClones, chain->populations,  seed,  chain->gammaParam, chain->totalPopSize, YES);
    if (programOptions->populationSampleSizesKnown == NO)
    {
        InitPopulationSampleSizes(chain->populations, programOptions->numCells, programOptions->numClones, chain->proportionsVector, seed);
    }
    //else fill the sample  sizes
    else{
        setChainPopulationSampleSizes(chain, sampleSizes, programOptions);
    }
    
    if (programOptions->populationSampleSizesKnown ==YES)
    {
        for( i = 0 ; i < chain->numClones; i++)
        {
            popI=*(chain->populations + i);
            do {
                popI->growthRate =popI->delta  / popI->effectPopSize;
                popI->popSize=popI->effectPopSize * popI->birthRate;
                popI->deathRate= popI->birthRate - popI->growthRate;
            }
            while(popI->popSize < popI->sampleSize);
        }
    }
    ListClonesAccordingTimeToOrigin2(chain->populations, chain->numClones);
    
//    InitListPossibleMigrations(chain->populations, chain->numClones);
    
//    MakeCoalescenceTree2 (seed, chain->populations,
//                          &(chain->numNodes),
//                          chain->numClones,
//                          programOptions,
//                          cumNumCA,
//                          meanNumCA,
//                          cumNumMIG,
//                          meanNumMIG,
//                          &numMIG,
//                          &numCA,
//                          &numEventsTot,
//                          &(chain->oldnodes),
//                          chain->oldtreeTips,
//                          &(chain->oldroot),
//                          ObservedCellNames,
//                          sampleSizes
//                          ) ;
    
    double newLogConditionalLikelihoodTree= LogConditionalLikelihoodTree(chain->oldroot, chain->oldnodes, chain->populations,  chain->numClones);
    
    fprintf (stderr, "\n>> New log conditional Likelihood tree after the new total effective population size,  of the chain %d is = %lf  \n", chain->chainNumber,newLogConditionalLikelihoodTree );
    
    char *newickString2;
    newickString2=NULL;
    newickString2 = toNewickString2 ( chain->oldroot, programOptions->mutationRate,     programOptions->doUseObservedCellNames);
    printf("\n newick after move= %s  \n", newickString2);
    newickString2 = toNewickString4 ( chain->oldroot, programOptions->mutationRate,     programOptions->doUseObservedCellNames);
    printf("\n newick after move= %s  \n", newickString2);
    free(newickString2);
    newickString2=NULL;
    

    double priorDensityNewProportionsVector =  DirichletDensity( chain->proportionsVector, chain->oldproportionsVector , chain->numClones);
    
    double priorDensityCurrentProportionsVector = DirichletDensity(chain->oldproportionsVector, chain->proportionsVector, chain->numClones);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewProportionsVector;
    
    double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensityCurrentProportionsVector;
    
    double randomNumber= RandomUniform (seed);
    
    double LogAcceptanceRate = (sumLogNumerators - sumLogDenominators) >0? (sumLogNumerators - sumLogDenominators) :0;
    
    if (log(randomNumber) < LogAcceptanceRate )
    {//accept the move
        printf("\n Accepted new proportions vector move\n");
        //free nodes, treeTips, treeRootInit and make them point to the new ones
        // free (chain->nodes);
        // chain->nodes=NULL;
        //free (chain->root);
        // chain->root=NULL;
        chain->nodes = chain->oldnodes;
        chain->root = chain->oldroot;
        chain->currentlogConditionalLikelihoodTree =newLogConditionalLikelihoodTree;
    }
    else
    {
        //reject the move
        printf("\n Rejected new proportions vector move\n");
        chain->totalPopSize = chain->oldtotalPopSize;
        for (j = 0; j < chain->numClones; j++)
        {
            popJ = *(chain->populations + j);
            popJ->effectPopSize = popJ->oldeffectPopSize;
        }
        //free oldnodes, oldtreeTips, oldtreeRootInit
        // free (chain->oldnodes);
        // chain->oldnodes=NULL;
        //free (chain->oldroot);
        // chain->oldroot=NULL;
    }
}
/***************************** newScaledGrowthRateMoveforPopulation*******************************/
/* newScaledGrowthRateMoveforPopulation*/
void newScaledGrowthRateMoveforPopulation(Chain *chain, Population *popI, long int *seed,  ProgramOptions *programOptions, char *ObservedCellNames[], pll_msa_t * msa, MCMCoptions * mcmcOptions, int *sampleSizes)
{

    ListClonesAccordingTimeToOrigin2(chain->populations, chain->numClones);
    double totalTreeLength;
    int        numCA, numMIG;
    double      cumNumCA, meanNumCA, cumNumMIG, meanNumMIG, numEventsTot;
    double newlogLikelihoodTree;
   
    //popI=*(populations + index);
    double randomDelta = RandomLogUniform(mcmcOptions->Deltafrom, mcmcOptions->Deltato, seed);
    //popI->delta = chain->proportionsVector[i] * randomDelta;
    popI->olddelta= popI->delta;
    popI->delta =  randomDelta;
    popI->growthRate =popI->delta  / popI->effectPopSize;
    popI->deathRate= popI->birthRate - popI->growthRate;

    ListClonesAccordingTimeToOrigin2(chain->populations, chain->numClones);
//
//    InitListPossibleMigrations(chain->populations, chain->numClones);
//
//    MakeCoalescenceTree2 (seed, chain->populations,
//                          &(chain->numNodes),
//                          chain->numClones,
//                          programOptions,
//                          cumNumCA,
//                          meanNumCA,
//                          cumNumMIG,
//                          meanNumMIG,
//                          &numMIG,
//                          &numCA,
//                          &numEventsTot,
//                          &(chain->oldnodes),
//                          chain->oldtreeTips,
//                          &(chain->oldroot),
//                          ObservedCellNames,
//                          sampleSizes
//                          ) ;
    
   // double newLogConditionalLikelihoodTree= LogConditionalLikelihoodTree(chain->oldroot, chain->oldnodes, chain->populations,  chain->numClones);
    
    double newLogConditionalLikelihoodTree= LogConditionalLikelihoodTree(chain->root, chain->nodes, chain->populations,  chain->numClones);
    
    fprintf (stderr, "\n>> New log conditional Likelihood tree after the new total growth rate for population %d,  of the chain %d is = %lf  \n", popI->index, chain->chainNumber,newLogConditionalLikelihoodTree );
    
    char *newickString2;
    newickString2=NULL;
    newickString2 = toNewickString2 ( chain->oldroot, programOptions->mutationRate,     programOptions->doUseObservedCellNames);
    printf("\n newick after move= %s  \n", newickString2);
    
    newickString2 = toNewickString4 ( chain->oldroot, programOptions->mutationRate,     programOptions->doUseObservedCellNames);
    printf("\n newick after move= %s  \n", newickString2);

    free(newickString2);
    newickString2=NULL;
    
    
    double priorDensityNewScaledGrowthRate =  LogUniformDensity(randomDelta, mcmcOptions->Deltafrom, mcmcOptions->Deltato);
 
    double priorDensitScaledGrowthRate =  LogUniformDensity(popI->olddelta, mcmcOptions->Deltafrom, mcmcOptions->Deltato);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewScaledGrowthRate;
    
    double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensitScaledGrowthRate;
    
    double randomNumber= RandomUniform (seed);
    
    double LogAcceptanceRate = (sumLogNumerators - sumLogDenominators) >0? (sumLogNumerators - sumLogDenominators) :0;
    
    if (log(randomNumber) < LogAcceptanceRate )
    {//accept the move
        printf("\n Accepted new growth rate move\n");
        //free nodes, treeTips, treeRootInit and make them point to the new ones
        // free (chain->nodes);
        // chain->nodes=NULL;
        //free (chain->root);
        // chain->root=NULL;
        chain->nodes = chain->oldnodes;
        chain->root = chain->oldroot;
        chain->currentlogConditionalLikelihoodTree =newLogConditionalLikelihoodTree;
    }
    else
    {
        //reject the move
        printf("\n Rejected new growth rate move\n");
         popI->olddelta= popI->delta;
        //free oldnodes, oldtreeTips, oldtreeRootInit
        // free (chain->oldnodes);
        // chain->oldnodes=NULL;
        //free (chain->oldroot);
        // chain->oldroot=NULL;
    }
}
void InitNumberNodes(double *TotalBirthRate, double *TotalDeathRate, int *TotalN,  Population **populations, ProgramOptions *programOptions) {
    programOptions->TotalNumSequences = 0;
    *TotalN = 0;
    *TotalBirthRate = 0.0;
    *TotalDeathRate = 0.0;
    Population* popI;
    int j;
    for (j = 0; j < programOptions->numClones; j++)
    {   popI = *(populations + j);
        popI->FatherPop =NULL;
        programOptions->TotalNumSequences = programOptions->TotalNumSequences + popI->sampleSize;
        *TotalN = *TotalN + popI->popSize;
        *TotalBirthRate = *TotalBirthRate + popI->birthRate;
        *TotalDeathRate = *TotalDeathRate + popI->deathRate;
    }
    programOptions->numNodes = 2 * programOptions->TotalNumSequences + programOptions->numClones+ 10;
    
    programOptions->numCells =programOptions->TotalNumSequences;
}
