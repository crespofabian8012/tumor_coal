//
//  fileFuntions.h
//  TumorCellCoal
//
//  Created by Fausto Fabian Crespo Fernandez on 7/8/19.
//  Copyright © 2019 Fausto Fabian Crespo Fernandez. All rights reserved.
//

#ifndef fileFuntions_h
#define fileFuntions_h

static void PrintInfo(int noisy, int thereisOutgroup, int outgroupSelection, double outgroupBranchLength_Root1Root2, double outgroupBranchLength_RootSample,  int numClones, int *TotalNumSequences, int *TotalN, double  *TotalBirthRate, double mutationRate, double *TotalDeathRate, int *CloneNameBegin, double *CloneTimeOriginInput, double *CloneTimeOriginInputSTDPOP, double *CloneDelta, double *expectedPopSize, int *ClonePopSizeBegin, int *CloneSampleSizeBegin, double *CloneBirthRateBegin, double *CloneDeathRateBegin);


void PrintTimes(int replicate, FILE            *fpTimes, double mutationRate, TreeNode *nodes,  int thereisOutgroup);
static void   ListTimes (int j, double mutationRate, TreeNode *nodes,  FILE *fpTimes, int thereisOutgroup);

static void PrintTimes2(int replicate, FILE  *fpTimes2, double mutationRate, TreeNode *nodes, int thereisOutgroup);
static void   ListTimes2 (int j, double mutationRate, TreeNode *nodes,  FILE *fpTimes, int thereisOutgroup);




static void PrintRunSettings (ProgramOptions programoptions, long int seed,  int TotalNumSequences,
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
                              double  *expectedPopSize,
                              double cumNumCA,
                              double meanNumCA,
                              double cumNumMIG,
                              double meanNumMIG,
                              double *numEventsTot
                              );
static void PrepareSeparateFiles(int ChainNumber, int paramSetNumber, int replicate,const FilePaths *filePaths, const ProgramOptions *programOptions,Files *files);
//static void PrepareSeparateFiles(int paramSetNumber, int replicate, FilePaths filePaths, ProgramOptions programOptions,Files *files
//                          ,int doPrintTrees, FILE **fpTrees, FILE **fpTrees2,char   resultsDir[MAX_NAME] ,  char        treeFile[MAX_NAME], char        treeDir[MAX_NAME],char        timesDir[MAX_NAME],
//                          char        SNVgenotypesDir[MAX_NAME],char        SNVhaplotypesDir[MAX_NAME],
//                          char        SNVGenotypesDir[MAX_NAME],
//                          char MLhaplotypesDir[MAX_NAME],
//                          char        FullhaplotypesDir[MAX_NAME], char fullGenotypesDir[MAX_NAME],char VCFdir[MAX_NAME],
//                          char CATGdir[MAX_NAME],
//                          char        TruehaplotypesDir[MAX_NAME],char fullHaplotypesDir[MAX_NAME],
//                          char        timesFile[MAX_NAME], char        SNVgenotypesFile[MAX_NAME],char SNVhaplotypesFile[MAX_NAME],char trueHaplotypesDir[MAX_NAME], char fullGenotypesFile[MAX_NAME],
//                          char trueHaplotypesFile[MAX_NAME],char MLhaplotypesFile[MAX_NAME], char fullHaplotypesFile[MAX_NAME],char VCFfile[MAX_NAME],
//                          char CATGfile[MAX_NAME] ,FILE            **fpSNVgenotypes, int doPrintTimes, int doSimulateData, FILE  **fpTimes,
//                          FILE  **fpTimes2,
//                          FILE **fpSNVhaplotypes, FILE **fpTrueHaplotypes,
//                          FILE **fpFullGenotypes,FILE **fpFullHaplotypes, FILE **fpVCF, FILE **fpCATG, FILE **fpMLhaplotypes, int doPrintSNVgenotypes, int  doPrintSNVhaplotypes, int doPrintTrueHaplotypes, int doPrintFullGenotypes, int doPrintFullHaplotypes, int doNGS,  int doPrintCATG,  int doPrintMLhaplotypes,int doSimulateReadCounts
//);
static void PrepareSeparateFilesGenotypes(int paramSetNumber, int TreeNum,int MutationAssignNum,
                                   const FilePaths *filePaths, const ProgramOptions *programOptions,Files *files);
//void PrepareSeparateFilesGenotypes(int paramSetNumber, int TreeNum,int MutationAssignNum,
 //                                  FilePaths filePaths, ProgramOptions programOptions,Files *files
//                                   int doPrintTrees, FILE **fpTrees, FILE **fpTrees2,char   resultsDir[MAX_NAME] ,  char        treeFile[MAX_NAME], char        treeDir[MAX_NAME], char        timesDir[MAX_NAME],
//                                  char        SNVgenotypesDir[MAX_NAME],char        SNVhaplotypesDir[MAX_NAME],
//                                  char        SNVGenotypesDir[MAX_NAME],
//                                  char MLhaplotypesDir[MAX_NAME],
//                                  char        FullhaplotypesDir[MAX_NAME], char fullGenotypesDir[MAX_NAME],char VCFdir[MAX_NAME],
//                                  char CATGdir[MAX_NAME],
//                                  char        TruehaplotypesDir[MAX_NAME],char fullHaplotypesDir[MAX_NAME],
//                                  char        timesFile[MAX_NAME], char        SNVgenotypesFile[MAX_NAME],char SNVhaplotypesFile[MAX_NAME],char trueHaplotypesDir[MAX_NAME], char fullGenotypesFile[MAX_NAME],
//                                  char trueHaplotypesFile[MAX_NAME],char MLhaplotypesFile[MAX_NAME], char fullHaplotypesFile[MAX_NAME],char VCFfile[MAX_NAME],
//                                  char CATGfile[MAX_NAME] ,FILE            **fpSNVgenotypes, int doPrintTimes, int doSimulateData, FILE  **fpTimes,
//                                  FILE  **fpTimes2,
//                                  FILE **fpSNVhaplotypes, FILE **fpTrueHaplotypes,
//                                  FILE **fpFullGenotypes,FILE **fpFullHaplotypes, FILE **fpVCF, FILE **fpCATG, FILE **fpMLhaplotypes, int doPrintSNVgenotypes, int  doPrintSNVhaplotypes, int doPrintTrueHaplotypes, int doPrintFullGenotypes, int doPrintFullHaplotypes, int doNGS,  int doPrintCATG, int numDataSets, int doPrintMLhaplotypes,int doSimulateReadCounts
//);
static void ReadParametersFromFile(ProgramOptions *programOptions, FilePaths *filePaths,
                            int **CloneNameBegin,
                            int **CloneSampleSizeBegin,
                            int **ClonePopSizeBegin,
                            double **CloneBirthRateBegin,
                            double **CloneDeathRateBegin,
                            double **CloneTimeOriginInput,
                            //                            int    *Nscaling,
                            //                            int    *numClones,
                            //                            int  *numDataSets,
                            //                            long int *userSeed,
                            //                            double  *mutationRate,
                            //                            int *noisy,
                            //                            int *outgroupSelection,
                            //                            int *thereisOutgroup,
                            //                            double *outgroupBranchLength_Root1Root2,
                            //                            double *outgroupBranchLength_RootSample,
                            //char treeFile[20],
                            //char timesFile[20],
                            //                            int   *doPrintTrees,
                            //                            int      *doPrintTimes,
                            //                            int       *alphabet,
                            //                            double  *propAltModelSites,
                            //                            int *altModel,
                            //                            int *doJC,
                            //                            int *doHKY,
                            //                            int *doGTR,
                            //                            int* doGTnR,
                            //                            int *thereIsMij,
                            double Mij[4][4],
                            //                            double *ADOrate,
                            //                            double *allelicImbalance,
                            //                            double * meanAmplificationError,
                            //                            double *varAmplificationError,
                            //                            int *simulateOnlyTwoTemplates,
                            //                            double *sequencingError,
                            //                            int *numCells, int *doSimulateFixedNumMutations,int *numFixedMutations,
                            double freq[4]
//                            , int *equalBaseFreq, double*alphaCoverage, int *rateVarCoverage,
//                            double *SNPrate, int *coverage, int *doSimulateReadCounts, double *genotypingError
);


void CreateFolder(FilePaths filePaths, ProgramOptions programOptions,Files *files//,
//                  FILE  **fpTrees, FILE **fpTimes, FILE **fpTrees2, FILE **fpTimes2,
//                  FILE            **fpSNVgenotypes,   FILE **fpSNVhaplotypes, FILE **fpTrueHaplotypes,
//                  FILE **fpFullGenotypes,FILE **fpFullHaplotypes, FILE **fpVCF, FILE **fpCATG
) ;
int openFile(FILE **file, char path[MAX_NAME] );

void PrepareGlobalFiles(int argc, char **argv, int doPrintTree, FILE *fpTrees, char   resultsDir[MAX_NAME] ,  char        treeFile[MAX_NAME],  char        timesFile[MAX_NAME], char        SNVgenotypesFile[MAX_NAME],char SNVhaplotypesFile[MAX_NAME],char trueHaplotypesFile[MAX_NAME], char fullGenotypesFile[MAX_NAME], char fullHaplotypesFile[MAX_NAME],char VCFfile[MAX_NAME],
                        char CATGfile[MAX_NAME] ,FILE            *fpSNVgenotypes, int doPrintTimes, int doSimulateData, FILE  *fpTimes, FILE *fpSNVhaplotypes, FILE *fpTrueHaplotypes,
                        FILE *fpFullGenotypes,FILE *fpFullHaplotypes, FILE *fpVCF, FILE *fpCATG, int doPrintSNVgenotypes, int  doPrintSNVhaplotypes, int doPrintTrueHaplotypes, int doPrintFullGenotypes, int doPrintFullHaplotypes, int doNGS,  int doPrintCATG, int numDataSets);


void PrintTrueFullHaplotypes (FILE *fp,TreeNode* nodes,TreeNode* treeRoot,int numNodes, int doPrintIUPAChaplotypes, int doPrintAncestors, int numSites, int numCells, int alphabet, int doUserTree , int ***data,   int doNGS,   char **cellNames, CellStr            *cell, int        HEALTHY_ROOT, int TUMOR_ROOT, char *cellnames[], int doUseObservedCellName);

void ReadParametersFromFastaFile(char *fileName, ProgramOptions *programOptions
                                // , int *numCells, int *numSites
                                 );

void ReadFastaFile(char *fileName, int** dataFromFile, char **cellNames, ProgramOptions *programOptions);

void PrintSNVGenotypes (FILE *fp,int numClones, Population **populations,TreeNode* nodes,TreeNode* treeRoot, int doPrintAncestors, int doNGS, int numCells, int numSNVs, int* SNVsites, int alphabet, int doUserTree, int ***data, char **cellNames, int TUMOR_ROOT, int HEALTHY_ROOT, CellStr* cell, int* SFS,int *numberDifferences, int TajimaD, double mutationRate);

static void   ReadParametersFromCommandLine (int argc, char **argv);
static void   ReadUntil (FILE *fv, char stopChar, char *what);

static void   PrintTitle ();
static void   PrintDate ();
static void   PrintUsage();
static void InitFilesPathsOptions( FilePaths *filePaths, ProgramOptions *programOptions)
{
    strcpy(filePaths->resultsDir, "Results");
    strcpy(filePaths->treeDir, "trees_dir");
    strcpy(filePaths->timesDir, "times_dir");
    strcpy(filePaths->SNVgenotypesDir, "snv_genotypes_dir");
    strcpy(filePaths->SNVhaplotypesDir, "snv_haplotypes_dir");
    strcpy(filePaths->trueHaplotypesDir, "true_haplotypes_dir");
    strcpy(filePaths->fullHaplotypesDir, "full_haplotypes_dir");
    strcpy(filePaths->MLhaplotypesDir, "ML_haplotypes_dir");
    strcpy(filePaths->fullGenotypesDir, "full_genotypes_dir");
    strcpy(filePaths->CATGdir, "catg_dir");
    strcpy(filePaths->VCFdir, "vcf_dir");
    
    strcpy(filePaths->SNVgenotypesFile, "snv_gen");
    strcpy(filePaths->SNVhaplotypesFile, "snv_hap");
    strcpy(filePaths->trueHaplotypesFile, "true_hap");
    strcpy(filePaths->fullHaplotypesFile, "full_hap");
    strcpy(filePaths->MLhaplotypesFile, "ML_hap");
    strcpy(filePaths->fullGenotypesFile, "full_gen");
    strcpy(filePaths->treeFile, "trees");
    strcpy(filePaths->timesFile, "times");
    if (strlen(filePaths->userTreeFile) == 0)
        strcpy(filePaths->userTreeFile, "usertree");
    if (strlen(filePaths->userGenomeFile) == 0)
        strcpy(filePaths->userGenomeFile, "usergenome");
    strcpy(filePaths->CATGfile, "catg");
    strcpy(filePaths->VCFfile, "vcf");
    strcpy(filePaths->logFile, "log");
    strcpy(filePaths->settingsFile, "log");
#ifdef MYDEBUG
    strcpy(filePaths.mutationsFile, "mutations");
#endif
    programOptions->doPrintSNVgenotypes=1;
    programOptions->doPrintSNVhaplotypes=1;
    programOptions->doPrintTrueHaplotypes=1;
    programOptions->doPrintFullGenotypes=1;
    programOptions->doPrintFullHaplotypes=1;
    programOptions->doNGS=1;
    programOptions->doPrintCATG=1;
    if (programOptions->doSimulateData == NO)
    {
        programOptions->doPrintSNVgenotypes = NO;
        programOptions->doPrintSNVhaplotypes = NO;
        programOptions->doPrintTrueHaplotypes = NO;
        programOptions->doPrintFullHaplotypes = NO;
        programOptions->doPrintMLhaplotypes = NO;
        programOptions->doPrintFullGenotypes = NO;
        programOptions->doPrintAncestors = NO;
        programOptions->doSimulateReadCounts = NO;
        programOptions->doPrintCATG = NO;
    }
}
#endif /* fileFuntions_h */
