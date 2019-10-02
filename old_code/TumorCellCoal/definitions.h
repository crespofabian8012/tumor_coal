
#include "pll_optimize.h"
#include "pll_tree.h"
#include "pllmod_algorithm.h"
#include "pllmod_common.h"

#ifndef macros_h
#define macros_h

/* Macros */
#define PROGRAM_NAME    "TumorCoal"
#define VERSION_NUMBER    "2.0"
#define NO          0
#define YES         1
#define INCREMENT_NODES   500

#define	infinity                999999
#define	MATERNAL                0
#define	PATERNAL                1
#define BINARY                  0
#define	DNA                     1
//#define    MAX_NAME             120
#define    MAX_NAME             120
#define    MAX_LINE                3500
#define    MAX_GENOME              1000000

#define ISMhap					0
#define Mk						1
#define finiteDNA				2

#define A						0
#define C						1
#define G						2
#define T						3
#define N                       4
#define R                       5
#define ADO						9
#define DELETION				7
#define M                       6
#define S                       8
#define Y                       10
#define K                       11


//#define genotypes
#define AA                      0
#define CC                      1
#define GG                      2
#define TT                      3
#define AC                      4
#define CA                      4
#define AG                      5
#define GA                      5
#define AT                      6
#define TA                      6
#define A_                      7
#define _A                      7
#define CG                      8
#define GC                      8
#define CT                      9
#define TC                      9
#define C_                      10
#define _C                      10
#define GT                      11
#define TG                      11
#define G_                      12
#define _G                      12
#define T_                      13
#define _T                      13
#define __                      14
//A/C => M
//A/G => R
//A/T => W
//A/_ => a
//
//C/A => M
//C/C => C
//C/G => S
//C/T => Y
//C/_ => c
//
//G/A => R
//G/C => S
//G/G => G
//G/T => K
//G/_ => g
//
//T/A => W
//T/C => Y
//T/G => K
//T/T => T
//T/_ => t
//
//_/A => a
//_/C => c
//_/G => g
//_/T => t
//_/_ => -


//#define M                       6 //A/C or C/A
//#define W                       8 //A/T or T/A
//#define S                        //C/G or G/C
//#define Y                       14 //C/T or T/C
//#define K                       15 //G/T or T/G
//#define a                       16 //A/_ or _/A
//#define c                       17 //C/_ or _/C
//#define g                       18 //G/_ or _/G
//#define t                       19 //T/_ or _/T

#define CG_AT					0	// C>A or G>T
#define CG_GC					1	// C>G or G>C
#define CG_TA					2	// C>T or G>A
#define TA_AT					3	// T>A or A>T
#define TA_CG					4	// T>C or A>G
#define TA_GC					5	// T>G or A>C

//#define M                       10 //A/C or C/A
////#define R                       11 //A/G or G/A
//#define W                       12 //A/T or T/A
//#define S                       13 //C/G or G/C
//#define Y                       14 //C/T or T/C
//#define K                       15 //G/T or T/G
//#define a                       16 //A/_ or _/A
//#define c                       17 //C/_ or _/C
//#define g                       18 //G/_ or _/G
//#define t                       19 //T/_ or _/T
//#define _                       20 //_/_


#define ISMhap					0
#define Mk2						1
#define finiteDNA				2

#define NUM_SIGNATURES			30

#define trinuc(i,j,k)   (i*16 + j*4 + k)
#define trimut(i,j,k)   (i*24 + j*4 + k)


#define CHECK_MUT_DISTRIBUTION
#undef CHECK_MUT_DISTRIBUTION

#define MYDEBUG
#undef MYDEBUG

#define LOAD_INTERNAL_SIGNATURES

#define PRINT_TRIMUTATIONS
#undef PRINT_TRIMUTATIONS

#define PRINT_TRINUC_GENOME
#undef PRINT_TRINUC_GENOME

#define PRINT_TRIMUTCOUNTER
#undef PRINT_TRIMUTCOUNTER


#define GT_MODEL "GTGTR4"
#define JC_MODEL "GTJC"
#define RATE_CATS 1
#define BRLEN_MIN 1e-6
#define BRLEN_MAX 1e+2

/* ascertainment bias correction */
#define PLL_ATTRIB_AB_LEWIS        (1 << 5)
#define PLL_ATTRIB_AB_FELSENSTEIN  (2 << 5)
#define PLL_ATTRIB_AB_STAMATAKIS   (3 << 5)
#define PLL_ATTRIB_AB_MASK         (7 << 5)
#define PLL_ATTRIB_AB_FLAG         (1 << 8)

#define PLL_ATTRIB_RATE_SCALERS    (1 << 9)

#define HOMO(state)   (state<4)
#define HETERO(state) (state>3)

/*                                        AA CC GG TT AC AG AT CG CT GT            */
static const char mut_dist[10][10] = {
                                        {  0, 2, 2, 2, 1, 1, 1, 2, 2, 2 },   /* AA */
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

/* Structures */
typedef struct node
{
  pll_unode_t     *nodeLeft;
  pll_unode_t     *nodeRight;
  pll_unode_t     *nodeBack;
  pll_tree_edge_t  *edgeLeft;
  pll_tree_edge_t  *edgeRight;
  pll_tree_edge_t  *edgeBack;
    
  struct node *left, *right, *anc1, *outgroup;
  int         index, label, isOutgroup;
  double      length, time,lengthModelUnits, timePUnits;
  int         class;
    int         indexOldClone, indexCurrentClone,orderCurrentClone;
    //indexCoalClone;
  double      effectPopSize;
  char        cellName[MAX_NAME];
  char        observedCellName[MAX_NAME];
  int*        maternalSequence;
  int*        paternalSequence;
  int*        numbersMutationsUnderSubtreePerSite;
  int*        numbersMaternalMutationsPerSite;
  int*        numbersPaternalMutationsPerSite;
  int         isLeaf;
    
}
TreeNode;

typedef TreeNode* pTreeNode;

typedef struct Population
{
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

} Population;
typedef struct
    {
    int		isSNV, isSNP, isVariant;
	int		numMutations, numMutationsMaternal, numMutationsPaternal;
    int     numDeletions, numDeletionsMaternal, numDeletionsPaternal;
    int     numCNLOH, numCNLOHmaternal, numCNLOHpaternal;
	int		hasADO;
	int		hasGenotypeError;
	int		countA, countC, countG, countT, countACGT, countCellswithData, countDropped;
	int		referenceAllele;
	int		*alternateAlleles;
	int		numAltAlleles;
    double  rateMultiplier;
    int    numberDiffReference;
	}
    SiteStr;

typedef struct
    {
	int 	tempLength;
	int 	numAvailablePositions;
	int 	*position;
	}
    TriNucStr;

typedef struct csite
    {
	int		numReads;
	int		*readCount;
	int		trueMaternalAllele, truePaternalAllele;
	int		maternalAllele, paternalAllele;
	int		thereIsMaternalAllele, thereIsPaternalAllele;
	double  maternalSiteAmplificationError, paternalSiteAmplificationError;
	double	**genLike;
	double	**scaledGenLike;
	int		MLmatAllele;
	int		MLpatAllele;
	int		numReadsDoublet;
	int		*readCountDoublet;
	int		maternalAlleleDoublet, paternalAlleleDoublet;
	double	**genLikeDoublet;
	double	**scaledGenLikeDoublet;
	int		MLmatAlleleDoublet;
	int		MLpatAlleleDoublet;
	}
    CellSiteStr;

typedef struct
    {
	struct csite	*site;
	int				hasDoublet;
	}
    CellStr;

typedef struct {
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
}
MCMCoptions;

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
}
FilePaths;

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

typedef struct {
    int doUseGenotypes;
    int doUseFixedTree;
    int doPrintTrees;
    int doPrintSNVgenotypes;
    int doPrintSNVhaplotypes;
    int doPrintTrueHaplotypes;
    int doPrintFullGenotypes;
    int doPrintFullHaplotypes;
    int doNGS;
    int doPrintTree;
    int doPrintCATG;
    int doPrintAncestors;
    int doPrintMLhaplotypes;
    int doSimulateReadCounts;
    int doPrintTimes;
    int doSimulateData;
    int doPrintIUPAChaplotypes;
    int doGeneticSignatures;
    int doJC;
    int doHKY;
    int doGTR;
    int doGTnR;
    int thereIsMij;
    int doSimulateFixedNumMutations;
    int    thereisOutgroup;
    int outgroupSelection;
    int doEstimateTimesOriginClones;
    int doAcceptTimes;
    int doPrintSeparateReplicates;
    double outgroupBranchLength_RootSample;
    double outgroupBranchLength_Root1Root2;
    int doUseObservedCellNames;
    double ADOrate;
    double Eij[4];
    double Mij[4];
    double SNPrate;
    double allelicImbalance;
    double alphaBranches;
    double alphaCoverage;
    double alphaSites;
    int alphabet;
    int altModel;
    int coverage;
    int doDemographics;
    int doExponential;
    int doUserGenome;
    int doUserTree;
    double doubletRate;
    int equalBaseFreq;
    double freq;
    double genotypingError;
    double growthRate;
    double haploidCoverageReduction;
    double healthyTipBranchLength;
    double meanAmplificationError;
    double mutationRate;
    int noisy;
    double nonISMRelMutRate;
    int numCells;
    int numDataSets;
    int numNodes;
    int Nscaling;
    int numPeriods;
    int numSites;
    int numUserSignatures;
    int ploidy;
    long int seed;
    long int userSeed;
    int simulateOnlyTwoTemplates;
    double varAmplificationError;
    double sequencingError;
    int rateVarCoverage;
    double transformingBranchLength;
    int rateVarAmongSites;
    int rateVarAmongLineages;
    double propAltModelSites;
    double titv;
    double seqErrorRate;
    double dropoutRate;
    int thereIsEij;
    int numClones;
    double altModelMutationRate;
    int numFixedMutations;
    int TotalNumSequences;
    int doSimulateFromPriors;
    int populationSampleSizesKnown;
    int assignationKnown;
    int numberClonesKnown;
    int MutationAssignNum;
}ProgramOptions;

typedef struct {
    int chainNumber;
    int numClones;
    int oldnumClones;
    int gammaParam;
    int currentNumberIerations;
    
  
    TreeNode *root;
    TreeNode *nodes;
    TreeNode  **treeTips;
    
    TreeNode *oldroot;
    TreeNode *oldnodes;
    TreeNode  **oldtreeTips;
    
    int numNodes;
    int oldnumNodes;
    
    double mutationRate;
    double seqErrorRate;
    double dropoutRate;
    double oldmutationRate;
    Population **populations;
    
    double *proportionsVector;
    double *oldproportionsVector;
    int totalPopSize;
    int oldtotalPopSize;
    double lambda ;
    double currentlogConditionalLikelihoodTree;
    double currentlogConditionalLikelihoodSequences;

    TreeNode    **treeRootInit;
    TreeNode    **oldtreeRootInit;
}Chain;
//extern double   Qij[16], mr;
#endif /* macros_h */
