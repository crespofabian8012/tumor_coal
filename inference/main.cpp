#include <iostream>
#include <string>
#include <stdarg.h>
#include <search.h>
#include <time.h>

#include "data_utils.hpp"
#include "population.hpp"
#include "Chain.hpp"
#include "definitions.hpp"
#include "utils.hpp"
#include "constants.hpp"

extern "C"
{
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
#include "libpll/pll_optimize.h"
#include "libpll/pllmod_algorithm.h"
#include "libpll/pllmod_common.h"
}

using namespace std;

int main(int argc, char* argv[] )
{
    // default evolution model: JC

    FILE *input_file;
    int i;

    char* input_path;

//    char *fileName; //="Houetal_HS.fasta";
//    char *fileNamePhylip ;//="Houetal_HS1.phylips";
    char *fileNameFasta = argv[2];
    char *fileNamePhylip = argv[3];
    char *treefileName = argv[4];

    ProgramOptions programOptions;
    Files files;
    FilePaths filePaths;
    MCMCoptions mcmcOptions;
    // 1. call function to parse the input file
//
    if (argc <= 2)
    {
        if ((input_file = freopen("input_path", "r", stdin)) != NULL)
       {
           // ReadMCMCParametersFromFile();
       }
       else
        {
            fprintf (stderr, "\nERROR: No parameters specified (use command line or parameter file)");
        }
    }
    programOptions.numberClonesKnown=YES;
    programOptions.populationSampleSizesKnown = YES;
    mcmcOptions.slidingWindowSizeTotalEffectPopSize = 20000;
    programOptions.doUseGenotypes = YES;
    programOptions.doUseFixedTree =NO;
    programOptions.seed = 1248697;
    programOptions.numClones = 3;


    mcmcOptions.tuningParameter = 1;
    mcmcOptions.thinning  = 1000;

    programOptions.seqErrorRate=programOptions.sequencingError;
    programOptions.dropoutRate=programOptions.ADOrate;

    vector<int> sampleSizes(programOptions.numClones);

    pllmod_treeinfo_t * treeinfo = pllmod_treeinfo_create(0, 0, 1, PLLMOD_COMMON_BRLEN_LINKED);

    //2. initialize data structures
//    /* set file dirs and names */
    InitFilesPathsOptions(filePaths, programOptions);
//
//    //3. do inference
//
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
//    //char *fileName ="Nietal_HS.fasta";
    ReadParametersFromFastaFile(fileNameFasta,  programOptions);
    int** ObservedData = (int**) malloc (programOptions.numCells * sizeof(int*));
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
  char *ObservedCellNames[programOptions.numCells];
   ReadFastaFile(fileNameFasta, ObservedData,  ObservedCellNames, programOptions);
//
//    for( i = 0 ; i < programOptions.numCells; i++)
//        fprintf (stderr, "observed cell name %s\n", ObservedCellNames[i]);
//    //fprintf (stderr, "observed data %d \n", *ObservedData[0]);
//
   pll_msa_t *msa = pll_phylip_load(fileNamePhylip, PLL_FALSE);
//
//    int chainNumber=0;
//
    mcmcOptions.numChains=2;

//    vector<Chain*> chains;
//
////    Chain *chains = (Chain *) malloc(mcmcOptions.numChains* ( sizeof(chain) +  4 * sizeof(pll_unode_t) + 4 * sizeof(TreeNode) + programOptions.numClones *sizeof(Population)));
////    if (!chains)
////    {
////        fprintf (stderr, "Could not allocate chains");
////        exit (1);
////    }
//    pll_utree_t * initialTree = pll_utree_parse_newick(treefileName);
//    if (!initialTree)
//       fprintf (stderr, "Error reading newick representation of initial tree \n");
//  //  pllmod_utree_set_length_recursive(initialTree, BRLEN_MIN, 1);
//
//    Chain::InitializeChains(chains, programOptions, mcmcOptions, sampleSizes, &programOptions.seed, ObservedCellNames, msa,  initialTree);
//    //////////////////////////
//    FILE    *outputShell;
//    char script[80];
//    char buf[1000];
//
//    unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
//    pll_partition_t * partition;
//
//    /* compute node count information */
//    tip_nodes_count = initialTree->tip_count;
//    inner_nodes_count = initialTree->inner_count;
//    nodes_count = inner_nodes_count + tip_nodes_count;
//    branch_count = initialTree->edge_count;


    //pllmod_treeinfo_t * treeinfo = pllmod_treeinfo_create(initialTree->vroot,
//                                                          tip_nodes_count,
//                                                          1, PLLMOD_COMMON_BRLEN_LINKED);


//    pll_unode_t ** tipnodes = initialTree->nodes;
//
//    /* create a libc hash table of size tip_nodes_count */
//    hcreate(tip_nodes_count);
//
//    /* populate a libc hash table with tree tip labels */
//    unsigned int * data = (unsigned int *)malloc(tip_nodes_count *
//                                                 sizeof(unsigned int));
//    char *label;
//    for (i = 0; i < tip_nodes_count; ++i)
//    {
//        data[i] = tipnodes[i]->clv_index;
//        ENTRY entry;
//        label= tipnodes[i]->label;
//        entry.key = tipnodes[i]->label;
//        entry.data = (void *)(data+i);
//        hsearch(entry, ENTER);
//    }
//
//    if (msa !=NULL)
//        printf("Original sequence (alignment) length : %d\n", msa->length);
//    else{
//        fprintf (stderr, "\nERROR: The multiple sequence alignment is empty\n\n");
//        PrintUsage();
//        return 0;
//    }
//    //pllmod_subst_model_t * model = pllmod_util_model_info_genotype(GT_MODEL);
//
////    partition = pll_partition_create(tip_nodes_count,
////                                     inner_nodes_count,
////                                     model->states,
////                                     (unsigned int)(msa->length),
////                                     1,
////                                     branch_count,
////                                     RATE_CATS,
////                                     inner_nodes_count,
////                                     PLL_ATTRIB_ARCH_AVX);
//    ///////////////
//    Chain *currentChain;
//    int currentIteration;
//    int sampleEvery = mcmcOptions.thinning;
//
//
//     for(int chainNumber=0; chainNumber< mcmcOptions.numChains;chainNumber++)
//   {
//        for (currentIteration = 0; currentIteration < mcmcOptions.Niterations; currentIteration++)
//       {
////            runChain(&(chains[chainNumber]),    &mcmcOptions,  &(programOptions.seed),  &filePaths, &files, &programOptions,
////                     varTimeGMRCA, ObservedCellNames, msa, sampleSizes);
////
//           if (currentIteration % sampleEvery == 0 )
//           {
//
////                PrintTrees(currentIteration, &(chains[chainNumber].root), files.fpTrees, programOptions.mutationRate, programOptions.doUseObservedCellNames);
////                PrintTrees2(currentIteration, &(chains[chainNumber].root), files.fpTrees2, programOptions.mutationRate,  ObservedCellNames, programOptions.doUseObservedCellNames);
////                PrintTimes(currentIteration, files.fpTimes, programOptions.mutationRate, chains[chainNumber].nodes, programOptions.thereisOutgroup);
////                PrintTimes2(currentIteration, files.fpTimes2, programOptions.mutationRate, chains[chainNumber].nodes, programOptions.thereisOutgroup);
//           }
//        }
//  }
//   pll_msa_destroy(msa);

  return 0;

}

