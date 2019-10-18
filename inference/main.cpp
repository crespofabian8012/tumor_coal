#include <iostream>
#include <string>
#include <stdarg.h>
#include <search.h>
#include <time.h>
#include <map>
#include <unordered_map>
#include <string>
#include "data_utils.hpp"
#include "population.hpp"
#include "mcmc_chain.hpp"
#include "definitions.hpp"
#include "utils.hpp"
#include "constants.hpp"

extern "C"
{
#include <libpll/pll_msa.h> //for pllmod_msa_empirical_frequencies
#include <libpll/pll_tree.h>
#include <libpll/pllmod_util.h>
#include <libpll/pllmod_common.h>
#include <libpll/pllmod_algorithm.h>
}
#include <stdarg.h>
#include <search.h>
#include <time.h>
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
    programOptions.numClones=3;
    
    std::map<std::string,int> tipsLabelling;
    std::map<pll_unode_t, Population> tipsAssign;
    
    mcmcOptions.tuningParameter = 1;
    mcmcOptions.thinning  = 1000;
    
    programOptions.seqErrorRate=programOptions.sequencingError;
    programOptions.dropoutRate=programOptions.ADOrate;
    
    vector<int> sampleSizes(programOptions.numClones);
    
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
    vector<vector<int> > ObservedData;
    
    char *ObservedCellNames[programOptions.numCells];
    ReadFastaFile(fileNameFasta, ObservedData,  ObservedCellNames, programOptions);
    //
    //    for( i = 0 ; i < programOptions.numCells; i++)
    //        fprintf (stderr, "observed cell name %s\n", ObservedCellNames[i]);
    //    //fprintf (stderr, "observed data %d \n", *ObservedData[0]);
    //
    pll_msa_t *msa = pll_phylip_load(fileNamePhylip, PLL_FALSE);
    if (!msa)
        fprintf (stderr, "Error reading phylip file \n");
    //    int chainNumber=0;
    //
    mcmcOptions.numChains=2;
    
    vector<Chain*> chains;
    

    pll_rtree_t * initialRootedTree = pll_rtree_parse_newick(treefileName);
    pll_utree_t * initialUnrootedTree = pll_utree_parse_newick(treefileName);
   

    if (!initialRootedTree)
    {
        fprintf (stderr, "Error reading newick representation of initial rooted tree \n");
        exit(1);
    }
    if (initialRootedTree->inner_count > 1)
        initialUnrootedTree = pll_rtree_unroot(initialRootedTree);
    else{
        fprintf (stderr, "The rooted tree has only one inner node \n");
        exit(1);
    }
    if (!initialRootedTree)
    {
        fprintf (stderr, "Error unrooting the tree \n");
        exit(1);
    }
     //pllmod_utree_set_length_recursive(initialUnrootedTree, BRLEN_MIN, 1);
    pll_unode_t *root = initialUnrootedTree->nodes[initialUnrootedTree->tip_count + initialUnrootedTree->inner_count - 1];
    pll_utree_reset_template_indices(root, initialUnrootedTree->tip_count);
       char * newick = pll_utree_export_newick(initialUnrootedTree->vroot,NULL);
       char * rootedNewick = pll_utree_export_newick_rooted(initialUnrootedTree->vroot, 6.13);
      printf("%s\n", newick);
      printf("%s\n", rootedNewick);
    
    Chain *chain;
    double loglh =  chain->LogConditionalLikelihoodSequences( msa,  newick, programOptions, 0, 0);
    
    free(newick);
    
      Chain::initializeChains(chains, programOptions, mcmcOptions, sampleSizes, &programOptions.seed, ObservedCellNames, msa,  initialUnrootedTree, initialRootedTree);

    Chain *currentChain;
    int currentIteration;
    int sampleEvery = mcmcOptions.thinning;


     for(int chainNumber=0; chainNumber< mcmcOptions.numChains;chainNumber++)
   {
       for (currentIteration = 0; currentIteration < mcmcOptions.Niterations; currentIteration++)
     {
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
      }
  }
   pll_msa_destroy(msa);

  return 0;

}

