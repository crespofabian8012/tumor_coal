#include <iostream>
#include <string>

#include <libpll/pll_msa.h>

#include "data_utils.hpp"
#include "population.hpp"
#include "chain.hpp"
#include "definitions.hpp"
#include "utils.hpp"
using namespace std;


int main(int argc, char* argv[] )
{
    // default evolution model: JC
    
    FILE *input_file;
    int i,j;
    
    char* input_path;
    if (argc <= 2)
        input_path = argv[1];
    else{
        fprintf (stderr, "\nERROR: No parameters specified (use command line or parameter file)");
    }


//    char *fileName; //="Houetal_HS.fasta";
//    char *fileNamePhylip ;//="Houetal_HS1.phylips";
    char *fileName = argv[2];
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
    ReadParametersFromFastaFile(fileName,  programOptions);
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
   ReadFastaFile(fileName, ObservedData,  ObservedCellNames, programOptions);
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

    chain *chains = (chain *) malloc(mcmcOptions.numChains* ( sizeof(chain) +  4 * sizeof(pll_unode_t) + 4 * sizeof(TreeNode) + programOptions.numClones *sizeof(Population)));
    if (!chains)
    {
        fprintf (stderr, "Could not allocate chains");
        exit (1);
    }

 //  InitializeChains(&chains, &programOptions, &mcmcOptions, sampleSizes, &(programOptions.seed), ObservedCellNames, msa);
//
//    Chain currentChain;
//    int currentIteration;
//    int sampleEvery = mcmcOptions.thinning;
//    for(chainNumber=0; chainNumber< mcmcOptions.numChains;chainNumber++)
//    {
//        for (currentIteration = 0; currentIteration < mcmcOptions.Niterations; currentIteration++)
//        {
//            runChain(&(chains[chainNumber]),    &mcmcOptions,  &(programOptions.seed),  &filePaths, &files, &programOptions,
//                     varTimeGMRCA, ObservedCellNames, msa, sampleSizes);
//
//            if (currentIteration % sampleEvery == 0 )
//            {
//
//                PrintTrees(currentIteration, &(chains[chainNumber].root), files.fpTrees, programOptions.mutationRate, programOptions.doUseObservedCellNames);
//                PrintTrees2(currentIteration, &(chains[chainNumber].root), files.fpTrees2, programOptions.mutationRate,  ObservedCellNames, programOptions.doUseObservedCellNames);
//
//                PrintTimes(currentIteration, files.fpTimes, programOptions.mutationRate, chains[chainNumber].nodes, programOptions.thereisOutgroup);
//                PrintTimes2(currentIteration, files.fpTimes2, programOptions.mutationRate, chains[chainNumber].nodes, programOptions.thereisOutgroup);
//
//            }
//
//        }
//
//    }
   pll_msa_destroy(msa);
    
  return 0;

}
