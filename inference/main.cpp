#include <iostream>
#include <string>

#include "data_utils.hpp"
#include "Population.hpp"
#include "definitions.h"
using namespace std;

int main(int argc, char* argv[] )
{
    // take input
    // 1. input fasta file and phylip files
    // 2. input tre file if doFixedTree ==YES
    // 3. read the number of clones if numberClonesKnown=YES
    // 4. number of chains
    // 5. number of  chains iterations
    // 6. number of warmup iterations
    // output directory
    
    // default evolution model: JC
    
    FILE    *input_file;
    int i,j;
    
    char* input_path;
    if (argc <= 2)
        input_path = argv[1];
    else{
        fprintf (stderr, "\nERROR: No parameters specified (use command line or parameter file)");
    }
    

//    char *fileName; //="Houetal_HS.fasta";
//    char *fileNamePhylip ;//="Houetal_HS1.phylips";
//    fileName ="SimulatedDataJC.fasta";
//    fileNamePhylip ="SimulatedDataJC.phylip";
//
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
//    programOptions.numberClonesKnown=YES;
//    programOptions.populationSampleSizesKnown = YES;
//    mcmcOptions.slidingWindowSizeTotalEffectPopSize = 20000;
//    programOptions.doUseGenotypes = YES;
//    programOptions.doUseFixedTree =NO;
//
//    mcmcOptions.tuningParameter = 1;
//    mcmcOptions.thinning  = 1000;
//
//    programOptions.seqErrorRate=programOptions.sequencingError;
//    programOptions.dropoutRate=programOptions.ADOrate;
//
//    int *sampleSizes =(int *) malloc(programOptions.numClones* (long) sizeof(int));
//    if (!sampleSizes)
//    {
//        fprintf (stderr, "Could not allocate samplesSizes");
//        exit (1);
//    }
//
//    //2. initialize data structures
//
//    Population **populations = (Population**)malloc (sizeof(struct Population*)  * programOptions.numClones);
//    if (!populations)
//    {
//        fprintf (stderr, "Could not allocate populations (%lu bytes)\n", (programOptions.numClones)  * (long) sizeof(Population*));
//        exit (1);
//    }
//    InitListClones(populations, programOptions.numClones, programOptions.noisy, CloneNameBegin, CloneSampleSizeBegin, CloneBirthRateBegin,  CloneDeathRateBegin, ClonePopSizeBegin, programOptions.TotalNumSequences);
//    InitListClonesTimes(populations, programOptions.numClones,  &programOptions.doEstimateTimesOriginClones,  CloneTimeOriginInput  );
//    InitNumberNodes(&TotalBirthRate, &TotalDeathRate, &TotalN, populations, &programOptions);
//    ListClonesAccordingTimeToOrigin(populations, programOptions.numClones);
//
//    /* set file dirs and names */
//    InitFilesPathsOptions(&filePaths, &programOptions);
//
//    //3. do inference
//
//    sampleSizes[0]=5;
//    sampleSizes[1]=3;
//    sampleSizes[2]=4;
//    programOptions.numClones = 3;
//
//    programOptions.numNodes = 2 * programOptions.TotalNumSequences + programOptions.numClones+ 10;
//    programOptions.numCells = programOptions.TotalNumSequences;
//
//    if (programOptions.numberClonesKnown==YES)
//    {
//        mcmcOptions.totalEffectPopSizefrom = round(log(programOptions.numCells +2));
//        mcmcOptions.totalEffectPopSizefrom = round(log(100* programOptions.numCells +200 ));
//
//    }
//
//    //char *fileName ="Nietal_HS.fasta";
//    ReadParametersFromFastaFile(fileName,  &programOptions);
//    ObservedData = (int**) malloc (programOptions.numCells * sizeof(int*));
//    if (!ObservedData)
//    {
//        fprintf (stderr, "Could not allocate the ObservedData memory\n");
//        exit (-1);
//    }
//    for (i=0; i< programOptions.numCells; i++)
//    {
//        ObservedData[i] = (int *) calloc (programOptions.numSites +1, sizeof(int));
//        if (!ObservedData[i])
//        {
//            fprintf (stderr, "Could not allocate the observed data structure\n");
//            exit (-1);
//        }
//    }
//    ReadFastaFile(fileName, ObservedData,  ObservedCellNames, &programOptions);
//
//    for( i = 0 ; i < programOptions.numCells; i++)
//        fprintf (stderr, "observed cell name %s\n", ObservedCellNames[i]);
//    //fprintf (stderr, "observed data %d \n", *ObservedData[0]);
//
//    msa = pll_phylip_load(fileNamePhylip, PLL_FALSE);
//
//    int chainNumber=0;
//
//    mcmcOptions.numChains=2;
//
//    Chain * chains = (Chain *) malloc(mcmcOptions.numChains* ( sizeof(Chain) +  4 * sizeof(TreeNode) + 4 * sizeof(TreeNode) + programOptions.numClones *sizeof(Population)));
//    if (!chains)
//    {
//        fprintf (stderr, "Could not allocate chains");
//        exit (1);
//    }
//
//    InitializeChains(&chains, &programOptions, &mcmcOptions, sampleSizes, &(programOptions.seed), ObservedCellNames, msa);
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
//    pll_msa_destroy(msa);
    
  return 0;

}
