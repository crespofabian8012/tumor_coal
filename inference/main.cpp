#include <iostream>
#include <string>
#include <stdarg.h>
#include <search.h>
#include <time.h>
#include <map>
#include <unordered_map>
#include <string>
#include <omp.h>
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
#include <libpll/pll.h>
}
#include <stdarg.h>
#include <search.h>
#include <time.h>
using namespace std;

int main(int argc, char* argv[] )
{

    FILE *input_file;
    char* input_path;
    //    char *fileName; //="Houetal_HS.fasta";
    //    char *fileNamePhylip ;//="Houetal_HS1.phylips";
    char *fileNameFasta ;//= argv[2];
    char *fileNamePhylip ;//= argv[3];
    char *treefileName; //= argv[4];
    
    ProgramOptions programOptions;
    Files files;
    FilePaths filePaths;
    MCMCoptions mcmcOptions;
    vector<double> varTimeGMRCA;
    if (argc <= 2)
        input_path = argv[1];
    else{
        fprintf (stderr, "\nERROR: No parameters specified (use command  parameter file)");
        PrintUsage();
    }
    
    // 1. call function to parse the input file
    if (argc <= 2)
    {
        if ((input_file = freopen(input_path, "r", stdin)) != NULL)
        {
            ReadMCMCParametersFromFile(programOptions, filePaths, mcmcOptions);
        }
        else
        {
            fprintf (stderr, "\nERROR: No parameters specified (use command line or parameter file)");
            exit(-1);
        }
    }
    programOptions.numberClonesKnown=YES;
    programOptions.populationSampleSizesKnown = YES;
    mcmcOptions.slidingWindowSizeTotalEffectPopSize = 20000;
    mcmcOptions.slidingWindowSizeGrowtRate=0.01;
    programOptions.doUseGenotypes = YES;
    programOptions.doUseFixedTree =NO;
    programOptions.seed = 1248697;
    
    programOptions.mutationRate= 9.1e-8;
    programOptions.doUsefixedMutationRate=0;
    programOptions.doSimulateData=NO;
    
    std::map<std::string,int> tipsLabelling;
    std::map<pll_unode_t, Population> tipsAssign;
    
    mcmcOptions.tuningParameter = 1;
    
    mcmcOptions.totalEffectPopSizefrom = 7;
    mcmcOptions.totalEffectPopSizeto = 11;
    mcmcOptions.MutRatefrom = -20;
    mcmcOptions.MutRateto = -13;
    mcmcOptions.Deltafrom = -4;
    mcmcOptions.Deltato = 1;
    mcmcOptions.fixedLambda=1;
    
    mcmcOptions.GrowthRatefrom = -20;
    mcmcOptions.GrowthRateto = log(0.2);
    //mcmcOptions.GrowthRateto = log(log(2));//maximum growth rate of log(2) = 0.6931472
    
    programOptions.seqErrorRate=programOptions.sequencingError=0;
    programOptions.dropoutRate=programOptions.ADOrate=0;
    
    vector<int> sampleSizes(programOptions.numClones);
    
    //2. initialize data structures
    //    /* set file dirs and names */
    InitFilesPathsOptions(filePaths, programOptions);
    //
    //    //3. do inference
    
    // programOptions.numClones = 3;
    
    //    //char *fileName ="Nietal_HS.fasta";
    
    fileNameFasta = filePaths.inputGenotypeFileFasta;
    ReadParametersFromFastaFile(fileNameFasta,  programOptions);
    vector<vector<int> > ObservedData;
    
    char *ObservedCellNames[programOptions.numCells];
    ReadFastaFile(fileNameFasta, ObservedData,  ObservedCellNames, programOptions);
    
    programOptions.numNodes = 2 * programOptions.TotalNumSequences + programOptions.numClones+ 10;
    programOptions.numCells = programOptions.TotalNumSequences;
    programOptions.doUseFixedTree = YES;
    
    if (programOptions.numberClonesKnown==YES)
    {
       // mcmcOptions.totalEffectPopSizefrom = round(log(10*programOptions.numCells));
       // mcmcOptions.totalEffectPopSizeto = round(log(200*programOptions.numCells ));
    }
    
    fileNamePhylip =filePaths.inputGenotypeFilePhylip;
    pll_msa_t *msa = pll_phylip_load(fileNamePhylip, PLL_FALSE);
    if (!msa)
        fprintf (stderr, "Error reading phylip file \n");
    
    vector<Chain*> chains(mcmcOptions.numChains);
    
    treefileName = filePaths.inputTreeFile;

    string healthyTipLabel = "healthycell";
    programOptions.healthyTipLabel ="healthycell";
    
    int currentIteration;
    //mcmcOptions.thinning  = 5000;
    //mcmcOptions.maxNumberProposalAttempts=10;
    int sampleEvery = mcmcOptions.thinning;
    mcmcOptions.paramMultiplierMoveTheta = 3;
    mcmcOptions.paramMultiplierEffectPopSize = 2;
    //mcmcOptions.Niterations = 10000000;
    mcmcOptions.numberWarmUpIterations = mcmcOptions.Niterations / 2.0;
    
    mcmcOptions.doInferenceWithoutData =0;
    //candidate for parallelizing
#pragma parallel for default(shared) private(chainNumber, currentIteration) firstprivate(files)
    for(int chainNumber=0; chainNumber< mcmcOptions.numChains;chainNumber++)
    {
        
        chains.at(chainNumber) = Chain::initializeChain( programOptions, mcmcOptions, sampleSizes, &programOptions.seed, ObservedCellNames, msa,  treefileName, healthyTipLabel);
        chains.at(chainNumber)->chainNumber =chainNumber;
        
        chains.at(chainNumber)->PrepareFiles(filePaths, programOptions, chains.at(chainNumber)->files, chainNumber);
        chains.at(chainNumber)->writeHeaderOutputChain(filePaths, programOptions,
                                                       chains.at(chainNumber)->files );
      
        for (currentIteration = 0; currentIteration < mcmcOptions.Niterations; currentIteration++)
        {
            fprintf (stderr, "\n Chain #%d, Iteration %d \n", chains.at(chainNumber)->chainNumber ,currentIteration );
            chains.at(chainNumber)->currentNumberIerations =currentIteration;
            chains.at(chainNumber)->runChain(mcmcOptions,  &(programOptions.seed),  filePaths, chains.at(chainNumber)->files, programOptions,ObservedCellNames, msa, sampleSizes, currentIteration );
            
            if (currentIteration % sampleEvery == 0 && currentIteration >= mcmcOptions.numberWarmUpIterations)
            {
                chains.at(chainNumber)->writeMCMCState(  currentIteration, filePaths, programOptions,chains.at(chainNumber)->files, mcmcOptions);
                //in the future i will write trees to a nexus file
                ////PrintTrees(currentIteration, &(chains[chainNumber].root), 
            }
        }
        chains.at(chainNumber)->currentNumberIerations =currentIteration;
        if (currentIteration >= mcmcOptions.Niterations -1)//last iteration
            {
                chains.at(chainNumber)->writeMCMCState(  currentIteration, filePaths, programOptions,chains.at(chainNumber)->files, mcmcOptions);
                
                fprintf (stderr, "\n Number accepted moves %d, number of rejected moves %d \n", chains.at(chainNumber)->totalAccepted,chains.at(chainNumber)->totalRejected );
            
            fclose(chains.at(chainNumber)->files.fplog);
            //chains.at(chainNumber)->closeFiles(filePaths,programOptions, files, chainNumber);
        }
    }
    //close files
    //fclose(files.fplog);
    pll_msa_destroy(msa);
    return 0;
}

