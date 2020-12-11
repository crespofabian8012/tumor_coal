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
 * MCMC inference of scaled growth rates
 */

#include <iostream>
#include <string>
#include <stdarg.h>
#include <search.h>
#include <time.h>
#include <map>
#include <unordered_map>
#include <string>
#include <omp.h>
#include "random.h"
#include "data_utils.hpp"
#include "population.hpp"
#include "mcmc_chain.hpp"
#include "chain_manager.hpp"
#include "definitions.hpp"
#include "utils.hpp"
#include "constants.hpp"
#include "output_functions.hpp"

#include <gsl/gsl_randist.h>

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
    vector<StructuredCoalescentTree *> structuredCoalTrees;
    vector<pll_rtree_t *> trueTrees;
    vector<long double> trueThetas;
    vector<vector<long double>> trueDeltaTs;
    vector<vector<long double>> trueTs;

    char *fileNameFasta ;//= argv[2];
    char *fileNamePhylip ;//= argv[3];
    char *treefileName; //= argv[4];
    
    ProgramOptions programOptions;
    //Files files;
    FilePaths filePaths;
    MCMCoptions mcmcOptions;
    
    mcmcOptions.noData = false;
    mcmcOptions.useGSLRandomGenerator = true;
    mcmcOptions.splitThetaDeltaTmoves=true;
    
    
    //test
    long double result = Chain::LogDensityCoalTimes({ 1.0, 1.5}, {1.0}, {1.5}, 1.5, 1.0, 2);
    // result must be -1.032054
    
    long double result2= log(Population::DensityTimeSTD(1.5, 1.0, 0.0));
    //result2 must be -1.282252
    
    
    printProgramHeader();
    setDefaultOptions(programOptions, mcmcOptions);
    if (argc == 1 )
    {
        
        std::cout << "\nERROR: No parameters specified (use command  parameter file)"<< endl;
        Output::PrintUsage();
    }
    if (argc <= 2 )//only a path to a MCMC configuration file
        input_path = argv[1];
    else
    {//many options from console
        
        if (!mcmcOptions.noData)
        {
             std::cout << "\nERROR: No parameters specified (use command  parameter file)"<<endl;
            Output::PrintUsage();
        }
    }
    
    // 1. call function to parse the MCMC configuration file
    if (argc <= 2)//only a path to a MCMC configuration file
    {
        if ((input_file = freopen(input_path, "r", stdin)) != NULL)
        {
            ReadMCMCParametersFromFile(programOptions, filePaths, mcmcOptions);
        }
        else
        {
            if (!mcmcOptions.noData)
            {
                
                 std::cout << "\nERROR: No parameters specified (use command line or parameter file)"<<endl;
                exit(-1);
            }
        }
    }
    // std::map<std::string,int> tipsLabelling;
    //std::map<pll_unode_t, Population> tipsAssign;
    
    std::vector<int> sampleSizes;
    
    //2. initialize data structures
    //    /* set file dirs and names */
    InitFilesPathsOptions(filePaths, programOptions);
    //
    //    //3. do inference

    pll_msa_t *msa;
    
    //    if (mcmcOptions.useSequencesLikelihood ==1)
    //   {
    fileNameFasta = filePaths.inputGenotypeFileFasta;
    Utils::ReadParametersFromFastaFile(fileNameFasta,  programOptions.numCells, programOptions.TotalNumSequences, programOptions.numSites);
    std::vector<std::vector<int> > ObservedData;
    
    char *ObservedCellNames[programOptions.numCells];
    Utils::ReadFastaFile(fileNameFasta, ObservedData,  ObservedCellNames, programOptions);
    programOptions.numNodes = 2 * programOptions.TotalNumSequences + programOptions.numClones+ 10;
    programOptions.numCells = programOptions.TotalNumSequences;
    
    programOptions.TotalTumorSequences=programOptions.TotalNumSequences-1;
    
    if (programOptions.numClones==1)
        sampleSizes.push_back(programOptions.numCells-1); //minus the healthycell
    
    fileNamePhylip =filePaths.inputGenotypeFilePhylip;
    msa = pll_phylip_load(fileNamePhylip, PLL_FALSE);
    if (!msa)
         std::cout << "Error reading phylip file \n"<< endl;
    // }
    treefileName = filePaths.inputTreeFile;
    if ((treefileName != NULL) && (treefileName[0] == '\0'))
    {
        std::cout << "\nERROR: No tree file  specified (use command line or parameter file)" <<endl;
        exit(-1);
    }
    std::string healthyTipLabel = "healthycell";
    programOptions.healthyTipLabel ="healthycell";

    //int sampleEvery = mcmcOptions.thinning;
    vector<Chain*> chains(mcmcOptions.numChains);
    vector<gsl_rng * > randomGenerators(mcmcOptions.numChains);
    
    int maxNumberThreads = omp_get_max_threads();
    Random::allocateListRandomNumbersGenerators(randomGenerators);
    
    pll_rtree_t * initialRootedTree = pll_rtree_parse_newick(treefileName);
    
    programOptions.doUseFixedTree=NO;
    mcmcOptions.maxNumberIndependentPosteriorValues=500;
   
    //int iterationToComputeThinnig= floor(mcmcOptions.percentIterationsToComputeThinnig * mcmcOptions.Niterations);
  
    
    vector<Chain*> chainsFortheSameTree(mcmcOptions.numberChainsPerTree);
    int numberTrees = floor((double)mcmcOptions.numChains / mcmcOptions.numberChainsPerTree);
    
    mcmcOptions.fixedValuesForSimulation = true;
    if (programOptions.doUseFixedTree == NO)
             simulateTrees( numberTrees,structuredCoalTrees,  trueTrees,         trueThetas,
                            trueDeltaTs,
                            trueTs,
                          sampleSizes, programOptions,mcmcOptions, randomGenerators,ObservedData,ObservedCellNames, msa,  healthyTipLabel  );
    
    ChainManager * chainManager= new ChainManager(mcmcOptions.numChains);
    omp_set_num_threads(maxNumberThreads);

    mcmcOptions.doMCMCMoveTimeOriginInputOldestPop=true;
    mcmcOptions.iterationToComputeThinnig= floor(0.25 *mcmcOptions.Niterations);
    mcmcOptions.numberWarmUpIterations = 10000;//floor(0.20 *mcmcOptions.Niterations);
    
    
    mcmcOptions.iterationsToMonitorChain= floor(0.05 *mcmcOptions.iterationToComputeThinnig);
    mcmcOptions.iterationsToMonitorChain =   floor(mcmcOptions.numberWarmUpIterations / 20);
    
    mcmcOptions.printChainStateEvery=mcmcOptions.iterationsToMonitorChain;
    chainManager->initializeChains(randomGenerators, programOptions, mcmcOptions, filePaths, sampleSizes, ObservedData, ObservedCellNames, msa, initialRootedTree, structuredCoalTrees, healthyTipLabel, trueTrees, trueThetas, trueDeltaTs, trueTs);
    
    int burnInIterations= 10000;
    mcmcOptions.Niterations= 1000000;
    mcmcOptions.maxNumberIndependentPosteriorValues=400;
    float start = clock();
    clock_t begin = omp_get_wtime();
    bool converged =false;
    //////////////////////
    chainManager->performWarmUp(programOptions, mcmcOptions, randomGenerators, filePaths);
    chainManager->resetChainValues(programOptions, mcmcOptions);
    
    for (size_t currentIteration = 0; currentIteration < mcmcOptions.Niterations; ++currentIteration)
    {
         //reset Chains values
           chainManager->stepAllChainsNoSave( currentIteration,programOptions, mcmcOptions, randomGenerators  );
        
        if (currentIteration  >=  burnInIterations){
              chainManager->saveChainState( currentIteration, programOptions, mcmcOptions );
              chainManager->writeChainsOutput(currentIteration, programOptions, mcmcOptions, filePaths);
        }
        
           if ( (currentIteration ) >=  burnInIterations +mcmcOptions.maxNumberIndependentPosteriorValues  )
               converged = chainManager->checkConvergence(mcmcOptions.numberChainsPerTree, programOptions, mcmcOptions, 0.01);
        
        if (converged){
            break;
            
        }
        
    }
    Random::freeListRandomNumbersGenerators(randomGenerators);
    trueTrees.clear();
    structuredCoalTrees.clear();
    trueThetas.clear();
    trueDeltaTs.clear();
    trueTs.clear();
    sampleSizes.clear();

    chains.clear();

    if (!mcmcOptions.noData && mcmcOptions.useSequencesLikelihood ==1)
        pll_msa_destroy(msa);
    clock_t end = omp_get_wtime();
    double elapsed_time = double(end - begin);
    std::cout << "\n\n*** Program finished ***" <<endl;
    /* execution time */
    double secs = (double)(clock() - start) / CLOCKS_PER_SEC;
    std::cout << "\n\n_________________________________________________________________" <<endl;
    std::cout << "\n sSerial Time processing: " << secs << " seconds and parallel Time processing " << elapsed_time  <<endl;
    std::cout << "\nIf you need help type '-?' in the command line of the program\n"<<endl;
    return 0;
 
}

