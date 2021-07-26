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
#include <boost/program_options.hpp>

//#include <spf.hpp>
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



int main(int argc, char* argv[] )
{
    FILE *input_file;
    const char* input_path;
    std::string config_file;
    
    namespace po = boost::program_options;
    po::options_description desc("Program options");
    desc.add_options()
    ("help", "Put a help message here.")
    ("config_file,c", po::value<std::string>(&config_file)->required(), "path to configuration file.")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    input_path = config_file.c_str();
    
    ProgramOptions programOptions;
    FilePaths filePaths;
    MCMCoptions mcmcOptions;
    bool doMCMC = true;

    if (config_file.empty())
     {
         std::cout << "\nERROR: No parameters specified (use command  parameter file)"<< std::endl;
         Output::PrintUsage();
     }
     else
     {
         input_path = config_file.c_str();
         // 1. call function to parse the MCMC configuration file
         if ((input_file = freopen(input_path, "r", stdin)) != NULL)
               {
                   ReadMCMCParametersFromFile(programOptions, filePaths, mcmcOptions);
               }
               else
               {
                   if (!mcmcOptions.noData)
                   {
                        std::cout << "\nERROR: No parameters specified (use command line or parameter file)"<< std::endl;
                       exit(-1);
                   }
               }
     }
    printProgramHeader();
    setDefaultOptions(programOptions, mcmcOptions);
    
    std::vector<StructuredCoalescentTree *> structuredCoalTrees;
    std::vector<pll_rtree_t *> trueTrees;
    std::vector<long double> trueThetas;
    std::vector<std::vector<long double>> trueDeltaTs;
    std::vector<std::vector<long double>> trueTs;

    char *fileNameFasta ;//= argv[2];
    char *fileNamePhylip ;//= argv[3];
    char *treefileName; //= argv[4];
    
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
         std::cout << "Error reading phylip file \n"<< std::endl;
    // }
    treefileName = filePaths.inputTreeFile;
    if ((treefileName != NULL) && (treefileName[0] == '\0'))
    {
        std::cout << "\nERROR: No tree file  specified (use command line or parameter file)" << std::endl;
        exit(-1);
    }
    std::string healthyTipLabel = "healthycell";
    programOptions.healthyTipLabel ="healthycell";

    //int sampleEvery = mcmcOptions.thinning;
    std::vector<Chain*> chains(mcmcOptions.numChains);
    std::vector<gsl_rng * > randomGenerators(mcmcOptions.numChains);
    
    int maxNumberThreads = omp_get_max_threads();
    Random::allocateListRandomNumbersGenerators(randomGenerators);
    std::vector<boost::mt19937 *> randomGeneratorsBoost = Random::allocateListRandomNumbersGeneratorsBoost(mcmcOptions.numChains);
    
    pll_rtree_t * initialRootedTree = pll_rtree_parse_newick(treefileName);
    
    programOptions.seqErrorRate = 0.01;
    programOptions.dropoutRate =0.2;
    mcmcOptions.noData = false;
    mcmcOptions.useGSLRandomGenerator = true;
    mcmcOptions.splitThetaDeltaTmoves=true;
    programOptions.doUseFixedTree=NO;
    //int iterationToComputeThinnig= floor(mcmcOptions.percentIterationsToComputeThinnig * mcmcOptions.Niterations);
    mcmcOptions.doMCMCMoveTimeOriginInputOldestPop=false;
    mcmcOptions.iterationToComputeThinnig= floor(0.25 *mcmcOptions.Niterations);
    mcmcOptions.numberWarmUpIterations = 10000;//floor(0.20 *mcmcOptions.Niterations);
    mcmcOptions.iterationsToMonitorChain= floor(0.05*mcmcOptions.iterationToComputeThinnig);
    mcmcOptions.iterationsToMonitorChain =   floor(mcmcOptions.numberWarmUpIterations / 20);
    mcmcOptions.verbose=0;
    mcmcOptions.burnInIterations= 5000;
    mcmcOptions.Niterations= 500000;
    mcmcOptions.maxNumberIndependentPosteriorValues=200;
    mcmcOptions.maxNumberIndependentPosteriorValues=500;
    mcmcOptions.numberIterationsAfterConvergence = 10000;
    
    
    std::vector<Chain*> chainsFortheSameTree(mcmcOptions.numberChainsPerTree);
    int numberTrees = floor((double)mcmcOptions.numChains / mcmcOptions.numberChainsPerTree);
    
    mcmcOptions.fixedValuesForSimulation = true;
    programOptions.doUsefixedMutationRate=true;
    programOptions.K=0.0;
    
    if (programOptions.doUseFixedTree == NO){
        std::cout << "\nSimulating " << numberTrees << " trees .... \n"<<std::endl;
         simulateTrees( numberTrees,structuredCoalTrees,  trueTrees,         trueThetas,
           trueDeltaTs,
           trueTs,
         sampleSizes, programOptions,mcmcOptions, randomGenerators,randomGeneratorsBoost, ObservedData,ObservedCellNames, msa,  healthyTipLabel,   programOptions.seqErrorRate,
             programOptions.dropoutRate  );
        std::cout << "\nSimulation trees ended. \n"<< std::endl;
        
    }
           
    programOptions.doUsefixedMutationRate = false;
    programOptions.K=0.8;

    if (doMCMC){
        
          ChainManager * chainManager= new ChainManager(mcmcOptions.numChains);
            omp_set_num_threads(maxNumberThreads);

           
           mcmcOptions.printChainStateEvery=mcmcOptions.iterationsToMonitorChain;
           chainManager->initializeChains(randomGenerators, programOptions, mcmcOptions, filePaths,
                                           sampleSizes, ObservedData, ObservedCellNames, msa, initialRootedTree,
                                           structuredCoalTrees, healthyTipLabel, trueTrees, trueThetas, trueDeltaTs, trueTs);
            
            float start = clock();
            clock_t begin = omp_get_wtime();
            //bool converged =false;
            //////////////////////
            std::cout << "\n Tuning MCMC parameters... \n"<<std::endl;
            chainManager->performWarmUp(programOptions, mcmcOptions, randomGenerators, filePaths);
            std::cout << "\n MCMC parameters tuned! \n"<<std::endl;
            chainManager->resetChainValues(programOptions, mcmcOptions);
            //chainManager->resizeStoredMCMCparametersAfterWarmUp(mcmcOptions);
            
            for (size_t currentIteration = 0; currentIteration < mcmcOptions.Niterations; ++currentIteration)
            {
                chainManager->stepAllChainsNoSave( currentIteration,programOptions, mcmcOptions, randomGenerators  );
                
               if (currentIteration >=  mcmcOptions.burnInIterations ){
                      chainManager->saveChainState( currentIteration, programOptions, mcmcOptions );
        //              //chainManager->saveChainState( currentIteration, programOptions, mcmcOptions );
        //              chainManager->writeChainsOutput(currentIteration, programOptions, mcmcOptions, filePaths);
                }
                 chainManager->writeChainsOutput(currentIteration, programOptions, mcmcOptions, filePaths);
                if ( (currentIteration >= (mcmcOptions.burnInIterations +mcmcOptions.maxNumberIndependentPosteriorValues)) &&  currentIteration % 1000 ==0 ){
                      chainManager->checkConvergence(mcmcOptions.numberChainsPerTree, programOptions, mcmcOptions, 0.01) ;
                }
                if (chainManager->numberFinishedChains >= 0.98*mcmcOptions.numChains)
                    break;
        //        if (converged){
        //            std::cout << "\nMCMC chains converged .... \n"<<std::endl;
        //            //break;
        //           // mcmcOptions.Niterations = mcmcOptions.Niterations + 10000;
        //        }
            }
        
        chains.clear();

        if (!mcmcOptions.noData && mcmcOptions.useSequencesLikelihood ==1)
             pll_msa_destroy(msa);
         clock_t end = omp_get_wtime();
         
         double elapsed_time = double(end - begin);
         std::cout << "\n\n*** Program finished ***" <<std::endl;
         /* execution time */
         double secs = (double)(clock() - start) / CLOCKS_PER_SEC;
         std::cout << "\n\n_________________________________________________________________" <<std::endl;
         std::cout << "\n sSerial Time processing: " << secs << " seconds and parallel Time processing " << elapsed_time  <<std::endl;
        
    }
    else/* SMC */
    {
        
        
        
        
        
    }
    
    /* clean memory*/
 
    Random::freeListRandomNumbersGenerators(randomGenerators);
    trueTrees.clear();
    structuredCoalTrees.clear();
    trueThetas.clear();
    trueDeltaTs.clear();
    trueTs.clear();
    sampleSizes.clear();

 
    std::cout << "\nIf you need help type '-?' in the command line of the program\n"<<std::endl;
    
    return 0;
 
}

