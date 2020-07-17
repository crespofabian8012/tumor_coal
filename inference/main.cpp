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
#include "definitions.hpp"
#include "utils.hpp"
#include "constants.hpp"

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
    
    char *fileNameFasta ;//= argv[2];
    char *fileNamePhylip ;//= argv[3];
    char *treefileName; //= argv[4];
    
    ProgramOptions programOptions;
    //Files files;
    FilePaths filePaths;
    MCMCoptions mcmcOptions;
    
    mcmcOptions.noData = false;
    mcmcOptions.useGSLRandomGenerator = true;
    mcmcOptions.splitThetaDeltaTmoves=false;
    
    setDefaultOptions(programOptions, mcmcOptions);
    if (argc == 1 )
    {
        
        fprintf (stderr, "\nERROR: No parameters specified (use command  parameter file)");
        PrintUsage();
    }
    if (argc <= 2 )//only a path to a MCMC configuration file
        input_path = argv[1];
    else
    {//many options from console
        
        if (!mcmcOptions.noData)
        {
            fprintf (stderr, "\nERROR: No parameters specified (use command  parameter file)");
            PrintUsage();
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
                
                fprintf (stderr, "\nERROR: No parameters specified (use command line or parameter file)");
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
    char *ObservedCellNames[programOptions.numCells];
    pll_msa_t *msa;
    
    
    //    if (mcmcOptions.useSequencesLikelihood ==1)
    //   {
    fileNameFasta = filePaths.inputGenotypeFileFasta;
    Utils::ReadParametersFromFastaFile(fileNameFasta,  programOptions);
    std::vector<std::vector<int> > ObservedData;
    
    Utils::ReadFastaFile(fileNameFasta, ObservedData,  ObservedCellNames, programOptions);
    programOptions.numNodes = 2 * programOptions.TotalNumSequences + programOptions.numClones+ 10;
    programOptions.numCells = programOptions.TotalNumSequences;
    
    if (programOptions.numClones==1)
        sampleSizes.push_back(programOptions.numCells-1); //minus the healthycell
    
    fileNamePhylip =filePaths.inputGenotypeFilePhylip;
    msa = pll_phylip_load(fileNamePhylip, PLL_FALSE);
    if (!msa)
        fprintf (stderr, "Error reading phylip file \n");
    // }
    treefileName = filePaths.inputTreeFile;
    if ((treefileName != NULL) && (treefileName[0] == '\0'))
    {
        fprintf (stderr, "\nERROR: No tree file  specified (use command line or parameter file)");
        exit(-1);
    }
    std::string healthyTipLabel = "healthycell";
    programOptions.healthyTipLabel ="healthycell";
    
    int currentIteration;
    int sampleEvery = mcmcOptions.thinning;
    vector<Chain*> chains(mcmcOptions.numChains);
    //parallelizing
    
    
    vector<gsl_rng * > randomGenerators(mcmcOptions.numChains);
    gsl_rng * currentRandomGenerator;
    Random::allocateListRandomNumbersGenerators(randomGenerators);
    
    pll_rtree_t * initialRootedTree = pll_rtree_parse_newick(treefileName);
    
    float start = clock();
    unsigned int chainNumber;
    unsigned int currentSimulation;
    programOptions.doUseFixedTree=NO;
    
    //#pragma omp parallel default(shared) shared(treefileName) private(chainNumber, currentIteration)
    //omp_set_num_threads(mcmcOptions.numChains);
    {
        //#pragma omp for
        for(chainNumber=0; chainNumber< mcmcOptions.numChains;chainNumber++)
        {
            currentRandomGenerator =randomGenerators.at(chainNumber);
            chains.at(chainNumber) = Chain::initializeChain( programOptions, mcmcOptions, sampleSizes, currentRandomGenerator, ObservedData,ObservedCellNames, msa,  initialRootedTree, healthyTipLabel);
            chains.at(chainNumber)->chainNumber =chainNumber;
            chains.at(chainNumber)->PrepareFiles(filePaths, programOptions, chains.at(chainNumber)->files, chainNumber);
            chains.at(chainNumber)->writeHeaderOutputChain(filePaths, programOptions,
                                                           chains.at(chainNumber)->files );
            chains.at(chainNumber)->initListMoves();
            for (currentIteration = 0; currentIteration < mcmcOptions.Niterations; currentIteration++)
            {
                if (currentIteration % mcmcOptions.printChainStateEvery == 0 || currentIteration >= (mcmcOptions.Niterations-1))
                    fprintf (stderr, "\n Chain #%d, Iteration %d \n", chains.at(chainNumber)->chainNumber ,currentIteration +1);
                chains.at(chainNumber)->currentNumberIerations =currentIteration;
                chains.at(chainNumber)->runChain(mcmcOptions,  currentRandomGenerator,  filePaths, chains.at(chainNumber)->files, programOptions,ObservedCellNames, msa, sampleSizes, currentIteration );
                if (currentIteration % sampleEvery == 0 && currentIteration >= mcmcOptions.numberWarmUpIterations)
                {
                    chains.at(chainNumber)->writeMCMCState(  currentIteration+1, filePaths, programOptions,chains.at(chainNumber)->files, mcmcOptions);
                    //in the future i will write trees to a nexus file
                    ////PrintTrees(currentIteration, &(chains[chainNumber].root),
                }
            }
            chains.at(chainNumber)->currentNumberIerations =currentIteration;
            if (currentIteration >= mcmcOptions.Niterations -1)//last iteration
            {
                chains.at(chainNumber)->printMovesSummary();
                string priorsType;
                if (mcmcOptions.priorsType==0)
                    priorsType="log uniform";
                else if(mcmcOptions.priorsType==1)
                    priorsType="exponential";
                else
                    priorsType="power law";
                string kernelType;
                if (mcmcOptions.kernelType==0)
                    kernelType="multiplier";
                else //(mcmcOptions.priorsType==1)
                    kernelType="normal";
                printf ("\n Use of %s priors and %s kernels \n", priorsType.c_str(), kernelType.c_str());
                fclose(chains.at(chainNumber)->files.fplog);
                //chains.at(chainNumber)->closeFiles(filePaths,programOptions, files, chainNumber);
            }
        }
    }
    Random::freeListRandomNumbersGenerators(randomGenerators);
    //for default(shared) private(chainNumber, currentIteration) firstprivate(files)
    //close files
    if (!mcmcOptions.noData && mcmcOptions.useSequencesLikelihood ==1)
        pll_msa_destroy(msa);
    fprintf(stderr, "\n\n*** Program finished ***");
    /* execution time */
    double secs = (double)(clock() - start) / CLOCKS_PER_SEC;
    fprintf(stderr, "\n\n_________________________________________________________________");
    fprintf(stderr, "\nTime processing: %G seconds\n", secs);
    fprintf(stderr, "\nIf you need help type '-?' in the command line of the program\n");
    return 0;
    return 0;
}

