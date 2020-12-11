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
 * chain_manager
 */
#ifndef chain_manager_hpp
#define chain_manager_hpp

#include "mcmc_chain.hpp"
using namespace std;

class Chain;
class ChainManager
{
public:
    int numberChains;
    std::vector<Chain *>chains;
public:
    ChainManager(int numberChains);
    void initializeChains(vector<gsl_rng * > &randomGenerators, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, FilePaths &filePaths,
                                        vector<int> &sampleSizes,  std::vector<std::vector<int> > &ObservedData,char* ObservedCellNames[], pll_msa_t *msa, pll_rtree_t * initialRootedTree, vector<StructuredCoalescentTree *> structuredCoalTrees,  string& healthyTipLabel,
                                        const vector<pll_rtree_t *> &trueTrees,        const  vector<long double> &trueThetas,
                                        const vector<vector<long double>> &trueDeltaTs,
                                        const  vector<vector<long double>> &trueTs
                          );
    void stepAllChains(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, vector<gsl_rng * > &randomGenerators );
    void monitorChains(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, vector<gsl_rng * > &randomGenerators, FilePaths &filePaths );
    vector<Chain *> getChainsWiththeSameTrueTreeThatAnotherChain(int chainNumber, int numberChainsPerTree);
    bool checkConvergenceGroupChains(int chainNumber, int numberChainsPerTree, ProgramOptions &programOptions, long double threshold);
    vector<long double> listPotentialScaleReduction(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions);
    void printGroupChainsSummary(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions);
    void deleteGroupChains(int chainNumber, int numberChainsPerTree);
    bool checkGroupChainsHaveCompletedIndependentIterations(int numberIterations, int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions);
    void closeFilesGroupChain(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions);
    bool checkConvergence(int numberChainsPerTree,ProgramOptions &programOptions,MCMCoptions &mcmcOptions, long double threshold );
    //void checkConvergence(int numberChainsPerTree,ProgramOptions &programOptions,MCMCoptions &mcmcOptions, long double threshold);
    void writeChainsOutput(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, FilePaths &filePaths);
    void performWarmUp(ProgramOptions &programOptions,MCMCoptions &mcmcOptions, vector<gsl_rng * > &randomGenerators,FilePaths &filePaths);
    void resetChainValues(ProgramOptions &programOptions,MCMCoptions &mcmcOptions);
    void stepAllChainsNoSave(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, vector<gsl_rng * > &randomGenerators);
    void saveChainState(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions );
    
};
#endif /* chain_manager_hpp */
