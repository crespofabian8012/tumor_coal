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
#include <Eigen/Dense>

class Chain;
class ChainManager
{
public:
    int numberChains;
    int numMCMCParametersChain;
    std::vector<Chain *>chains;
    Eigen::Matrix<Eigen::MatrixXd, Eigen::Dynamic, 1> samples_;
    Eigen::VectorXi warmup_;
    int numberConvergedChains;
    int numberFinishedChains;
    
public:
    ChainManager(int numberChains);
    void initializeChains(std::vector< gsl_rng * > &randomGenerators, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, FilePaths &filePaths,
                          std::vector<int> &sampleSizes,  std::vector<std::vector<int> > &ObservedData,char* ObservedCellNames[], pll_msa_t *msa, pll_rtree_t * initialRootedTree, std::vector<StructuredCoalescentTree *> structuredCoalTrees,  std::string& healthyTipLabel,
                          const std::vector<pll_rtree_t *> &trueTrees,        const  std::vector<long double> &trueThetas,
                          const std::vector<std::vector<long double>> &trueDeltaTs,
                          const  std::vector<std::vector<long double>> &trueTs
                          );
    void stepAllChains(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, std::vector<const gsl_rng * > &randomGenerators );
    void monitorChains(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, std::vector<const gsl_rng * > &randomGenerators, FilePaths &filePaths );
    std::vector<Chain *> getChainsWiththeSameTrueTreeThatAnotherChain(int chainNumber, int numberChainsPerTree);
    bool checkPSRFGroupChains(int chainNumber, int numberChainsPerTree, ProgramOptions &programOptions, long double threshold);
    std::vector<long double> listPotentialScaleReduction(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions);
    void printGroupChainsSummary(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions);
    void deleteGroupChains(int chainNumber, int numberChainsPerTree);
    bool checkGroupChainsHaveCompletedIndependentIterations(int numberIterations, int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions);
    void closeFilesGroupChain(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions);
    bool checkConvergence(int numberChainsPerTree,ProgramOptions &programOptions,MCMCoptions &mcmcOptions, long double threshold );
    //void checkConvergence(int numberChainsPerTree,ProgramOptions &programOptions,MCMCoptions &mcmcOptions, long double threshold);
    void writeChainsOutput(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, FilePaths &filePaths);
    void performWarmUp(ProgramOptions &programOptions,MCMCoptions &mcmcOptions, std::vector<gsl_rng * > &randomGenerators,FilePaths &filePaths);
    void resetChainValues(ProgramOptions &programOptions,MCMCoptions &mcmcOptions);
    void stepAllChainsNoSave(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, std::vector<gsl_rng * > &randomGenerators);
    void saveChainState(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions );
    void addSample(const int chain, const Eigen::MatrixXd& sample) ;
    void addSample(const Eigen::MatrixXd& sample);
    void addSample(const std::vector<std::vector<double> >& sample);
    Eigen::VectorXd samples(const int chain, const int index, int numberSamples) const;
    int numSamples(const int chain, const int idxMCMCParameter) const {
       assert(chain >=0 && chain <chains.size());
        int  sampleSize= chains[chain]->getSampleSizeByIndex(idxMCMCParameter);
        return sampleSize;
    }
   int numSamplesEigen(const int chain) const {
       
       return samples_(chain).rows() ;
       
   }
    
    double splitPotentialScaleReductionFactor( const int index, std::vector<Chain *> chainsWithSameTree, int numberChainsPerTree) const;
    double splitPotentialScaleReductionFactor(const std::string& name,  std::vector<Chain *> chainsWithSameTree, int numberChainsPerTree) const;
    double effectiveSampleSize(const int index,  std::vector<Chain *> chainsWithSameTree, int numberChainsPerTree) const ;
    double effectiveSampleSize( const std::string& name,  std::vector<Chain *> chainsWithSameTree, int numberChainsPerTree) const;
    double splitEffectiveSampleSize( const int index,  std::vector<Chain *> chainsWithSameTree, int numberChainsPerTree) const;
    double splitEffectiveSampleSize( const std::string& name,  std::vector<Chain *> chainsWithSameTree, int numberChainsPerTree) const;
    std::vector<long double> listPotentialScaleReductionFactors(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions);
    std::vector<long double> listEffectiveSampleSizes(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions);
    bool checkEffectiveSampleSizeGroupChains( int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions);
    void setWarmup(const int chain, const int warmup) {
      warmup_(chain) = warmup;
     //   std::cout << "\n Chain " << chain <<  " updated number  Warmup " << warmup << std::endl;
    }
    inline int warmup(const int chain) const { return warmup_(chain); }
    inline int numberSamplesAfterWarmUp(const int chain, const int idxMCMCParameter) const {
        int numberSamples = numSamples(chain, idxMCMCParameter);
        int numberWarmupSamples= warmup(chain);
        int result;
        if(numberSamples >=numberWarmupSamples)
          result= numberSamples;
        else//if (numberSamples==0)(unsaved parameters)
          result=  numberSamples;
        
//        std::cout << "\n Chain " << chain <<  " numberSamplesAfterWarmup samples " << numberSamples << "  numberWarmup "<< numberWarmupSamples << " result "<< result << std::endl;
        return result;
    }
    void resizeStoredMCMCparametersAfterWarmUp(MCMCoptions &mcmcOptions);
    double* getLastChainSamples(const int chainId, const int idxMCMCParameter, int numberSamples) const;
    void  finish(int numberChainsPerTree, ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void setConvergedGroupChains(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions);
    bool groupChainsFinished(std::vector<Chain *> groupChains);
    bool checksplitPSRFGroupChains(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions, long double threshold);
};
#endif /* chain_manager_hpp */
