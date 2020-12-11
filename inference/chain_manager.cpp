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
 * ChainManager class
 */
#include <numeric>

#include "chain_manager.hpp"
#include "mcmc_chain.hpp"
#include "utils.hpp"
#include<omp.h>

ChainManager::ChainManager(int numberChains){
    if (numberChains >0)
        this->numberChains= numberChains;
    else
        fprintf (stderr, "\n ERROR: The number of chains cannot be negative \n");
    
    for(size_t chainNumber=0; chainNumber <numberChains;++chainNumber){
        chains.push_back(nullptr);
    }
}
void ChainManager::initializeChains(vector<gsl_rng * > &randomGenerators, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, FilePaths &filePaths,
                                    vector<int> &sampleSizes,  std::vector<std::vector<int> > &ObservedData,char* ObservedCellNames[], pll_msa_t *msa, pll_rtree_t * initialRootedTree, vector<StructuredCoalescentTree *> structuredCoalTrees,  string& healthyTipLabel,
                                    const vector<pll_rtree_t *> &trueTrees,        const  vector<long double> &trueThetas,
                                    const vector<vector<long double>> &trueDeltaTs,
                                    const  vector<vector<long double>> &trueTs
                                    )
{

    gsl_rng * currentRandomGenerator;
    StructuredCoalescentTree * currenStructuredCoalTreeToInfer;
    for(size_t chainNumber=0; chainNumber< numberChains;++chainNumber)
    {
        currentRandomGenerator =randomGenerators.at(chainNumber);
        
        int idx = chainNumber / mcmcOptions.numberChainsPerTree;
        currenStructuredCoalTreeToInfer =structuredCoalTrees.at(idx);
        
        chains.at(chainNumber) = Chain::initializeChain( chainNumber, programOptions, mcmcOptions, sampleSizes, currentRandomGenerator, ObservedData,ObservedCellNames, msa,  initialRootedTree,currenStructuredCoalTreeToInfer, healthyTipLabel, filePaths);
        //chains.at(chainNumber)->chainNumber =chainNumber;
        //chains.at(chainNumber)->PrepareFiles(filePaths, programOptions, chains.at(chainNumber)->files, chainNumber);
        
        chains.at(chainNumber)->writeHeaderOutputChain(filePaths, programOptions,
                                                       chains.at(chainNumber)->files,trueTrees.at(idx) ,trueThetas.at(idx), trueDeltaTs.at(idx), trueTs.at(idx), sampleSizes );
        chains.at(chainNumber)->initListMoves(programOptions, mcmcOptions);
    }
}
void ChainManager::stepAllChains(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, vector<gsl_rng * > &randomGenerators  )
{
    gsl_rng * currentRandomGenerator;
    for(size_t chainNumber=0; chainNumber< numberChains;++chainNumber)
    {
        currentRandomGenerator =randomGenerators.at(chainNumber);
        if (currentIteration % mcmcOptions.printChainStateEvery == 0 || currentIteration >= (mcmcOptions.Niterations-1))
            std::cout << "\n Chain " << chains.at(chainNumber)->chainNumber << " Starting Iteration " << currentIteration +1 << endl;
        chains.at(chainNumber)->currentNumberIerations =currentIteration;
        
        chains.at(chainNumber)->stepAllMoves(mcmcOptions,  currentRandomGenerator,  programOptions);
        
        chains.at(chainNumber)->saveMCMCState( currentIteration, programOptions,mcmcOptions);
    }
    
    
}
void ChainManager::stepAllChainsNoSave(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, vector<gsl_rng * > &randomGenerators  )
{
    gsl_rng * currentRandomGenerator;
    for(size_t chainNumber=0; chainNumber< numberChains;++chainNumber)
    {
        currentRandomGenerator =randomGenerators.at(chainNumber);
        if (currentIteration % mcmcOptions.printChainStateEvery == 0 || currentIteration >= (mcmcOptions.Niterations-1))
            std::cout << "\n Chain " << chains.at(chainNumber)->chainNumber << " Starting Iteration " << currentIteration +1 << endl;
        chains.at(chainNumber)->currentNumberIerations =currentIteration;
        
        chains.at(chainNumber)->stepAllMoves(mcmcOptions,  currentRandomGenerator,  programOptions);
        
    }
    
    
}
void ChainManager::saveChainState(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions )
{
    gsl_rng * currentRandomGenerator;
    for(size_t chainNumber=0; chainNumber< numberChains;++chainNumber)
    {
      chains.at(chainNumber)->saveMCMCState( currentIteration, programOptions,mcmcOptions);
        
    }
    
    
}
void ChainManager::monitorChains(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, vector<gsl_rng * > &randomGenerators, FilePaths &filePaths )
{
    
    for(size_t chainNumber=0; chainNumber <numberChains;++chainNumber)
    {
        if (mcmcOptions.doThinning)
        {  chains.at(chainNumber)->computeThinnig( mcmcOptions);
            chains.at(chainNumber)->resetPosteriorValues();
            std::cout << "\n Chain "<< chains.at(chainNumber)->chainNumber << " thinning: "<< chains.at(chainNumber)->thinning;
            //adjust if acceptance rate is  too low (<0.15) or too high(>0.5)
            chains.at(chainNumber)->printMovesSummary(programOptions, mcmcOptions);
        }
      
        
    }
    
}
void ChainManager::writeChainsOutput(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, FilePaths &filePaths)
{
   
    for(size_t chainNumber=0; chainNumber <numberChains;++chainNumber)
    {
     if (chains.at(chainNumber)->thinning >1 && currentIteration % chains.at(chainNumber)->thinning == 0 )
   //if (chains.at(chainNumber)->thinning >=1 )
     {
       
            chains.at(chainNumber)->numberIndependentLongUpdates = chains.at(chainNumber)->numberIndependentLongUpdates +1;
         
             // fprintf (stderr, "\n Updated number of independent samples for chain %lu  is %d \n",chainNumber, chains.at(chainNumber)->numberIndependentLongUpdates );
         
            chains.at(chainNumber)->writeMCMCState(  currentIteration+1, filePaths, programOptions,chains.at(chainNumber)->files, mcmcOptions);
        
     }
    }
    
}
vector<Chain *> ChainManager::getChainsWiththeSameTrueTreeThatAnotherChain(int chainNumber, int numberChainsPerTree)
{
    vector<Chain *> chainsFortheSameTree;
    int idxGroup=0;
    for(size_t i=0; i <numberChainsPerTree;++i)
    {
        idxGroup= chainNumber / numberChainsPerTree;
        chainsFortheSameTree.push_back(chains.at(idxGroup + i));
    }
    return chainsFortheSameTree;
}
bool ChainManager::checkConvergenceGroupChains(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions, long double threshold )
{
    bool result=false;
    vector<long double> potentScaleReduction = listPotentialScaleReduction(chainNumber, numberChainsPerTree,programOptions);
    long double meanPotentialScaleReduction = std::accumulate( potentScaleReduction.begin(), potentScaleReduction.end(), 0.0) / potentScaleReduction.size();
    std::cout << "\n The current mean Potential Scale Reduction Factor is  " << meanPotentialScaleReduction <<  endl;
    if (abs(meanPotentialScaleReduction -1.0) < threshold ){
         std::cout << "\n The final mean Potential Scale Reduction Factor is  " << meanPotentialScaleReduction <<  endl;
        result=true;
    }
    return result;
    
}
    
vector<long double> ChainManager::listPotentialScaleReduction(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions)
{
    vector<Chain *> groupChains= getChainsWiththeSameTrueTreeThatAnotherChain(chainNumber, numberChainsPerTree);
    long double average;
    int numClones = programOptions.numClones;
    vector<long double> elements;
    vector<long double> variances;
    vector<long double> listPotentialScaleReductions;
    vector<long double> meanChainThetas;
    vector<vector<long double>> meanChainDeltaTs;
    vector<long double> meanDeltaTs;
    vector<vector<long double>> meanChainTimeOriginSTDs;
    vector<long double> meanTimeOriginSTDs;
    vector<vector<long double>> meanChainProportions;
    vector<long double> meanProportions;
    
    long double meanOfMeansTheta;
    long double BTheta;
    long double WTheta;
    long double potentialScaleReductionTheta;
    vector<long double> varianceChainTheta;
    

    vector<long double> varianceDeltaTs;
    vector<long double> varianceTs;
    vector<long double> varianceProportions;
    
    vector<vector<long double>> varianceChainDeltaT;
    vector<vector<long double>> varianceChainT;
    vector<vector<long double>> varianceChainProportion;
    
    int m=groupChains.size();
    int n;
    for(size_t i=0; i <groupChains.size();++i)
    {
        meanDeltaTs.clear();
        meanTimeOriginSTDs.clear();
        meanProportions.clear();
        varianceDeltaTs.clear();
        varianceTs.clear();
        varianceProportions.clear();
        auto chain = chains.at(i);
        average = Utils::mean( chain->posteriorTheta);
        n= chain->posteriorTheta.size();
        meanChainThetas.push_back(average  );
        varianceChainTheta.push_back(Utils::variance(chain->posteriorTheta) );
        
       
        for(size_t j=0; j <programOptions.numClones;++j)
        {
            auto popI =  chain->getPopulationbyIndex(j);
             assert(popI->posteriorDeltaT.size()==n);
            average = Utils::mean( popI->posteriorDeltaT);
            meanDeltaTs.push_back(average);
            varianceDeltaTs.push_back(Utils::variance(popI->posteriorDeltaT));
    
            assert(popI->posteriorTimeOriginSTD.size()==n);
            average = Utils::mean( popI->posteriorTimeOriginSTD);
            meanTimeOriginSTDs.push_back(average);
            varianceTs.push_back(Utils::variance(popI->posteriorTimeOriginSTD));
            if (programOptions.numClones >1)
            {
                assert(popI->posteriorProportion.size()==n);
                average = Utils::mean( popI->posteriorProportion);
                meanProportions.push_back(average);
                varianceProportions.push_back(Utils::variance(popI->posteriorProportion));
            }
            
        }
        assert(meanDeltaTs.size()==programOptions.numClones);
        assert(meanTimeOriginSTDs.size()==programOptions.numClones);
        assert(varianceDeltaTs.size()==programOptions.numClones);
        assert(varianceTs.size()==programOptions.numClones);
        
        meanChainDeltaTs.push_back(meanDeltaTs);
        meanChainTimeOriginSTDs.push_back(meanTimeOriginSTDs);
        
        varianceChainDeltaT.push_back(varianceDeltaTs);
        varianceChainT.push_back(varianceTs);
        
        if (programOptions.numClones >1)
        {
            assert(meanProportions.size()==programOptions.numClones);
            assert(varianceProportions.size()==programOptions.numClones);
            meanChainProportions.push_back(meanProportions);
            varianceChainProportion.push_back(varianceProportions);
        }
    }
    meanOfMeansTheta = Utils::mean( meanChainThetas);
    BTheta = n * Utils::variance(meanChainThetas) ;
    WTheta = (1.0 /m)*Utils::mean(varianceChainTheta) ;
    potentialScaleReductionTheta= (1.0-1.0/n)* WTheta+ (1.0/n)*BTheta;
    potentialScaleReductionTheta = sqrt(potentialScaleReductionTheta / WTheta);
    listPotentialScaleReductions.push_back(potentialScaleReductionTheta);
    
  
    vector<long double > potentialScaleReductionsDeltaT =Utils::potentialScaleReductionArray(numClones,n,  m, meanChainDeltaTs, varianceChainDeltaT );
    
    listPotentialScaleReductions.insert(listPotentialScaleReductions.end(), potentialScaleReductionsDeltaT.begin(), potentialScaleReductionsDeltaT.end());
    
    vector<long double > potentialScaleReductionsT = Utils::potentialScaleReductionArray(numClones,n,  m, meanChainTimeOriginSTDs, varianceChainT );
    
    listPotentialScaleReductions.insert(listPotentialScaleReductions.end(), potentialScaleReductionsT.begin(), potentialScaleReductionsT.end());
    
    if (programOptions.numClones >1)
    {
          vector<long double > potentialScaleReductionsProportion =Utils::potentialScaleReductionArray(numClones,n,  m, meanChainProportions, varianceChainProportion );
        listPotentialScaleReductions.insert(listPotentialScaleReductions.end(), potentialScaleReductionsProportion.begin(), potentialScaleReductionsProportion.end());
    }
    return listPotentialScaleReductions;
}
void ChainManager::printGroupChainsSummary(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions)
{
     vector<Chain *> groupChains= getChainsWiththeSameTrueTreeThatAnotherChain(chainNumber, numberChainsPerTree);
    for(size_t i=0; i <groupChains.size();++i)
    {
        auto chain = groupChains.at(i);
        chain->printLastMovesSummary();
    }
}
void ChainManager::deleteGroupChains(int chainNumber, int numberChainsPerTree)
{
    vector<Chain *> groupChains= getChainsWiththeSameTrueTreeThatAnotherChain(chainNumber, numberChainsPerTree);
    for(size_t i=0; i <groupChains.size();++i){
        auto chain = groupChains.at(i);
        delete chain;
    }
}
bool ChainManager::checkGroupChainsHaveCompletedIndependentIterations(int numberIterations, int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions)
{
    bool result=true;
    vector<Chain *> chainsFortheSameTree;
    int idxGroup=0;
    for(size_t i=0; i <numberChainsPerTree;++i)
    {
        idxGroup= chainNumber / numberChainsPerTree;
        result = chains.at(idxGroup + i)->numberIndependentLongUpdates >= numberIterations;
        
        fprintf (stderr, "\n The number of independent samples for chain %lu  is %d \n",idxGroup + i, chains.at(idxGroup + i)->numberIndependentLongUpdates );
       
        if (!result)
            break;
    }
    return result;
}
void ChainManager::closeFilesGroupChain(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions)
{
    
    vector<Chain *> groupChains= getChainsWiththeSameTrueTreeThatAnotherChain(chainNumber, numberChainsPerTree);
    for(size_t i=0; i <groupChains.size();++i)
    {
        auto chain = groupChains.at(i);
        chain->closeFiles(programOptions);
    }
}
bool ChainManager::checkConvergence(int numberChainsPerTree,ProgramOptions &programOptions,MCMCoptions &mcmcOptions, long double threshold )
{
    bool converged=false;
    for(size_t chainNumber=0; chainNumber< numberChains; chainNumber += numberChainsPerTree)
    {
        bool completedIterations=checkGroupChainsHaveCompletedIndependentIterations(mcmcOptions.maxNumberIndependentPosteriorValues, chainNumber,numberChainsPerTree,programOptions);
        if (completedIterations)
        {
            converged= checkConvergenceGroupChains(chainNumber, numberChainsPerTree,programOptions, threshold );
            if (converged)
            {
                printGroupChainsSummary(chainNumber, mcmcOptions.numberChainsPerTree,programOptions);
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
                std::cout << "\n Use of " << priorsType.c_str() << "priors and " << kernelType.c_str() << endl;
                closeFilesGroupChain(chainNumber, mcmcOptions.numberChainsPerTree, programOptions);
                deleteGroupChains(chainNumber, mcmcOptions.numberChainsPerTree);
                break;
            }
        }
    }
    return(converged);
}
void ChainManager::performWarmUp(ProgramOptions &programOptions,MCMCoptions &mcmcOptions, vector<gsl_rng * > &randomGenerators,FilePaths &filePaths)
{
    //gsl_rng * currentRandomGenerator;
//#pragma omp parallel   shared(programOptions,mcmcOptions,filePaths,   randomGenerators )
       {
           gsl_rng * currentRandomGenerator;
  //  #pragma omp for  nowait
    for(int chainNumber=0; chainNumber< numberChains;++chainNumber)
           {
               currentRandomGenerator =randomGenerators.at(chainNumber);
               //gsl_rng_set(currentRandomGenerator, omp_get_thread_num() * 12567);
            
               for (int currentIteration = 0; currentIteration < (mcmcOptions.numberWarmUpIterations); ++currentIteration)
               {
                   
                   if (currentIteration % mcmcOptions.printChainStateEvery == 0 || currentIteration >= (mcmcOptions.Niterations-1))
                       std::cout << "\n Chain " << chains.at(chainNumber)->chainNumber << " Starting Iteration " << currentIteration +1 << endl;
                   chains.at(chainNumber)->currentNumberIerations =currentIteration;
                   
                   chains.at(chainNumber)->stepAllMoves(mcmcOptions,  currentRandomGenerator,  programOptions);
                   
                   chains.at(chainNumber)->saveMCMCState( currentIteration, programOptions,mcmcOptions);
                   
                   if (currentIteration>0 && currentIteration %  mcmcOptions.iterationsToMonitorChain ==0 && currentIteration < mcmcOptions.numberWarmUpIterations )
                   {
                       if (mcmcOptions.doThinning)
                       {  chains.at(chainNumber)->computeThinnig( mcmcOptions);
                           chains.at(chainNumber)->resetPosteriorValues();
                           std::cout << "\n Chain "<< chains.at(chainNumber)->chainNumber << " thinning: "<< chains.at(chainNumber)->thinning;
                           chains.at(chainNumber)->printMovesSummary(programOptions, mcmcOptions);
                       }
                   }
                   
               }
           }
   }
}
void ChainManager::resetChainValues(ProgramOptions &programOptions,MCMCoptions &mcmcOptions){
    
    
    
    for(int chainNumber=0; chainNumber< numberChains;++chainNumber)
    {
        
        chains.at(chainNumber)->resetPosteriorValues();
        
    }
    
    
    
}
