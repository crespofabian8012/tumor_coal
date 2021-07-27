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

#include "autocovariance.hpp"

#include "chain_manager.hpp"

#include "utils.hpp"
//#include<omp.h>

using namespace std;
ChainManager::ChainManager(int numberChains){
    if (numberChains >0)
        this->numberChains= numberChains;
    else
        fprintf (stderr, "\n ERROR: The number of chains cannot be negative \n");
    
    for(size_t chainNumber=0; chainNumber <numberChains;++chainNumber)
    {
        chains.push_back(nullptr);
    }
    numberConvergedChains = 0;
    numberFinishedChains = 0;
    warmup_.resize(numberChains);
}
void ChainManager::initializeChains(std::vector< gsl_rng * > &randomGenerators, ProgramOptions &programOptions,          MCMCoptions &mcmcOptions, FilePaths &filePaths,
                                    std::vector<int> &sampleSizes,  std::vector<std::vector<int> > &ObservedData,char* ObservedCellNames[], pll_msa_t *msa, pll_rtree_t * initialRootedTree, std::vector<StructuredCoalescentTree *> structuredCoalTrees,  std::string& healthyTipLabel,
                                    const std::vector<pll_rtree_t *> &trueTrees,        const  std::vector<long double> &trueThetas,
                                    const std::vector<std::vector<long double>> &trueDeltaTs,
                                    const  std::vector<std::vector<long double>> &trueTs
                                    )
{
    StructuredCoalescentTree * currenStructuredCoalTreeToInfer;
    for(size_t chainNumber=0; chainNumber< numberChains;++chainNumber)
    {
        int idx = chainNumber / mcmcOptions.numberChainsPerTree;
        currenStructuredCoalTreeToInfer =structuredCoalTrees.at(idx);
        chains.at(chainNumber) = Chain::initializeChain( chainNumber, programOptions, mcmcOptions, sampleSizes, randomGenerators.at(chainNumber), ObservedData,ObservedCellNames,                                              msa,  initialRootedTree,currenStructuredCoalTreeToInfer, healthyTipLabel, filePaths);
        setWarmup(chainNumber, mcmcOptions.burnInIterations);
        if (chainNumber==0)
        {
            numMCMCParametersChain =  chains.at(chainNumber)->numParams();
        }
        else
        {
            assert(numMCMCParametersChain ==  chains.at(chainNumber)->numParams());
        }
        chains.at(chainNumber)->writeHeaderOutputChain(filePaths, programOptions,
                                                       chains.at(chainNumber)->files,trueTrees.at(idx) ,
                                                       trueThetas.at(idx), trueDeltaTs.at(idx),
                                                       trueTs.at(idx), sampleSizes, currenStructuredCoalTreeToInfer );
        chains.at(chainNumber)->initListMoves(programOptions, mcmcOptions);
    }
}
void ChainManager::stepAllChains(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, std::vector<const gsl_rng * > &randomGenerators  )
{
    for(size_t chainNumber=0; chainNumber< numberChains;++chainNumber)
    {
        if (currentIteration % mcmcOptions.printChainStateEvery == 0 || currentIteration >= (mcmcOptions.Niterations-1))
            std::cout << "\n Chain " << chains.at(chainNumber)->chainNumber << " Starting Iteration " << currentIteration +1 << std::endl;
        chains.at(chainNumber)->currentNumberIerations =currentIteration;
        chains.at(chainNumber)->stepAllMoves(mcmcOptions,  randomGenerators.at(chainNumber),  programOptions);
        chains.at(chainNumber)->saveMCMCState( currentIteration, programOptions,mcmcOptions);
        //addSample(chainNumber, chains.at(chainNumber)->getStoredSamples());
    }
}
void ChainManager::stepAllChainsNoSave(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, std::vector<gsl_rng * > &randomGenerators  )
{
    string state;
    for(size_t chainNumber=0; chainNumber< numberChains;++chainNumber)
    {
        if (chains.at(chainNumber)->numberIndependentLongUpdates < mcmcOptions.numberIterationsAfterConvergence ){
            state = (currentIteration<= mcmcOptions.numberWarmUpIterations)? "(Warmup)": "(Sampling)";
            if (currentIteration % mcmcOptions.printChainStateEvery == 0 || currentIteration >= (mcmcOptions.Niterations-1))
                std::cout << "\n Chain " << chains.at(chainNumber)->chainNumber << ":  Iteration " << currentIteration +1 << " / " << mcmcOptions.Niterations << " [" <<      (currentIteration +1)*100/mcmcOptions.Niterations<< "] " << state<< std::endl;
            chains.at(chainNumber)->currentNumberIerations =currentIteration;
            
            chains.at(chainNumber)->stepAllMoves(mcmcOptions,  randomGenerators.at(chainNumber),  programOptions);
        }
    }
}
void ChainManager::saveChainState(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions )
{
    for(size_t chainNumber=0; chainNumber< numberChains;++chainNumber)
    {
        if (chains.at(chainNumber)->numberIndependentLongUpdates < mcmcOptions.numberIterationsAfterConvergence )
        {
            chains.at(chainNumber)->saveMCMCState( currentIteration, programOptions,mcmcOptions);
        }
        //addSample(chainNumber, chains.at(chainNumber)->getStoredSamples());
    }
}
void ChainManager::monitorChains(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, std::vector<const gsl_rng * > &randomGenerators, FilePaths &filePaths )
{
    for(size_t chainNumber=0; chainNumber <numberChains;++chainNumber)
    {
        if (mcmcOptions.doThinning)
        {
            chains.at(chainNumber)->computeThinnig( mcmcOptions);
            chains.at(chainNumber)->resetPosteriorValues();
            std::cout << "\n Chain "<< chains.at(chainNumber)->chainNumber << " thinning: "<< chains.at(chainNumber)->thinning;
            chains.at(chainNumber)->printMovesSummary(programOptions, mcmcOptions);
        }
    }
}
void ChainManager::writeChainsOutput(int currentIteration, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, FilePaths &filePaths)
{
    for(size_t chainNumber=0; chainNumber <numberChains;++chainNumber)
    {
        //if (chains.at(chainNumber)->thinning >1 && currentIteration % chains.at(chainNumber)->thinning == 0 )
        //  if (chains.at(chainNumber)->thinning >=1 )
        //   {
        //if (chains.at(chainNumber)->converged && chains.at(chainNumber)->numberIndependentLongUpdates < mcmcOptions.numberIterationsAfterConvergence && !(chains.at(chainNumber)->endReached))
        if (currentIteration % chains.at(chainNumber)->thinning == 0 && chains.at(chainNumber)->numberIndependentLongUpdates < mcmcOptions.numberIterationsAfterConvergence && !(chains.at(chainNumber)->endReached))
        {
            chains.at(chainNumber)->numberIndependentLongUpdates = chains.at(chainNumber)->numberIndependentLongUpdates +1;
            // fprintf (stderr, "\n Updated number of independent samples for chain %lu  is %d \n",chainNumber, chains.at(chainNumber)->numberIndependentLongUpdates );
            chains.at(chainNumber)->writeMCMCState(  currentIteration+1, filePaths, programOptions,chains.at(chainNumber)->files, mcmcOptions);
        }
        if(chains.at(chainNumber)->converged && chains.at(chainNumber)->numberIndependentLongUpdates >= mcmcOptions.numberIterationsAfterConvergence && !(chains.at(chainNumber)->endReached)){
            
            chains.at(chainNumber)->endReached=true;
            
            std::cout << "\nMCMC chain: "<< chainNumber << " reached "<< mcmcOptions.numberIterationsAfterConvergence << " samples "<<  " after convergence. \n"<< std::endl;
            
            closeFilesGroupChain(chainNumber, mcmcOptions.numberChainsPerTree, programOptions);
            
            numberFinishedChains++;
            std::cout << "\nNumber of MCMC chains that finished: "<< numberFinishedChains << " with  "<< mcmcOptions.numberIterationsAfterConvergence << " saved samples after convergence "<< std::endl;
            
        }
        // }
    }
}
std::vector<Chain *> ChainManager::getChainsWiththeSameTrueTreeThatAnotherChain(int chainNumber, int numberChainsPerTree)
{
    std::vector<Chain *> chainsFortheSameTree;
    int idxGroup=0;
    for(size_t i=0; i <numberChainsPerTree;++i)
    {
        idxGroup= chainNumber / numberChainsPerTree;
        chainsFortheSameTree.push_back(chains.at(idxGroup + i));
    }
    return chainsFortheSameTree;
}
bool ChainManager::checksplitPSRFGroupChains(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions, long double threshold )
{
    bool result=false;
    //std::vector<long double> potentScaleReduction = listPotentialScaleReduction(chainNumber, numberChainsPerTree,programOptions);
    std::vector<long double> potentScaleReduction = listPotentialScaleReductionFactors( chainNumber,  numberChainsPerTree, programOptions);
    
    long double meanPotentialScaleReduction = std::accumulate( potentScaleReduction.begin(), potentScaleReduction.end(), 0.0) / potentScaleReduction.size();
    std::cout << "\n The current mean  split Potential Scale Reduction Factor for chains " << chainNumber << " to " << (chainNumber+numberChainsPerTree-1) << " is  " << meanPotentialScaleReduction <<  std::endl;
    if (abs(meanPotentialScaleReduction -1.0) < threshold ){
        std::cout << "\n The final mean split Potential Scale Reduction Factor for chains " << chainNumber << " to " << (chainNumber+numberChainsPerTree-1) << " is  " << meanPotentialScaleReduction <<  std::endl;
        result=true;
    }
    return result;
}
bool ChainManager::checkPSRFGroupChains(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions, long double threshold )
{
    bool result=false;
    std::vector<long double> potentScaleReduction = listPotentialScaleReduction(chainNumber, numberChainsPerTree,programOptions);
    
    long double meanPotentialScaleReduction = std::accumulate( potentScaleReduction.begin(), potentScaleReduction.end(), 0.0) / potentScaleReduction.size();
    std::cout << "\n The current mean Potential Scale Reduction Factor for chains " << chainNumber << " to " << (chainNumber+numberChainsPerTree-1) << " is  " << meanPotentialScaleReduction <<  std::endl;
    if (abs(meanPotentialScaleReduction -1.0) < threshold ){
        std::cout << "\n The final mean Potential Scale Reduction Factor for chains " << chainNumber << " to " << (chainNumber+numberChainsPerTree-1) << " is  " << meanPotentialScaleReduction <<  std::endl;
        result=true;
    }
    return result;
}

std::vector<long double> ChainManager::listPotentialScaleReduction(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions)
{
    std::vector<Chain *> groupChains= getChainsWiththeSameTrueTreeThatAnotherChain(chainNumber, numberChainsPerTree);
    long double average;
    int numClones = programOptions.numClones;
    std::vector<long double> elements;
    std::vector<long double> variances;
    std::vector<long double> listPotentialScaleReductions;
    std::vector<long double> meanChainThetas;
    std::vector<std::vector<long double>> meanChainDeltaTs;
    std::vector<long double> meanDeltaTs;
    std::vector<std::vector<long double>> meanChainTimeOriginSTDs;
    std::vector<long double> meanTimeOriginSTDs;
    std::vector<std::vector<long double>> meanChainProportions;
    std::vector<long double> meanProportions;
    
    long double meanOfMeansTheta;
    long double BTheta;
    long double WTheta;
    long double potentialScaleReductionTheta;
    std::vector<long double> varianceChainTheta;
    
    
    std::vector<long double> varianceDeltaTs;
    std::vector<long double> varianceTs;
    std::vector<long double> varianceProportions;
    
    std::vector<std::vector<long double>> varianceChainDeltaT;
    std::vector<std::vector<long double>> varianceChainT;
    std::vector<std::vector<long double>> varianceChainProportion;
    
    int m=groupChains.size();
    int n;
    long double currentVariance;
    
    for(size_t i=0; i <groupChains.size();++i)
    {
        meanDeltaTs.clear();
        meanTimeOriginSTDs.clear();
        meanProportions.clear();
        varianceDeltaTs.clear();
        varianceTs.clear();
        varianceProportions.clear();
        auto chain = chains.at(i);
        
        std::vector<long double > vec =chain->thetaPar->getParameter()->getChainValues();
        //assert(std::equal(chain->posteriorTheta.begin(), chain->posteriorTheta.end(), vec.begin()));
        
        average = chain->thetaPar->getParameter()->getCurrentMean(); // must be equal Utils::mean(vec);
        //assert(average== chain->thetaPar->getParameter()->getCurrentMean());
        
        n= vec.size();
        meanChainThetas.push_back(average  );
        
        currentVariance= chain->thetaPar->getParameter()->getCurrentVariance(); //Utils::variance(vec);
        varianceChainTheta.push_back(currentVariance );
        //assert(currentVariance== chain->thetaPar->getParameter()->getCurrentVariance());
        
        for(size_t j=0; j <programOptions.numClones;++j)
        {
            auto popI =  chain->getPopulationbyIndex(j);
            
            //assert(popI->posteriorDeltaT.size()==n);
            assert(popI->DeltaT->getParameter()->getChainValues().size()==n);
            average = popI->DeltaT->getParameter()->getCurrentMean();//Utils::mean( popI->DeltaT->getParameter()->getChainValues());
            meanDeltaTs.push_back(average);
            varianceDeltaTs.push_back(Utils::variance(popI->DeltaT->getParameter()->getChainValues()));
            
            assert(popI->TimeOriginSTD->getParameter()->getChainValues().size()==n);
            average = popI->TimeOriginSTD->getParameter()->getCurrentMean();//Utils::mean( popI->TimeOriginSTD->getParameter()->getChainValues());
            meanTimeOriginSTDs.push_back(average);
            //varianceTs.push_back(Utils::variance(popI->TimeOriginSTD->getParameter()->getChainValues()));
            varianceTs.push_back(popI->TimeOriginSTD->getParameter()->getCurrentVariance());
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
    
    
    std::vector<long double > potentialScaleReductionsDeltaT =Utils::potentialScaleReductionArray(numClones,n,  m, meanChainDeltaTs, varianceChainDeltaT );
    
    listPotentialScaleReductions.insert(listPotentialScaleReductions.end(), potentialScaleReductionsDeltaT.begin(), potentialScaleReductionsDeltaT.end());
    
    std::vector<long double > potentialScaleReductionsT = Utils::potentialScaleReductionArray(numClones,n,  m, meanChainTimeOriginSTDs, varianceChainT );
    
    listPotentialScaleReductions.insert(listPotentialScaleReductions.end(), potentialScaleReductionsT.begin(), potentialScaleReductionsT.end());
    
    if (programOptions.numClones >1)
    {
        std::vector<long double > potentialScaleReductionsProportion =Utils::potentialScaleReductionArray(numClones,n,  m, meanChainProportions, varianceChainProportion );
        listPotentialScaleReductions.insert(listPotentialScaleReductions.end(), potentialScaleReductionsProportion.begin(), potentialScaleReductionsProportion.end());
    }
    return listPotentialScaleReductions;
}

void ChainManager::printGroupChainsSummary(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions)
{
    std::vector<Chain *> groupChains= getChainsWiththeSameTrueTreeThatAnotherChain(chainNumber, numberChainsPerTree);
    for(size_t i=0; i <groupChains.size();++i)
    {
        auto chain = groupChains.at(i);
        chain->printLastMovesSummary();
    }
}
void ChainManager::deleteGroupChains(int chainNumber, int numberChainsPerTree)
{
    std::vector<Chain *> groupChains= getChainsWiththeSameTrueTreeThatAnotherChain(chainNumber, numberChainsPerTree);
    for(size_t i=0; i <groupChains.size();++i){
        auto chain = groupChains.at(i);
        delete chain;
    }
}
bool ChainManager::checkGroupChainsHaveCompletedIndependentIterations(int numberIterations, int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions)
{
    bool result=false;
    std::vector<Chain *> chainsFortheSameTree;
    int idxGroup=0;
    long double ESS;
    long double ESSThreshold = 200;
    for(size_t i=0; i <numberChainsPerTree;++i)
    {
        idxGroup= chainNumber / numberChainsPerTree;
        //result = chains.at(idxGroup + i)->numberIndependentLongUpdates >= numberIterations;
        
        ESS= chains.at(idxGroup + i)->thetaPar->ESS();
        result = (ESS > ESSThreshold) ;
        if (!result){
            std::cout << "\n The chain "
            << idxGroup + i <<  " has not have reached ESS = 200: " << ESS << " for theta" <<  std::endl;
            break;
        }
        //std::cout << "\n The number of independent theta samples for chain " <<  idxGroup + i <<  " is " << chains.at(idxGroup + i)->numberIndependentLongUpdates<< " and ESS " <<  ESS << std::endl;
        
        for(size_t j=0; j <programOptions.numClones;++j)
        {
            auto popI =  chains.at(idxGroup + i)->getPopulationbyIndex(j);
            ESS= popI->DeltaT->ESS();
            result = result & (ESS > ESSThreshold );
            //std::cout << "\n The number of independent Delta samples for population with order "<< popI->order << " in chain " <<  idxGroup + i <<  " is " << chains.at(idxGroup + i)->numberIndependentLongUpdates<< " and ESS " <<  ESS << std::endl;
            
            if (!result){
                std::cout << "\n The chain "
                << idxGroup + i <<  " has not have reached ESS =200: " << ESS << "   for Delta in population " << popI->order <<  std::endl;
                break;
            }
            
            ESS= popI->TimeOriginSTD->ESS();
            result = result & (ESS > ESSThreshold );
            
            // std::cout << "\n The number of independent time origin STD samples for  population with order "<< popI->order << " in chain " <<  idxGroup + i <<  " is " << chains.at(idxGroup + i)->numberIndependentLongUpdates<< " and ESS " <<  ESS << std::endl;
            
            if (!result){
                std::cout << "\n The chain "
                << idxGroup + i <<  " has not have reached ESS =200: " << ESS << "   for timeorigin STD in population " << popI->order <<  std::endl;
                break;
            }
        }
        if (!result){
            std::cout << "\n The chain "
            << idxGroup + i <<  " has not have reached ESS = 200" <<  std::endl;
            break;
        }
        
    }
    return result;
}
bool ChainManager::checkEffectiveSampleSizeGroupChains( int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions)
{
    bool result=true;
    long double ESS;
    
    std::vector<Chain *> groupChains= getChainsWiththeSameTrueTreeThatAnotherChain(chainNumber, numberChainsPerTree);
    
    std::vector<long double> ESSs = listEffectiveSampleSizes( chainNumber,  numberChainsPerTree,programOptions);
    
    for(size_t i=0; i < ESSs.size();++i)
    {
        ESS= ESSs.at(i);
        if (ESS > 0){
            // groupChains[0]->get
            std::cout << "\n The current combined ESS for parameter " << groupChains[0]->getParameterNameByIndex(i) << "  in chains " << chainNumber << " to " << (chainNumber+numberChainsPerTree-1) << " is  " << ESS <<  std::endl;
        }
        
        if (ESS >0 && ESS < 400 ){
            result=false;
            return result;
        }
    }
    return result;
}
void ChainManager::setConvergedGroupChains(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions){
    
    std::vector<Chain *> groupChains= getChainsWiththeSameTrueTreeThatAnotherChain(chainNumber, numberChainsPerTree);
    for(size_t i=0; i <groupChains.size();++i)
    {
        auto chain = groupChains.at(i);
        chain->converged = true;
    }
}
void ChainManager::closeFilesGroupChain(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions)
{
    
    std::vector<Chain *> groupChains= getChainsWiththeSameTrueTreeThatAnotherChain(chainNumber, numberChainsPerTree);
    for(size_t i=0; i <groupChains.size();++i)
    {
        auto chain = groupChains.at(i);
        chain->closeFiles(programOptions);
    }
}
bool ChainManager::checkConvergence(int numberChainsPerTree,ProgramOptions &programOptions,MCMCoptions &mcmcOptions, long double threshold )
{
    bool checkPSRF=false;
    bool checkESS=false;
    bool checksplitPSRF=false;
    for(size_t chainNumber=0; chainNumber< numberChains; chainNumber += numberChainsPerTree)
    {
        checkPSRF=false;
        checkESS=false;
        checksplitPSRF=false;
        if (chains[chainNumber]->converged)
            continue;
        
        std::vector<Chain *> groupChains= getChainsWiththeSameTrueTreeThatAnotherChain(chainNumber, numberChainsPerTree);
        
        if (groupChainsFinished(groupChains))
            continue;
        
        bool checkIndividualChainESS=checkGroupChainsHaveCompletedIndependentIterations(mcmcOptions.maxNumberIndependentPosteriorValues, chainNumber,numberChainsPerTree,programOptions);
        
        //if (checkESS)
        if (checkIndividualChainESS)
        {
            std::cout << "\nMCMC  chains: "<< chainNumber << " to "<< (chainNumber +numberChainsPerTree-1) << " have individual ESS greater than 200 .... \n"<< std::endl;
            
            checksplitPSRF= checksplitPSRFGroupChains(chainNumber, numberChainsPerTree,programOptions, threshold );
            
            checkESS = checkEffectiveSampleSizeGroupChains(chainNumber,  numberChainsPerTree, programOptions);
            
            
            if (checkPSRF && checkESS )
            {
                //checkESS = checkEffectiveSampleSizeGroupChains(chainNumber,  numberChainsPerTree, programOptions);
                
                checkPSRF= checkPSRFGroupChains(chainNumber, numberChainsPerTree,programOptions, threshold );
                if (checkPSRF && !(groupChains[0]->converged)){
                    setConvergedGroupChains(chainNumber, numberChainsPerTree,  programOptions);
                    numberConvergedChains += mcmcOptions.numberChainsPerTree;
                    //                     for(size_t j=0; j< numberChainsPerTree; j++)
                    //                           {
                    //                                       setWarmup(chainNumber+j, chains[chainNumber+j]->currentNumberIerations);
                    //                               
                    //                            }
                    printGroupChainsSummary(chainNumber, mcmcOptions.numberChainsPerTree,programOptions);
                    std::string priorsType;
                    if (mcmcOptions.priorsType==0)
                        priorsType="log uniform";
                    else if(mcmcOptions.priorsType==1)
                        priorsType="exponential";
                    else
                        priorsType="power law";
                    std::string kernelType;
                    if (mcmcOptions.kernelType==0)
                        kernelType="multiplier";
                    else //(mcmcOptions.priorsType==1)
                        kernelType="normal";
                    std::cout << "\n Use of " << priorsType.c_str() << "priors and " << kernelType.c_str() << std::endl;
                    std::cout << "\nMCMC  chains: "<< chainNumber << " to "<< (chainNumber +numberChainsPerTree-1) << " started to sample from the posterior .... \n"<< std::endl;
                }
                
            }// if (checkPSRF)
            //closeFilesGroupChain(chainNumber, mcmcOptions.numberChainsPerTree, programOptions);
            // deleteGroupChains(chainNumber, mcmcOptions.numberChainsPerTree);
            //break;
        }// if (checkIndividualChainESS)
    }// for
    return(checkPSRF);
}
bool ChainManager::groupChainsFinished(std::vector<Chain *> groupChains){
    
    bool result = true;
    for (int idx = 0; idx < groupChains.size(); idx++) {
        if (!(groupChains[idx]->endReached))
            return false;
        
        
    }
    return result;
}
//void ChainManager::addSamples(){
//    std::vector<const double*> draws(chains.size());
//    std::vector<size_t> sizes(chains.size());
//    Chain* currentChain;
//   // for (int index = 4; index < chains.num_params(); index++) {
//      for (int chain = 0; chain < chains.size(); ++chain) {
//          currentChain = chains[chain];
//
//        samples_(chain) = currentChain->getSample();
//        draws[chain] = &samples_(chain)(0);
//        sizes[chain] = samples_(chain).size();
//      }
//}
void ChainManager::performWarmUp(ProgramOptions &programOptions,MCMCoptions &mcmcOptions, std::vector<gsl_rng * > &randomGenerators,FilePaths &filePaths)
{
    //gsl_rng * currentRandomGenerator;
#pragma omp parallel   shared(programOptions,mcmcOptions,filePaths,   randomGenerators )
    {
#pragma omp for // nowait
        for(int chainNumber=0; chainNumber< numberChains;++chainNumber)
        {
            // gsl_rng_set(currentRandomGenerator, omp_get_thread_num() * 12567);
            for (int currentIteration = 0; currentIteration < (mcmcOptions.numberWarmUpIterations); ++currentIteration)
            {
                if (currentIteration % mcmcOptions.printChainStateEvery == 0 || currentIteration >= (mcmcOptions.Niterations-1))
                    std::cout << "\n Chain " << chains.at(chainNumber)->chainNumber << " Starting Iteration " << currentIteration +1 << std::endl;
                chains.at(chainNumber)->currentNumberIerations =currentIteration;
                
                chains.at(chainNumber)->stepAllMoves(mcmcOptions,  randomGenerators.at(chainNumber),  programOptions);
                
                chains.at(chainNumber)->saveMCMCState( currentIteration, programOptions,mcmcOptions);
                
                if (currentIteration>0 && currentIteration %  mcmcOptions.iterationsToMonitorChain ==0 && currentIteration < mcmcOptions.numberWarmUpIterations )
                {
                    if (mcmcOptions.doThinning)
                    {
                        //  addSample(chainNumber, chains.at(chainNumber)->getStoredSamples());
                        chains.at(chainNumber)->computeThinnig( mcmcOptions);
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
void ChainManager::addSample(const int chain, const Eigen::MatrixXd& sample) {
    
    int numChains = samples_.size();
    
    int num_params = numMCMCParametersChain;
    if (numChains == 0 || chain >= numChains) {
        int n = numChains;
        
        Eigen::Matrix<Eigen::MatrixXd, Eigen::Dynamic, 1> samples_copy(numChains);
        Eigen::VectorXi warmup_copy(numChains);
        for (int i = 0; i < n; i++) {
            samples_copy(i) = samples_(i);
            warmup_copy(i) = warmup_(i);
        }
        samples_.resize(chain + 1);
        warmup_.resize(chain + 1);
        for (int i = 0; i < n; i++) {
            samples_(i) = samples_copy(i);
            warmup_(i) = warmup_copy(i);
        }
        for (int i = n; i < chain + 1; i++) {
            samples_(i) = Eigen::MatrixXd(0, num_params);
            warmup_(i) = 0;
        }
    }
    int row = samples_(chain).rows();
    Eigen::MatrixXd new_samples(row + sample.rows(), num_params);
    new_samples << samples_(chain), sample;
    samples_(chain) = new_samples;
}
void ChainManager::addSample(const std::vector<std::vector< double> >& sample) {
    int n_row = sample.size();
    if (n_row == 0)
        return;
    int n_col = sample[0].size();
    Eigen::MatrixXd sample_copy(n_row, n_col);
    for (int i = 0; i < n_row; i++) {
        sample_copy.row(i)
        = Eigen::VectorXd::Map(&sample[i][0], sample[0].size());
    }
    addSample(sample_copy);
}
void ChainManager::addSample(const Eigen::MatrixXd& sample) {
    if (sample.rows() == 0)
        return;
    if (sample.cols() != numMCMCParametersChain)
        throw std::invalid_argument(
                                    "addSample(sample): number of columns in"
                                    " sample does not match chains");
    addSample(samples_.size(), sample);
}
Eigen::VectorXd ChainManager::samples(const int chain, const int index, int numberSamples) const{
    return samples_(chain).col(index).bottomRows(numberSamples);
}

double ChainManager::splitPotentialScaleReductionFactor(const int index, std::vector<Chain *> chainsWithSameTree, int numberChainsPerTree) const {
    assert(chainsWithSameTree.size() ==numberChainsPerTree);
    int n_chains = chainsWithSameTree.size();
    
    std::vector<const double*> draws(n_chains);
    std::vector<size_t> sizes(n_chains);
    int n_kept_samples = 0;
    for (int chain = 0; chain < n_chains; ++chain) {
        n_kept_samples = numberSamplesAfterWarmUp(chain, index);
        if (n_kept_samples <= 50)
            return 1.0;
        draws[chain]= getLastChainSamples(chain, index, n_kept_samples);
        //= samples_(chain).col(index).bottomRows(n_kept_samples).data();
        sizes[chain] = n_kept_samples;
    }
    
    return mcmc_utils::compute_split_potential_scale_reduction(draws, sizes);
}

double ChainManager::splitPotentialScaleReductionFactor(const std::string& name, std::vector<Chain *> chainsWithSameTree, int numberChainsPerTree) const {
    
    assert(chainsWithSameTree.size() ==numberChainsPerTree);
    int index =-1;
    if (chainsWithSameTree.size() >0)
        index = chainsWithSameTree[0]->index(name);
    assert (index !=-1);
    return splitPotentialScaleReductionFactor(index,  chainsWithSameTree,  numberChainsPerTree);
    
}
double ChainManager::effectiveSampleSize(const int index, std::vector<Chain *> chainsWithSameTree, int numberChainsPerTree) const {
    
    assert(chainsWithSameTree.size() ==numberChainsPerTree);
    int n_chains = chainsWithSameTree.size();
    std::vector<const double*> draws(n_chains);
    std::vector<size_t> sizes(n_chains);
    int n_kept_samples = 0;
    for (int chain = 0; chain < n_chains; ++chain) {
        n_kept_samples = numberSamplesAfterWarmUp(chain, index);
        if (n_kept_samples <= 50)
            return 0.0;
        draws[chain] = getLastChainSamples(chain, index, n_kept_samples);
        // = samples_(chain).col(index).bottomRows(n_kept_samples).data();
        sizes[chain] = n_kept_samples;
    }
    return mcmc_utils::compute_effective_sample_size(draws, sizes);
}

double ChainManager::effectiveSampleSize(const std::string& name, std::vector<Chain *> chainsWithSameTree, int numberChainsPerTree) const {
    
    assert(chainsWithSameTree.size() ==numberChainsPerTree);
    int index =-1;
    if (chainsWithSameTree.size() >0)
        index = chainsWithSameTree[0]->index(name);
    assert (index !=-1);
    return effectiveSampleSize(index, chainsWithSameTree, numberChainsPerTree);
}

double ChainManager::splitEffectiveSampleSize(const int index, std::vector<Chain *> chainsWithSameTree, int numberChainsPerTree) const {
    
    assert(chainsWithSameTree.size() ==numberChainsPerTree);
    int n_chains = chainsWithSameTree.size();
    std::vector<const  double*> draws(n_chains);
    std::vector<size_t> sizes(n_chains);
    int n_kept_samples = 0;
    for (int chain = 0; chain < n_chains; ++chain) {
        n_kept_samples = numberSamplesAfterWarmUp(chain, index);
        if (n_kept_samples <= 50)
            return 0.0;
        draws[chain]= getLastChainSamples(chain, index, n_kept_samples);
        //  = samples_(chain).col(index).bottomRows(n_kept_samples).data();
        sizes[chain] = n_kept_samples;
    }
    return mcmc_utils::compute_split_effective_sample_size(draws, sizes);
}

double ChainManager::splitEffectiveSampleSize( const std::string& name,  std::vector<Chain *> chainsWithSameTree, int numberChainsPerTree) const {
    
    assert(chainsWithSameTree.size() ==numberChainsPerTree);
    int index =-1;
    if (chainsWithSameTree.size() >0)
        index = chainsWithSameTree[0]->index(name);
    assert (index !=-1);
    return splitEffectiveSampleSize(index, chainsWithSameTree,numberChainsPerTree );
}
std::vector<long double> ChainManager::listPotentialScaleReductionFactors(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions)
{
    std::vector<Chain *> groupChains= getChainsWiththeSameTrueTreeThatAnotherChain(chainNumber, numberChainsPerTree);
    
    std::vector<long double> listPotentialScaleReductions;
    double PSRF= 0.0;
    
    for (int i = 0; i < numMCMCParametersChain; ++i) {
        PSRF = splitPotentialScaleReductionFactor(i,groupChains,  numberChainsPerTree);
        listPotentialScaleReductions.push_back(PSRF);
    }
    return listPotentialScaleReductions;
    
}
std::vector<long double> ChainManager::listEffectiveSampleSizes(int chainNumber, int numberChainsPerTree,ProgramOptions &programOptions)
{
    std::vector<Chain *> groupChains= getChainsWiththeSameTrueTreeThatAnotherChain(chainNumber, numberChainsPerTree);
    
    std::vector<long double> listPotentialScaleReductions;
    double ESS= 0.0;
    for (int i = 0; i < numMCMCParametersChain; ++i) {
        ESS = effectiveSampleSize(i,groupChains,  numberChainsPerTree);
        listPotentialScaleReductions.push_back(ESS);
    }
    return listPotentialScaleReductions;
    
}
void ChainManager::resizeStoredMCMCparametersAfterWarmUp(MCMCoptions &mcmcOptions){
    
    for(int chainNumber=0; chainNumber< numberChains;++chainNumber){
        chains.at(chainNumber)->resizeStoredMCMCparametersAfterWarmUp(mcmcOptions);
    }
}
double* ChainManager::getLastChainSamples(const int chainId, const int idxMCMCParameter, int numberSamples) const{
    
    assert(chainId >=0 && chainId <chains.size());
    std::vector<long double> samples= chains[chainId]->getSampleByIndex(idxMCMCParameter);
    assert(numberSamples <=samples.size());
    std::vector<double> result(samples.end() - numberSamples, samples.end());
    return &result[0];
    
}
void  ChainManager::finish(int numberChainsPerTree, ProgramOptions &programOptions, MCMCoptions &mcmcOptions){
    
    for(size_t chainNumber=0; chainNumber< numberChains; chainNumber += numberChainsPerTree){
        //closeFilesGroupChain(chainNumber, mcmcOptions.numberChainsPerTree, programOptions);
        deleteGroupChains(chainNumber, mcmcOptions.numberChainsPerTree);
    }
}
