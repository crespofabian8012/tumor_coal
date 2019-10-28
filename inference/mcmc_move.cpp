//
//  mcmc_move.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 28/10/19.
//

#include "mcmc_move.hpp"
#include "data_utils.hpp"
#include "random.h"
MCMCmove::MCMCmove(Chain *chain)
{
    this->chain = chain;
}
Chain * MCMCmove::getChain()
{
    return chain;
}

void NewTotalEffectPopSizeMove::safeCurrentValue()
{
      Chain *chain=getChain();
    int i=0;
     Population *popI;
    chain->oldtotalPopSize = chain->totalPopSize;
    for( i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        popI->oldeffectPopSize=popI->effectPopSize;
        popI->oldPopSize = popI->popSize;
        popI->oldDeathRate = popI->deathRate;
        popI->oldGrowthRate = popI->growthRate;

    }
}
void NewTotalEffectPopSizeMove::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    int i=0;
    Population *popI;
    Chain *chain=getChain();
    safeCurrentValue();
    double newTotalPopulationSize;

    bool allPopulationPopSizesSet=false;
    do {
        newTotalPopulationSize= chain->proposalSlidingWindow(chain->oldtotalPopSize,  mcmcOptions.slidingWindowSizeTotalEffectPopSize);
        chain->totalPopSize= newTotalPopulationSize;
        chain->updateEffectPopSizesCurrentProportionsVector();
        for( i = 0 ; i < chain->numClones; i++)
        {
          popI=chain->populations[i];
          popI->growthRate =popI->delta  / popI->effectPopSize;
          popI->popSize=popI->effectPopSize * popI->birthRate;
          popI->deathRate= popI->birthRate - popI->growthRate;
            if (popI->popSize < popI->sampleSize)
            {
                allPopulationPopSizesSet=false;
                break;
            }
        }
        allPopulationPopSizesSet=true;
    }
    while(!allPopulationPopSizesSet);
}
void NewTotalEffectPopSizeMove::rollbackMove()
{
    Chain *chain=getChain();
    chain->totalPopSize= chain->oldtotalPopSize ;
     int i=0;
    
    Population *popI;
    chain->oldtotalPopSize = chain->totalPopSize;
    for( i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        popI->effectPopSize=popI->oldeffectPopSize;
        popI->popSize = popI->oldPopSize;
        popI->deathRate = popI->oldDeathRate;
        popI->growthRate = popI->oldGrowthRate;
    }
}
double NewTotalEffectPopSizeMove::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();

    double newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions);
    
    fprintf (stderr, "\n>> New log conditional Likelihood tree after the new total effective population size,  of the chain %d is = %lf  \n", chain->chainNumber,newLogConditionalLikelihoodTree );
  
    double priorDensityNewTotalEffectivePopulationSize= LogUniformDensity(chain->totalPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    double priorDensityCurrentTotalEffectivePopulationSize= LogUniformDensity(chain->oldtotalPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewTotalEffectivePopulationSize;
    
    double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensityCurrentTotalEffectivePopulationSize;
    
    double randomNumber= randomUniformFromGsl();
    
    double LogAcceptanceRate = (sumLogNumerators - sumLogDenominators) >0? (sumLogNumerators - sumLogDenominators) :0;
  
    return LogAcceptanceRate;
}
void MCMCmove::move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    
    safeCurrentValue();
    makeProposal(programOptions, mcmcOptions);
    double logAcceptanceRate =computeLogAcceptanceProb(programOptions, mcmcOptions);
    double randomNumber= randomUniformFromGsl();
    if (log(randomNumber) < logAcceptanceRate )
    {//accept the move
        printf("\n Accepted new move");
    }
    else
    {
        rollbackMove();
    }
    
    
}
//void NewProportionsVectorMove::makeProposal()
//{
//
//
//}
//void NewProportionsVectorMove::rollbackMove()
//{
//
//
//}
//void NewProportionsVectorMove::computeAcceptanceProb()
//{
//
//}
//void NewGrowthRateMove::makeProposal()
//{
//
//
//}
//void NewGrowthRateMove::rollbackMove()
//{
//
//
//}
//void NewGrowthRateMove::computeAcceptanceProb()
//{
//
//
//}

