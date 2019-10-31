//
//  mcmc_move.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 28/10/19.
//

#include "mcmc_move.hpp"
#include "data_utils.hpp"
#include "random.h"
MCMCmove::MCMCmove(Chain *chain,  string nameMove)
{
    this->chain = chain;
    this->nameMove = nameMove;
}
Chain * MCMCmove::getChain()
{
    return chain;
}
NewTotalEffectPopSizeMove::NewTotalEffectPopSizeMove(Chain * chain, string nameMove ):MCMCmove(chain, nameMove)
{

}
NewProportionsVectorMove::NewProportionsVectorMove(Chain * chain, string nameMove ):MCMCmove(chain, nameMove)
{
    
}
NewEffectPopSizeMoveForPopulation::NewEffectPopSizeMoveForPopulation(Chain *chain, string nameMove, Population *pop):MCMCmove(chain, nameMove)
{
    
}
void NewTotalEffectPopSizeMove::safeCurrentValue()
{
      Chain *chain=getChain();
    int i=0;
     Population *popI;
    chain->oldTotalEffectPopSize = chain->totalEffectPopSize;
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

  
    mcmcOptions.slidingWindowSizeTotalEffectPopSize= 2* (chain->totalEffectPopSize - mcmcOptions.fixedLambda * chain->totalSampleSize());
    
    bool allPopulationPopSizesSet=false;
    do {
        newTotalPopulationSize= chain->proposalSlidingWindow(chain->oldTotalEffectPopSize,  mcmcOptions.slidingWindowSizeTotalEffectPopSize);
        chain->totalEffectPopSize= newTotalPopulationSize;
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
    chain->totalEffectPopSize= chain->oldTotalEffectPopSize ;
     int i=0;
    
    Population *popI;
    chain->oldTotalEffectPopSize = chain->totalEffectPopSize;
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
    
    fprintf (stderr, "\n>> log conditional Likelihood tree for proposal new total effective population size,  of the chain %d is = %lf  \n", chain->chainNumber,newLogConditionalLikelihoodTree );
  
    double priorDensityNewTotalEffectivePopulationSize= LogUniformDensity(chain->totalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    double priorDensityCurrentTotalEffectivePopulationSize= LogUniformDensity(chain->oldTotalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewTotalEffectivePopulationSize;
    
    double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensityCurrentTotalEffectivePopulationSize;
    
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
        printf("\n Accepted new move %s", nameMove.c_str());
    }
    else
    {
        rollbackMove();
    }
}
//////////////////////////////////////////////
void NewProportionsVectorMove::safeCurrentValue()
{
    Chain *chain=getChain();
    int i=0;
    Population *popI;
     chain->oldproportionsVector= chain->proportionsVector;
    for( i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        popI->oldeffectPopSize=popI->effectPopSize;
        popI->oldPopSize = popI->popSize;
        popI->oldDeathRate = popI->deathRate;
        popI->oldGrowthRate = popI->growthRate;
    }
}
void NewProportionsVectorMove::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    int i=0;
    Population *popI;

    safeCurrentValue();
    double *proportionsVectorArray;
    
    bool allPopulationPopSizesSet=false;
    do {
        proportionsVectorArray=&(chain->oldproportionsVector[0]);
         //init array with the current proportions vector
        randomDirichletFromGsl(chain->numClones, proportionsVectorArray, &(chain->proportionsVector[0]));
        
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
void NewProportionsVectorMove::rollbackMove()
{
    Chain *chain=getChain();
    chain->proportionsVector= chain->oldproportionsVector ;
    int i=0;
    Population *popI;
    for( i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        popI->effectPopSize=popI->oldeffectPopSize;
        popI->popSize = popI->oldPopSize;
        popI->deathRate = popI->oldDeathRate;
        popI->growthRate = popI->oldGrowthRate;
    }
}
double NewProportionsVectorMove::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
   // ListClonesAccordingTimeToOrigin(chain->populations);
    double newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions);
    double priorDensityNewProportionsVector =  chain->DirichletDensity( chain->proportionsVector, chain->oldproportionsVector , chain->numClones);
    
    double priorDensityCurrentProportionsVector = chain->DirichletDensity(chain->oldproportionsVector, chain->proportionsVector, chain->numClones);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewProportionsVector;
    
    double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensityCurrentProportionsVector;
    
    double logAcceptanceRate = (sumLogNumerators - sumLogDenominators) >0? (sumLogNumerators - sumLogDenominators) :0;
    
    return(logAcceptanceRate);
}
NewGrowthRateMoveForPopulation::NewGrowthRateMoveForPopulation(Chain *chain,string nameMove,  Population *pop):MCMCmove(chain, nameMove)
{
    this->pop =pop;
}
void NewGrowthRateMoveForPopulation::safeCurrentValue()
{
   pop->olddelta= pop->delta;
}
void NewGrowthRateMoveForPopulation::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    safeCurrentValue();
    double randomDelta = RandomLogUniform(mcmcOptions.Deltafrom, mcmcOptions.Deltato);
    pop->delta =  randomDelta;
    pop->growthRate =pop->delta  / pop->effectPopSize;
    pop->deathRate= pop->birthRate - pop->growthRate;
}
void NewGrowthRateMoveForPopulation::rollbackMove()
{
  pop->delta= pop->olddelta;
}
double NewGrowthRateMoveForPopulation::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    double newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions);
    
    fprintf (stderr, "\n>> New log conditional Likelihood tree after the new total growth rate for population %d,  of the chain %d is = %lf  \n", pop->index, chain->chainNumber, newLogConditionalLikelihoodTree );
    
    double priorDensityNewScaledGrowthRate =  LogUniformDensity(pop->delta, mcmcOptions.Deltafrom, mcmcOptions.Deltato);
    
    double priorDensitScaledGrowthRate =  LogUniformDensity(pop->olddelta, mcmcOptions.Deltafrom, mcmcOptions.Deltato);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewScaledGrowthRate;
    
    double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensitScaledGrowthRate;
    
    
    double logAcceptanceRate = (sumLogNumerators - sumLogDenominators) >0? (sumLogNumerators - sumLogDenominators) :0;
   
    return logAcceptanceRate;
}
void NewEffectPopSizeMoveForPopulation::safeCurrentValue()
{    Chain *chain=getChain();
    pop->oldeffectPopSize= pop->effectPopSize;
    chain->oldTotalEffectPopSize =chain->totalEffectPopSize;
}
void NewEffectPopSizeMoveForPopulation::rollbackMove()
{  Chain *chain=getChain();
   pop->effectPopSize= pop->oldeffectPopSize;
   chain->totalEffectPopSize =chain->oldTotalEffectPopSize;
}
void NewEffectPopSizeMoveForPopulation::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    int i;
    Population *popI;
    Chain *chain=getChain();
    safeCurrentValue();
    bool populationPopSizesSet=false;
    double newEffectPopulationSize;
    int totalSampleSize = chain->totalSampleSize();
    
    double slidingWindowSize= 2* (pop->oldeffectPopSize - pop->birthRate * pop->sampleSize);
    do {
        newEffectPopulationSize= chain->proposalSlidingWindow(pop->oldeffectPopSize,  slidingWindowSize);
        
        pop->effectPopSize= newEffectPopulationSize;
        
        chain->totalEffectPopSize =  chain->totalEffectPopSize + pop->effectPopSize - pop->oldeffectPopSize;
       
        chain->updateEffectPopSizesCurrentProportionsVector();
        pop=chain->populations[i];
        pop->growthRate =pop->delta  / pop->effectPopSize;
        pop->popSize=pop->effectPopSize * pop->birthRate;
        pop->deathRate= pop->birthRate - pop->growthRate;
        if (pop->popSize < pop->sampleSize)
        {
                populationPopSizesSet=false;
        }
    }
    while(!populationPopSizesSet);
}
double NewEffectPopSizeMoveForPopulation::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions){
    
    Chain *chain=getChain();
    double newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions);
    
     double slidingWindowSize= 2* (pop->oldeffectPopSize - pop->birthRate * pop->sampleSize);
    
    double priorDensityNewEffectivePopulationSize= LogUniformDensity(pop->effectPopSize,  mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    
    double priorDensityCurrentTotalEffectivePopulationSize= LogUniformDensity(chain->oldTotalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewEffectivePopulationSize;
    
    double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensityCurrentTotalEffectivePopulationSize;
    
    double LogAcceptanceRate = (sumLogNumerators - sumLogDenominators) >0? (sumLogNumerators - sumLogDenominators) :0;
    
    return LogAcceptanceRate;
}
void NewTotalEffectPopSizeMove::move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    MCMCmove::move(programOptions,mcmcOptions);
}
void NewProportionsVectorMove::move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    MCMCmove::move(programOptions,mcmcOptions);
}
void NewGrowthRateMoveForPopulation::move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    MCMCmove::move(programOptions,mcmcOptions);
}
void NewEffectPopSizeMoveForPopulation::move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    MCMCmove::move(programOptions,mcmcOptions);
}
void NewTimeOriginOnTreeforPopulationMove::move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    MCMCmove::move(programOptions,mcmcOptions);
}
NewTimeOriginOnTreeforPopulationMove::NewTimeOriginOnTreeforPopulationMove(Chain *chain,string nameMove,  Population *pop):MCMCmove(chain, nameMove)
{
    this->pop =pop;
}
void NewTimeOriginOnTreeforPopulationMove::safeCurrentValue()
{
    Chain *chain=getChain();
    int i=0;
    Population *popI;
    //save information for the population
    pop->oldrMRCA =  pop->rMRCA ;
    pop->oldFatherPop =  pop->FatherPop;
    pop->oldCoalescentEventTimes = pop->CoalescentEventTimes;
    pop->oldimmigrantsPopOrderedByModelTime= pop->immigrantsPopOrderedByModelTime;
     pop->oldSampleSize = pop->sampleSize;
    
     //save information for the other populations
    for( i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        
        pop->oldFatherPop =  pop->FatherPop;
        pop->oldCoalescentEventTimes = pop->CoalescentEventTimes;
        pop->oldimmigrantsPopOrderedByModelTime= pop->immigrantsPopOrderedByModelTime;
    }
}
void NewTimeOriginOnTreeforPopulationMove::rollbackMove()
{
    Chain *chain=getChain();
    int i=0;
    Population *popI;
    //save information for the population
    pop->rMRCA = pop->oldrMRCA ;
    pop->FatherPop = pop->oldFatherPop ;
    pop->CoalescentEventTimes = pop->oldCoalescentEventTimes ;
    pop->immigrantsPopOrderedByModelTime = pop->oldimmigrantsPopOrderedByModelTime;
    pop->sampleSize = pop->oldSampleSize;
    
    //save information for the other populations
    for( i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
         pop->FatherPop = pop->oldFatherPop ;
        pop->CoalescentEventTimes = pop->oldCoalescentEventTimes ;
        pop->immigrantsPopOrderedByModelTime = pop->oldimmigrantsPopOrderedByModelTime;
    }
}
void NewTimeOriginOnTreeforPopulationMove::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
     Chain *chain=getChain();
    //first compute the list of  edges without events
    //then select at random an edge
    //recompute the sample sizes
    //recompute the proportions vector
    // recompute the coalescent and migration events
    bool existsZeroSampleSizePop=false;
    double alpha[chain->numClones];
    int numberPoints = chain->numClones -1;
    std::map<pll_rnode_t*, Population*>  rmrcaOfPopulation;

    rmrcaOfPopulation=  chain->chooseAvailableEdgeOnRootedTreeForPopulation(pop, rmrcaOfPopulation, programOptions.healthyTipLabel);
    
    for (unsigned int i = 0; i < chain->numClones; ++i){
        auto pop =  chain->populations[i];
        alpha[i]= pop->sampleSize;
    }
    
    
    int totalSampleSize=chain->initialRootedTree->tip_count-1;//not the healthytip
    std::transform(alpha, alpha + chain->numClones , alpha,std::bind2nd(std::divides<double>(),totalSampleSize));
    chain->initProportionsVector();
    chain->generateProportionsVectorFromDirichlet(alpha);
    chain->initEffectPopulationSizesFromProportionsVector();
    chain->initTimeOriginSTD();
    chain->initPopulationMigration();//after setting the timeSTD
    chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(rmrcaOfPopulation, programOptions.healthyTipLabel);
    chain->filterSortPopulationsCoalescentEvents();
}
double  NewTimeOriginOnTreeforPopulationMove::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    
    double newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions);
    
    double priorDensityNewTotalEffectivePopulationSize=0;
//    priorDensityNewTotalEffectivePopulationSize = LogUniformDensity(chain->totalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    double priorDensityCurrentTotalEffectivePopulationSize=0;
//    priorDensityCurrentTotalEffectivePopulationSize= LogUniformDensity(chain->oldTotalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewTotalEffectivePopulationSize;
    
    double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensityCurrentTotalEffectivePopulationSize;
    
    double LogAcceptanceRate = (sumLogNumerators - sumLogDenominators) >0? (sumLogNumerators - sumLogDenominators) :0;
    
    return LogAcceptanceRate;
    
    
}
