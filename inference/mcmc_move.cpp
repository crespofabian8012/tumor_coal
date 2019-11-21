//
//  mcmc_move.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 28/10/19.
//

#include "mcmc_move.hpp"
#include "data_utils.hpp"
#include "random.h"
using namespace std;
MCMCmove::MCMCmove(Chain *chain,  string nameMove)
{
    this->chain = chain;
    this->nameMove = nameMove;
    this->numberAccept=0;
    this->numberReject=0;
    this->newLogConditionalLikelihoodTree =chain->currentlogConditionalLikelihoodTree;
    this->newLogConditionalLikelihoodSequences=chain->currentlogConditionalLikelihoodSequences;
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
        popI->oldCoalescentEventTimes = popI->CoalescentEventTimes;
        popI->CoalescentEventTimes.clear();
        popI->oldimmigrantsPopOrderedByModelTime= popI->immigrantsPopOrderedByModelTime;
        popI->immigrantsPopOrderedByModelTime.clear();
    }
  
}
void NewTotalEffectPopSizeMove::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    int i=0;
    Population *popI;
    Chain *chain=getChain();

    double newTotalPopulationSize;

    mcmcOptions.slidingWindowSizeTotalEffectPopSize= 2 * (chain->totalEffectPopSize - (chain->totalSampleSize() / mcmcOptions.fixedLambda));
    
    bool allPopulationPopSizesSet=false;
    do {
        allPopulationPopSizesSet=true;
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
        
    }
    while(!allPopulationPopSizesSet);
    chain->initPopulationMigration();//after setting the timeSTD
    chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
    chain->filterSortPopulationsCoalescentEvents();
}
void NewTotalEffectPopSizeMove::rollbackMove()
{
    Chain *chain=getChain();
    chain->totalEffectPopSize= chain->oldTotalEffectPopSize ;
     int i=0;
    
    Population *popI;
   
    for( i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        popI->effectPopSize=popI->oldeffectPopSize;
        popI->popSize = popI->oldPopSize;
        popI->deathRate = popI->oldDeathRate;
        popI->growthRate = popI->oldGrowthRate;
        popI->CoalescentEventTimes = popI->oldCoalescentEventTimes;
        popI->immigrantsPopOrderedByModelTime =  popI->oldimmigrantsPopOrderedByModelTime;
    }
    
}
double NewTotalEffectPopSizeMove::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();

 newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions);
    
    fprintf (stderr, "\n>> log conditional Likelihood tree for the new total effective population size %d(old %d),  of the chain %d is = %lf  \n",chain->totalEffectPopSize, chain->oldTotalEffectPopSize,chain->chainNumber,newLogConditionalLikelihoodTree );
  
    double priorDensityNewTotalEffectivePopulationSize= LogUniformDensity(chain->totalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    double priorDensityCurrentTotalEffectivePopulationSize= LogUniformDensity(chain->oldTotalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewTotalEffectivePopulationSize;
    
    double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensityCurrentTotalEffectivePopulationSize;
    
    double LogAcceptanceRate = std::min(sumLogNumerators - sumLogDenominators,0.0);
  
    return LogAcceptanceRate;
}

void MCMCmove::move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    safeCurrentValue();
    makeProposal(programOptions, mcmcOptions);
    double logAcceptanceRate =computeLogAcceptanceProb(programOptions, mcmcOptions);
    double randomNumber= randomUniformFromGsl();
    double logRandom=log(randomNumber);
    if (logRandom < logAcceptanceRate )
    {//accept the move
        printf("\n Accepted new move %s \n", nameMove.c_str());
        this->numberAccept++;
        Chain *chain=getChain();
        if (chain->currentlogConditionalLikelihoodTree !=newLogConditionalLikelihoodTree)
             chain->currentlogConditionalLikelihoodTree = newLogConditionalLikelihoodTree;
        if (chain->currentlogConditionalLikelihoodSequences !=newLogConditionalLikelihoodSequences)
               chain->currentlogConditionalLikelihoodSequences = newLogConditionalLikelihoodSequences;
    }
    else
    {
        rollbackMove();
        this->numberReject++;
         printf("\n Rejected new move %s \n", nameMove.c_str());
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
        popI->oldCoalescentEventTimes = popI->CoalescentEventTimes;
        popI->CoalescentEventTimes.clear();
        popI->oldimmigrantsPopOrderedByModelTime= popI->immigrantsPopOrderedByModelTime;
        popI->immigrantsPopOrderedByModelTime.clear();
    }
}
void NewProportionsVectorMove::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    int i=0;
    Population *popI;

    double *proportionsVectorArray;
    
    bool allPopulationPopSizesSet=false;
    do {
        allPopulationPopSizesSet=true;
        //proportionsVectorArray=&(chain->oldproportionsVector[0]);
         //init array with the current proportions vector
        //print old proportions vector
         for( i = 0 ; i < chain->numClones; i++)
             fprintf (stderr, "\n old proportions vector at %d: %lf \n",i,chain->oldproportionsVector.at(i) );
        randomDirichletFromVector (chain->oldproportionsVector, chain->proportionsVector);
//        randomDirichletFromGsl(chain->numClones, proportionsVectorArray, &(chain->proportionsVector[0]));
        
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
    }
    while(!allPopulationPopSizesSet);
    chain->initPopulationMigration();//after setting the timeSTD
    chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
    chain->filterSortPopulationsCoalescentEvents();
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
        popI->CoalescentEventTimes = popI->oldCoalescentEventTimes;
        popI->immigrantsPopOrderedByModelTime =  popI->oldimmigrantsPopOrderedByModelTime;
    }
}
double NewProportionsVectorMove::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
   // ListClonesAccordingTimeToOrigin(chain->populations);
  newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions);
    double priorDensityNewProportionsVector =  chain->DirichletDensity( chain->proportionsVector, chain->oldproportionsVector , chain->numClones);
    
    double priorDensityCurrentProportionsVector = chain->DirichletDensity(chain->oldproportionsVector, chain->proportionsVector, chain->numClones);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewProportionsVector;
    
    double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensityCurrentProportionsVector;
    
    double logAcceptanceRate = std::min(sumLogNumerators - sumLogDenominators,0.0);
    
    return(logAcceptanceRate);
}
NewGrowthRateMoveForPopulation::NewGrowthRateMoveForPopulation(Chain *chain,string nameMove,  Population *pop):MCMCmove(chain, nameMove)
{
    this->pop =pop;
}
void NewGrowthRateMoveForPopulation::safeCurrentValue()
{
     Chain *chain=getChain();
   pop->olddelta= pop->delta;
    Population *popI;
    for( int i = 0 ; i < chain->numClones; i++)
    {
        popI= chain->populations[i];
        popI->oldCoalescentEventTimes = popI->CoalescentEventTimes;
        popI->CoalescentEventTimes.clear();
        popI->oldimmigrantsPopOrderedByModelTime= popI->immigrantsPopOrderedByModelTime;
        popI->immigrantsPopOrderedByModelTime.clear();
    }
}
void NewGrowthRateMoveForPopulation::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();

    double randomDelta = RandomLogUniform(mcmcOptions.Deltafrom, mcmcOptions.Deltato);
    pop->delta =  randomDelta;
    pop->growthRate =pop->delta  / pop->effectPopSize;
    pop->deathRate= pop->birthRate - pop->growthRate;
    
    chain->initPopulationMigration();//after setting the timeSTD
    chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
    chain->filterSortPopulationsCoalescentEvents();
}
void NewGrowthRateMoveForPopulation::rollbackMove()
{
    Chain *chain=getChain();
    Population *popI;
  pop->delta= pop->olddelta;
    for(int i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        popI->CoalescentEventTimes = popI->oldCoalescentEventTimes;
        popI->immigrantsPopOrderedByModelTime =  popI->oldimmigrantsPopOrderedByModelTime;
    }
}
double NewGrowthRateMoveForPopulation::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions);
    
    fprintf (stderr, "\n>> New log conditional Likelihood tree after the new total growth rate for population %d,  of the chain %d is = %lf  \n", pop->index, chain->chainNumber, newLogConditionalLikelihoodTree );
    
    double priorDensityNewScaledGrowthRate =  LogUniformDensity(pop->delta, mcmcOptions.Deltafrom, mcmcOptions.Deltato);
    
    double priorDensitScaledGrowthRate =  LogUniformDensity(pop->olddelta, mcmcOptions.Deltafrom, mcmcOptions.Deltato);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewScaledGrowthRate;
    
    double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensitScaledGrowthRate;
    
    
    double logAcceptanceRate = std::min(sumLogNumerators - sumLogDenominators,0.0);
   
    return logAcceptanceRate;
}
void NewEffectPopSizeMoveForPopulation::safeCurrentValue()
{    Chain *chain=getChain();
    pop->oldeffectPopSize= pop->effectPopSize;
    chain->oldTotalEffectPopSize =chain->totalEffectPopSize;
    Population *popI;
    for( int i = 0 ; i < chain->numClones; i++)
    {
        popI= chain->populations[i];
        popI->oldCoalescentEventTimes = popI->CoalescentEventTimes;
        popI->CoalescentEventTimes.clear();
        popI->oldimmigrantsPopOrderedByModelTime= popI->immigrantsPopOrderedByModelTime;
        popI->immigrantsPopOrderedByModelTime.clear();
    }
    
}
void NewEffectPopSizeMoveForPopulation::rollbackMove()
{  Chain *chain=getChain();
   pop->effectPopSize= pop->oldeffectPopSize;
   chain->totalEffectPopSize =chain->oldTotalEffectPopSize;
    Population *popI;
    for(int i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        popI->CoalescentEventTimes = popI->oldCoalescentEventTimes;
        popI->immigrantsPopOrderedByModelTime =  popI->oldimmigrantsPopOrderedByModelTime;
    }
}
void NewEffectPopSizeMoveForPopulation::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    int i;
 
    Chain *chain=getChain();

    bool populationPopSizesSet=false;
    double newEffectPopulationSize;
    //int totalSampleSize = chain->totalSampleSize();
    
    double slidingWindowSize= 2* (pop->oldeffectPopSize - pop->birthRate * pop->sampleSize);
    do {
        newEffectPopulationSize= chain->proposalSlidingWindow(pop->oldeffectPopSize,  slidingWindowSize);
        
        pop->effectPopSize= newEffectPopulationSize;
        
        chain->totalEffectPopSize =  chain->totalEffectPopSize + pop->effectPopSize - pop->oldeffectPopSize;
       
        chain->updateEffectPopSizesCurrentProportionsVector();
     
        pop->growthRate =pop->delta  / pop->effectPopSize;
        pop->popSize=pop->effectPopSize * pop->birthRate;
        pop->deathRate= pop->birthRate - pop->growthRate;
        if (pop->popSize < pop->sampleSize)
        {
                populationPopSizesSet=false;
        }
    }
    while(!populationPopSizesSet);
    chain->initPopulationMigration();//after setting the timeSTD
    chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
    chain->filterSortPopulationsCoalescentEvents();
}
double NewEffectPopSizeMoveForPopulation::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions){
    
    Chain *chain=getChain();
    newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions);
    
    // double slidingWindowSize= 2* (pop->oldeffectPopSize - pop->birthRate * pop->sampleSize);
    
    double priorDensityNewEffectivePopulationSize= LogUniformDensity(pop->effectPopSize,  mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    
    double priorDensityCurrentTotalEffectivePopulationSize= LogUniformDensity(chain->oldTotalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewEffectivePopulationSize;
    
    double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensityCurrentTotalEffectivePopulationSize;
    
    double LogAcceptanceRate = std::min(sumLogNumerators - sumLogDenominators,0.0);
    
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
void NewTimeOriginOnEdgeforPopulationMove::move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    MCMCmove::move(programOptions,mcmcOptions);
}
NewTimeOriginOnTreeforPopulationMove::NewTimeOriginOnTreeforPopulationMove(Chain *chain,string nameMove,  Population *pop):MCMCmove(chain, nameMove)
{
    this->pop =pop;
}
NewTimeOriginOnEdgeforPopulationMove::NewTimeOriginOnEdgeforPopulationMove(Chain *chain,string nameMove,  Population *pop):MCMCmove(chain, nameMove)
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
         popI->oldSampleSize = popI->sampleSize;
        pop->oldFatherPop =  pop->FatherPop;
        popI->oldCoalescentEventTimes = popI->CoalescentEventTimes;
        popI->CoalescentEventTimes.clear();
        popI->oldimmigrantsPopOrderedByModelTime= popI->immigrantsPopOrderedByModelTime;
        popI->immigrantsPopOrderedByModelTime.clear();
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
        popI->sampleSize = popI->oldSampleSize;
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
    std::map<pll_rnode_t*, Population*>  newrmrcaOfPopulation;

    chain->proposedrMRCAPopulation=  chain->chooseAvailableEdgeOnRootedTreeForPopulation(pop, chain->rMRCAPopulation, programOptions.healthyTipLabel);
    
    for (unsigned int i = 0; i < chain->numClones; ++i){
        auto pop =  chain->populations[i];
        alpha[i]= pop->sampleSize;
    }
    
    int totalSampleSize=chain->initialRootedTree->tip_count-1;//not the healthytip
    std::transform(alpha, alpha + chain->numClones , alpha,[totalSampleSize](double a) {return a /totalSampleSize; } );
    chain->initProportionsVector();
    chain->generateProportionsVectorFromDirichlet(alpha);
    chain->initEffectPopulationSizesFromProportionsVector();
    chain->initTimeOriginSTD();
    chain->initPopulationMigration();//after setting the timeSTD
    chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->proposedrMRCAPopulation, programOptions.healthyTipLabel);
    chain->filterSortPopulationsCoalescentEvents();
}
double  NewTimeOriginOnTreeforPopulationMove::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    
   newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions);
    double currentSumAvailBranchLengths = chain->sumAvailableBranchLengths(chain->rMRCAPopulation);
      double newSumAvailBranchLengths = chain->sumAvailableBranchLengths(chain->proposedrMRCAPopulation);
    double numeratorQ=log(pop->oldrMRCA->length) + log(newSumAvailBranchLengths);
//    priorDensityNewTotalEffectivePopulationSize = LogUniformDensity(chain->totalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    double denominatorQ= log(pop->rMRCA->length) + log(currentSumAvailBranchLengths);
//    priorDensityCurrentTotalEffectivePopulationSize= LogUniformDensity(chain->oldTotalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +numeratorQ;
    
    double sumLogDenominators=chain->currentlogConditionalLikelihoodTree + denominatorQ;
    
    double LogAcceptanceRate = std::min(sumLogNumerators - sumLogDenominators,0.0);
    return LogAcceptanceRate;
}
void NewTimeOriginOnEdgeforPopulationMove::safeCurrentValue()
{
    Chain *chain=getChain();
    int i=0;
    Population *popI;
    //save information for the population
    
    pop->oldFatherPop =  pop->FatherPop;
    pop->oldCoalescentEventTimes = pop->CoalescentEventTimes;
    pop->oldimmigrantsPopOrderedByModelTime= pop->immigrantsPopOrderedByModelTime;
    pop->oldSampleSize = pop->sampleSize;
    
    //save information for the other populations
    for( i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        
        popI->oldFatherPop =  pop->FatherPop;
        popI->oldCoalescentEventTimes = popI->CoalescentEventTimes;
        popI->CoalescentEventTimes.clear();
        popI->oldimmigrantsPopOrderedByModelTime= popI->immigrantsPopOrderedByModelTime;
        popI->immigrantsPopOrderedByModelTime.clear();
    }
}
void NewTimeOriginOnEdgeforPopulationMove::rollbackMove()
{
    Chain *chain=getChain();
    int i=0;
    Population *popI;
    //save information for the population
    //pop->rMRCA = pop->oldrMRCA ;
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
void NewTimeOriginOnEdgeforPopulationMove::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    
    std::map<pll_rnode_t*, Population*>  newrmrcaOfPopulation;
    
    chain->chooseNewTimeofOriginOnEdge(pop);

    chain->initPopulationMigration();//after setting the timeSTD
    chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
    chain->filterSortPopulationsCoalescentEvents();
}
double  NewTimeOriginOnEdgeforPopulationMove::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    
   newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions);
    double currentSumAvailBranchLengths;// chain->sumAvailableBranchLengths(chain->rMRCAPopulation);
    double newSumAvailBranchLengths ;// chain->sumAvailableBranchLengths(chain->proposedrMRCAPopulation);
   // double numeratorQ=log(pop->oldrMRCA->length) + log(newSumAvailBranchLengths);
    double numeratorQ =0;
    //    priorDensityNewTotalEffectivePopulationSize = LogUniformDensity(chain->totalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    //double denominatorQ= log(pop->rMRCA->length) + log(currentSumAvailBranchLengths);
    double denominatorQ= 0;
    //    priorDensityCurrentTotalEffectivePopulationSize= LogUniformDensity(chain->oldTotalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +numeratorQ;
    
    double sumLogDenominators=chain->currentlogConditionalLikelihoodTree + denominatorQ;
    
   double LogAcceptanceRate = std::min(sumLogNumerators - sumLogDenominators,0.0);
    return LogAcceptanceRate;
}
