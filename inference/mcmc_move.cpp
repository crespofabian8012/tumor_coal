//
//  mcmc_move.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 28/10/19.
//

#include "mcmc_move.hpp"
#include "data_utils.hpp"
#include "random.h"
#include <algorithm>
using namespace std;
MCMCmove::MCMCmove(Chain *chain,  string nameMove)
{
    this->chain = chain;
    this->nameMove = nameMove;
    this->numberAccept=0;
    this->numberReject=0;
    this->numberAttemps=0;
    this->newLogConditionalLikelihoodTree =chain->currentlogConditionalLikelihoodTree;
    this->newLogConditionalLikelihoodSequences=chain->currentlogConditionalLikelihoodSequences;
}
double MCMCmove::numberAccepted()
{
    return this->numberAccept;
}
double MCMCmove::numberRejected()
{
    
      return this->numberReject;
    
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
        popI->oldFatherPop = popI->FatherPop;
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
    mcmcOptions.slidingWindowSizeTotalEffectPopSize= 2 * 0.01* (chain->totalEffectPopSize - (chain->totalSampleSize() / mcmcOptions.fixedLambda));
    bool allPopulationPopSizesSet=false;
   
    do {
        numberAttemps++;
        allPopulationPopSizesSet=true;
        newTotalPopulationSize= chain->proposalSlidingWindow(chain->oldTotalEffectPopSize,  mcmcOptions.slidingWindowSizeTotalEffectPopSize);
      
        chain->totalEffectPopSize= newTotalPopulationSize;
          fprintf (stderr, "\n proposed total effect pop size %d \n",chain->totalEffectPopSize );
        chain->updateEffectPopSizesCurrentProportionsVector();
        for( i = 0 ; i < chain->numClones; i++)
        {
          popI=chain->populations[i];
          popI->growthRate =popI->delta  / popI->effectPopSize;
          popI->popSize=popI->effectPopSize * popI->birthRate;
          popI->deathRate= popI->birthRate - popI->growthRate;
          popI->numCompletedCoalescences=0;
          popI->CoalescentEventTimes.clear();
          popI->immigrantsPopOrderedByModelTime.clear();
          popI->numIncomingMigrations=0;
            
            if (popI->popSize < popI->sampleSize)
            {
                allPopulationPopSizesSet=false;
                break;
            }
        }
        if (!allPopulationPopSizesSet)
             fprintf (stderr, "\n The  popSize < sampleSize\n");
        if (!allPopulationPopSizesSet)
            break;
        chain->initPopulationMigration();//after setting the timeSTD
        chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
        chain->filterSortPopulationsCoalescentEvents();
        allPopulationPopSizesSet = chain->checkMigrationsOrder();
        if (!allPopulationPopSizesSet)
            fprintf (stderr, "\n The order of migrations if not correct, attempts %d\n", numberAttemps);
        if (numberAttemps >mcmcOptions.maxNumberProposalAttempts)
            break;
    }
    while(!allPopulationPopSizesSet);
    
    if (allPopulationPopSizesSet)
       fprintf (stderr, "\n new total effect pop size %d \n",chain->totalEffectPopSize );
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
        popI->FatherPop = popI->oldFatherPop;
        popI->CoalescentEventTimes = popI->oldCoalescentEventTimes;
        popI->numCompletedCoalescences=popI->oldCoalescentEventTimes.size();
        popI->oldCoalescentEventTimes.clear();
        popI->immigrantsPopOrderedByModelTime =  popI->oldimmigrantsPopOrderedByModelTime;
        popI->numIncomingMigrations=popI->oldimmigrantsPopOrderedByModelTime.size();
        popI->oldimmigrantsPopOrderedByModelTime.clear();
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
    numberAttemps=0;
    makeProposal(programOptions, mcmcOptions);
    if (numberAttemps > mcmcOptions.maxNumberProposalAttempts)
    {
        rollbackMove();
        this->numberReject++;
        printf("\n Rejected new move %s because too many  proposals \n", nameMove.c_str());
        return;
    }
        
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
        popI->oldFatherPop = popI->FatherPop;
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
    vector<double> sampleSizesvector(chain->numClones);
    bool allPopulationPopSizesSet=false;
    //double totalSampleSize= chain->totalSampleSize();
    for( i = 0 ; i < chain->numClones; i++)
    {
        sampleSizesvector.at(i) = chain->oldproportionsVector.at(i) * chain->totalEffectPopSize;
        fprintf (stderr, "\n old proportions vector at %d: %lf \n",i,sampleSizesvector.at(i) );
    }
    do {
        allPopulationPopSizesSet=true;
          numberAttemps++;
        randomDirichletFromVector (sampleSizesvector, chain->proportionsVector);
//        randomDirichletFromGsl(chain->numClones, proportionsVectorArray, &(chain->proportionsVector[0]));
        chain->updateEffectPopSizesCurrentProportionsVector();
        for( i = 0 ; i < chain->numClones; i++)
        {
            popI=chain->populations[i];
            popI->growthRate =popI->delta  / popI->effectPopSize;
            popI->popSize=popI->effectPopSize * popI->birthRate;
            popI->deathRate= popI->birthRate - popI->growthRate;
            popI->numCompletedCoalescences=0;
            popI->CoalescentEventTimes.clear();
            popI->immigrantsPopOrderedByModelTime.clear();
            popI->numIncomingMigrations=0;
            if (popI->popSize < popI->sampleSize)
            {
                allPopulationPopSizesSet=false;
            }
         }
        if (!allPopulationPopSizesSet)
            break;
        chain->initPopulationMigration();//after setting the timeSTD
        chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
            chain->filterSortPopulationsCoalescentEvents();
        allPopulationPopSizesSet = chain->checkMigrationsOrder();
        if (numberAttemps >mcmcOptions.maxNumberProposalAttempts)
            break;
    }
    while(!allPopulationPopSizesSet);
    if (allPopulationPopSizesSet)
    {
//        for( i = 0 ; i < chain->numClones; i++)
//         fprintf (stderr, "\n new proportions vector at %d: %lf \n",i,chain->proportionsVector.at(i) );
    }
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
        popI->FatherPop = popI->oldFatherPop;
        popI->deathRate = popI->oldDeathRate;
        popI->growthRate = popI->oldGrowthRate;
        popI->CoalescentEventTimes = popI->oldCoalescentEventTimes;
           popI->numCompletedCoalescences=popI->oldCoalescentEventTimes.size();
        popI->oldCoalescentEventTimes.clear();
        popI->immigrantsPopOrderedByModelTime =  popI->oldimmigrantsPopOrderedByModelTime;
           popI->numIncomingMigrations=popI->oldimmigrantsPopOrderedByModelTime.size();
        popI->oldimmigrantsPopOrderedByModelTime.clear();
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
        popI->oldFatherPop = popI->FatherPop;
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
        popI->FatherPop = popI->oldFatherPop;
        popI->CoalescentEventTimes = popI->oldCoalescentEventTimes;
          popI->numCompletedCoalescences=popI->oldCoalescentEventTimes.size();
        popI->oldCoalescentEventTimes.clear();
        popI->immigrantsPopOrderedByModelTime =  popI->oldimmigrantsPopOrderedByModelTime;
          popI->numIncomingMigrations=popI->oldimmigrantsPopOrderedByModelTime.size();
        popI->oldimmigrantsPopOrderedByModelTime.clear();
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
        popI->oldFatherPop = popI->FatherPop;
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
        popI->FatherPop = popI->oldFatherPop;
        popI->CoalescentEventTimes = popI->oldCoalescentEventTimes;
          popI->numCompletedCoalescences=popI->oldCoalescentEventTimes.size();
        popI->oldCoalescentEventTimes.clear();
        popI->immigrantsPopOrderedByModelTime =  popI->oldimmigrantsPopOrderedByModelTime;
          popI->numIncomingMigrations=popI->oldimmigrantsPopOrderedByModelTime.size();
        popI->oldimmigrantsPopOrderedByModelTime.clear();
    }
}
void NewEffectPopSizeMoveForPopulation::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    int i;
 
    Chain *chain=getChain();

    bool populationPopSizesSet=false;
    double newEffectPopulationSize;
    //int totalSampleSize = chain->totalSampleSize();
    
    double slidingWindowSize= 2* 0.1* (pop->oldeffectPopSize - pop->birthRate * pop->sampleSize);
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
        if (!populationPopSizesSet)
            break;
        chain->initPopulationMigration();//after setting the timeSTD
        chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
        chain->filterSortPopulationsCoalescentEvents();
        populationPopSizesSet = chain->checkMigrationsOrder();
    }
    while(!populationPopSizesSet);
   
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
    pop->oldrMRCA =  pop->rMRCA ;//this only changes to population pop
    pop->oldTimeOriginInput = pop->timeOriginInput;
   
     //save information for the other populations
    for( i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        popI->oldSampleSize = popI->sampleSize;
        popI->oldeffectPopSize = popI->effectPopSize;
         popI->oldTimeOriginSTD = popI->timeOriginSTD;
        popI->oldFatherPop =  popI->FatherPop;
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
    
    pop->timeOriginInput = pop->oldTimeOriginInput;
    
    //save information for the other populations
    for( i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        popI->sampleSize = popI->oldSampleSize;
         popI->effectPopSize = popI->oldeffectPopSize;
        popI->timeOriginSTD = popI->oldTimeOriginSTD;
        popI->FatherPop = popI->oldFatherPop ; //the father population can change also for other populations
        popI->CoalescentEventTimes = popI->oldCoalescentEventTimes ;
          popI->numCompletedCoalescences=popI->oldCoalescentEventTimes.size();
        popI->oldCoalescentEventTimes.clear();
        popI->immigrantsPopOrderedByModelTime = popI->oldimmigrantsPopOrderedByModelTime;
          popI->numIncomingMigrations=popI->oldimmigrantsPopOrderedByModelTime.size();
        popI->oldimmigrantsPopOrderedByModelTime.clear();
    }
}
void NewTimeOriginOnTreeforPopulationMove::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
     Chain *chain=getChain();

    //first compute the list of  edges without events
    //then select at random an edge(proportional to the length)
    //recompute the sample sizes
    //recompute the proportions vector
    // recompute the coalescent and migration events
    bool existsZeroSampleSizePop=false;
    double alpha[chain->numClones];
    int numberPoints = chain->numClones -1;
    std::map<pll_rnode_t*, Population*>  newrmrcaOfPopulation;
      map<pll_rnode_t*, Population*>::iterator it;
    for ( it = chain->rMRCAPopulation.begin(); it != chain->rMRCAPopulation.end(); it++ )
    {
//        fprintf (stderr, "\n Before. The population %d with order %d has MRCA node %d, time origin input %lf, std %lf, and sample size %d\n", it->second->index,  it->second->order,  it->first->node_index, it->second->timeOriginInput, it->second->timeOriginSTD, it->second->sampleSize);
    }
        chain->proposedrMRCAPopulation=  chain->chooseAvailableEdgeOnRootedTreeForPopulation(pop, chain->rMRCAPopulation, programOptions.healthyTipLabel);
    
//       for ( it = chain->proposedrMRCAPopulation.begin(); it != chain->proposedrMRCAPopulation.end(); it++ )
//        {
//        fprintf (stderr, "\n After. The population %d with order %d has MRCA node %d \n", it->second->index,  it->second->order,  it->first->node_index );
//         }

    for (unsigned int i = 0; i < chain->numClones; ++i)
    {
            auto pop =  chain->populations[i];
            alpha[i]= pop->sampleSize;
           // fprintf (stderr, "\n New sample size %d for population order %d \n", pop->sampleSize, pop->order );
    }
    int totalSampleSize=chain->initialRootedTree->tip_count-1;//not the healthytip
    std::transform(alpha, alpha + chain->numClones , alpha,[totalSampleSize](double a) {return a /totalSampleSize; } );
    //chain->initProportionsVector();
    chain->copyProportionsVector(alpha);
   // chain->generateProportionsVectorFromDirichlet(alpha);
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
    //pop->oldFatherPop =  pop->FatherPop;
    //pop->oldCoalescentEventTimes = pop->CoalescentEventTimes;
    //pop->oldimmigrantsPopOrderedByModelTime= pop->immigrantsPopOrderedByModelTime;
    //pop->oldSampleSize = pop->sampleSize;
//    for(unsigned i=0; i< pop->FatherPop->immigrantsPopOrderedByModelTime.size();i++){
//
//         fprintf (stderr, "\n inmigrant %d, time %lf to father population order %d of population of order  %d \n", i , pop->FatherPop->immigrantsPopOrderedByModelTime.at(i).first, pop->FatherPop->order, pop->order);
//
//    }
    
//    for(unsigned i=0; i< pop->CoalescentEventTimes.size();i++){
//
//        fprintf (stderr, "\n coalescent %d, time %lf  of population order %d and sample size %d  \n", i , pop->CoalescentEventTimes.at(i), pop->order, pop->sampleSize);
//
//    }
   
    pop->FatherPop->oldimmigrantsPopOrderedByModelTime =  pop->FatherPop->immigrantsPopOrderedByModelTime;
   // pop->FatherPop->immigrantsPopOrderedByModelTime.clear();
    
    pop->oldTimeOriginInput = pop->timeOriginInput;
    pop->oldTimeOriginSTD = pop->timeOriginSTD;
    //save information for the other populations
//    for( i = 0 ; i < chain->numClones; i++)
//    {
//        popI=chain->populations[i];
//       // popI->oldFatherPop =  pop->FatherPop;
//        popI->oldCoalescentEventTimes = popI->CoalescentEventTimes;
//        popI->CoalescentEventTimes.clear();
//        popI->numCompletedCoalescences=0;
//        popI->oldimmigrantsPopOrderedByModelTime= popI->immigrantsPopOrderedByModelTime;
//        popI->immigrantsPopOrderedByModelTime.clear();
//        popI->numIncomingMigrations=0;
//    }
}
void NewTimeOriginOnEdgeforPopulationMove::rollbackMove()
{
    Chain *chain=getChain();
    int i=0;
    Population *popI;
    //save information for the population
    
    pop->FatherPop->immigrantsPopOrderedByModelTime = pop->FatherPop->oldimmigrantsPopOrderedByModelTime;
    
    //pop->FatherPop->oldimmigrantsPopOrderedByModelTime.clear();
    pop->timeOriginInput = pop->oldTimeOriginInput;
    pop->timeOriginSTD = pop->oldTimeOriginSTD;
    //save information for the other populations
//    for( i = 0 ; i < chain->numClones; i++)
//    {
//        popI=chain->populations[i];
//        //pop->FatherPop = pop->oldFatherPop ;
//        pop->CoalescentEventTimes = pop->oldCoalescentEventTimes ;
//        pop->oldCoalescentEventTimes.clear();
//        pop->immigrantsPopOrderedByModelTime = pop->oldimmigrantsPopOrderedByModelTime;
//        pop->oldimmigrantsPopOrderedByModelTime.clear();
//    }
}
void NewTimeOriginOnEdgeforPopulationMove::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    
    std::map<pll_rnode_t*, Population*>  newrmrcaOfPopulation;
    if (pop->timeOriginInput <=0)
        fprintf(stderr, "the time origin input is 0");
    
    chain->chooseNewTimeofOriginOnEdge(pop);
    //chain->initTimeOriginSTD();
//    chain->initPopulationMigration();//after setting the timeSTD
//
//    chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
//    chain->filterSortPopulationsCoalescentEvents();
    sort(pop->FatherPop->immigrantsPopOrderedByModelTime.begin(), pop->FatherPop->immigrantsPopOrderedByModelTime.end(), Population::comparePopulationsPairByTimeOrigin);
 
//    for (int i = 0; i < pop->FatherPop->immigrantsPopOrderedByModelTime.size(); ++i)
//        printf("\n ordered migrations: time(father pop units) : %lf, pop order: %d, time of origin input%lf \n", pop->FatherPop->immigrantsPopOrderedByModelTime[i].first,  pop->FatherPop->immigrantsPopOrderedByModelTime[i].second->order , pop->FatherPop->immigrantsPopOrderedByModelTime[i].second->timeOriginInput);
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
