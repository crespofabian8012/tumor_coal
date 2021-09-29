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
 * tumor population class
 */

#include <boost/multiprecision/cpp_dec_float.hpp>
#include "population.hpp"
#include "definitions.hpp"
#include "poset_smc.hpp"
#include "partial_tree_node.hpp"
#include "state.hpp"
#include <algorithm>
#include <vector>
#include <iterator>
#include "data_types.hpp"
#include "data_utils.hpp"
#include "mcmc_parameter.hpp"
#include "utils.hpp"
#include <gsl/gsl_randist.h>
#include <random>

extern "C"
    {
#include "random.h"
    }

#include <boost/multiprecision/cpp_dec_float.hpp>

//#include <algorithm>
//#include <vector>
//#include <iterator>

//#include "data_types.hpp"
//#include "data_utils.hpp"
//#include "mcmc_parameter.hpp"

//#include "data_utils.hpp"

//Parametrized Constructor
Population::Population(int ind, int ord, long double timeOriginInput,
                       int sampleSize,  int popSize, long double birthRate,
                       long  double deathRate, bool estimateTOR, MCMCoptions &mcmcOptions)
{
    if (ind >=0)
        index=ind;
    else
        fprintf (stderr, "\n ERROR: Population  index cannot be negative \n");
    if (ord >= 0)
        order=ord;
    else
        fprintf (stderr, "\n ERROR: Population order cannot be negative \n");
    
    //    if (popSize >= 0)
    //    {
    //        this->popSize=popSize;
    //    }
    //    else
    //        fprintf (stderr, "\n ERROR: Population size cannot be negative \n");
    //
    if (sampleSize >= 0)
    {
        this->sampleSize=sampleSize;
    }
    else
        fprintf (stderr, "\n ERROR: Population sample size cannot be negative \n");
    
    if (birthRate >0){
        this->birthRate =birthRate;
    }
    else
        fprintf (stderr, "\n ERROR: Population birth rate cannot be negative \n");
    
    if (deathRate >0){
        this->deathRate = deathRate;
    }
    else
        fprintf (stderr, "\n ERROR: Population death rate cannot be negative \n");
    
    if (birthRate > deathRate)
    {
        growthRate = birthRate - deathRate;
    }
    else
        fprintf (stderr, "\n ERROR: Population growth rate cannot be negative \n");
    
    if (timeOriginInput >= 0)
    {
        this->timeOriginInput = timeOriginInput;
    }
    else
        fprintf (stderr, "\n ERROR: Population time of origin   cannot be negative \n");
    //    if (birthRate >0)
    //        effectPopSize = popSize / birthRate;
    
    //delta= growthRate * effectPopSize;
    delta = 0.0;
    
    isAlive = YES;
    
    //    if (effectPopSize>0)
    //      timeOriginSTD = timeOriginInput /effectPopSize;
    //    else
    timeOriginSTD=0;
    
    numActiveGametes = sampleSize;
    
    numGametes = sampleSize;
    
    nodeIdAncestorMRCA = 0;
    numCompletedCoalescences = 0;
    nextAvailableIdInmigrant = 0;
    numIncomingMigrations = 0;
    numPossibleMigrations = 0;
    doEstimateTimeOrigin = estimateTOR;
    lowerBoundTimeOriginInput=0;
    
    MRCA=0;
    
    deltaT = 0.0;
    olddeltaT = 0.0;
    
    
    x= 1;
    
    oldx = 0.0;
    theta = 0.0;
    oldTheta = 0.0;
    oldOrder=0;
    oldScaledTimeOriginInput=0.0;
    olddelta=0.0;
    //oldeffectPopSize =0.0;
    oldDeathRate=0.0;
    //oldPopSize=0.0;
    oldGrowthRate =0.0;
    oldTimeOriginSTD=0.0;
    oldTimeOriginInput=0.0;
    indexFirstObservedCellName=0;
    scaledtimeOriginInput=0.0;
    currentModelTime = 0.0;
    
    long double zero =0.0;
    DeltaT = new MCMCParameterWithKernel("DeltaT_pop"+std::to_string(index), 0.0,mcmcOptions.paramMultiplierGrowthRate, 0.0);
    
    if(mcmcOptions.priorsType==0)
    {
        DeltaT->getParameter()->setPrior(&Distributions::LogUniformDensity);
        
        DeltaT->getParameter()->setSecondParLogPrior(mcmcOptions.GrowthRatefrom); //from parameter
        DeltaT->getParameter()->setThirdParLogPrior(mcmcOptions.GrowthRateto); //to parameter
    }
    else if(mcmcOptions.priorsType==2)
    {
        DeltaT->getParameter()->setPrior(&Distributions::LogPowerLawDistibutionDensity);
        DeltaT->getParameter()->setSecondParLogPrior(mcmcOptions. parameterPowerLawDistributionGrowthRate);
        DeltaT->getParameter()->setThirdParLogPrior(zero); //from parameter
    }
    
    else if(mcmcOptions.priorsType==1) {
        DeltaT->getParameter()->setPrior(&Distributions::LogExponentialDensity);
        DeltaT->getParameter()->setThirdParLogPrior(mcmcOptions.lambdaExponentialPriorGrowthRate);//lambda parameter
        DeltaT->getParameter()->setSecondParLogPrior(zero);//from paremeter
    }
    
    
    TimeOriginSTD = new MCMCParameterWithKernel("TimeOriginSTD_pop"+std::to_string(index),0.0,mcmcOptions.paramMultiplierTimeOriginOldestPop, 0.0);
    
    Theta=new MCMCParameterWithKernel("Theta_pop"+std::to_string(index),0.0,mcmcOptions.paramMultiplierMoveTheta,0.0);
    long double doubleSampleSize= (long double) sampleSize;
    
    X = new MCMCParameter<long double>("Proportion_pop"+std::to_string(index),doubleSampleSize, 0.0);
    sampleSizePar = new MCMCParameter<long double>("SampleSize_pop"+std::to_string(index),doubleSampleSize, 0.0);
}

Population::Population(int ind, int ord, long double timeOriginInput,
                       int sampleSize,int totalSampleSize,  int popSize,int totalPopSize,  long double birthRate,
                       long  double deathRate, bool estimateTOR)
{
    if (ind >=0)
        index=ind;
    else
        fprintf (stderr, "\n ERROR: Population  index cannot be negative \n");
    if (ord >= 0)
        order=ord;
    else
        fprintf (stderr, "\n ERROR: Population order cannot be negative \n");
    
    //    if (popSize >= 0)
    //    {
    //        this->popSize=popSize;
    //    }
    //    else
    //        fprintf (stderr, "\n ERROR: Population size cannot be negative \n");
    
    if (sampleSize >= 0)
    {
        this->sampleSize=sampleSize;
    }
    else
        fprintf (stderr, "\n ERROR: Population sample size cannot be negative \n");
    
    if (birthRate >0){
        this->birthRate =birthRate;
    }
    else
        fprintf (stderr, "\n ERROR: Population birth rate cannot be negative \n");
    
    if (deathRate >0){
        this->deathRate = deathRate;
    }
    else
        fprintf (stderr, "\n ERROR: Population death rate cannot be negative \n");
    
    if (birthRate > deathRate)
    {
        growthRate = birthRate - deathRate;
    }
    else
        fprintf (stderr, "\n ERROR: Population growth rate cannot be negative \n");
    
    if (timeOriginInput >= 0)
    {
        this->timeOriginInput = timeOriginInput;
    }
    else
        fprintf (stderr, "\n ERROR: Population time of origin   cannot be negative \n");
    
    long double effectPopSize=1.0;
    long double totalEffectPopSize=1.0;
    if (birthRate >0){
        effectPopSize = popSize / birthRate;
        totalEffectPopSize = totalPopSize / birthRate;
        if (effectPopSize>0 && totalEffectPopSize >0){
            
            delta= growthRate * (effectPopSize);
            deltaT = growthRate * totalEffectPopSize;
        }
        else{
            delta=growthRate;
            deltaT = growthRate;
        }
        
    }
    else
        delta =0.0;
    
    isAlive = YES;
    
    if (!estimateTOR){
        
        if ( effectPopSize>0){
            // timeOriginSTD = timeOriginInput /effectPopSize;
            timeOriginSTD = timeOriginInput ;
        }
        else
            timeOriginSTD=timeOriginInput;
        
        
    }
    
    numActiveGametes = sampleSize;
    
    numGametes = sampleSize;
    
    nodeIdAncestorMRCA = 0;
    numCompletedCoalescences = 0;
    nextAvailableIdInmigrant = 0;
    indexNextInmigrant = 0;
    numIncomingMigrations = 0;
    numPossibleMigrations = 0;
    doEstimateTimeOrigin = estimateTOR;
    lowerBoundTimeOriginInput=0;
    
    MRCA=0;
    currentModelTime = 0.0;
    deltaT = 0.0;
    olddeltaT = 0.0;
    
    if (totalSampleSize >0)
        x = (1.0 *popSize) /totalPopSize;
    else
        x = 1.0;
    
    oldx = 0.0;
    theta = 0.0;
    oldTheta = 0.0;
    oldOrder=0;
    oldScaledTimeOriginInput=0.0;
    olddelta=0.0;
    oldDeathRate=0.0;
    
    oldGrowthRate =0.0;
    oldTimeOriginSTD=0.0;
    oldTimeOriginInput=0.0;
    indexFirstObservedCellName=0;
    scaledtimeOriginInput=0.0;
    
    DeltaT = new MCMCParameterWithKernel("DeltaT",0.0,1, 0.0);
    
    TimeOriginSTD = new MCMCParameterWithKernel("TimeOriginSTD",0.0,1,0.0);
    TimeOriginInput = new MCMCParameterWithKernel("TimeOriginInput",0.0,1,0.0);
    Theta=new MCMCParameterWithKernel("Theta",0.0,1, 0.0);
    
    long double doubleSampleSize= (long double) sampleSize;
    X = new MCMCParameter<long double>("Proportion",doubleSampleSize, 0.0);
    sampleSizePar = new MCMCParameter<long double>("SampleSize",doubleSampleSize, 0.0);
}
Population::Population(int ind, int ord, int sampleSize, long double delta, long double theta,
                       ProgramOptions &programOptions)
{
    if (ind >=0)
        index=ind;
    else
        fprintf (stderr, "\n ERROR: Population  index cannot be negative \n");
    if (ord >= 0)
        order=ord;
    else
        fprintf (stderr, "\n ERROR: Population order cannot be negative \n");
    
    
    if (sampleSize >= 0)
    {
        this->sampleSize=sampleSize;
    }
    else
        fprintf (stderr, "\n ERROR: Population sample size cannot be negative \n");
    
    
    if (delta > 0)
    {
        this->delta = delta;
    }
    else
        fprintf (stderr, "\n ERROR: Population scaled growth rate cannot be negative \n");
    
    if (theta > 0)
    {
        this->theta = theta;
    }
    else
        fprintf (stderr, "\n ERROR: Population scaled growth rate cannot be negative \n");
    
    
    
    
    isAlive = YES;
    
    numActiveGametes = sampleSize;
    
    numGametes = sampleSize;
    
    nodeIdAncestorMRCA = 0;
    numCompletedCoalescences = 0;
    nextAvailableIdInmigrant = 0;
    numIncomingMigrations = 0;
    numPossibleMigrations = 0;
    indexNextInmigrant = 0;
    
    lowerBoundTimeOriginInput=0;
    
    MRCA=0;
    
    deltaT = 0.0;
    olddeltaT = 0.0;
    
    x = 1.0;
    oldx = 0.0;
    oldTheta = 0.0;
    oldOrder=0;
    oldScaledTimeOriginInput=0.0;
    olddelta=0.0;
    oldDeathRate=0.0;
    
    oldGrowthRate =0.0;
    oldTimeOriginSTD=0.0;
    oldTimeOriginInput=0.0;
    indexFirstObservedCellName=0;
    scaledtimeOriginInput=0.0;
    currentModelTime = 0.0;
    
    DeltaT = new MCMCParameterWithKernel("DeltaT",0.0,1, 0.0);
    
    TimeOriginSTD = new MCMCParameterWithKernel("TimeOriginSTD",0.0,1,0.0);
    TimeOriginInput = new MCMCParameterWithKernel("TimeOriginInput",0.0,1,0.0);
    Theta=new MCMCParameterWithKernel("Theta",0.0,1, 0.0);
    
    long double doubleSampleSize= (long double) sampleSize;
    X = new MCMCParameter<long double>("Proportion",doubleSampleSize, 0.0);
    sampleSizePar = new MCMCParameter<long double>("SampleSize",doubleSampleSize, 0.0);
}
Population::Population(const Population &original){
    
    index = original.index;
    order = original.order;
    sampleSize = original.sampleSize;
    birthRate = original.birthRate;
    deathRate = original.deathRate;
    delta = original.delta;
    deltaT = original.deltaT;
    theta = original.theta;
    x = original.x;
    growthRate = original.growthRate;
    idsGametes = original.idsGametes;
    numCompletedCoalescences = original.numCompletedCoalescences;
    numIncomingMigrations = original.numIncomingMigrations;
    numPossibleMigrations = original.numPossibleMigrations;
    indexNextInmigrant = original.indexNextInmigrant ;
    currentModelTime = original.currentModelTime;
    
    timeOriginSTD = original.timeOriginSTD;
    scaledtimeOriginInput = original.scaledtimeOriginInput;
    isAlive = original.isAlive;
    numActiveGametes = original.numActiveGametes;
    numGametes = original.numGametes;
    
    FatherPop = original.FatherPop;
    
    CoalescentEventTimes = original.CoalescentEventTimes;
    immigrantsPopOrderedByModelTime = original.immigrantsPopOrderedByModelTime;
    idsActiveGametes = original.idsActiveGametes;
    pairCurrentProposals = original.pairCurrentProposals;
    
}
long double Population::ProbabilityComeFromPopulation(Population *PopJ, std::vector<Population*> &populations, int numClones, long double K)
{
    long double ProbabilityIJ, AboveTerm, BelowTerm;
    int     l, j;
    long double h, a, b, c, d, e, t, cum;
    Population *p;
    ProbabilityIJ = 0.0;
    AboveTerm = 0.0;
    BelowTerm = 0.0;
    l = 0;
    h = 0.0;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    d = 0.0;
    e = 0.0;
    t = 0.0;
    cum = 0.0;
    // calculate h
    
    t = (timeOriginSTD ) * (x) / ( PopJ->x);
    h = exp(LogCalculateH(t, PopJ->timeOriginSTD, PopJ->delta, K));
    //AboveTerm = ( PopJ->popSize) * h;
    AboveTerm = h;
    j=0;
    for (l = order + 1; l < numClones; l++)
    {
        p = populations[l];
        t = (timeOriginSTD ) * (x) / ( p->x);
        h = exp(LogCalculateH(t, p->timeOriginSTD, p->delta, K));
        
        // cum = cum + ( ( p->popSize) * h);
        cum = cum +  h;
    }
    
    BelowTerm = cum;
    
    ProbabilityIJ = AboveTerm / BelowTerm;
    
    return ProbabilityIJ;
}

long double Population::CalculateH (long double t, long double TOrigin, long double delta, long double K)
{
    long double H, AboveTerm, BelowTerm, firstTerm, secondTerm;
    long double a, b;
    long double extraTerm;
    assert(delta>0.0);
    assert(TOrigin>t);
    
    
    
    a = 1.0 - exp(-1.0 * delta * (TOrigin - t));
    firstTerm = a * a;
    secondTerm = exp(-1.0 * delta * t);
    AboveTerm = firstTerm * secondTerm;
    b = 1.0 - exp(-1.0 * delta * TOrigin);
    BelowTerm = b * b;
    if (BelowTerm == 0.0)
        printf ("\n Error!: BelowTerm = 0.0 \n");
    if (a == 0.0)
        printf ("\n Error!: a = 0.0 \n");
    
    extraTerm = 1+ (K /delta)*b*(secondTerm -1.0) / a ;
    
    H = (AboveTerm / BelowTerm)*extraTerm;
    
    return H;
}
long double Population::LogCalculateH (long double t, long double TOrigin, long double delta, long double K)
{
    long double logH, AboveTerm, BelowTerm, firstTerm, secondTerm, thirdTerm;
    long double a, b;
    long double extraTerm;
    assert(delta> 0.0);
    assert(TOrigin>=t);
    if (t ==TOrigin)
        return DOUBLE_NEG_INF;
    
    
    a = 1.0 - exp(-1.0 * delta * (TOrigin - t));
    firstTerm = 2.0 * log(a);
    secondTerm =  - delta * t;
    thirdTerm = exp(delta * t);
    AboveTerm = firstTerm + secondTerm;
    
    
    b = 1.0 - exp(-1.0 * delta * TOrigin);
    BelowTerm = 2.0 * log(b);
    
    if (a == 0.0)
        std::cout << "\n Error!: a = 0.0 in LogCalculateH \n" << std::endl;
    
    extraTerm = log(1+ (K /delta)*b*(thirdTerm -1.0) / a ) ;
    
    if (K==0){
        
        logH = AboveTerm - BelowTerm;
    }
    else {
        logH =  AboveTerm - BelowTerm + extraTerm;
    }
    
    
    return logH;
}
long double Population::LogLambda (long double t, long double TOrigin, long double delta, long double K)
{
    long double loglambda;
    loglambda = log(2.0)-Population::LogCalculateH ( t,  TOrigin,  delta, K);
    return loglambda;
}
long double Population::FmodelTstandard (long double t, long double TOrigin, long double delta,   long double K)
{
    long double ModelTimeF, firstTerm, secondTerm, thirdTerm;
    long double a, b, c, d;
    
    if (t==0)
        return 0.0;
    
    ModelTimeF = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    thirdTerm = 0.0;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    d=0.0;
    
    ModelTimeF = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    thirdTerm = 0.0;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    
    if (K==0){
        a = exp(delta * t) - 1.0;
        b = 1.0 - exp(-1.0 * delta * TOrigin);
        c = 1.0 - exp(-1.0 * delta * (TOrigin - t));
        
        if ( c == 0.0)
            std::cout << "\n  c = 0.0, TOrigin = " << TOrigin << " t = " <<  t << std::endl;
        if ( delta == 0.0)
            std::cout << "\n delta  = 0.0 in FmodelTstandard "<< std::endl;
        
        ModelTimeF = a * b / (delta * c);
        
        //std::cout << "ModelTimeF " <<  ModelTimeF<< std::endl;
        
    }
    else{
        
        c = exp(-1.0 * delta * TOrigin);
        d = 1.0 -c;
        a = ((K / delta) * d) - c;
        b = 1.0 - (K/ delta)*d;
        
        long double numerator =(a*exp(delta*t)+b);
        
        if ( numerator == 0)
            std::cout << "\n t= b => log(0) in FmodelTstandard \n"<< std::endl;
        
        long double denominator = 1-exp(-delta*(TOrigin-t));
        
        if ( denominator == 0)
            std::cout << "\n  log(Inf) in FmodelTstandard \n" << std::endl;
        
        ModelTimeF= (2.0 / K )*log(numerator/denominator);
    }
    return ModelTimeF;
}

long double Population::GstandardTmodel (long double V, long double TOrigin, long double delta, long double K)
{
    long double StandardTimeG, firstTerm, secondTerm, thirdTerm;
    long double a, b, c, d, e, x;
    
    StandardTimeG = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    thirdTerm = 0.0;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    d = 0.0;
    e = 0.0;
    
    if (V==0)
        return 0.0;
    
    firstTerm = TOrigin;
    
    secondTerm = 1.0 / delta;
    
    if (K==0){
        
        a =  exp(-1.0 * delta * TOrigin);
        
        b = (1 - a) * (1 - a) * (1.0 / a);
        
        c = 1 - a;
        
        d = (1 - a) * (1.0 / a);
        e = V + d;
        
        thirdTerm = log(1 - b / e);
        
        thirdTerm = log(1 - ((1 - a) * (1 - a) * (1.0 / a)) / (V * delta + (1 - a) * (1.0 / a)));
        
        thirdTerm = log(1 + delta * V - a) - log(1 + (delta * V - 1) * a);
        
        StandardTimeG = secondTerm * thirdTerm;
        
        if ( (1 + delta * V - a) <= 0 ||   (1 + (delta * V - 1)*a ) <= 0 ) // do approximation if required
        {
            std::cout << "\nApplying approximation of math formula to avoid log(0)\n" << std::endl;
            StandardTimeG = 0.0;
            firstTerm = 0.0;
            secondTerm = 0.0;
            thirdTerm = 0.0;
            a = 0.0;
            b = 0.0;
            c = 0.0;
            d = 0.0;
            e = 0.0;
            a = 1 / delta;
            b = log(1 + delta * V);
            firstTerm = a * b;
            d = (V * V * delta * exp(-1.0 * delta * TOrigin)) / (1 + V);
            secondTerm =  d;
            StandardTimeG = firstTerm - secondTerm;
        }
        
    }
    else {
        x = exp(K*V /2.0);
        c = exp(-1.0 * delta * TOrigin);
        d= 1.0 -c;
        a = ((K/delta) * d ) - c;
        b = 1.0 - (K/ delta)*d;
        
        long double numerator = x -b;
        if ( numerator == 0.0)
            fprintf (stderr, "\n x= b => log(0) \n");
        
        long double denominator = (x*c) +a;
        
        if ( denominator == 0.0)
            fprintf (stderr, "\n (x*c +a) == 0 => log(Inf) \n");
        
        thirdTerm= log(numerator/denominator);
        
        StandardTimeG = secondTerm * thirdTerm;
        
    }
    return StandardTimeG;
}

void Population::InitListPossibleMigrations(int order)
{
    this->order = order;
    numPossibleMigrations = order + 1;
    numIncomingMigrations = 1; //the time of origin counts as one migration
    immigrantsPopOrderedByModelTime.push_back(std::make_pair(timeOriginSTD, this));
}
int Population::resetMigrationsList(){
    
    immigrantsPopOrderedByModelTime.clear();
    InitListPossibleMigrations(this->order);
    return 0;
}
void Population::setCurrentModelTime(long double newModelTime )
{
    assert(newModelTime>currentModelTime );
    currentModelTime= newModelTime;
}
void Population::UpdateListMigrants( int numClones, Population *PopChild, Population *PopFather  )
{
    assert(PopChild!=NULL);
    assert(PopFather!=NULL);
    if (PopChild->index == PopFather->index) {
        fprintf (stderr, "\nError. The target Population %d for  migration must be different than the Population of origin %d \n", PopFather->index, PopChild->index);
        exit (-1);
    }
    if (PopFather->order <= PopChild->order )
    {
        fprintf (stderr, "\nError. The father Population of order   %d for  migration must be older than the Population of origin order %d \n", PopFather->order, PopChild->order);
        
        fprintf (stderr, "\nError. The father Population  %d for  migration must be older than the Population of origin %d \n", PopFather->index, PopChild->index);
        fprintf (stderr, "\nError. The MRCA node for the father is %d and the MRCA for the child is %d \n", PopFather->rMRCA->node_index, PopChild->rMRCA->node_index);
        fprintf (stderr, "\nError. The previous  MRCA for the child is %d \n",
                 PopChild->oldrMRCA->node_index);
        //exit (-1);
    }
    int updatedNumIncomingMigrations = PopFather->numIncomingMigrations;
    
    long double updatedMigrationTime = (PopChild->timeOriginSTD) * (PopChild->x) / (PopFather->x);
    if(updatedNumIncomingMigrations + 1 <= PopFather->numPossibleMigrations){
        updatedNumIncomingMigrations = updatedNumIncomingMigrations + 1;
        PopFather->numIncomingMigrations = updatedNumIncomingMigrations;
        PopFather->immigrantsPopOrderedByModelTime.push_back(std::make_pair(updatedMigrationTime, PopChild));
    }
    sort(PopFather->immigrantsPopOrderedByModelTime.begin(), PopFather->immigrantsPopOrderedByModelTime.end(), comparePopulationsPairByTimeOrigin);
    //     printf("\n pop order  %d choose pop father of order %d \n", PopChild->order, PopFather->order);
    //    for (int i = 0; i < PopFather->immigrantsPopOrderedByModelTime.size(); ++i)
    //        printf("\n ordered migrations: time(father pop units) : %lf, pop order: %d, time of origin %lf \n", PopFather->immigrantsPopOrderedByModelTime[i].first,  PopFather->immigrantsPopOrderedByModelTime[i].second->order , PopFather->immigrantsPopOrderedByModelTime[i].second->timeOriginSTD);
}

bool Population::comparePopulationsPairByTimeOrigin(const std::pair<long double, Population *> s1, const std::pair<long double, Population *> s2)
{
    return (s1.first < s2.first);
}

int Population::compare (const void * a, const void * b)
{
    long double*p1 = (long double*)a;
    
    long double*p2 = (long double*)b;
    
    if (*p1 > *p2) return 1;
    
    else if (*p2 > *p1 )return -1;
    
    else return 0;
    
}
void Population::InitCoalescentEvents(int numClones)
{
    int j = 0;
    assert(sampleSize >0);
    numCompletedCoalescences =0;
    CoalescentEventTimes.clear();
    for (j = 0; j < sampleSize + numClones - 2; j++) {
        CoalescentEventTimes.push_back(0);
    }
}
void Population::InitIdsActiveGametes()
{
    int j = 0;
    assert(sampleSize >0);
    idsActiveGametes.clear();
    for (j = 0; j < sampleSize ; j++) {
        idsActiveGametes.push_back(0);
    }
}

void Population::InitIdsGametes(int numClones)
{
    int j = 0;
    
    assert(sampleSize >0);
    idsGametes.clear();
    for (j = 0; j < (2* sampleSize + numClones - 2) ; j++) {
        idsGametes.push_back(0);
    }
}
void Population::InitRTips()
{
    int j = 0;
    assert(sampleSize >0);
    rtips.clear();
    for (j = 0; j < sampleSize  ; j++) {
        rtips.push_back(NULL);
    }
}
void Population::resetGametesCounters() {
    numCompletedCoalescences = 0;
    nodeIdAncestorMRCA = 0;
    numActiveGametes=0;
    numGametes=0;
}

void Population::resetActiveGametes()
{
    resetGametesCounters();
    idsActiveGametes.clear();
    idsGametes.clear();
}
void Population::ChooseRandomIndividual(int *firstInd,   int numClones,   int *secondInd, long *seed, int choosePairIndividuals)
{
    long double random;
    int k, w;
    long double*cumPopulPart = (long double*) malloc((numActiveGametes + 1)* (long) sizeof(double));
    if (!cumPopulPart)
    {
        fprintf (stderr, "Could not allocate cumPopulPart (%lu bytes)\n", (numActiveGametes + 1) * (long) sizeof(double));
        exit (-1);
    }
    cumPopulPart[0] = 0;
    for (k = 1; k <= numActiveGametes; k++)
        cumPopulPart[k] = 0;
    for (k = 1; k <= numActiveGametes; k++)
        cumPopulPart[k] = cumPopulPart[k - 1] + 1.0 / (numActiveGametes);
    
    w = 0;
    
    random = Random::RandomUniform(seed);
    *firstInd = Population::bbinClones(random, cumPopulPart, numActiveGametes)-1;
    w = 0;
    
    if (*firstInd >= numActiveGametes || *firstInd < 0 ) /* checking */
    {
        fprintf (stderr, "\n\nERROR: firstInd out of range!\n");
        exit (-1);
    }
    if (choosePairIndividuals== YES && numActiveGametes > 1) {
        
        do//choose randomly another individual to coalesce
        {
            random = Random::RandomUniform(seed);
            *secondInd = Population::bbinClones(random, cumPopulPart, numActiveGametes)-1;
            
        } while (*firstInd == *secondInd  );
    }
    free (cumPopulPart);
    cumPopulPart=NULL;
}
void Population::ChooseRandomIndividual(int *firstInd,   int numClones,   int *secondInd, const gsl_rng *randomGenerator, int choosePairIndividuals)
{
    
    if (choosePairIndividuals== YES && numActiveGametes > 1){
        *firstInd = gsl_rng_uniform_int(randomGenerator, numActiveGametes);
        do {
            *secondInd = gsl_rng_uniform_int(randomGenerator, numActiveGametes);
        } while (*firstInd == *secondInd);
        
        
    }
    else{
        *firstInd = gsl_rng_uniform_int(randomGenerator, numActiveGametes);
        *secondInd = -1;
        
    }
}

/***************** bbinClones *****************/
/* binary search in the probabilities with clones */
int Population::bbinClones (long double dat, long double *v, int n)
{
    int init, end, middle;
    
    if (dat >= 0 && dat <= v[1])
        return (1); /* first population */
    
    init = 1;
    end = n;
    
    while (init <= end)
    {
        middle = (init + end) / 2;
        
        if (dat > v[middle - 1] && dat <= v[middle])
            return (middle);
        else if (dat > v[middle])
            init = middle + 1;
        else
            end = middle - 1;
    }
    fprintf (stderr, "\n Warning in bbinClones function");
    exit (-1);
    return -1;
}
long double Population::LogDensityTime2(long double u){
    boost::multiprecision::cpp_dec_float_50  term1 = exp(-1.0*delta*u);
    double doubleTerm= term1.convert_to<double>();
    boost::multiprecision::cpp_dec_float_50 term2 = delta * term1;
    doubleTerm= term2.convert_to<double>();
    boost::multiprecision::cpp_dec_float_50 term3 = 1.0-term1;
    boost::multiprecision::cpp_dec_float_50 result = log(delta * term2 /(term3 * term3));
    doubleTerm= result.convert_to<double>();
    result = result - term2/term3;
    doubleTerm= result.convert_to<double>();
    
    return (long double)result;
}
long double Population::LogDensityTime(long double u){
    long double term1 = exp(-1.0*delta*u);
    long double term2 = delta * term1;
    long double term3 = 1.0-term1;
    long double result = 2*log(delta) -1.0*delta*u -2*log(term3);
    result = result - term2/term3;
    if (term3 ==0 || term2 == 0)
        fprintf (stderr, "\n LogDensityTime -inf, result %Lf  \n ", result);
    return result;
}
long double Population::DensityTimeSTD(long double u, long double from){
    long double  result=1.0;
    long double term1 = exp(-1.0*delta*u);
    long double term2 = delta * term1;
    long double term3 = 1.0-term1;
    if (u >=from){
        
        result = delta * term2 /(term3 * term3);
        result = result * exp( -1 * term2/term3);
    }
    else {
        result=0;
    }
    return result;
}
long double Population::DensityTimeSTD(long double u, long double deltaPar, long double from){
    long double  result=1.0;
    long double term1 = exp(-1.0*deltaPar*u);
    long double term2 = deltaPar * term1;
    long double term3 = 1.0-term1;
    if (u >=from){
        
        result = deltaPar * term2 /(term3 * term3);
        result = result * exp( -1 * term2/term3);
    }
    else {
        result=0;
    }
    return result;
}
long double Population::LogDensityTimeSTDFrom(long double u, long double from){
    long double  densityFrom, result;
    densityFrom = DensityTimeSTD(u, from );
    if (densityFrom ==0)
        result=log(0);
    else
        result=log(densityFrom);
    return result;
}
long double Population::LogProbNoCoalescentEventBetweenTimes(long double from, long double to, int numberActiveInd,  long double K)
{
    int j=numberActiveInd;
    long double result=0.0;
    if (j==0 || j==1)
        return 0;
    result= log(j *(j-1) /2.0) -1.0 * j* (j-1)*(Population::FmodelTstandard(to,timeOriginSTD, delta,   K)-Population::FmodelTstandard(from, timeOriginSTD, delta, K))/ 2.0;
    return result;
}
long double Population::LogProbNoCoalescentEventBetweenTimes(long double from, long double to, int numberActiveInd, long double TOrigin, long double deltaPar, long double K)
{
    int j=numberActiveInd;
    long double result=0.0;
    if (j==0 || j==1)
        return 0;
    result= log(j *(j-1) /2.0) -1.0 * j* (j-1)*(Population::FmodelTstandard(to,TOrigin, deltaPar, K)-Population::FmodelTstandard(from, TOrigin, deltaPar, K))/ 2.0;
    return result;
}
void Population::filterAndSortCoalescentEvents(){
    
    auto it = remove_if(CoalescentEventTimes.begin(), CoalescentEventTimes.end(), Population::isNotPositive);
    
    if (it != CoalescentEventTimes.end())
        CoalescentEventTimes.erase(it);
    
    sort(CoalescentEventTimes.begin(), CoalescentEventTimes.end());
}
bool Population::isNotPositive(long double d) //test condition for remove_if algo.
{
    return d <= 0;
}
void Population::multiplyCoalescentsEventByFactor(long double factor)
{
    int j = 0;
    for (j = 0; j < CoalescentEventTimes.size(); j++) {
        CoalescentEventTimes.at(j) = CoalescentEventTimes.at(j) * factor;
    }
}
void Population::multiplyMigrationsTimesByFactor(long double factor)
{
    int j = 0;
    for (j = 0; j < immigrantsPopOrderedByModelTime.size(); j++) {
        immigrantsPopOrderedByModelTime.at(j).first = immigrantsPopOrderedByModelTime.at(j).first *factor;
    }
}
void Population::setLowerBoundTimeOriginInput(long double from)
{
    lowerBoundTimeOriginInput = from;
    
}
void Population::restoreOldCoalescentTimes(){
    
    CoalescentEventTimes = oldCoalescentEventTimes ;
    numCompletedCoalescences= oldCoalescentEventTimes.size();
}
void Population::restoreOldImmigrationTimes(){
    
    immigrantsPopOrderedByModelTime = oldimmigrantsPopOrderedByModelTime;
    numIncomingMigrations= oldimmigrantsPopOrderedByModelTime.size();
}

void Population::savePosteriorValues(){
    
    
    posteriorDeltaT.push_back(deltaT);
    DeltaT->saveCurrentValue();
    DeltaT->saveCurrentValueToList();
    
    posteriorTimeOriginSTD.push_back(timeOriginSTD);
    TimeOriginSTD->saveCurrentValue();
    TimeOriginSTD->saveCurrentValueToList();
    std::vector<long double> vec = TimeOriginSTD->getParameter()->getChainValues();
    assert(std::equal(posteriorTimeOriginSTD.begin(), posteriorTimeOriginSTD.end(), vec.begin()));
    
    Theta->saveCurrentValue();
    Theta->saveCurrentValueToList();
    
    posteriorProportion.push_back(x);
    X->saveCurrentValue();
    
}
std::pair<Pair, std::pair<double, double>> Population::initPairProposalsNextEventTime(gsl_rng *rngGsl, double K){
    
    assert(numActiveGametes > 1);

    std::vector<Pair > pairs= PosetSMC::getCombinationsRandomOrder(numActiveGametes);
    
    auto rd = std::random_device {};
    auto rng = std::default_random_engine { rd() };
 
    std::shuffle(std::begin(pairs), std::end(pairs), rng);
    
    double minModeltime = DOUBLE_INF;
    double minKingmanTime = DOUBLE_INF;
    double proposalKingmanTime, proposalModelTime;
    Pair currPair;
    std::pair<Pair, std::pair<double, double>> currEntry;
    std::pair<Pair, std::pair<double, double>> winnerModelTime;
    std::pair<Pair, std::pair<double, double>> winnerKingmanTime;
    for(std::vector<Pair>::iterator it = pairs.begin(), end = pairs.end(); it != end; ++it)
    {
        currPair = *it;
        proposalKingmanTime =  Random::RandomExponential(1.0,0,  true, rngGsl,NULL );
        proposalModelTime = Population::GstandardTmodel(proposalKingmanTime, timeOriginSTD, delta, K);
        currEntry = std::make_pair(currPair,std::make_pair(0.0, proposalModelTime));
        if (proposalModelTime<minModeltime){
            
            winnerModelTime = currEntry;
            minModeltime =proposalModelTime;
        }
        if (proposalKingmanTime<minKingmanTime){
                   
                   winnerKingmanTime = currEntry;
            minKingmanTime = proposalKingmanTime;
               }
        
        pairCurrentProposals.insert(currEntry);
    }
    
   // std::cout <<" pair min Kingman time " << winnerKingmanTime.first.first << " , "<< winnerKingmanTime.first.second  << std::endl;
   // std::cout <<" pair min Model time " << winnerModelTime.first.first << " , "<< winnerModelTime.first.second << std::endl;
   // std::cout <<" min Model time " << winnerModelTime.second.second << std::endl;
    return winnerModelTime;

}
std::pair<Pair, std::pair<double, double>> Population::updatePairCurrentProposalsMap(int posNewNodeIdsGametes, double modelTimeNewNode,double K, double &logWeightDiff, std::pair<Pair, std::pair<double, double>> copyPairMinModelTime, State *currState,  gsl_rng *rngGsl, bool normalizeClv){
    
    double waitingTimeKingman, proposalKingmanTime, proposalModelTime;
    double minModeltime = DOUBLE_INF;
    double currKingmanCoalTime = Population::FmodelTstandard(modelTimeNewNode, timeOriginSTD, delta, K);
    std::pair<Pair, std::pair<double, double>>  pairMinModelTime;
    int idxFirst, idxSecond, idxFirstRoot, idxSecondRoot;
    double logLikNewHeight = 0.0;
    
    assert(pairCurrentProposals.size()== (numActiveGametes)*(numActiveGametes-1)/ 2.0 -1);
    PairDoubleMap::iterator it = pairCurrentProposals.begin();
    double lik_factor_newNode;
    while ( it != pairCurrentProposals.end()) {
        
        auto currPair = it->first;
        //assert(it.second.second >modelTimeNewNode);
        
        if (Utils::pairsIntersected(currPair, copyPairMinModelTime.first)){
            
            double dropoutNodeProposalModelTime =   it->second.second;
            double dropoutNodeModelTimeEntry =   it->second.first;
            //numerator
            logWeightDiff +=   -log(theta)+log(numActiveGametes*(numActiveGametes-1)/ 2.0) +
                                Population::LogLambda(dropoutNodeProposalModelTime, timeOriginSTD, delta,K)-
                                 (numActiveGametes*(numActiveGametes-1)/ 2.0)*Population::FmodelTstandard(dropoutNodeProposalModelTime, timeOriginSTD, delta,K) ;
        
            //denominator
             logWeightDiff -= -Population::FmodelTstandard(dropoutNodeProposalModelTime, timeOriginSTD, delta,K)+ Population::FmodelTstandard(dropoutNodeModelTimeEntry, timeOriginSTD, delta,K);
            
            idxFirst = currPair.first;
            idxSecond = currPair.second;
            
            idxFirstRoot = idsActiveGametes[idxFirst];
            idxSecondRoot = idsActiveGametes[idxSecond];
            
            assert(idxFirstRoot != idxSecondRoot);
            
            logLikNewHeight = 0.0;
            //std::unique_ptr<PartialTreeNode>
             auto newNode =  currState->uniquePtrNewNode( idxFirstRoot,  idxSecondRoot, index, dropoutNodeProposalModelTime, logLikNewHeight, normalizeClv );
           
            lik_factor_newNode = newNode->likelihood_factor();
            
            if (lik_factor_newNode > -1e10){
              //logWeightDiff -= lik_factor_newNode;
            }
            else
                std::cout << "log lik factor -inf for pair "<< currPair.first << " , " << currPair.second << std::endl;
            //only erase the pairs related with the last position of idsActiveGametes(the other ones will be updated)
            if (currPair.first == copyPairMinModelTime.first.second  || currPair.second ==copyPairMinModelTime.first.second)
                 it = pairCurrentProposals.erase(it);
            else{
                waitingTimeKingman =  Random::RandomExponential(1.0,0,  true, rngGsl,NULL );
                      
                proposalKingmanTime = waitingTimeKingman + currKingmanCoalTime;
                
                proposalModelTime = Population::GstandardTmodel(proposalKingmanTime, timeOriginSTD, delta, K);
                it->second.first = modelTimeNewNode;//creation time
                it->second.second  = proposalModelTime;//new proposal time
                if (proposalModelTime<minModeltime){
                                    
                                    pairMinModelTime = *it;
                                    minModeltime = proposalModelTime;
                                }
                it++;
            }
        }
        else{
            //find the minimum proposal model time in the pairs that didnt dropout
            if (it->second.second < minModeltime){
                
                minModeltime = it->second.second;
                pairMinModelTime = *it;
            }
             it++;
        }
    }

     assert(pairCurrentProposals.size()== (numActiveGametes-1)*(numActiveGametes-2)/ 2.0);
     return pairMinModelTime;
}
void Population::resetPosteriorValues(){
    
    posteriorDeltaT.clear();
    DeltaT->resetParameterValues();
    
    posteriorTimeOriginSTD.clear();
    TimeOriginSTD->resetParameterValues();
    
    Theta->resetParameterValues();
    
    posteriorProportion.clear();
    X->resetValues();
    
}
void Population::setTimeOriginInputTree(long double timeOriginInputTree){
    
    assert(timeOriginInputTree>0);
    timeOriginInput = timeOriginInputTree;
    if(TimeOriginInput)
        TimeOriginInput->setParameterValue(timeOriginInputTree);
    scaledtimeOriginInput = timeOriginInput / theta;
    
    timeOriginSTD = scaledtimeOriginInput / x;
    if(TimeOriginSTD)
        TimeOriginSTD->setParameterValue(timeOriginSTD);
    
}
void Population::setTimeOriginSTD(long double parTimeOriginSTD){
    
    assert(parTimeOriginSTD>0);
    timeOriginSTD = parTimeOriginSTD;
    if (TimeOriginSTD)
        TimeOriginSTD->setParameterValue(parTimeOriginSTD);
    
    scaledtimeOriginInput = timeOriginSTD * x;
    timeOriginInput = scaledtimeOriginInput * theta;
    if (TimeOriginInput)
        TimeOriginInput->setParameterValue(timeOriginInput);
}
void Population::setScaledTimeOriginInputTree(long double parScaledTimeOriginInputTree){
    assert(parScaledTimeOriginInputTree>0);
    scaledtimeOriginInput = parScaledTimeOriginInputTree;
}
void Population::setTheta(long double parTheta){
    assert(parTheta>=0);
    theta = parTheta;
    if (Theta)
        Theta->setParameterValue(parTheta);
}
void Population::setDelta(long double parDelta){
    assert(parDelta>=0);
    delta = parDelta;
}
void Population::setProportion(long double parX){
    x=parX;
    X->setValue(parX);
}
void  Population::setPopulationToriginConditionalDelta( const gsl_rng *rngGsl ){
    
    assert(delta >0);
    timeOriginSTD= Random::RandomDensityModelTimeOrigin(delta, true, 0.0, rngGsl, NULL);
    //timeOriginSTD = Random::RandomDensityModelTimeOriginLambda(birthRate,true,delta, sampleSize*10,  rngGsl, NULL);
    timeOriginInput = timeOriginSTD *x;
    
    //std::cout << "Population with sample size "<< sampleSize << " torigin std: "<<timeOriginSTD<< " and scaled physical torigin: "  << timeOriginInput << " ,  order "<< order << " and delta "<< delta << std::endl;
    
}
long double Population::proposeTimeNextCoalEvent(gsl_rng* rngGsl, int numActiveLineages,  double K){
    
    long double  RateCA = (long double)  numActiveLineages * ((long double) numActiveLineages - 1) / 2.0;
    long double  ThisTimeCA_W = Random::RandomExponentialStartingFrom (RateCA,0,  true, rngGsl,NULL ) ;
    
    return ThisTimeCA_W;
    
}
long double Population::proposeWaitingKingmanTimeNextCoalEvent(gsl_rng* rngGsl   ){
    assert(numActiveGametes>1);
    
    long double  RateCA = (long double)  numActiveGametes * ((long double) numActiveGametes - 1) / 2.0;
    long double  ThisTimeCA_W = Random::RandomExponential(RateCA,0,  true, rngGsl,NULL ) ;
    
    return ThisTimeCA_W;
    
}
long double Population::proposeWaitingKingmanTimeNextCoalEvent(gsl_rng* rngGsl, size_t numActiveGametesPar   ){
    assert(numActiveGametesPar>1);
    
    long double  RateCA = (long double)  numActiveGametesPar * ((long double) numActiveGametesPar - 1) / 2.0;
    long double  ThisTimeCA_W = Random::RandomExponential(RateCA,0,  true, rngGsl,NULL ) ;
    
    return ThisTimeCA_W;
    
}
long double Population::logConditionalLikelihoodNextCoalescentTime(long double timeNextEvent,   double K){
    
    long double result=0.0;
    long double temp, termOnlyAfterFirstCoalEvent;
    if (numActiveGametes > 1)
    {
        temp=log(numActiveGametes * (numActiveGametes-1.0)/2.0);
        result= result + temp;
        temp =  Population::LogLambda(timeNextEvent,timeOriginSTD, delta, K);
        result= result + temp;
        //        termOnlyAfterFirstCoalEvent =(currentCoalescentEventInThisEpoch == 0)?Population::FmodelTstandard(currentTime, timeOriginSTD, delta, K):Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEventInThisEpoch-1], timeOriginSTD, delta, K);//if no coalescent
        termOnlyAfterFirstCoalEvent =(currentModelTime>0)? Population::FmodelTstandard(currentModelTime, timeOriginSTD, delta, K):0.0;//if no coalescent
        temp =  (numActiveGametes / 2.0)* (numActiveGametes - 1.0)*(Population::FmodelTstandard(timeNextEvent, timeOriginSTD, delta, K)-termOnlyAfterFirstCoalEvent);
        result= result - temp;
        
    }
    return result;
}
long double Population::logConditionalDensityWaitingTimeScaledByTheta(long double waitingTimeScaledbyTheta,  double K){
    
    long double result=0.0;
    long double temp;
    
    assert(theta >0);
    temp = -log(theta);
    if (numActiveGametes > 1)
    {
        temp=log(numActiveGametes * (numActiveGametes-1.0)/2.0);
        result= result + temp;
        temp =  Population::LogLambda(currentModelTime +waitingTimeScaledbyTheta/theta,timeOriginSTD, delta, K);
        result= result + temp;
        //        termOnlyAfterFirstCoalEvent =(currentCoalescentEventInThisEpoch == 0)?Population::FmodelTstandard(currentTime, timeOriginSTD, delta, K):Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEventInThisEpoch-1], timeOriginSTD, delta, K);//if no coalescent
        temp =  ((numActiveGametes )* (numActiveGametes - 1.0)/2.0)*(Population::FmodelTstandard(currentModelTime +waitingTimeScaledbyTheta/theta, timeOriginSTD, delta, K)-Population::FmodelTstandard(currentModelTime, timeOriginSTD, delta, K));
        result= result - temp;
        
    }
    return result;
}
long double Population::logConditionalDensityWaitingTimeScaledByTheta(long double waitingTimeScaledbyTheta, size_t numActiveGametesP, double currentModelTimeP,  double K){
    
    long double result=0.0;
    long double temp;
    
    assert(theta >0);
    temp = -log(theta);
    if (numActiveGametesP > 1)
    {
        temp=log(numActiveGametesP * (numActiveGametesP-1.0)/2.0);
        result= result + temp;
        temp =  Population::LogLambda(currentModelTimeP +waitingTimeScaledbyTheta/theta,timeOriginSTD, delta, K);
        result= result + temp;
        //        termOnlyAfterFirstCoalEvent =(currentCoalescentEventInThisEpoch == 0)?Population::FmodelTstandard(currentTime, timeOriginSTD, delta, K):Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEventInThisEpoch-1], timeOriginSTD, delta, K);//if no coalescent
        temp =  ((numActiveGametesP )* (numActiveGametesP - 1.0)/2.0)*(Population::FmodelTstandard(currentModelTimeP +waitingTimeScaledbyTheta/theta, timeOriginSTD, delta, K)-Population::FmodelTstandard(currentModelTimeP, timeOriginSTD, delta, K));
        result= result - temp;
        
    }
    return result;
}
void Population::sampleEventTimesScaledByProportion(gsl_rng *random, double K){
    
    int numMigrations = numIncomingMigrations;
    int indexNextMigration = 0;
    double currentTime = 0.;
    double timeNextMigration = 0.;
    double currentTimeKingman = 0., waitingTimeKingman = 0., timeNextEvent = 0.;
    Population *incomingPop;
    double currentModelTime = 0.;
    numCompletedCoalescences =0;
    while (indexNextMigration < numMigrations) {
        timeNextMigration = immigrantsPopOrderedByModelTime[indexNextMigration].first;
        
        if (numActiveGametes >= 2) {
            
            currentTimeKingman = Population::FmodelTstandard (currentTime ,timeOriginSTD, delta,   K);//this is current time in Kingman coalescent
            currentModelTime = Population::GstandardTmodel(currentTimeKingman, timeOriginSTD, delta,  K);//this is current time in Kingman coal time
            
            assert(abs(currentTime-currentModelTime)<0.0000001);
            waitingTimeKingman=proposeWaitingKingmanTimeNextCoalEvent(random) ;
            currentTimeKingman = currentTimeKingman + waitingTimeKingman;
            //now from Kingman coal to  current population model time
            timeNextEvent =   Population::GstandardTmodel(currentTimeKingman, timeOriginSTD, delta, K);
            
        }
        else
        {
            timeNextEvent = timeNextMigration + 1.0; // it reached a "provisional" MRCA
        }
        if ( timeNextEvent < timeNextMigration)
        {
            //coalescent event
            currentTime = timeNextEvent; // update current time in model time
            assert(currentTime>0);
            assert(x>0);
            CoalescentEventTimes[numCompletedCoalescences]=  currentTime *x;
            numCompletedCoalescences= numCompletedCoalescences+1;
            numActiveGametes = numActiveGametes - 1;
            
        }
        else
        {
            if (indexNextMigration < numMigrations - 1) //indexNextMigration corresponds to one of the true migrations
            {
                currentTime = timeNextMigration;
                incomingPop = immigrantsPopOrderedByModelTime[indexNextMigration].second;
                indexNextMigration = indexNextMigration + 1;
                
                incomingPop->numActiveGametes = incomingPop->numActiveGametes - 1; /* now the other clone has 1 less node */
                
                numActiveGametes = numActiveGametes + 1; /* now this clone has 1 more node */
                
            }
            else {
                //origin reached
                currentTime = timeNextMigration;
                indexNextMigration = indexNextMigration + 1;
                
            }
        }
        
    }
}
double  Population::nextCoalEventTime(int idxNextCoal, int indexNextMigration,  double currentTime, bool& isThereInmigration, Population* inmigrantPop ){
    
    assert(idxNextCoal< CoalescentEventTimes.size());
    assert(indexNextMigration< immigrantsPopOrderedByModelTime.size());
    
    double nextCoalTime = CoalescentEventTimes[idxNextCoal];
    double nextInmigrantTime = immigrantsPopOrderedByModelTime[indexNextMigration].first;
    
    if (nextCoalTime < nextInmigrantTime){
        isThereInmigration = false;
        inmigrantPop = nullptr;
    }
    else{
        isThereInmigration = true;
        inmigrantPop = immigrantsPopOrderedByModelTime[indexNextMigration].second;
    }
    return nextCoalTime;
}

void Population::printCoalEventTimes(std::ostream &stream){
    
    for (std::vector<long double>::const_iterator i = CoalescentEventTimes.begin(); i != CoalescentEventTimes.end(); ++i){
        stream << *i << ' ';
    }
    stream <<  "\n";
}
void Population::printInmigrationEventTimes(std::ostream &stream){
    std::pair<long double, Population *> temp;
    for (std::vector<std::pair<long double, Population *>>::const_iterator i = immigrantsPopOrderedByModelTime.begin(); i != immigrantsPopOrderedByModelTime.end(); ++i){
        temp = *i;
        stream << temp.first << ' ';
    }
    stream <<  "\n";
}
// Class PopulationSet
PopulationSet::PopulationSet(int numClones){
    
    if (numClones >=0)
        this->numClones=numClones;
    else
        fprintf (stderr, "\n ERROR: The number of clones cannot be negative \n");
    
    initPopulation();
    initProportionsVector();
    
    setPopulationsBirthRate( 1);
    
}
std::vector<Population *>& PopulationSet::getPopulations(){
    
    return populations;
    
}
void PopulationSet::initPopulation()
{
    
    for (unsigned z = 0; z <= (numClones - 1); z++)
    {
        int ind = z;
        int ord = 0;
        long double  timeOriginInput = 0.0;
        int sampleSize = 0;
        int totalSampleSize =0;
        int totalPopSize =0;
        int popSize =0;
        long double  birthRate = 1.0;
        long double  deathRate = 0.99;
        
        
        auto pop = new Population(ind, ord, timeOriginInput, sampleSize,totalSampleSize, popSize,totalPopSize, birthRate, deathRate, NO);
        populations.push_back(pop);
    }
}
std::vector<long double> PopulationSet::samplePopulationGrowthRateFromPriors(MCMCoptions &mcmcOptions, const gsl_rng * randomGenerator )
{
    int i;
    Population *popI;
    long double  randomGrowthRate;
    std::vector<long double> result(numClones);
    for( i = 0 ; i < numClones; i++)
    {
        popI=populations[i];
        if (mcmcOptions.priorsType ==0)
        {
            randomGrowthRate = Random::RandomLogUniform(mcmcOptions.GrowthRatefrom, mcmcOptions.GrowthRateto, mcmcOptions.useGSLRandomGenerator, randomGenerator, NULL);
        }
        if (mcmcOptions.priorsType ==2)
        {
            randomGrowthRate =Random::RandomPowerLawDistribution(mcmcOptions.parameterPowerLawDistributionGrowthRate, 0, mcmcOptions.useGSLRandomGenerator, randomGenerator, NULL);
            
        }
        if (mcmcOptions.priorsType ==1)
        {
            if (mcmcOptions.fixedValuesForSimulation)
                randomGrowthRate= 1.0 / mcmcOptions.lambdaExponentialGrowthRateSimulation;
            else
                randomGrowthRate =Random::RandomExponentialStartingFrom(mcmcOptions.lambdaExponentialGrowthRateSimulation,0, mcmcOptions.useGSLRandomGenerator, randomGenerator, NULL);
        }
        // randomGrowthRate=100;
        popI->deltaT =  randomGrowthRate;
        popI->DeltaT->setParameterValue(randomGrowthRate);
        result.at(i)= randomGrowthRate;
        
        std::cout << "\n True scaled growth rate: " << randomGrowthRate << std::endl;
        
        if (mcmcOptions.verbose>=2)
            std::cout << "\n initial  deltaT"<< popI->deltaT << std::endl;
    }
    return(result);
}
void  PopulationSet::setPopulationsToriginConditionalDelta( const gsl_rng *rngGsl )
{
    int i;
    Population *popI;
    for( i = 0 ; i < numClones; i++)
    {
        popI=populations[ i];
        if (popI != NULL)
            popI->setPopulationToriginConditionalDelta(rngGsl);
    }
}
void  PopulationSet::setPopulationsBirthRate( long double  lambda )
{
    int i;
    Population *popI;
    for( i = 0 ; i < numClones; i++)
    {
        popI=populations[ i];
        if (popI != NULL)
            popI->birthRate = lambda;
    }
}
void PopulationSet::initPopulationSampleSizes(std::vector<int> &sampleSizes)
{
    for (size_t i = 0; i <numClones; ++i)
    {
        auto popI =  populations[i];
        popI->sampleSize = sampleSizes.at(i);
    }
}
Population * PopulationSet::getPopulationbyIndex(int indexPopulation)
{
    Population *pop;
    for (size_t i = 0; i < numClones; ++i){
        pop = populations[i];
        if (pop->index ==indexPopulation)
            return pop;
    }
    return NULL;
}
void PopulationSet::initPopulationGametes()
{
    for (size_t i = 0; i < numClones; ++i)
    {
        auto popI =  populations[i];
        
        popI->InitIdsGametes(numClones);
        popI->InitIdsActiveGametes();
        
    }
}
void PopulationSet::initPopulationRTips()
{
    for (size_t i = 0; i < numClones; ++i)
    {
        auto popI =  populations[i];
        
        popI->InitRTips();
        
    }
}
void PopulationSet::resetNumActiveGametesCounter(){
    
    for (size_t j = 0; j <numClones; j++)
    {
        auto pop =  populations[j];
        pop->numActiveGametes = pop->idsActiveGametes.size();
        pop->numGametes = pop->idsGametes.size();
    }
    
}
void PopulationSet::initProportionsVectorFromSampleSizes(std::vector<int> &sampleSizes )
{
    int sum =0;
    for (size_t i = 0; i < numClones; i++)
    {
        sum+=sampleSizes.at(i);
    }
    for (size_t i = 0; i < numClones; i++)
    {
        proportionsVector.at(i) =(long double)sampleSizes.at(i) / sum ;
    }
}
void PopulationSet::initPopulationsThetaDelta(long double theta)
{
    for (size_t i = 0; i < numClones; ++i)
    {
        auto popI =  populations[i];
        popI->x = proportionsVector[i];
        if(popI->X)
            popI->X->setValue(proportionsVector[i]);
        popI->theta = popI->x * theta;
        if (popI->Theta)
            popI->Theta->setParameterValue(popI->theta);
        
        popI->delta = popI->deltaT * popI->x;
    }
}
// this method assumes that populations are sorted in ascending order  by time of origin
double PopulationSet::proposeNextCoalEventTime( gsl_rng *random, int& idxReceiverPop, int& idxIncomingPop, int& idxIncomingNode, double &logLik, double K){
    Population  *incomingPop;
    double currentTimeKingman = 0.;
    double timeNextCoalEvent = 0.;
    double waitingTimeKingman = 0.;
    double timeNextMigration = 0.;
    bool isCoalescentEvent;
    int numberActiveGametesCurrentPop = 0;
    double modelTimeCurrentPop = 0.0;
    double nextEventTimeScaledByProportion = 0.0;
    logLik = 0.0;
    idxIncomingNode = -1;
    
    Population*  currentPop = getCurrentPopulation();
    int currentIdxNextInmigrant = currentPop->indexNextInmigrant;
    numberActiveGametesCurrentPop = currentPop->numActiveGametes;
    modelTimeCurrentPop = currentPop->getCurrentModelTime() ;
    bool nextCoalEventFound= false;
    while(!nextCoalEventFound){
        
        
        timeNextMigration = currentPop->immigrantsPopOrderedByModelTime[currentIdxNextInmigrant].first;
        
        if (numberActiveGametesCurrentPop >=2){
            
            currentTimeKingman = Population::FmodelTstandard (modelTimeCurrentPop , currentPop->timeOriginSTD, currentPop->delta,   K);//this is current time in Kingman coal time
            
            waitingTimeKingman= currentPop->proposeWaitingKingmanTimeNextCoalEvent(random,numberActiveGametesCurrentPop );
            
            currentTimeKingman = currentTimeKingman+waitingTimeKingman;
            
            timeNextCoalEvent =   Population::GstandardTmodel(currentTimeKingman, currentPop->timeOriginSTD, currentPop->delta, K);;
            isCoalescentEvent=true;
            
        }
        else{
            timeNextCoalEvent = timeNextMigration + 1.0;
        }
        if ( timeNextCoalEvent < timeNextMigration)
        {
            
            
            idxReceiverPop = currentPop->index;
            idxIncomingPop = currentPop->index;
            
            double waitingTimeScaledByTheta= (timeNextCoalEvent- modelTimeCurrentPop) *currentPop->x*currentPop->theta;
            
            logLik += currentPop->logConditionalDensityWaitingTimeScaledByTheta(waitingTimeScaledByTheta,  numberActiveGametesCurrentPop, modelTimeCurrentPop,   K);
            //            pop->setCurrentModelTime( timeNextCoalEvent);
            //            pop->CoalescentEventTimes[pop->numCompletedCoalescences]=  pop->getCurrentModelTime();
            //            pop->numCompletedCoalescences= pop->numCompletedCoalescences+1;
            nextCoalEventFound = true;
            nextEventTimeScaledByProportion += timeNextCoalEvent *currentPop->x;
            return(nextEventTimeScaledByProportion);
            
            
        }
        else{
            // logLik += probab no coalescent before the migration
            
            if (currentIdxNextInmigrant < currentPop->numIncomingMigrations - 1) //indexNextMigration corresponds to one of the true migrations
            {
                
                logLik += currentPop->LogProbNoCoalescentEventBetweenTimes(modelTimeCurrentPop, timeNextMigration, numberActiveGametesCurrentPop, K );
                // pop->setCurrentModelTime( timeNextMigration);
                
                modelTimeCurrentPop = timeNextMigration;
                
                incomingPop = currentPop->immigrantsPopOrderedByModelTime[currentPop->indexNextInmigrant].second;
                //  pop->indexNextInmigrant = pop->indexNextInmigrant + 1;
                
                currentIdxNextInmigrant = currentIdxNextInmigrant +1;
                assert(incomingPop->numActiveGametes<=1);
                idxIncomingNode = incomingPop->idsActiveGametes[incomingPop->numActiveGametes];
                //incomingPop->numActiveGametes = incomingPop->numActiveGametes - 1; /* now the other clone has 1 less node */
                
                idxIncomingPop = incomingPop->index;
                idxReceiverPop = currentPop->index;
                
                nextEventTimeScaledByProportion += timeNextMigration *currentPop->x;
                //  pop->idsActiveGametes[pop->numActiveGametes] = idxIncomingNode;
                //   pop->numActiveGametes = pop->numActiveGametes + 1; /* now this clone has 1 more node */
                numberActiveGametesCurrentPop = numberActiveGametesCurrentPop+1;
                /* now this clone has 1 more node */
            }
            else {
                //origin reached
                modelTimeCurrentPop = timeNextMigration;
                //pop->setCurrentModelTime( timeNextMigration);
                //  pop->indexNextInmigrant = pop->indexNextInmigrant + 1;
                currentIdxNextInmigrant = currentIdxNextInmigrant +1;
                //                if (i< numClones-1)   //if it is not the last one
                //                {
                //                    //choose the father population from which the population i came
                //                    fatherPop= ChooseFatherPopulation(  pop, random,  0, K);
                //                    pop->FatherPop = fatherPop;
                //                    //update list of migrant times
                //                    Population::UpdateListMigrants( numClones, pop, fatherPop);
                //                }
                std::cout << "weird"<< std::endl;
                return (timeNextMigration*currentPop->x);
            }
        }
        
    }
    //            if ( params.doFixedEventimes){
    //                timeNextCoalEvent = params.getPopulationEvent(i, result->getIdNextCoalEventForPopulation(i));
    //            }
    //            else{
    ////                timeNextCoalEvent = pop->nextCoalEventTime(result->getIdNextCoalEventForPopulation(i), result->getIdNextInmigrationEventForPopulation(i),   result->getHeightModelTime() ,isThereInmigration,  inmigrantPop );
    //            }
    //this line would never be reached
    std::cout << "super weird"<< std::endl;
    return -1.0;
}
void  PopulationSet::acceptNextCoalEventTime(gsl_rng *random, double receiverModelTime, int idxReceiverPop, int idxIncomingPop,  double K){
    
    Population *receiverPop = getPopulationbyIndex( idxReceiverPop);
    
    
    receiverPop->setCurrentModelTime(receiverModelTime);
    if (idxReceiverPop == idxIncomingPop ){//coalescent
        
       
        receiverPop->CoalescentEventTimes[receiverPop->numCompletedCoalescences]=  receiverPop->getCurrentModelTime();
        receiverPop->numCompletedCoalescences= receiverPop->numCompletedCoalescences +1;
        
        if ( receiverPop->numActiveGametes <= 1 && idxReceiverPop< numClones-1)   //if it is not the last one
        {
            Population *fatherPop;
            //choose the father population from which the population i came
            fatherPop= ChooseFatherPopulation(  receiverPop, random,  0, K);
            receiverPop->FatherPop = fatherPop;
            //update list of migrant times
            Population::UpdateListMigrants( numClones, receiverPop, fatherPop);
        }
    }
    else{
        Population *incomingPop = getPopulationbyIndex( idxIncomingPop);
      
        receiverPop->indexNextInmigrant = receiverPop->indexNextInmigrant+1;
        incomingPop->numActiveGametes =  incomingPop->numActiveGametes-1;
    }
    
}
std::pair<Pair, std::pair<double, double>> PopulationSet::initPairProposalsNextEventTimePerPopulation(gsl_rng *rngGsl, double K){
    
    Pair pairMinModelTime;
    Population *popI = getCurrentPopulation();
    
    std::pair<Pair, std::pair<double, double>> pairMinModelTimeCurrPop;
    
    if (popI->numActiveGametes > 1)
        pairMinModelTimeCurrPop = popI->initPairProposalsNextEventTime(rngGsl, K);
    
    
//    for (unsigned int i = 0; i < numClones; ++i){
//        auto popI =  populations[i];
//        if (popI->numActiveGametes > 1){
//
//            Pair pairMinModelTimeCurrPop = popI->initPairProposalsNextEventTime(rngGsl, K);
//            if (pairMinModelTimeCurrPop
//
//        }
//
//
//    }

    return pairMinModelTimeCurrPop;
}
void PopulationSet::initDeltaThetaFromPriors( const gsl_rng *rngGsl, long double &theta){
    long double delta;
    
    theta =  Random::RandomExponential(1, NULL, true, rngGsl, NULL);
    
    for (unsigned int i = 0; i < numClones; ++i){
        auto popI =  populations[i];
        delta = Random::RandomExponential(0.01, NULL, true, rngGsl, NULL);
        
        popI->theta = theta * popI->x;
        popI->delta = delta;
    }
}
void PopulationSet::initTheta( long double &theta){
    assert(theta>0);
    for (unsigned int i = 0; i < numClones; ++i){
        auto popI =  populations[i];
        assert(popI->x>0);
        popI->theta = theta * popI->x;
        
    }
}
PopulationSet::PopulationSet(const PopulationSet &original)
//: populations(original.populations.size())
{
    numClones = original.numClones;
    proportionsVector = original.proportionsVector;
    populations.reserve(original.populations.size());
    for (auto p : original.populations)
        populations.push_back(new Population(*p));
    //    try {
    //          // auto is preferable here
    //          std::vector<Population*>::iterator thisit = populations.begin();
    //           std::vector<Population*>::const_iterator thatit = original.populations.cbegin();
    //
    //          for (; thatit != original.populations.cend(); ++thisit, ++thatit)
    //              *thisit = new Population(**thatit);
    //      } catch (...) {
    //          for (std::vector<Population*>::iterator i = populations.begin(); i != populations.end(); ++i)
    //              if (!*i)
    //                  break;
    //              else
    //                  delete *i;
    //
    //          throw;
    //      }
    
}
void  PopulationSet::initListPossibleMigrations()
{
    sort(populations.begin(), populations.end(), comparePopulationsByTimeOrigin);
    Population *p;
    int i;
    for (i = 0; i < numClones; i++) {
        p = populations[i];
        p->InitListPossibleMigrations(i);
    }
}
void PopulationSet::initPopulationsCoalescentEvents( )
{
    int i;
    Population *popI;
    for( i = 0 ; i < numClones; i++)
    {
        popI = populations[i];
        popI->InitCoalescentEvents(numClones);
    }
}
void PopulationSet::initProportionsVector(){
    for (size_t i = 0; i < numClones; i++) {
        proportionsVector.push_back(0);
        oldproportionsVector.push_back(0);
    }
    
}
void PopulationSet::sortPopulationsByTorigin(){
    
    sort(populations.begin(), populations.end(), comparePopulationsByTimeOrigin);
    //qsort(populations, numClones, sizeof(Population), comparePopulationsByTimeOrigin);
    
}
Population* PopulationSet::ChooseFatherPopulation( Population  *PopChild, const gsl_rng *randomGenerator, int noisy, long double K)
{
    
    Population  *p;
    long double  pij, ran;
    int  j, k;
    long double  cumProb[numClones  - (int)(PopChild->order)];
    cumProb[0] = 0.0;
    for (j = PopChild->order + 1; j < numClones; j++)
    {
        cumProb[j - PopChild->order] = 0.0;
        p = populations[ j];
        pij= PopChild->ProbabilityComeFromPopulation ( p, populations,  numClones, K);
        //pij = ProbabilityCloneiFromClonej2(PopChild, p, populations, numClones);
        cumProb[j - PopChild->order] = cumProb[j - 1 - PopChild->order] + pij;
        
    }
    // now selecting the ancestral clone
    ran =  Random::randomUniformFromGsl2(randomGenerator);
    //fprintf (stderr, "\n ran = %lf ", ran);
    int w = -1;
    for (k = 1; k < numClones ; k++)
    {
        if (ran <= cumProb[k])
        {
            w = k;
            break;
        }
    }
    Population *result =  populations[ PopChild->order + w];
    if (noisy > 3)
        fprintf (stderr, "\nClone %d derived from clone %d\n", PopChild->index, result->index); // clone ThisCloneNumber (i) is originated from clone ThisOriginCloneNumber (j)
    /* Update list of migation times considering that clone i comes from clone j */
    if (noisy > 1)
        fprintf (stderr, "\n*** Updating list of migration times (considering that clone %d comes from clone %d) ..\n", PopChild->index, result->index);
    return (result); //the father population has  order  (PopChild->order) + w
}
void PopulationSet::AssignSequencesToPopulations(std::vector<pll_rnode_t*> rnodes,
                                                 ProgramOptions &programOptions,
                                                 int noisy,  int TotalTumorSequences,
                                                 int &numActiveGametes, int &nextAvailable,
                                                 int &labelNodes, char* ObservedCellNames[], int doUseObservedCellNames, std::vector<int> &sampleSizes)

{
    Population *pop;
    pll_rnode_t *p;
    TreeNode *u;
    int i, j;
    long double  currentSampleSize;
    numActiveGametes = 0;
    int indexFirstObservedCellName;
    
    std::vector<int> CumSumNodes(numClones+1);
    CumSumNodes[0] = 0;
    
    for (i = 1; i <= numClones; i++)
    {
        pop = populations[i-1];
        pop->resetGametesCounters();
        currentSampleSize = pop->sampleSize;
        CumSumNodes[i] = CumSumNodes[i - 1] + currentSampleSize;
    }
    if (noisy > 1)
        fprintf (stderr, "\n Initial relation nodes-clones:");
    int  cumIndivid = 0;
    int currentPopIndex;
    numActiveGametes = 0;
    
    for (j = 0; j < TotalTumorSequences; j++)
    {
        cumIndivid++;
        p = rnodes[j];
        p->node_index = j;
        
        labelNodes = j;
        
        if (p->data == NULL)
            p->data = new TreeNode(0);
        u= (TreeNode *)(p->data);
        u->nodeClass = 1;
        
        for (i = 1; i <= numClones; i++)
        {
            pop = populations[ i - 1];
            currentPopIndex = pop->index;
            indexFirstObservedCellName= pop->indexFirstObservedCellName;
            // Identify to which clone belongs this node (sample)
            if (cumIndivid <= CumSumNodes[i] && cumIndivid > CumSumNodes[i - 1])
            {
                pop->idsActiveGametes[pop->numActiveGametes]=j;
                pop->rtips[pop->numActiveGametes]=p;
                pop->numActiveGametes=pop->numActiveGametes+1;
                pop->idsGametes[pop->numGametes]=j;
                pop->numGametes=pop->numGametes+1;
                
                u->indexOldClone = currentPopIndex;
                u->indexCurrentClone = currentPopIndex;
                u->orderCurrentClone = pop->order;
                break;
            }
        }
        if (noisy > 1)
            fprintf (stderr,"\n > The node %d(%d) belongs to clone %d", u->index, cumIndivid, u->indexOldClone);
        numActiveGametes = numActiveGametes + 1;
    }
    //AssignObservedCellNamestoTips(nodes, populations, sampleSizes,ObservedCellNames,  programOptions);
    // AssignObservedCellNamestoTips2(nodes, populations, sampleSizes,ObservedCellNames,  programOptions);
    
    CumSumNodes.clear();
    nextAvailable = numActiveGametes;
    labelNodes = labelNodes + 1;
}
std::vector <long double> PopulationSet::getDeltaTs()
{
    std::vector <long double> result(numClones);
    Population *p;
    int i;
    for (i = 0; i < numClones; i++) {
        p = populations[i];
        result.at(i)=p->deltaT;
    }
    return result;
    
}
std::vector <long double> PopulationSet::getDeltas()
{
    std::vector <long double> result(numClones);
    Population *p;
    int i;
    for (i = 0; i < numClones; i++) {
        p = populations[i];
        result.at(i)=p->delta;
    }
    return result;
    
}
std::vector <long double> PopulationSet::getTs()
{
    std::vector <long double> result(numClones);
    Population *p;
    int i;
    for (i = 0; i < numClones; i++) {
        p = populations[i];
        result.at(i)=p->timeOriginSTD;
    }
    return result;
    
}
std::vector <long double> PopulationSet::getSampleSizes()
{
    std::vector <long double> result(numClones);
    Population *p;
    int i;
    for (i = 0; i < numClones; i++) {
        p = populations[i];
        result[i] = p->sampleSize;
    }
    return result;
    
}
void PopulationSet::initPopulationDeltas(std::vector<double> &deltas){
    
    assert(numClones==deltas.size());
    Population *p;
    int i;
    for (i = 0; i < numClones; i++) {
        p = populations[i];
        assert(deltas[i] >0);
        p->setDelta(deltas[i]);
    }
    
}
void PopulationSet::initPopulationTOriginSTD(std::vector<double> &TOriginSTDs){
    assert(numClones==TOriginSTDs.size());
    Population *p;
    int i;
    for (i = 0; i < numClones; i++) {
        p = populations[i];
        assert(TOriginSTDs[i] >0);
        p->setTimeOriginSTD(TOriginSTDs[i]);
    }
    
}
void PopulationSet::sampleEventTimesScaledByProportion(gsl_rng * random, double K, int noisy){
    
    Population *currentPop;
    Population *fatherPop;
    int i=0;
    
    while (i < numClones) {
        //currentPop = *(populations + i);
        currentPop = populations[i];
        assert(currentPop->numCompletedCoalescences==0);
        currentPop->sampleEventTimesScaledByProportion(random, K);
        
        currentPop->printCoalEventTimes(std::cout);
        if (i< numClones-1)   //if it is not the last one
        {
            //choose the father population from which the population i came
            fatherPop= ChooseFatherPopulation(  currentPop, random,  noisy, K);
            currentPop->FatherPop = fatherPop;
            //update list of migrant times
            Population::UpdateListMigrants( numClones, currentPop, fatherPop);
        }
        i = i + 1;
    }
}
Population* PopulationSet::getCurrentPopulation(){
    Population * popI;
    for(size_t i = 0; i != numClones; ++i)
    {
        popI = populations[i];
        if ( popI->currentModelTime < popI->timeOriginSTD)
        {
            return populations[i];
            
        }
    }
    return popI;
}
/////////////////////////////////////////////////////////////////////////////////////////////////
StructuredCoalescentTree::StructuredCoalescentTree(int numClones, std::vector<int> &sampleSizes, long double theta,
                                                   MCMCoptions &mcmcOptions,ProgramOptions &programOptions, const gsl_rng *  rngGsl, boost::mt19937 *rngBoost,
                                                   std::vector<std::vector<int> > &ObservedData,char* ObservedCellNames[], pll_msa_t *msa, std::string& healthyTipLabel, long double seqErrorRate,
                                                   long double dropoutRate){
    if (seqErrorRate >=0)
        this->seqErrorRate = seqErrorRate;
    else
        fprintf (stderr, "\n ERROR: The sequencing error rate cannot be negative \n");
    
    if (dropoutRate >=0)
        this->dropoutRate = dropoutRate;
    else
        fprintf (stderr, "\n ERROR: The dropout error rate cannot be negative \n");
    
    
    if (numClones >0)
        this->numClones = numClones;
    else
        fprintf (stderr, "\n ERROR: The number of populations cannot be negative \n");
    
    this->populationSet = new PopulationSet(numClones);
    assert(sampleSizes.size()==numClones);
    copy(sampleSizes.begin(), sampleSizes.end(), back_inserter(this->sampleSizes));
    
    auto deltaTs  = populationSet->samplePopulationGrowthRateFromPriors(mcmcOptions, rngGsl );
    
    if (theta >0)
        this->theta = theta;
    else
        fprintf (stderr, "\n ERROR: Theta cannot be negative \n");
    
    populationSet->initPopulationSampleSizes(sampleSizes);
    
    populationSet->initPopulationGametes();
    populationSet->initPopulationRTips();
    populationSet->initProportionsVectorFromSampleSizes(sampleSizes);
    populationSet->initPopulationsThetaDelta( theta);
    int numberTips =0;
    for(std::vector<int>::iterator it = sampleSizes.begin(); it != sampleSizes.end(); ++it)
        numberTips += *it;
    
    Population* oldestPop=populationSet->getPopulationbyIndex(numClones -1);
    
    oldestPop->timeOriginSTD  =   Random::RandomDensityModelTimeOrigin (oldestPop->delta, mcmcOptions.useGSLRandomGenerator, 0, rngGsl, rngBoost );
    oldestPop->TimeOriginSTD->setParameterValue(oldestPop->timeOriginSTD);
    std::cout << "\n True time of origin oldest pop: " << oldestPop->timeOriginSTD << std::endl;
    
    initEdgesRootedTree(numberTips);
    currentNumberEdgesRootedTree=0;
    
    populationSet->initListPossibleMigrations();
    populationSet->initPopulationsCoalescentEvents();
    
    healthyTip=NULL;
    
    root =   MakeCoalescenceTree (rngGsl,rngBoost,
                                  msa,
                                  programOptions.numNodes,
                                  programOptions,
                                  ObservedData,
                                  ObservedCellNames,
                                  sampleSizes
                                  ) ;
    
    
    rtree = pll_rtree_wraptree(root,programOptions.numCells);
    initLogLikelihoodSequences(programOptions, msa);
    
}
StructuredCoalescentTree::StructuredCoalescentTree(PopulationSet *populationSet, std::vector<int> &sampleSizes,
                                                   long double theta, MCMCoptions &mcmcOptions,ProgramOptions &programOptions,const  gsl_rng *rngGsl,
                                                   boost::mt19937 *rngBoost,
                                                   std::vector<std::vector<int> > &ObservedData,char* ObservedCellNames[], pll_msa_t *msa, std::string& healthyTipLabel, long double seqErrorRate,
                                                   long double dropoutRate )
{
    if (seqErrorRate >=0)
        this->seqErrorRate = seqErrorRate;
    else
        fprintf (stderr, "\n ERROR: The sequencing error rate cannot be negative \n");
    
    if (dropoutRate >=0)
        this->dropoutRate = dropoutRate;
    else
        fprintf (stderr, "\n ERROR: The dropout error rate cannot be negative \n");
    
    if (populationSet->numClones >0)
    {
        this->numClones=populationSet->numClones;
        this->populationSet = populationSet;
        
    }
    else
        fprintf (stderr, "\n ERROR: The number of clones cannot be negative \n");
    
    assert(sampleSizes.size()==numClones);
    
    auto deltaTs  = populationSet->samplePopulationGrowthRateFromPriors(mcmcOptions, rngGsl );
    
    if (theta >0)
        this->theta = theta;
    else
        fprintf (stderr, "\n ERROR: Theta cannot be negative \n");
    
    
    copy(sampleSizes.begin(), sampleSizes.end(), back_inserter(this->sampleSizes));
    populationSet->initPopulationSampleSizes(sampleSizes);
    
    populationSet->initPopulationGametes();
    populationSet->initPopulationRTips();
    populationSet->initProportionsVectorFromSampleSizes(sampleSizes);
    populationSet->initPopulationsThetaDelta( theta);
    int numberTips =0;
    for(std::vector<int>::iterator it = sampleSizes.begin(); it != sampleSizes.end(); ++it)
        numberTips += *it;
    Population* oldestPop=populationSet->getPopulationbyIndex(numClones -1);
    
    oldestPop->timeOriginSTD  =   Random::RandomDensityModelTimeOrigin (oldestPop->delta, mcmcOptions.useGSLRandomGenerator, 0, rngGsl, NULL );
    
    
    std::cout << "\n True time of origin oldest pop: " << oldestPop->timeOriginSTD << std::endl;
    initEdgesRootedTree(numberTips);
    currentNumberEdgesRootedTree=0;
    
    populationSet->initListPossibleMigrations();
    populationSet->initPopulationsCoalescentEvents();
    healthyTip=NULL;
    
    root =   MakeCoalescenceTree (rngGsl, rngBoost,
                                  msa,
                                  programOptions.numNodes,
                                  programOptions,
                                  ObservedData,
                                  ObservedCellNames,
                                  sampleSizes) ;
    
    rtree = pll_rtree_wraptree(root,programOptions.numCells);
}
pll_rnode_t * StructuredCoalescentTree::MakeCoalescenceTree (
                                                             const gsl_rng *rngGsl,
                                                             boost::mt19937 *rngBoost,
                                                             pll_msa_t *msa,
                                                             int &numNodes,
                                                             ProgramOptions &programOptions,
                                                             std::vector<std::vector<int> > &ObservedData,
                                                             char* ObservedCellNames[],
                                                             std::vector<int> &sampleSizes
                                                             ) {
    
    int      c,  i,  m, cumIndivid, isCoalescence, whichInd,
    newInd, eventNum, numActiveGametes,
    isMigration, whichClone, currentNumberAliveClones;
    int     labelNodes;
    long double     currentTime, eventTime;
    pll_rnode_t  *p;
    TreeNode *t;
    
    int     *numParcialActiveGametes;
    int         *fromClone;
    int         ThisCloneNumber, ThisOriginCloneNumber;
    long double       minValue;
    long double       ThisRateCA;
    long double       ThisTimeCA_W;
    long double       ThisTimeCA_V1;
    long double       ThisTimeCA_V2;
    
    int         doAmigration;
    int         ThisCloneNumberMigrations, ThisM;
    int nextAvailable;
    Population *pop;
    std::string nodeLabel;
    
    /* defaults */
    isCoalescence = NO;
    isMigration = NO;
    newInd = whichClone = labelNodes = cumIndivid = 0;
    numParcialActiveGametes = NULL;
    fromClone = NULL;
    eventTime = 0.0;
    c = m = 0;
    whichInd = 0;
    minValue = 0.0;
    ThisCloneNumber = 0;
    ThisOriginCloneNumber = 0;
    ThisRateCA = 0.0;
    ThisTimeCA_W = 0.0;
    ThisTimeCA_V1 = 0.0;
    ThisTimeCA_V2 = 0.0;
    doAmigration = -1;
    ThisCloneNumberMigrations = -1;
    ThisM = -1;
    std::vector<Population *>& populations =  populationSet->getPopulations();
    
    
    //*numNodes = 2 * TotalNumSequences * numClones+ 1; // (2 * TotalNumSequences) + numClones(superfluos) - 1, but let's allocate some more..
    currentNumberAliveClones = numClones;
    
    //InitListPossibleMigrations(populations, numClones);
    for (i = 0; i < programOptions.numClones; i++){
        pop =populations[i];
        pop->sampleSize= sampleSizes[i];
        pop->resetMigrationsList();
    }
    //resetMigrationsList( populations,  numClones);
    
    for (i=0; i< msa->count; i++){
        p=new pll_rnode_t();
        p->left =NULL;
        p->right=NULL;
        p->parent=NULL;
        p->label =new char[strlen(msa->label[i]) + 1]{};
        std::copy(msa->label[i], msa->label[i] + strlen(msa->label[i]), p->label);
        
        nodeLabel = p->label;
        
        t = new TreeNode(msa->length);
        t->initNumberTipsVector(programOptions.numClones);
        //strcpy( t->observedCellName,p->label);
        
        std::copy(p->label, p->label + strlen(p->label), t->observedCellName);
        
        t->genotypeSequence = ObservedData[i];
        p->data = t;
        if (std::strcmp(p->label, programOptions.healthyTipLabel.c_str())==0)
            healthyTip = p;
        else{
            rtreeTips.push_back(p);
            rnodes.push_back(p);
        }
        
    }
    
    for (i = 0; i < (numNodes - programOptions.TotalTumorSequences); i++)
    {
        p = new pll_rnode_t();
        p->left =NULL;
        p->right=NULL;
        p->parent=NULL;
        t = new TreeNode(0);
        t->initNumberTipsVector(programOptions.numClones);
        p->data = t;
        rnodes.push_back(p);
    }
    
    
    
    populationSet->AssignSequencesToPopulations( rnodes, programOptions,  programOptions.noisy, programOptions.TotalTumorSequences, numActiveGametes,  nextAvailable,
                                                labelNodes, ObservedCellNames, programOptions.doUseObservedCellNames, sampleSizes);
    Population *currentPop;
    Population *fatherPop;
    i=0;
    currentTime=0.0;
    while (i < numClones) {
        //currentPop = *(populations + i);
        currentPop = populations[i];
        assert(currentPop->numCompletedCoalescences==0);
        SimulatePopulation( *currentPop,  programOptions, rngGsl,rngBoost,
                           programOptions.numNodes,
                           numClones,
                           nextAvailable ,
                           numActiveGametes,
                           labelNodes,
                           currentTime,
                           eventNum,
                           programOptions.K);
        if (i< numClones-1)   //if it is not the last one
        {
            //choose the father population from which the population i came
            fatherPop= populationSet->ChooseFatherPopulation(  currentPop, rngGsl,  programOptions.noisy, programOptions.K);
            currentPop->FatherPop = fatherPop;
            //update list of migrant times
            Population::UpdateListMigrants( numClones, currentPop, fatherPop);
        }
        i = i + 1;
    }
    pll_rnode_t *root=NULL;
    
    root =  BuildTree( currentPop,rngGsl,
                      programOptions,
                      root,
                      nextAvailable,
                      newInd,
                      currentTime,
                      labelNodes
                      );
    
    
    // RelabelNodes(treeRootInit[0], treeRootInit, &intLabel );
    return root;
}
void StructuredCoalescentTree::initEdgesRootedTree(int numberTips)
{
    int i;
    assert(numberTips >0);
    edges.clear();
    edgeLengths.clear();
    
    
    //edgeLengths (2*numberTips-2, std::make_pair(0.0, NULL));
    for( i = 0 ; i < 2*numberTips-2; i++)
    {
        
        edges.push_back(NULL);
        edgeLengths.push_back(std::make_pair(0.0, nullptr));
    }
    
}
void StructuredCoalescentTree::initLogLikelihoodSequences(const ProgramOptions &programOptions, const pll_msa_t *msa){
    
    long double result=0;
    long double currentlogSumDensitiesTimeOriginSTDPopulations = SumLogDensitiesTimeOriginSTDPopulations();
    result=currentlogSumDensitiesTimeOriginSTDPopulations;
    long double currentlogProbFatherPopulations = SumLogProbFatherPopulations(programOptions.K);
    result= result + currentlogProbFatherPopulations;
    long double currentlogDensityCoalescentTimesForPopulation= LogDensityCoalescentTimesForPopulation(programOptions.K);
    
    result = result + currentlogDensityCoalescentTimesForPopulation;
    logLikelihoodTree= result;
    
    std::cout << "True log likelihood of the tree  is " << logLikelihoodTree << std::endl;
    // if (mcmcOptions.useSequencesLikelihood ==1)
    //  {
    char *  rootedNewick2 = pll_rtree_export_newick( root, NULL);
    //initialRootedTree = pll_rtree_parse_newick_string(rootedNewick3);
    
    //logLikelihoodSequences= LogConditionalLikelihoodSequences( msa,  rootedNewick2, programOptions, seqErrorRate, dropoutRate);
    
    // std::cout << "True log likelihood of the sequences  is " << logLikelihoodTree << std::endl;
    free(rootedNewick2);
    
}
long double StructuredCoalescentTree::SumLogDensitiesTimeOriginSTDPopulations() {
    Population* popI;
    
    long double temp;
    long double result = 0;
    unsigned int i;
    
    for ( i = 0; i < numClones; i++)
    {
        popI= populationSet->getPopulationbyIndex(i);
        temp=popI->LogDensityTimeSTDFrom( popI->timeOriginSTD, 0);
        //if (isnan(temp) || isinf(temp))
        //   fprintf (stderr, "\n isNan temp LogDensityTime, Delta: %Lf, T:%.20Lf \n",popI->delta ,popI->timeOriginSTD);
        result = result + temp;
        // fprintf (stderr, "\n Product log Density Time   = %lf after pop order %d,  popI->LogDensityTime( popI->timeOriginSTD) %lf,  popI->timeOriginSTD: %lf, popI->delta: %lf\n", product, popI->order, popI->LogDensityTime( popI->timeOriginSTD), popI->timeOriginSTD,popI->delta );
    }
    return result;
}

long double StructuredCoalescentTree::SumLogDensitiesTimeOriginSTDPopulations(std::vector<Population *> populations, int numClones, long double K) {
    Population* popI;
    
    long double temp;
    long double result = 0;
    unsigned int i;
    
    for ( i = 0; i < numClones; i++)
    {
        popI= populations[i];
        temp=popI->LogDensityTimeSTDFrom( popI->timeOriginSTD, 0);
        //if (isnan(temp) || isinf(temp))
        //   fprintf (stderr, "\n isNan temp LogDensityTime, Delta: %Lf, T:%.20Lf \n",popI->delta ,popI->timeOriginSTD);
        result = result + temp;
        // fprintf (stderr, "\n Product log Density Time   = %lf after pop order %d,  popI->LogDensityTime( popI->timeOriginSTD) %lf,  popI->timeOriginSTD: %lf, popI->delta: %lf\n", product, popI->order, popI->LogDensityTime( popI->timeOriginSTD), popI->timeOriginSTD,popI->delta );
    }
    return result;
}

long double StructuredCoalescentTree::SumLogProbFatherPopulations(long double K) {
    
    Population *fatherPop;
    Population *popJ;
    long double  result=0;
    long double temp;
    for ( unsigned int j = 0; j < numClones-1 ; j++)
    {
        popJ = populationSet->getPopulationbyIndex(j);
        if (popJ->FatherPop !=NULL)
        {
            fatherPop = popJ->FatherPop;
            //result = result + log( fatherPop->popSize);
            temp=popJ->timeOriginSTD * popJ->x / fatherPop->x;
            temp=Population::LogCalculateH(temp, fatherPop->timeOriginSTD, fatherPop->delta, K);
            result = result  + temp;
            // fprintf (stderr, "\n Product calculate H    = %lf after pop order %d \n", product, popJ->order);
        }
        if(popJ->order < numClones-1 && popJ->FatherPop == NULL)
            fprintf (stderr, "\n the pop with order %d has father pop null\n",popJ->order );
    }
    return(result);
}
long double StructuredCoalescentTree::SumLogProbFatherPopulations(std::vector<Population *> populations, int numClones, long double K) {
    
    Population *fatherPop;
    Population *popJ;
    long double  result=0;
    long double temp;
    for ( unsigned int j = 0; j < numClones-1 ; j++)
    {
        popJ = populations[j];
        if (popJ->FatherPop !=NULL)
        {
            fatherPop = popJ->FatherPop;
            //result = result + log( fatherPop->popSize);
            temp=popJ->timeOriginSTD * popJ->x / fatherPop->x;
            temp=Population::LogCalculateH(temp, fatherPop->timeOriginSTD, fatherPop->delta, K);
            result = result  + temp;
            // fprintf (stderr, "\n Product calculate H    = %lf after pop order %d \n", product, popJ->order);
        }
        if(popJ->order < numClones-1 && popJ->FatherPop == NULL)
            fprintf (stderr, "\n the pop with order %d has father pop null\n",popJ->order );
    }
    return(result);
}
long double  StructuredCoalescentTree::LogDensityCoalescentTimesForPopulation(long double K){
    Population *popI;
    long double result =0;
    long double temp;
    long double timeCurrentEvent;
    long int numberAliveCells;
    long int currentCoalescentEventInThisEpoch=0;
    long int currentMigrationEvent=0;
    long double termOnlyAfterFirstCoalEvent;
    long double lastEventTimeBeforeMigration=0;
    std::vector<long double> immigrantsTimes;
    
    // vector<double> allEventsSorted;
    for ( unsigned i = 0; i < numClones; i++)
    {
        popI = populationSet->getPopulationbyIndex(i);
        numberAliveCells = popI->sampleSize;
        currentCoalescentEventInThisEpoch=0;
        currentMigrationEvent=0;
        lastEventTimeBeforeMigration=0;
        
        immigrantsTimes.clear();
        std::transform(popI->immigrantsPopOrderedByModelTime.begin(), popI->immigrantsPopOrderedByModelTime.end(),
                       std::back_inserter(immigrantsTimes),
                       [](auto const& pair){ return pair.first; });
        std::vector<long double> allEventsSorted(popI->CoalescentEventTimes.size()+immigrantsTimes.size(), 0.0);
        merge(popI->CoalescentEventTimes.begin(), popI->CoalescentEventTimes.end(), immigrantsTimes.begin(), immigrantsTimes.end(), allEventsSorted.begin());
        bool allCoalTimesPositive = Utils::checkAllElementsGreaterThanZero(popI->CoalescentEventTimes, popI->sampleSize-1);
        if (!allCoalTimesPositive){
            fprintf (stderr, "\nERROR: Not all coalescent times are positive \n\n");
            fprintf (stderr, "\nERROR: Number of completed coalescent events %d \n", popI->numCompletedCoalescences);
            //exit(1);
        }
        else{
            
            //fprintf (stderr, "\n All coalescent times are positive \n\n");
            
        }
        //assert(std::find(popI->CoalescentEventTimes.begin(), popI->CoalescentEventTimes.end(), 0) == popI->CoalescentEventTimes.end());
        
        for ( unsigned j = 0; j < allEventsSorted.size(); j++)
        {
            timeCurrentEvent= allEventsSorted.at(j);
            if (timeCurrentEvent== popI->timeOriginSTD)
                break;
            if(std::binary_search(popI->CoalescentEventTimes.begin(), popI->CoalescentEventTimes.end(), timeCurrentEvent) )//is a coalescent event
            {
                if (numberAliveCells > 1)
                {
                    temp=log(numberAliveCells * (numberAliveCells-1.0)/2.0);
                    result= result + temp;
                    temp =  Population::LogLambda(timeCurrentEvent,popI->timeOriginSTD, popI->delta, K);
                    result= result + temp;
                    termOnlyAfterFirstCoalEvent =(currentCoalescentEventInThisEpoch == 0)?Population::FmodelTstandard(lastEventTimeBeforeMigration, popI->timeOriginSTD, popI->delta, K):Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEventInThisEpoch-1], popI->timeOriginSTD, popI->delta, K);//if no coalescent
                    temp =  (numberAliveCells / 2.0)* (numberAliveCells - 1.0)*(Population::FmodelTstandard(timeCurrentEvent,popI->timeOriginSTD, popI->delta, K)-termOnlyAfterFirstCoalEvent);
                    result= result - temp;
                    lastEventTimeBeforeMigration = timeCurrentEvent;
                    currentCoalescentEventInThisEpoch++;
                    numberAliveCells--;
                }
            }
            else if(std::binary_search(immigrantsTimes.begin(),
                                       immigrantsTimes.end(), timeCurrentEvent) )//is a migration event
            {
                if (numberAliveCells > 1)
                {
                    temp= popI->LogProbNoCoalescentEventBetweenTimes(lastEventTimeBeforeMigration,timeCurrentEvent, numberAliveCells, K );
                    result= result+ temp;
                }
                lastEventTimeBeforeMigration=timeCurrentEvent;
                currentMigrationEvent++;
                numberAliveCells++;
                currentCoalescentEventInThisEpoch=0;//new epoch is starting
            }
        }
    }
    // fprintf (stderr, "\n Result = %lf \n", result);
    return result;
}
long double  StructuredCoalescentTree::LogDensityCoalescentTimesForPopulation(std::vector<Population *> populations, int numClones, long double K){
    Population *popI;
    long double result =0;
    long double temp;
    long double timeCurrentEvent;
    long int numberAliveCells;
    long int currentCoalescentEventInThisEpoch=0;
    long int currentMigrationEvent=0;
    long double termOnlyAfterFirstCoalEvent;
    long double lastEventTimeBeforeMigration=0;
    std::vector<long double> immigrantsTimes;
    
    // vector<double> allEventsSorted;
    for ( unsigned i = 0; i < numClones; i++)
    {
        popI = populations[i];
        numberAliveCells = popI->sampleSize;
        currentCoalescentEventInThisEpoch=0;
        currentMigrationEvent=0;
        lastEventTimeBeforeMigration=0;
        
        immigrantsTimes.clear();
        std::transform(popI->immigrantsPopOrderedByModelTime.begin(), popI->immigrantsPopOrderedByModelTime.end(),
                       std::back_inserter(immigrantsTimes),
                       [](auto const& pair){ return pair.first; });
        std::vector<long double> allEventsSorted(popI->CoalescentEventTimes.size()+immigrantsTimes.size(), 0.0);
        merge(popI->CoalescentEventTimes.begin(), popI->CoalescentEventTimes.end(), immigrantsTimes.begin(), immigrantsTimes.end(), allEventsSorted.begin());
        bool allCoalTimesPositive = Utils::checkAllElementsGreaterThanZero(popI->CoalescentEventTimes, popI->sampleSize-1);
        if (!allCoalTimesPositive){
            fprintf (stderr, "\nERROR: Not all coalescent times are positive \n");
            fprintf (stderr, "\nERROR: Th sample size is %d\n",popI->sampleSize );
            fprintf (stderr, "\nERROR: Number of completed coalescent events %d \n", popI->numCompletedCoalescences);
            //exit(1);
        }
        else{
            
            //fprintf (stderr, "\n All coalescent times are positive \n\n");
            
        }
        //assert(std::find(popI->CoalescentEventTimes.begin(), popI->CoalescentEventTimes.end(), 0) == popI->CoalescentEventTimes.end());
        
        for ( unsigned j = 0; j < allEventsSorted.size(); j++)
        {
            timeCurrentEvent= allEventsSorted.at(j);
            if (timeCurrentEvent== popI->timeOriginSTD)
                break;
            if(std::binary_search(popI->CoalescentEventTimes.begin(), popI->CoalescentEventTimes.end(), timeCurrentEvent) )//is a coalescent event
            {
                if (numberAliveCells > 1)
                {
                    temp=log(numberAliveCells * (numberAliveCells-1.0)/2.0);
                    result= result + temp;
                    temp =  Population::LogLambda(timeCurrentEvent,popI->timeOriginSTD, popI->delta, K);
                    result= result + temp;
                    termOnlyAfterFirstCoalEvent =(currentCoalescentEventInThisEpoch == 0)?Population::FmodelTstandard(lastEventTimeBeforeMigration, popI->timeOriginSTD, popI->delta, K):Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEventInThisEpoch-1], popI->timeOriginSTD, popI->delta, K);//if no coalescent
                    temp =  (numberAliveCells / 2.0)* (numberAliveCells - 1.0)*(Population::FmodelTstandard(timeCurrentEvent,popI->timeOriginSTD, popI->delta, K)-termOnlyAfterFirstCoalEvent);
                    result= result - temp;
                    lastEventTimeBeforeMigration = timeCurrentEvent;
                    currentCoalescentEventInThisEpoch++;
                    numberAliveCells--;
                }
            }
            else if(std::binary_search(immigrantsTimes.begin(),
                                       immigrantsTimes.end(), timeCurrentEvent) )//is a migration event
            {
                if (numberAliveCells > 1)
                {
                    temp= popI->LogProbNoCoalescentEventBetweenTimes(lastEventTimeBeforeMigration,timeCurrentEvent, numberAliveCells, K );
                    result= result+ temp;
                }
                lastEventTimeBeforeMigration=timeCurrentEvent;
                currentMigrationEvent++;
                numberAliveCells++;
                currentCoalescentEventInThisEpoch=0;//new epoch is starting
            }
        }
    }
    // fprintf (stderr, "\n Result = %lf \n", result);
    return result;
}
void StructuredCoalescentTree::MakeCoalescenceEvent( Population &population,const gsl_rng *randomGenerator, int noisy,   int &numActiveGametes, int &nextAvailable,
                                                    int &labelNodes, long double  &currentTime, int &numNodes)
{
    
    // TreeNode  *p, *q, *r;
    pll_rnode_t  *p, *q,  *r ;
    TreeNode * tp, *tq, *tr ;
    int firstInd,  secondInd=0, newInd=0;
    int choosePairIndividuals = YES;
    
    population.ChooseRandomIndividual(&firstInd, numClones,   &secondInd, randomGenerator, choosePairIndividuals);
    
    newInd = nextAvailable;
    if (noisy > 1)
        fprintf (stderr, "Coalescence involving %d and %d to create node %d (in clone %d)", population.idsActiveGametes[firstInd], population.idsActiveGametes[secondInd], newInd, population.index);
    /*  set pointers between nodes */
    assert(firstInd< population.idsActiveGametes.size() );
    assert(secondInd< population.idsActiveGametes.size() );
    p = rnodes[population.idsActiveGametes[firstInd]];
    q = rnodes[population.idsActiveGametes[secondInd]];
    
    r = rnodes[newInd];
    
    r->node_index = nextAvailable;
    
    tr= (TreeNode *)(r->data);
    tp= (TreeNode *)(p->data);
    tq= (TreeNode *)(q->data);
    
    tr->index =   nextAvailable;
    tr->label = labelNodes;
    labelNodes=labelNodes+1;
    tr->indexOldClone = tr->indexCurrentClone = population.index;//here the clone number is updated
    
    tr->orderCurrentClone =population.order;
    tr->nodeClass  = 4;
    // link the nodes
    r->left = p;
    r->right = q;
    p->parent = r;
    q->parent = r;
    
    tr->time = currentTime;
    tr->scaledByThetaTimeInputTreeUnits =currentTime * theta;
    if (population.order == numClones-1)//oldest population
    {
        
        tr->timePUnits =currentTime * theta;
    }
    else{
        
        tr->timePUnits =currentTime * theta * population.x;
    }
    
    tr->timeInputTreeUnits = tr->timePUnits;
    
    p->length = tr->timePUnits -tp->timePUnits;
    q->length= tr->timePUnits -tq->timePUnits;
    
    if (noisy > 1)
        fprintf (stderr, "\t|\tCurrentTime (input units) = %Lf", tr->timePUnits);
    /* readjust active nodes */
    
    population.idsActiveGametes[firstInd] = r->node_index;// the r3 node
    population.idsActiveGametes[secondInd] = population.idsActiveGametes[population.numActiveGametes - 1];
    
    numActiveGametes = numActiveGametes - 1; /* less 1 active node */
    
    //update list ids nodes
    population.idsGametes[population.numGametes] = newInd;
    population.numGametes = population.numGametes +1;
    
    
    nextAvailable=nextAvailable+1; /* 1 node more is available */
    
    assert(tr->time>0);
    
    population.CoalescentEventTimes[ population.numCompletedCoalescences]=  tr->time;
    population.numCompletedCoalescences= population.numCompletedCoalescences+1;
    
    population.numActiveGametes = population.numActiveGametes - 1; /* now this clone
                                                                    has 1 less node */
    
    
    /* memory for number of nodes */
    if (nextAvailable >= numNodes)  /* if there aren't enough nodes it go into and it addition more */
    {
        /* ReallocNodes(&numNodes, activeGametes); */
        if (noisy == 4)
            fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
        numNodes += INCREMENT_NODES;
        
    }
    
}
void StructuredCoalescentTree::SimulatePopulation(Population &popI,
                                                  ProgramOptions &programOptions,
                                                  const gsl_rng *rngGsl,
                                                  boost::mt19937 *rngBoost,
                                                  int &numNodes,
                                                  int numClones,
                                                  int &nextAvailable,
                                                  int &numActiveGametes,
                                                  int &labelNodes,
                                                  long double  &currentTime,
                                                  int &eventNum,
                                                  long double K)
{
    int   i,  k, isCoalescence,
    firstInd,  newInd,
    isMigration, whichClone;
    long double      eventTime;
    pll_rnode_t  *p, *r, *r1 ;
    TreeNode  *u,*v, *u1;
    
    eventNum = 0;
    
    long double  ThisRateCA = 0.0;
    long double  ThisTimeCA_W = 0.0;
    long double  ThisTimeCA_V1 = 0.0;
    long double  ThisTimeCA_V2 = 0.0;
    int numParcialActiveGametes;
    
    
    numParcialActiveGametes = popI.numActiveGametes;
    
    int numMigrations = popI.numIncomingMigrations; //taking into account also the    time of origin as a migration
    
    long double  timeNextMigration;
    int indexNextMigration = 0;
    Population *incomingPop;
    
    if (programOptions.noisy > 1)
        fprintf (stderr, "\n\n>> Simulating evolutionary history of clone %d (number active gametes %d, original time to origin %Lf)\n", popI.index, popI.numActiveGametes, popI.timeOriginInput);
    
    if (programOptions.noisy > 1)
        fprintf (stderr, "\n\n> Simulating evolutionary history of clone  or order  %d ..\n", popI.order);
    
    currentTime=0;
    popI.numCompletedCoalescences =0;
    while (indexNextMigration < numMigrations) {
        timeNextMigration = popI.immigrantsPopOrderedByModelTime[indexNextMigration].first;
        if ( popI.numActiveGametes >= 2) {
            ThisRateCA = (long double)  popI.numActiveGametes * ((long double) popI.numActiveGametes - 1) / 2.0;
            ThisTimeCA_W = Random::RandomExponentialStartingFrom (ThisRateCA,0,  true, rngGsl,NULL ) ;
            ThisTimeCA_V1 = Population::FmodelTstandard (currentTime, popI.timeOriginSTD, popI.delta, K);
            ThisTimeCA_V1 = ThisTimeCA_V1 + ThisTimeCA_W;
            // from standard time to model time, GstandardTmodel(V, T, delta)
            ThisTimeCA_V2 = Population::GstandardTmodel(ThisTimeCA_V1, popI.timeOriginSTD, popI.delta, K);
        }
        else
        {
            ThisTimeCA_V2 = timeNextMigration + 1.0; // it reached a "provisional" MRCA
        }
        if ( ThisTimeCA_V2 < timeNextMigration)
        {
            //choose randomly two lineages to coalesce
            isCoalescence = YES;
            isMigration = NO;
            eventNum= eventNum +1;
            whichClone = popI.index;
            currentTime = ThisTimeCA_V2; // update current time in model time
            eventTime = currentTime;
            
            if (programOptions.noisy > 1)
            {
                fprintf (stderr, "\n\n*** Event %3d *** currentTime (model time) = %Lf, currentTime (standard time) = %Lf\n", eventNum, ThisTimeCA_V2, ThisTimeCA_V1 );
                fprintf (stderr, "\n\n*** Event %3d *** currentTime (input units) = %Lf\n", eventNum, ThisTimeCA_V2);
            }
            if (programOptions.noisy == 4)
                fprintf (stderr, "* Coalescence *\n");
            
            assert(currentTime>0);
            MakeCoalescenceEvent( popI,  rngGsl, programOptions.noisy, numActiveGametes, nextAvailable, labelNodes, currentTime,  programOptions.numNodes);
            
        }
        else
        {
            if (indexNextMigration < numMigrations - 1) //indexNextMigration corresponds to one of the true migrations
            {
                isCoalescence = NO;
                isMigration = YES;
                eventNum= eventNum +1;
                currentTime = timeNextMigration;
                eventTime = currentTime;
                //update migration times in model time
                if (programOptions.noisy > 1)
                {
                    fprintf (stderr, "\n\n*** Event %3d *** *currentTime (model units) = %Lf\n", eventNum, ThisTimeCA_V2);
                }
                if (programOptions.noisy == 4)
                {fprintf (stderr, "* Migration *\n");}
                
                incomingPop = popI.immigrantsPopOrderedByModelTime[indexNextMigration].second;
                
                p = rnodes[incomingPop->nodeIdAncestorMRCA]; // root of younger clone
                
                indexNextMigration = indexNextMigration + 1;
                
                v= (TreeNode *)(p->data);
                
                printf( "\n The incoming population %d to  population %d with node %d and nodelet %d time %Lf", incomingPop->order, popI.order, v->index, p->node_index, v->timePUnits );
                
                v->indexCurrentClone = popI.index;
                v->indexOldClone = incomingPop->index;
                v->orderCurrentClone = popI.order;
                
                k = v->indexCurrentClone;
                incomingPop->numActiveGametes = incomingPop->numActiveGametes - 1; /* now the other clone has 1 less node */
                // remove node from old clone in list of active gametes and add the new node of the current clone
                //popI->idsActiveGametes[popI->numActiveGametes]=r->index;//adding the superfluos node
                popI.idsActiveGametes[popI.numActiveGametes]=p->node_index;//adding the superfluos node
                
                popI.numActiveGametes = popI.numActiveGametes + 1; /* now this clone has 1 more node */
                //                if (noisy > 1)
                
                if (programOptions.noisy > 1)
                    fprintf (stderr, "\t|\tCurrentTime (input units) = %Lf", v->timePUnits);
                /* memory for number of nodes */
                if (nextAvailable >= numNodes)  /* if there aren't enough nodes it go into and it addition more */
                {
                    /* ReallocNodes(&numNodes, activeGametes); */
                    if (programOptions.noisy == 4)
                        fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
                    numNodes += INCREMENT_NODES;
                }
            }
            else {
                //origin reached
                currentTime = timeNextMigration;
                indexNextMigration = indexNextMigration + 1;
                isCoalescence = NO;
                isMigration = NO;
                eventTime = currentTime;
                if (programOptions.noisy > 1)
                    fprintf (stderr, "\n\n*** Clone origin ***\n");
                if (programOptions.noisy == 4)
                    fprintf (stderr, "Clone origin %d at time (model units) = %Lf\n", popI.index, currentTime);
                if (popI.order < numClones - 1) // do not do it for the last clone
                {
                    newInd = nextAvailable;
                    
                    r1 = rnodes[newInd];   /* new nodelet ancester */
                    
                    r1->node_index = newInd;
                    
                    if (r1->data ==NULL)
                        r1->data =  new TreeNode(0);
                    
                    u1= (TreeNode *)(r1->data);
                    
                    u1->index = nextAvailable;
                    u1->label = labelNodes;
                    labelNodes=labelNodes+1;
                    u1->indexOldClone =u1->indexCurrentClone = popI.index;
                    
                    u1->orderCurrentClone = popI.order;
                    //u1->effectPopSize= popI->effectPopSize;
                    popI.nodeIdAncestorMRCA=newInd;
                    
                    u1->nodeClass  =4;
                    firstInd = nextAvailable - 1;
                    //p = nodes + activeGametes[firstInd]; // descendant node (previously generated node, nextAvailable - 1)
                    p = rnodes[firstInd]; // descendant node (previously generated node, nextAvailable - 1)
                    // link the nodes
                    //r->left = p;
                    r1->left = p;
                    r1->right = NULL;
                    p->parent= r1;
                    
                    u1->time = currentTime;
                    u1->scaledByThetaTimeInputTreeUnits = currentTime * theta;
                    if (popI.order == numClones-1)//oldest population
                    {
                        u1->timePUnits =  currentTime * theta ;
                        p->length = currentTime * theta;
                    }
                    else{
                        
                        u1->timePUnits =  currentTime * theta * popI.x;
                        p->length = currentTime * theta * popI.x;
                    }
                    
                    u1->timeInputTreeUnits = u1->timePUnits;
                    popI.rMRCA = p;
                    /* readjust active nodes */
                    nextAvailable=nextAvailable+1; /* 1 node more is available */
                    popI.idsActiveGametes[0] = newInd;//always will be in the 0  position because there is only one left
                    
                    if (programOptions.noisy > 1)
                        fprintf (stderr, "Creating origin node, it creates node %d derived from node %d", newInd, firstInd);
                    if (programOptions.noisy > 1)
                        fprintf (stderr, "\t|\tCurrentTime (input units) = %Lf", u1->timePUnits);
                    /* memory for number of nodes */
                    if (nextAvailable >= numNodes)  /* if there aren't enough nodes it go into and it addition more */
                    {
                        /* ReallocNodes(&numNodes); */
                        if (programOptions.noisy == 4)
                            fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
                        numNodes += INCREMENT_NODES;
                    }
                }
                else  {//origin of oldest pop reached
                    popI.nodeIdAncestorMRCA=nextAvailable-1;//for the last population, nodeIdAncesterMRCA is the MRCA instead of ancester of MRCA
                    r = rnodes[nextAvailable-1];//popI->idsActiveGametes[0]
                    u= (TreeNode *)(r->data);
                    u->indexOldClone = u->indexCurrentClone = popI.index;
                    u->orderCurrentClone = popI.order;
                    popI.rMRCA= r;
                }
            }
        }
        if (programOptions.noisy > 3)
        {
            fprintf (stderr, "\nActive nodes (%d):",  popI.numActiveGametes);
            for (i = 0; i < popI.numActiveGametes; i++)
                fprintf (stderr, " %d", popI.idsActiveGametes[i]);
            fprintf (stderr, "\t|\tNext node available = %d", nextAvailable);
        }
    }
    if (programOptions.noisy > 1)
        fprintf (stderr, "\n\nEvolutionary history of clone %d is completed \n", popI.index);
    
}
PopulationSet& StructuredCoalescentTree::getPopulationSet()
{
    return *populationSet;
}
pll_rtree_t* StructuredCoalescentTree::getTree()
{
    
    return rtree;
}

pll_rnode_t* StructuredCoalescentTree::BuildTree(Population *CurrentPop,
                                                 const gsl_rng *randomGenerator,
                                                 ProgramOptions &programOptions,
                                                 pll_rnode_t *tumour_mrca,
                                                 int &nextAvailable,
                                                 int &newInd,
                                                 long double   &currentTime,
                                                 int &labelNodes)
{
    int i, j;
    int indexCurrentTip;
    int  foundSuperflousNode;
    pll_rnode_t *p, *q, *r;
    TreeNode *u, *anc, *healthyR, *u1, *u2;
    pll_rnode_t *treeRootInit;
    /********** BUILDING TREES ***********/
    if (programOptions.noisy > 1)
    {
        fprintf (stderr, "\n>> Building trees ..");
    }
    
    j=0;
    indexCurrentTip=0;
    /* get rid of superflous nodes */
    foundSuperflousNode = YES;
    while (foundSuperflousNode == YES)
    {
        foundSuperflousNode = NO;
        for (i = 0; i < nextAvailable; i++) // available all the nodes
        {
            p = rnodes[i];
            
            //fprintf (stderr, "\n\np->index = %d", p->index);
            
            if (p->left == NULL && p->right == NULL && p->parent == NULL)
            {
                // nothing to do with this node because it is not connected to anything
            }
            else if (p->left == NULL && p->right == NULL && p->parent != NULL)
            {
                // (*treeTips)[indexCurrentTip]=p;
                if(indexCurrentTip <  programOptions.TotalNumSequences)
                {
                    //treeTips[indexCurrentTip]=p;
                    if(std::find(rtreeTips.begin(), rtreeTips.end(), p) == rtreeTips.end())
                    {
                        rtreeTips.push_back(p);
                        
                    }
                    indexCurrentTip++;
                }
                
                // do not do anything with this node because it is a tip
            }
            else if (p->left != NULL && p->right == NULL && p->parent != NULL)
            {
                // this is a superflous node and can be removed(this superfluos nodes are the MRCA nodes of the demes
                foundSuperflousNode = YES;
                q = p->left;
                r = p->parent;
                
                if (p->parent->left == p)
                {
                    
                    r->left = q;
                    q->parent = r;
                    p->left = NULL;
                    p->parent = NULL;
                    
                }
                else
                {
                    r->right = q;
                    q->parent = r;
                    p->left = NULL;
                    p->parent = NULL;
                    
                }
                
                //fprintf (stderr, "\n - this is a superflous node and can be removed (1)");
            }
            else if (p->left == NULL && p->right != NULL && p->parent != NULL)
            {
                // this is a superflous node and can be removed
                foundSuperflousNode = YES;
                q = p->right;
                r = p->parent;
                if (p->parent->left == p)
                {
                    r->left = q;
                    q->parent = r;
                    p->right = NULL;
                    p->parent = NULL;
                    
                }
                else
                {
                    r->right = q;
                    q->parent = r;
                    p->right = NULL;
                    p->parent = NULL;
                    
                }
                
                //fprintf (stderr, "\n - this is a superflous node and can be removed (2)");
            }
            else if (p->left != NULL && p->right != NULL && p->parent != NULL)
            {
                
                // this is an internal node formed by a coalescence event, do not touch
                
            }
            else if (p->left != NULL && p->right != NULL && p->parent == NULL)
            {
                
                // this is the  MRCA
                
            }
            
            else if (p->left != NULL && p->right == NULL && p->parent == NULL)
            {
                // Seems to be the last coalescent event among sequences with non-ancestral material
                // it is not superfluous, we just remove it
                p->left->parent = NULL;
                
            }
            else if (p->left == NULL && p->right != NULL && p->parent == NULL)
            {
                // not clear what this node could be doing, but we will remove it anyway
                fprintf (stderr, "strange\n");
                p->left=NULL;
                p->right->parent =NULL;
                
            }
            else
            {
                fprintf (stderr, "You should not be here, I think\n");
                //fprintf (stderr, "%d %d-- %d %d %d\n", Index(p), j, Index(p->left), Index(p->right), Index(p->anc1));
            }
            if (p->parent != NULL)
            {//update length field
                
                // p->length = p->parent->time- p->time;
                //p->length = (p->anc1->timePUnits- p->timePUnits);
                u = (TreeNode *)(p->data);
                anc = (TreeNode *)(p->parent->data);
                u->length =(anc->timePUnits- u->timePUnits);
                assert(p->length-(anc->timePUnits- u->timePUnits)<=0.0000000001);
                p->length = (anc->timePUnits- u->timePUnits);
                //*mutationRate;
                u->lengthModelUnits = (anc->time- u->time);
                addEdgeFromNode(p );
                //*mutationRate;
                
            }
        }
        //fprintf (stderr, "\n");
    }//while
    
    /* about the MRCA */
    newInd=nextAvailable-1;
    p = rnodes[newInd]; /* because the last one event is the last coalescence */
    //fprintf (stderr, "\n\n\n>> newInd = %d\n", newInd);
    
    if (programOptions.thereisOutgroup == NO)
    {
        p = rnodes[newInd];
        u = (TreeNode *)(p->data);
        u->nodeClass = 5;
        tumour_mrca = p;
        //CurrentPop->rMRCA=p;
        //treeRootInit[0] = p;
        p->parent=NULL;
        u->anc1 = NULL;
        treeRootInit =tumour_mrca;
    }
    if (programOptions.thereisOutgroup == YES && programOptions.outgroupSelection > 0)  /*** Root and outgroup ***/
    {
        p = rnodes[newInd]; // MRCA
        u = (TreeNode *)(p->data);
        u->nodeClass = 4;
        
        if (programOptions.noisy > 1)
            std::cout << "\n\n>> Attaching outgroup .. "<< std::endl;
        
        
        if (programOptions.outgroupSelection == 1)  /*** Root 2 times and outgroup ***/
            currentTime = CurrentPop->timeOriginSTD;
        // origin of the clone; // currentTime + (outgroupBranchLength_Root1Root2 / mutationRate); // set time of the new root (from which the MRCA and outgroup nodes are derived)
        else if (programOptions.outgroupSelection == 2) { /*** Root 2 times and outgroup ***/
            currentTime = CurrentPop->timeOriginSTD + (programOptions.outgroupBranchLength_Root1Root2 / (theta * CurrentPop->x)); // origin of the clone + time given by the user
            
        }
        else
        {
            std::cout<< "\n\nError simulationg the outgroup. Check input settings\n"<<std::endl;
            
        }
        
        pll_rnode_t*       healthyRoot = rnodes[nextAvailable];
        
        if (healthyRoot->data ==NULL)
            healthyRoot->data= new TreeNode(0);
        u = (TreeNode *)(p->data);
        healthyR = (TreeNode *)(healthyRoot->data);
        healthyR->index = nextAvailable;
        healthyR->label = labelNodes;
        healthyR->indexOldClone = healthyR->indexCurrentClone = CurrentPop->index;
        healthyR->orderCurrentClone = CurrentPop->order;
        //healthyR->effectPopSize= u->effectPopSize;
        labelNodes= labelNodes+1;
        healthyRoot->left = p;//coalTreeMRCA;
        //        coalTreeMRCA->anc = healthyRoot;
        //p->anc1 = healthyRoot;
        p->parent = healthyRoot;
        
        
        healthyR->nodeClass = 5;
        //        coalTreeMRCA->length = transformingBranchLength/mutationRate;
        p->length = 0;
        u->length = 0;
        //        coalTreeMRCA->branchLength = transformingBranchLength;
        u->lengthModelUnits = 0;
        //        healthyRoot->time = currentTime +  transformingBranchLength/mutationRate;
        healthyR->time = currentTime  ;
        healthyR->scaledByThetaTimeInputTreeUnits = currentTime * theta;
        
        healthyR->timePUnits = currentTime * theta;
        //healthyR->timePUnits = currentTime * theta *CurrentPop->x;
        healthyR->timeInputTreeUnits = healthyR->timePUnits;
        //p->length = (p->anc1->timePUnits- p->timePUnits);
        anc = (TreeNode *)(p->parent->data);
        p->length = (anc->timePUnits- u->timePUnits);
        u->length = (anc->timePUnits- u->timePUnits);
        //*mutationRate;
        u->lengthModelUnits = (anc->time- u->time);
        addEdgeFromNode(p );
        //*mutationRate;
        
        healthyRoot->length = 0;
        healthyR->length = 0;
        nextAvailable++;
        
        pll_rnode_t* healthyTip1= rnodes[nextAvailable];
        if (healthyTip1->data==NULL)
            healthyTip1->data =  new TreeNode(0);
        
        
        u1= (TreeNode *)(healthyTip1->data);
        if (healthyTip!=NULL)
        {
            healthyTip1->label = healthyTip->label;
            if (healthyTip->data!=NULL)
            {
                u2= (TreeNode *)(healthyTip->data);
                u1->genotypeSequence = u2->genotypeSequence;
            }
        }
        else
        {
            healthyTip1->label = new char[((std::string)(programOptions.healthyTipLabel)).length() + 1];
            healthyTip1->label= strcpy(healthyTip1->label, programOptions.healthyTipLabel.c_str());
        }
        
        healthyRoot->right = healthyTip1;
        healthyTip1->parent=healthyRoot;
        u1->time  =0;
        u1->timePUnits =0;
        double  healthyTipBranchLengthRatio = Random::randomUniformFromGsl2(randomGenerator);
        
        u1->time = healthyTipBranchLengthRatio * healthyR->time;
        u1->timePUnits = healthyTipBranchLengthRatio * healthyR->timePUnits;
        u1->timeInputTreeUnits = u1->timePUnits;
        
        healthyTip1->length= healthyR->timePUnits- u1->timePUnits;
        
        u1->length = 0;
        u1->lengthModelUnits = 0;
        
        u1->isOutgroup= YES;
        
        treeRootInit=healthyRoot;
        
    }
    
    int intLabel = 0;
    if (programOptions.noisy > 1)
        std::cout << "\n\n>> Relabeling nodes on tree... \n\n" << std::endl;
    if (programOptions.thereisOutgroup == YES)
        intLabel = programOptions.TotalTumorSequences + 1;
    else
        intLabel = programOptions.TotalTumorSequences;
    
    RelabelNodes2( treeRootInit, intLabel );
    
    return treeRootInit;
}
void StructuredCoalescentTree::addEdgeFromNode(pll_rnode_t *node ){
    pll_tree_edge_t *edge;
    TreeNode *u, *v;
    edge = new pll_tree_edge_t();
    u= (TreeNode *)(node->data);
    if (node->parent!=NULL)
    {
        v= (TreeNode *)(node->parent->data);
        edge->edge.rtree.parent = u->timeInputTreeUnits >= v->timeInputTreeUnits ? node :  node->parent;
        edge->edge.rtree.child= u->timeInputTreeUnits < v->timeInputTreeUnits ? node :  node->parent;
        assert(abs(edge->edge.rtree.child->length -(v->timeInputTreeUnits- u->timeInputTreeUnits))<=0.000000001 );
        edge->length= edge->edge.rtree.child->length;
        //edgeLengths.push_back(make_pair(edge->length, edge));
        //edges.push_back(edge);
        edgeLengths[currentNumberEdgesRootedTree]=std::make_pair(edge->length, edge);
        edges[currentNumberEdgesRootedTree]=edge;
        currentNumberEdgesRootedTree++;
    }
    
}

void StructuredCoalescentTree::RelabelNodes2(pll_rnode_t *p, int &intLabel)
{
    if (p != NULL)
    {
        TreeNode *u;
        RelabelNodes2 (p->left, intLabel);
        RelabelNodes2 (p->right, intLabel);
        /*RelabelNodes (p->outgroup);*/
        if (p->left == NULL && p->right == NULL) /* is tip */
        {
            u=(TreeNode*) p->data;
            u->label = u->index ;
        }
        else                  /* all ancester */
        {
            u=(TreeNode*) p->data;
            u->label = intLabel;
            intLabel=intLabel+1;
        }
    }
}
std::vector <long double> StructuredCoalescentTree::getDeltaTs()
{
    return populationSet->getDeltaTs();
}
std::vector <long double> StructuredCoalescentTree::getDeltas(){
    
    return populationSet->getDeltas();
}
std::vector <long double> StructuredCoalescentTree::getTs(){
    
    return populationSet->getTs();
}
std::vector <long double> StructuredCoalescentTree::getSampleSizes(){
    
    return populationSet->getSampleSizes();
}
pll_rnode_t* StructuredCoalescentTree::getRoot(){
    
    return root;
}
