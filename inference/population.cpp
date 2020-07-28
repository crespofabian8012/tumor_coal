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
#include "random.h"
#include <algorithm>
#include <vector>
#include <iterator>

using namespace std;

//#include "data_utils.hpp"


//Parametrized Constructor
Population::Population(int ind, int ord, long double timeOriginInput,
                       int sampleSize, int popSize, long double birthRate,
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
    
    if (popSize >= 0)
    {
        this->popSize=popSize;
    }
    else
        fprintf (stderr, "\n ERROR: Population size cannot be negative \n");
    
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
    if (birthRate >0)
        effectPopSize = popSize / birthRate;
    
    delta= growthRate * effectPopSize;
    
    isAlive = YES;
    
    if (effectPopSize>0)
      timeOriginSTD = timeOriginInput /effectPopSize;
    else
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
    x = 1.0;
    oldx = 0.0;
    theta = 0.0;
    oldTheta = 0.0;
    oldOrder=0;
    oldScaledTimeOriginInput=0.0;
    olddelta=0.0;
    oldeffectPopSize =0.0;
    oldDeathRate=0.0;
    oldPopSize=0.0;
    oldGrowthRate =0.0;
    oldTimeOriginSTD=0.0;
    oldTimeOriginInput=0.0;
    indexFirstObservedCellName=0;
    scaledtimeOriginInput=0.0;
}
long double Population::ProbabilityComeFromPopulation(Population *PopJ, vector<Population*> &populations, int numClones)
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
    
    t = (timeOriginSTD ) * (effectPopSize) / ( PopJ->effectPopSize);
    h = CalculateH(t, PopJ->timeOriginSTD, PopJ->delta);
    AboveTerm = ( PopJ->popSize) * h;
    j=0;
    for (l = order + 1; l < numClones; l++)
    {
        p = populations[l];
        t = (timeOriginSTD ) * (effectPopSize) / ( p->effectPopSize);
        h = CalculateH(t, p->timeOriginSTD, p->delta);
        
        cum = cum + ( ( p->popSize) * h);
    }
    
    BelowTerm = cum;
    
    ProbabilityIJ = AboveTerm / BelowTerm;

    return ProbabilityIJ;
}

long double Population::CalculateH (long double t, long double TOrigin, long double delta)
{
    long double H, AboveTerm, BelowTerm, firstTerm, secondTerm;
    long double a, b;
    
    AboveTerm = 0.0;
    BelowTerm = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    a = 0.0;
    b = 0.0;
    H = 0.0;

    a = 1.0 - exp(-1.0 * delta * (TOrigin - t));
    firstTerm = a * a;
    secondTerm = exp(-1.0 * delta * t);
    AboveTerm = firstTerm * secondTerm;
    b = 1.0 - exp(-1.0 * delta * TOrigin);
    BelowTerm = b * b;
    if (BelowTerm == 0.0)
        printf ("\n BelowTerm = 0.0 \n");
    
    H = AboveTerm / BelowTerm;
 
    return H;
}
long double Population::LogCalculateH (long double t, long double TOrigin, long double delta)
{
    long double logH, AboveTerm, BelowTerm, firstTerm, secondTerm;
    long double a, b;
    
    AboveTerm = 0.0;
    BelowTerm = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    a = 0.0;
    b = 0.0;
    logH = 0.0;

    a = 1.0 - exp(-1.0 * delta * (TOrigin - t));
    firstTerm = 2.0 * log(a);
    secondTerm = -1.0 * delta * t;
    AboveTerm = firstTerm + secondTerm;

    
    b = 1.0 - exp(-1.0 * delta * TOrigin);
    BelowTerm = 2.0 * log(b);
    
    logH = AboveTerm - BelowTerm;

    return logH;
}

long double Population::FmodelTstandard (long double t, long double TOrigin, long double delta)
{
    long double ModelTimeF, firstTerm, secondTerm, thirdTerm;
    long double a, b, c;
    
    ModelTimeF = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    thirdTerm = 0.0;
    a = 0.0;
    b = 0.0;
    c = 0.0;

    // New formula from Carsten!
    ModelTimeF = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    thirdTerm = 0.0;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    
    
    a = exp(delta * t) - 1.0;
    b = 1.0 - exp(-1.0 * delta * TOrigin);
    c = 1.0 - exp(-1.0 * delta * (TOrigin - t));
    
    if ( c == 0.0)
        fprintf (stderr, "\n  c = 0.0, TOrigin = %Lf, t = %Lf \n", TOrigin, t);
    if ( delta == 0.0)
        fprintf (stderr, "\n delta  = 0.0 \n");
    
    ModelTimeF = a * b / (delta * c);
    //fprintf (stderr, "ModelTimeF = %lf\n", ModelTimeF);
    
    return ModelTimeF;
}

long double Population::GstandardTmodel (long double V, long double TOrigin, long double delta)
{
    long double StandardTimeG, firstTerm, secondTerm, thirdTerm;
    long double a, b, c, d, e;
    
    StandardTimeG = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    thirdTerm = 0.0;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    d = 0.0;
    e = 0.0;
 
    firstTerm = TOrigin;

    secondTerm = 1 / delta;
 
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
        fprintf (stderr, "\nApplying approximation of math formula to avoid log(0)\n");
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
    
    return StandardTimeG;
}

void Population::InitListPossibleMigrations(int order)
{
    this->order = order;
    numPossibleMigrations = order + 1;
    numIncomingMigrations = 1; //the time of origin counts as one migration
    immigrantsPopOrderedByModelTime.push_back(make_pair(timeOriginSTD, this));
    
}
int Population::resetMigrationsList(){
    
    immigrantsPopOrderedByModelTime.clear();
    InitListPossibleMigrations(this->order);
    return 0;
}
void Population::UpdateListMigrants( int numClones, Population *PopChild, Population *PopFather  )
{
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
        PopFather->immigrantsPopOrderedByModelTime.push_back(make_pair(updatedMigrationTime, PopChild));
    }
    sort(PopFather->immigrantsPopOrderedByModelTime.begin(), PopFather->immigrantsPopOrderedByModelTime.end(), comparePopulationsPairByTimeOrigin);
    //     printf("\n pop order  %d choose pop father of order %d \n", PopChild->order, PopFather->order);
    //    for (int i = 0; i < PopFather->immigrantsPopOrderedByModelTime.size(); ++i)
    //        printf("\n ordered migrations: time(father pop units) : %lf, pop order: %d, time of origin %lf \n", PopFather->immigrantsPopOrderedByModelTime[i].first,  PopFather->immigrantsPopOrderedByModelTime[i].second->order , PopFather->immigrantsPopOrderedByModelTime[i].second->timeOriginSTD);
}

bool Population::comparePopulationsPairByTimeOrigin(const pair<long double, Population *> s1, const pair<long double, Population *> s2)
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
    for (j = 0; j < sampleSize + numClones - 2; j++) {
        CoalescentEventTimes.push_back(0);
    }
}
void Population::InitIdsActiveGametes()
{
    int j = 0;
    assert(sampleSize >0);
    for (j = 0; j < sampleSize ; j++) {
        idsActiveGametes.push_back(0);
    }
}

void Population::InitIdsGametes(int numClones)
{
    int j = 0;
   
    assert(sampleSize >0);
    for (j = 0; j < (2* sampleSize + numClones - 2) ; j++) {
        idsGametes.push_back(0);
    }
}
void Population::InitRTips()
{
    int j = 0;
    assert(sampleSize >0);
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
void Population::ChooseRandomIndividual(int *firstInd,   int numClones,   int *secondInd, gsl_rng *randomGenerator, int choosePairIndividuals)
{
    long double random;
    int k, w;
    //long double*cumPopulPart = (long double*) malloc((numActiveGametes + 1)* (long) sizeof(double));
    vector<long double> cumPopulPart(numActiveGametes + 1);
//    if (!cumPopulPart)
//    {
//        fprintf (stderr, "Could not allocate cumPopulPart (%lu bytes)\n", (numActiveGametes + 1) * (long) sizeof(double));
//        exit (-1);
//    }
    cumPopulPart[0] = 0;
    for (k = 1; k <= numActiveGametes; k++)
        cumPopulPart[k] = 0;
    for (k = 1; k <= numActiveGametes; k++)
        cumPopulPart[k] = cumPopulPart[k - 1] + 1.0 / (numActiveGametes);
    
    w = 0;
    
    random = Random::randomUniformFromGsl2(randomGenerator);
    *firstInd = Population::bbinClones(random, &cumPopulPart[0], numActiveGametes)-1;
    w = 0;
    
    if (*firstInd >= numActiveGametes || *firstInd < 0 ) /* checking */
    {
        fprintf (stderr, "\n\nERROR: firstInd out of range!\n");
        exit (-1);
    }
    if (choosePairIndividuals== YES && numActiveGametes > 1) {
        
        do//choose randomly another individual to coalesce
        {
            random = Random::randomUniformFromGsl2(randomGenerator);
            *secondInd = Population::bbinClones(random, &cumPopulPart[0], numActiveGametes)-1;
            
        } while (*firstInd == *secondInd  );
    }
   // free (cumPopulPart);
   // cumPopulPart=NULL;
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
    long double result = log(delta * term2 /(term3 * term3));
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
long double Population::LogDensityTimeSTDFrom(long double u, long double from){
    long double  densityFrom, result;
    densityFrom = DensityTimeSTD(u, from );
    if (densityFrom ==0)
        result=log(0);
    else
        result=log(densityFrom);
    return result;
}
long double Population::LogProbNoCoalescentEventBetweenTimes(long double from, long double to, int numberActiveInd)
{
    int j=numberActiveInd;
    long double result=0.0;
    if (j==0 || j==1)
        return 0;
    result= log(j *(j-1) /2.0) -1.0 * j* (j-1)*(Population::FmodelTstandard(to,timeOriginSTD, delta)-Population::FmodelTstandard(from, timeOriginSTD, delta))/ 2.0;
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


