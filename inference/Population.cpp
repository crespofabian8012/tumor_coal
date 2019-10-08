//
//  Population.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 4/10/19.
//

#include "Population.hpp"

using namespace std;

//Parametrized Constructor
population::population(int ind, int ord,   double  timeOriginInputp,
               int           sampleSizep, int popSizep, double birthRatep,
               double deathRatep)
    {
        if (ind >=0)
           index=ind;
        else
            fprintf (stderr, "\n ERROR: Population  index cannot be negative \n");
        if (ord >= 0)
          order=ord;
        else
            fprintf (stderr, "\n ERROR: Population order cannot be negative \n");
        
        if (popSizep >= 0)
        {
            popSize=popSizep;
            oldpopSize = popSizep;
        }
        else
            fprintf (stderr, "\n ERROR: Population size cannot be negative \n");
        
        if (sampleSizep >= 0)
        {
           sampleSize=sampleSizep;
           oldsampleSize=sampleSizep;
        }
        else
            fprintf (stderr, "\n ERROR: Population sample size cannot be negative \n");
      
        
        if (birthRatep >0){
            birthRate =birthRatep;
            oldbirthRate = birthRatep;
        }
        else
            fprintf (stderr, "\n ERROR: Population birth rate cannot be negative \n");
      
        if (deathRatep >0){
            deathRate= deathRatep;
            olddeathRate =deathRatep;
        }
        else
            fprintf (stderr, "\n ERROR: Population death rate cannot be negative \n");
     
        if (birthRatep > deathRatep)
        {
            growthRate = birthRate - deathRate;
            oldgrowthRate = growthRate;
        }
        else
            fprintf (stderr, "\n ERROR: Population growth rate cannot be negative \n");
        
        if (timeOriginInputp > 0)
        {
            timeOriginInput = timeOriginInputp;
            oldtimeOriginInput= timeOriginInputp;
        }
        else
            fprintf (stderr, "\n ERROR: Population time of origin   cannot be negative \n");
        
        effectPopSize = popSize / birthRate;
        oldeffectPopSize=effectPopSize;
       
        
        delta= growthRate * effectPopSize;
        olddelta=delta;
        
        isAlive=YES;
        
        timeOriginSTD = timeOriginInput /effectPopSize;
        
        numActiveGametes = sampleSize;
        oldsampleSize = sampleSize;
        
        numGametes=sampleSize;
        
        numCompletedCoalescences=0;
        nextAvailableIdInmigrant =0;
        numIncomingMigrations = 0;
        numPossibleMigrations = 0;
        
        migrationTimes =NULL;
        MRCA=NULL;
        FatherPop=NULL;
        oldFatherPop=NULL;
        CoalescentEventTimes=NULL;
        oldCoalescentEventTimes=NULL;
        immigrantsPopOrderedModelTime=NULL;
        
        idsActiveGametes=NULL;
        idsGametes=NULL;
    
}
double population::ProbabilityComeFromPopulation ( population* PopJ, population **populations, int numClones)
{
    double  ProbabilityIJ, AboveTerm, BelowTerm;
    int     l, j;
    double  h, a, b, c, d, e, t, cum;
    population *p;
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
    //fprintf (stderr, "AboveTerm = %lf\n", AboveTerm);
    j=0;
    for (l = order + 1; l < numClones; l++)
    {    p = *(populations + l);
        //fprintf (stderr, "\ni = %d, j = %d, l = %d\n", i, j, l);
        t = (timeOriginSTD ) * (effectPopSize) / ( p->effectPopSize);
        h = CalculateH(t, p->timeOriginSTD, p->delta);
        
        cum = cum + ( ( p->popSize) * h);
    }
    
    BelowTerm = cum;
    //fprintf (stderr, "BelowTerm = %lf\n", BelowTerm);
    
    ProbabilityIJ = AboveTerm / BelowTerm;
    //fprintf (stderr, "ProbabilityIJ = %lf\n", ProbabilityIJ);
    
    return ProbabilityIJ;
}

double population::CalculateH (double t, double TOrigin, double delta)
{
    double  H, AboveTerm, BelowTerm, firstTerm, secondTerm;
    double  a, b;
    
    AboveTerm = 0.0;
    BelowTerm = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    a = 0.0;
    b = 0.0;
    H = 0.0;
    
    //printf ("\nInput(H), t=%lf T=%lf delta=%lf ", t, T, delta);
    
    a = 1.0 - exp(-1.0 * delta * (TOrigin - t));
    firstTerm = a * a;
    secondTerm = exp(-1.0 * delta * t);
    AboveTerm = firstTerm * secondTerm;
    //printf ("\nAboveTerm(H) = %lf (%lf %lf %lf) / delta = %lf, T = %lf, t = %lf", AboveTerm, a, firstTerm, secondTerm, delta, T, t);
    
    b = 1.0 - exp(-1.0 * delta * TOrigin);
    BelowTerm = b * b;
    //printf ("\nBelowTerm(H) = %lf", BelowTerm);
    
    H = AboveTerm / BelowTerm;
    //printf ("\nH = %lf", H);
    
    return H;
}

double population::FmodelTstandard (double t, double TOrigin, double delta)
{
    double  ModelTimeF, firstTerm, secondTerm, thirdTerm;
    double  a, b, c;
    
    ModelTimeF = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    thirdTerm = 0.0;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    
    /*fprintf (stderr, "t = %lf\n", t);
     fprintf (stderr, "T = %lf\n", T);
     fprintf (stderr, "delta = %lf\n", delta);*/
    
    a = 1.0 - exp(-1.0 * delta * TOrigin);
    firstTerm = a * a;
    //fprintf (stderr, "firstTerm = %lf\n", firstTerm);
    
    secondTerm = exp(delta * TOrigin);
    //fprintf (stderr, "secondTerm = %lf\n", secondTerm);
    
    b = 1.0 / (1.0 - exp(-1.0 * delta * (TOrigin - t)));
    c = 1.0 / (1.0 - exp(-1.0 * delta * TOrigin));
    thirdTerm = b - c;
    //fprintf (stderr, "thirdTerm = %lf\n", thirdTerm);
    
    
    ModelTimeF = firstTerm * secondTerm * thirdTerm;
    //fprintf (stderr, "ModelTimeF = %lf\n", ModelTimeF);
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
    
    ModelTimeF = a * b / (delta * c);
    //fprintf (stderr, "ModelTimeF = %lf\n", ModelTimeF);
    
    return ModelTimeF;
}
int   population::InitListPossibleMigrations( ){
    
    //qsort(populations, numClones, sizeof(Population*), comparePopulationsByTimeOrigin);
    //after qsort
    population *p;
    int i, j;
    double d;
   
        numPossibleMigrations = order + 1;
        numIncomingMigrations = 1; //the time of origin counts as one migration
        // if (!(p->migrationTimes)){
        migrationTimes =(double *)  malloc( (order + 1) * sizeof( double));
        if (!(migrationTimes)){
            
            fprintf (stderr, "Could not allocate migrationTimes (%lu bytes)\n", (order + 1) * (long) sizeof(double));
            exit (1);
            
        }
        // }
        
        d = (double)(timeOriginSTD);
        migrationTimes[0] = d;
        
        
        if(order >0 )//&& !(p->immigrantsPopOrderedModelTime) )
            immigrantsPopOrderedModelTime = (population**) malloc(order * sizeof(population*)); //besides the other possible immigrants we need to add this populations itself
    
        
        for (j = 1; j < (order +1); j++)
        {

            immigrantsPopOrderedModelTime[j-1]=NULL;
    
            d = 2 * ((double)(timeOriginSTD)); //  a value greater than time of origin  standarized by  the population
            migrationTimes[j] = d;
        }
    return 0;
}
int  population::resetMigrationsList(){
    
    population *p;
    int i, j;
    struct Population** pops ;
    double* migrationTimes ;
    double d;
    
    
        numIncomingMigrations = 1; //the time of origin counts as one migration
        d = (double)(timeOriginSTD);
        if (migrationTimes !=NULL)
           migrationTimes[0] = d;
        else
            return 1;
        for (j = 1; j < (order +1); j++)
        {
        
             immigrantsPopOrderedModelTime[j-1]=NULL;

            d = 2 * ((double)(timeOriginSTD)); //  a value greater than time of origin  standarized by  the population
            migrationTimes[j] = d;
       
        }
    return 0;
}
void  population::UpdateListMigrants( int numClones, population *PopChild, population *PopFather  )
{
    if (PopChild->index ==  PopFather->index) {
        fprintf (stderr, "\nError. The target population %d for  migration must be different than the population of origin %d \n", PopFather->index, PopChild->index);
        exit (-1);
    }
    if (PopFather->order <= PopChild->order ) {
        fprintf (stderr, "\nError. The target population %d for  migration must be older than the population of origin %d \n", PopFather->index, PopChild->index);
        exit (-1);
    }
    population *pOrigin, *pTarget;
    double *ptr;
    int  j;
    int lengthMigrationsArray = (int)(PopFather->order) + 1;
    int updatedNumIncomingMigrations = PopFather->numIncomingMigrations;
    // printf ( "\n lengthMigrationsArray= %d \n", lengthMigrationsArray );
    
    population **pops  = PopFather->immigrantsPopOrderedModelTime;
    ptr = PopFather->migrationTimes;
    for (j = 0; j < lengthMigrationsArray; j++) {
        if (ptr[j] >  (PopFather->timeOriginSTD) ) {
            break; //first null pointer
        }
    }
    double updatedMigrationTime = (PopChild->timeOriginSTD) * (PopChild->effectPopSize) / (PopFather->effectPopSize);
    PopFather-> migrationTimes[j] = updatedMigrationTime;
    if(updatedNumIncomingMigrations + 1 <= PopFather->numPossibleMigrations){
        updatedNumIncomingMigrations = updatedNumIncomingMigrations + 1;
        PopFather->numIncomingMigrations = updatedNumIncomingMigrations;
        PopFather->immigrantsPopOrderedModelTime[j-1] = PopChild;
    }
    //  fprintf (stderr ,"\n updatedNumIncomingMigrations %d \n",PopFather->numIncomingMigrations);
    //PopFather->immigrantsPopOrderedModelTime[j-1] = PopChild;
    //order immigrant population by time of origin
    if(PopFather->numIncomingMigrations > 1 )
        qsort(PopFather->migrationTimes, PopFather->numIncomingMigrations, sizeof(double), compare);
    if(PopFather->numIncomingMigrations -1 > 1 )
        qsort(PopFather->immigrantsPopOrderedModelTime, PopFather->numIncomingMigrations -1,  sizeof(population*), comparePopulationsByTimeOrigin);
    
}

int population::comparePopulationsByTimeOrigin(const void *s1, const void *s2)
{
    population *p1 = *(population **)s1;
    population *p2 = *(population **)s2;
    if (  p1->timeOriginInput  > p2 ->timeOriginInput)
        return 1;
    else if (p1->timeOriginInput < p2->timeOriginInput)
        return -1;
    else
        return 0;
}
int population::compare (const void * a, const void * b)

{
    
    double *p1 = (double *)a;
    
    double *p2 = (double *)b;
    
    
    
    if (*p1 > *p2) return 1;
    
    else if (*p2 > *p1 )return -1;
    
    else return 0;
    
}
void population::InitCoalescentEvents(int numClones){
    
    CoalescentEventTimes=(double *) calloc((sampleSize -1)+ numClones -1, (long) sizeof(double));
    if (!(CoalescentEventTimes))
    {
        fprintf (stderr, "CoalescentEventTimes (%lu bytes)\n", (sampleSize -1 +numClones -1) * (long) sizeof(double));
        exit (-1);
    }
    
    
    
    
}
