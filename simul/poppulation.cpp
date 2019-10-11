//
//  Population.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 4/10/19.
//

#include "poppulation.hpp"
#include "definitions.hpp"

using namespace std;

//Parametrized Constructor
Population::Population(int ind, int ord, double timeOriginInput,
                       int sampleSize, int popSize, double birthRate,
                       double deathRate, bool estimateTOR)
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
        
        if (timeOriginInput > 0)
        {
            this->timeOriginInput = timeOriginInput;
        }
        else
            fprintf (stderr, "\n ERROR: Population time of origin   cannot be negative \n");
        
        effectPopSize = popSize / birthRate;
        
        delta= growthRate * effectPopSize;
        
        isAlive=YES;
        
        timeOriginSTD = timeOriginInput /effectPopSize;
        
        numActiveGametes = sampleSize;

        numGametes=sampleSize;

        nodeIdAncestorMRCA = 0;
        numCompletedCoalescences=0;
        nextAvailableIdInmigrant =0;
        numIncomingMigrations = 0;
        numPossibleMigrations = 0;
        doEstimateTimeOrigin = estimateTOR;
        
        MRCA=0;
}
double Population::ProbabilityComeFromPopulation(Population *PopJ, vector<Population*> &populations, int numClones)
{
    double  ProbabilityIJ, AboveTerm, BelowTerm;
    int     l, j;
    double  h, a, b, c, d, e, t, cum;
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
    //fprintf (stderr, "AboveTerm = %lf\n", AboveTerm);
    j=0;
    for (l = order + 1; l < numClones; l++)
    {
        p = populations[l];
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

double Population::CalculateH (double t, double TOrigin, double delta)
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

double Population::FmodelTstandard (double t, double TOrigin, double delta)
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

double Population::GstandardTmodel (double V, double TOrigin, double delta)
{
    double  StandardTimeG, firstTerm, secondTerm, thirdTerm;
    double  a, b, c, d, e;
    
    StandardTimeG = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    thirdTerm = 0.0;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    d = 0.0;
    e = 0.0;
    
    //fprintf (stderr, "\nV = %lf, T = %lf, delta = %lf\n", V, T, delta);
    
    
    firstTerm = TOrigin;
    //fprintf (stderr, "\nfirstTerm = %lf\n", firstTerm);
    
    secondTerm = 1 / delta;
    // fprintf (stderr, "secondTerm = %lf\n", secondTerm);
    
    
    //a = 1.0 - exp(-1.0 * delta * T);
    a =  exp(-1.0 * delta * TOrigin);
    
    //fprintf (stderr, "a = %lf\n", a);
    //b = (a * a) * exp(delta * T);
    b = (1 - a) * (1 - a) * (1.0 / a);
    //fprintf (stderr, "b= %lf\n", b);
    
    //c = 1.0 - exp (-1.0 * delta * T);
    c = 1 - a;
    // fprintf (stderr, "c = %lf\n", c);
    //d = c * exp(delta * T);
    d = (1 - a) * (1.0 / a);
    // fprintf (stderr, "d = %lf\n", d);
    e = V + d;
    //fprintf (stderr, "e = %lf\n", e);
    thirdTerm = log(1 - b / e);
    //fprintf (stderr, "valueOfLog = %lf\n", thirdTerm);
    thirdTerm = log(1 - ((1 - a) * (1 - a) * (1.0 / a)) / (V * delta + (1 - a) * (1.0 / a)));
    //fprintf (stderr, "ArgumentOfLog = %lf\n", 1 - b/e);
    //fprintf (stderr, "valueOfLog = %lf\n", thirdTerm);
    
    
    thirdTerm = log(1 + delta * V - a) - log(1 + (delta * V - 1) * a);
    
    //StandardTimeG = firstTerm + (secondTerm * thirdTerm);
    StandardTimeG = secondTerm * thirdTerm;
    //fprintf (stderr, "StandardTimeG = %lf\n", StandardTimeG);
    
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
        //fprintf (stderr, "\nfirstTerm = %lf\n", firstTerm);
        
        
        d = (V * V * delta * exp(-1.0 * delta * TOrigin)) / (1 + V);
        secondTerm =  d;
        //fprintf (stderr, "secondTerm = %lf\n", secondTerm);
        
        StandardTimeG = firstTerm - secondTerm;
        //fprintf (stderr, "StandardTimeG = %lf\n", StandardTimeG);
    }
    
    return StandardTimeG;
}

void Population::InitListPossibleMigrations(int order)
{
//    int j;
   
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
    if (PopFather->order <= PopChild->order ) {
        fprintf (stderr, "\nError. The target Population %d for  migration must be older than the Population of origin %d \n", PopFather->index, PopChild->index);
        exit (-1);
    }
    int updatedNumIncomingMigrations = PopFather->numIncomingMigrations;
     int lengthMigrationsArray = (int)(PopFather->order) + 1;
    // printf ( "\n lengthMigrationsArray= %d \n", lengthMigrationsArray );

    double updatedMigrationTime = (PopChild->timeOriginSTD) * (PopChild->effectPopSize) / (PopFather->effectPopSize);
    if(updatedNumIncomingMigrations + 1 <= PopFather->numPossibleMigrations){
        updatedNumIncomingMigrations = updatedNumIncomingMigrations + 1;
        PopFather->numIncomingMigrations = updatedNumIncomingMigrations;
        PopFather->immigrantsPopOrderedByModelTime.push_back(make_pair(updatedMigrationTime, PopChild));
    }

      sort(PopFather->immigrantsPopOrderedByModelTime.begin(), PopFather->immigrantsPopOrderedByModelTime.end(), comparePopulationsPairByTimeOrigin);
 
     printf("\n pop order  %d choose pop father of order %d \n", PopChild->order, PopFather->order);
    for (int i = 0; i < PopFather->immigrantsPopOrderedByModelTime.size(); ++i)
        printf("\n ordered migrations: time : %lf, pop order: %d, time of origin %lf \n", PopFather->immigrantsPopOrderedByModelTime[i].first,  PopFather->immigrantsPopOrderedByModelTime[i].second->order , PopFather->immigrantsPopOrderedByModelTime[i].second->timeOriginSTD);
   

    //  fprintf (stderr ,"\n updatedNumIncomingMigrations %d \n",PopFather->numIncomingMigrations);
    //PopFather->immigrantsPopOrderedModelTime[j-1] = PopChild;
    //order immigrant Population by time of origin
//    if (PopFather->numIncomingMigrations > 1 )
//        sort(PopFather->migrationTimes.begin(), PopFather->migrationTimes.end(), compare);
//    if (PopFather->numIncomingMigrations -1 > 1 )
//        sort(PopFather->immigrantsPopOrderedModelTime.begin(), PopFather->immigrantsPopOrderedModelTime.end(), comparePopulationsByTimeOrigin);
}

bool Population::comparePopulationsPairByTimeOrigin(const pair<double, Population *> s1, const pair<double, Population *> s2)
{
    //printf("\n s1: %lf ", s1.first);
    // printf("\n s2: %lf \n", s2.first);
    return (s1.first < s2.first);
}


//int Population::comparePopulationsByTimeOrigin(const void *s1, const void *s2)
//{
//    Population *p1 = *(Population **)s1;
//    Population *p2 = *(Population **)s2;
//    if (  p1->timeOriginInput  > p2 ->timeOriginInput)
//        return 1;
//    else if (p1->timeOriginInput < p2->timeOriginInput)
//        return -1;
//    else
//        return 0;
//}
int Population::compare (const void * a, const void * b)
{
    double *p1 = (double *)a;

    double *p2 = (double *)b;
    
    if (*p1 > *p2) return 1;
    
    else if (*p2 > *p1 )return -1;
    
    else return 0;
    
}
void Population::InitCoalescentEvents(int numClones)
{
    int j = 0;
    for (j = 0; j < sampleSize + numClones - 2; j++) {
        CoalescentEventTimes.push_back(0);
    }
}

void Population::resetActiveGametes()
{
    numCompletedCoalescences = 0;
    nodeIdAncestorMRCA = 0;
    numActiveGametes=0;
    numGametes=0;
    idsActiveGametes.clear();
    idsGametes.clear();
}
