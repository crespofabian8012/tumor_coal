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
        migrationTimes =(double *)  malloc( (p->order + 1) * sizeof( double));
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
void population::SimulatePopulation(  Population** populations,
                        ProgramOptions *programOptions,
                        long int *seed,
                        int *numNodes,
                        int numClones,
                        double      cumNumCA,
                        double meanNumCA,
                        double cumNumMIG,
                        double meanNumMIG,
                        int  *numMIG,
                        int  *numCA,
                        double *numEventsTot,
                        pll_unode_t    **nodes,
                        int *nextAvailable,
                        int*  numActiveGametes,
                        int* labelNodes,
                        double *currentTime,
                        int* eventNum
                        )
{
    int  c, d, i, j, w, k, m, cumIndivid, isCoalescence, whichInd,
    firstInd, secondInd, newInd,  foundSuperflousNode,
    isMigration, whichClone, currentNumberAliveClones;
    double     eventTime;
    pll_unode_t  *p, *q, *r, *r1, *r2, *r3 ;
    TreeNode  *u,*v, *u1, *u2, *u3;
    double    ran;
    double    *cumPopulPart;
    *eventNum = 0;
    int numSimClones = 0;
    double ThisRateCA = 0.0;
    double ThisTimeCA_W = 0.0;
    double ThisTimeCA_V1 = 0.0;
    double ThisTimeCA_V2 = 0.0;
    int numParcialActiveGametes;
    int choosePairIndividuals;
    
      numParcialActiveGametes = this->numActiveGametes;
    
    int numMigrations = numIncomingMigrations; //taking into account also the    time of origin as a migration

    double timeNextMigration;
    int indexNextMigration = 0;
    population *incommingPop;

    if (programOptions->noisy > 1)
        fprintf (stderr, "\n\n>> Simulating evolutionary history of clone %d (number active gametes %d, original time to origin %lf)\n", index, this->numActiveGametes, timeOriginInput);
    
    if (programOptions->noisy > 1)
        fprintf (stderr, "\n\n> Simulating evolutionary history of clone  or order  %d ..\n", order);
    
    *currentTime=0;
    while (indexNextMigration < numMigrations) {
        timeNextMigration = (double)(migrationTimes)[indexNextMigration];
        //fprintf (stderr, "\n\n> numParcialActiveGametes= %d \n", numParcialActiveGametes);
        if ( this->numActiveGametes >= 2) {
            ThisRateCA = (double)  this->numActiveGametes * ((double) this->numActiveGametes - 1) / 2.0;
            ThisTimeCA_W = RandomExponential (ThisRateCA, seed) ;
            ThisTimeCA_V1 = FmodelTstandard (*currentTime, timeOriginSTD, delta);
            ThisTimeCA_V1 = ThisTimeCA_V1 + ThisTimeCA_W;
            // from standard time to model time, GstandardTmodel(V, T, delta)
            ThisTimeCA_V2 = GstandardTmodel(ThisTimeCA_V1, timeOriginSTD, delta);
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
            *numCA = *numCA + 1;
            *eventNum= *eventNum +1;
            whichClone = index;
            *currentTime = ThisTimeCA_V2; // update current time in model time
            eventTime = *currentTime;
            
            if (programOptions->noisy > 1)
            {
                fprintf (stderr, "\n\n*** Event %3d *** currentTime (model time) = %lf, currentTime (standard time) = %lf\n", *eventNum, ThisTimeCA_V2, ThisTimeCA_V1 );
                fprintf (stderr, "\n\n*** Event %3d *** currentTime (input units) = %lf\n", *eventNum, ThisTimeCA_V2);
            }
            if (programOptions->noisy == 4)
                fprintf (stderr, "* Coalescence *\n");
            MakeCoalescenceEvent(  nodes, numClones, seed, programOptions->noisy, numActiveGametes,  nextAvailable,
                                 labelNodes, currentTime,  &(programOptions->numNodes));
           
       
        }
        else
        {
            if (indexNextMigration < numMigrations - 1) //indexNextMigration corresponds to one of the true migrations
            {
                isCoalescence = NO;
                isMigration = YES;
                *numMIG = *numMIG + 1;
                *eventNum= *eventNum +1;
                *currentTime = timeNextMigration;
                eventTime = *currentTime;
                //update migration times in model time
                if (programOptions->noisy > 1)
                {
                    fprintf (stderr, "\n\n*** Event %3d *** *currentTime (model units) = %lf\n", *eventNum, ThisTimeCA_V2);
                }
                if (programOptions->noisy == 4)
                {fprintf (stderr, "* Migration *\n");}
                newInd = *nextAvailable;
                
                
                r1 = *nodes + newInd;   /* new nodelet ancester */
                r1 = *nodes + (newInd + 1);   /* new  nodelet ancester */
                r1 = *nodes + (newInd + 2);   /* new nodelet ancester */
                
              
                r1->node_index = newInd;
                r2->node_index = newInd+1;
                r3->node_index = newInd+2;
                
                r1->data =  malloc(sizeof(TreeNode));
                r2->data =  malloc(sizeof(TreeNode));
                r3->data =  malloc(sizeof(TreeNode));
                
                u1= (TreeNode *)(r1->data);
                u2= (TreeNode *)(r2->data);
                u3= (TreeNode *)(r3->data);
                
                u1->label = u2->label= u3->label= *labelNodes;
                *labelNodes=*labelNodes+1;
                
                u1->indexCurrentClone = u2->indexCurrentClone= u3->indexCurrentClone= index;
                u1->indexCurrentClone = u2->indexCurrentClone = u3->indexCurrentClone = index;
                u1->orderCurrentClone =u2->orderCurrentClone= u3->orderCurrentClone= order;
                u1->nodeClass = u2->nodeClass= u3->nodeClass= 4;
                
                
                //    p = *nodes + MatrixMigrationIDnodeMRCA[ThisCloneNumber][doAmigration]; // root of younger clone
                incommingPop = *((immigrantsPopOrderedModelTime) + indexNextMigration );
                
                if (!incommingPop){
                    fprintf (stderr, "\nError. The incoming population to  poulation %d is empty \n", index);
                    exit (-1);
                }
                //r->indexOldClone = incommingPop->index;
                p = *nodes + (incommingPop->nodeIdAncesterMRCA); // root of younger clone
                indexNextMigration = indexNextMigration + 1;
                
                v= (TreeNode *)p;
                
                v->indexCurrentClone = index;
                v->indexOldClone = incommingPop->index;
                v->orderCurrentClone = order;
                // link the nodes
                //r->left = p;
                r1->back = p;
                r1->next = r2;
                r2->next = r3;
                r3->next = r1;
                
                //r->right = NULL;
                r2->back = NULL;
                r3->back = NULL;
                //choosePairIndividuals = NO;
                
                //ChooseRandomIndividual(&firstInd, numClones, popI,  &secondInd, seed, choosePairIndividuals);
                //q=*nodes + firstInd;
                //r->right = q;//choose another random living individual of the population
                
                p->back = r1;
                //p->anc1 = r;
                //q->anc1 = r;
                
                //connectNodes(p, NULL, r);
                //p->time = *currentTime;
                // p->timePUnits = *currentTime * (popI->effectPopSize);
                
                
                u1->time =u2->time= u3->time= *currentTime;// this is equal to the time of the migration
                u1->timePUnits = u2->timePUnits= u3->timePUnits= *currentTime * (effectPopSize);
                *nextAvailable=*nextAvailable+3; /* 3 nodelets more are available */
                
                k = v->indexCurrentClone;
                incommingPop->numActiveGametes = incommingPop->numActiveGametes - 1; /* now the other clone has 1 less node */
                // remove node from old clone in list of active gametes and add the new node of the current clone
                //popI->idsActiveGametes[popI->numActiveGametes]=r->index;//adding the superfluos node
                idsActiveGametes[this->numActiveGametes]=p->node_index;//adding the superfluos node
                numActiveGametes = numActiveGametes + 1; /* now this clone has 1 more node */
                //                if (noisy > 1)
                //                    fprintf (stderr, "Migration, creating node %d (clone %d) derived from node %d (clone %d)", newInd, popI->index, incommingPop->nodeIdAncesterMRCA, k);
                //fprintf (stderr, "Migration, creating node %d (clone %d) derived from node %d (clone %d)", newInd, ThisCloneNumber, MatrixMigrationIDnodeMRCA[ThisCloneNumber][doAmigration], k);
                if (programOptions->noisy > 1)
                    fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", v->timePUnits);
                // fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", r->timePUnits);
                /* memory for number of nodes */
                if (*nextAvailable >= *numNodes)  /* if there aren't enough nodes it go into and it addition more */
                {
                    /* ReallocNodes(&numNodes, activeGametes); */
                    if (programOptions->noisy == 4)
                        fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
                    *numNodes += INCREMENT_NODES;
                    /* realloc */
                    // *nodes = (TreeNode *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode));
                    *nodes = (pll_unode_t *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode));
                    if (!(*nodes))
                    {
                        fprintf (stderr, "Could not reallocate nodes (%lu bytes)\n", *numNodes  * (long) sizeof(pll_unode_t));
                        exit (-1);
                    }
                    //                    activeGametes = (int *) realloc (activeGametes, *numNodes * (long) sizeof(int));
                    //                    if (!activeGametes)
                    //                    {
                    //                        fprintf (stderr, "Could not reallocate activeGametes (%lu bytes)\n", *numNodes * (long) sizeof(int));
                    //                        exit (-1);
                    //                    }
                }
                
            }
            else {
                //origin reached
                *currentTime = timeNextMigration;
                indexNextMigration = indexNextMigration + 1;
                isCoalescence = NO;
                isMigration = NO;
                eventTime = *currentTime;
                if (programOptions->noisy > 1)
                    fprintf (stderr, "\n\n*** Clone origin ***\n");
                if (programOptions->noisy == 4)
                    fprintf (stderr, "Clone origin %d at time (model units) = %lf\n", index, *currentTime);
                if (order < numClones - 1) // do not do it for the last clone
                {
                    newInd = *nextAvailable;
                    
                    r1 = *nodes + newInd;   /* new nodelet ancester */
                    r1 = *nodes + (newInd + 1);   /* new  nodelet ancester */
                    r1 = *nodes + (newInd + 2);   /* new nodelet ancester */
                    
                    r1->node_index = newInd;
                    r2->node_index = newInd +1;
                    r3->node_index = newInd +2;
                    
                    r1->data =  malloc(sizeof(TreeNode));
                    r2->data =  malloc(sizeof(TreeNode));
                    r3->data =  malloc(sizeof(TreeNode));
                    
                    u1= (TreeNode *)(r1->data);
                    u2= (TreeNode *)(r2->data);
                    u3= (TreeNode *)(r3->data);
                    
                    u1->index =u2->index= u3->index= *nextAvailable;
                    u1->label =u2->label= u3->label= *labelNodes;
                    *labelNodes=*labelNodes+1;
                    u1->indexOldClone =u1->indexCurrentClone = index;
                    u2->indexOldClone =u2->indexCurrentClone = index;
                    u3->indexOldClone =u3->indexCurrentClone = index;
                    
                    u1->orderCurrentClone = u2->orderCurrentClone= u3->orderCurrentClone=order;
                    u1->effectPopSize= u2->effectPopSize=  u3->effectPopSize=effectPopSize;
                    nodeIdAncesterMRCA=newInd;
                    
                    u1->nodeClass = u2->nodeClass =u3->nodeClass =4;
                    
                    //firstInd = *nextAvailable - 1;
                    firstInd = *nextAvailable - 1;
                    //p = nodes + activeGametes[firstInd]; // descendant node (previously generated node, nextAvailable - 1)
                    p = *nodes + firstInd; // descendant node (previously generated node, nextAvailable - 1)
                    // link the nodes
                    //r->left = p;
                    r1->back = p;
                    p->back= r1;
                    r1->next=r2;
                    r2->next=r3;
                    r3->next=r1;
                    r2->back=NULL;
                    r3->back=NULL;
                    //r->right = NULL;
                    //p->anc1 = r;
                    
                    u1->time = u2->time= u3->time= *currentTime;
                    u1->timePUnits =  u2->timePUnits =  u3->timePUnits = *currentTime * effectPopSize;
                    MRCA = p;
                    
                    //connectNodes(p, NULL, r);
                    //fprintf (stderr, "\n r->index = %d, r->time = %lf\n", r->index, r->time);
                    /* readjust active nodes */
                    *nextAvailable=*nextAvailable+1; /* 1 node more is available */
                    idsActiveGametes[0] = newInd;//always will be in the 0  position because there is only one left
                    
                    //popI->idsGametes[popI->numGametes] = newInd; r is a superflous node and it will be removed so no need to add it
                    //popI->numGametes = popI->numGametes +1;
                    
                    if (programOptions->noisy > 1)
                        fprintf (stderr, "Creating origin node, it creates node %d derived from node %d", newInd, firstInd);
                    if (programOptions->noisy > 1)
                        fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", u1->timePUnits);
                    /* memory for number of nodes */
                    if (*nextAvailable >= *numNodes)  /* if there aren't enough nodes it go into and it addition more */
                    {
                        /* ReallocNodes(&numNodes); */
                        if (programOptions->noisy == 4)
                            fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
                        *numNodes += INCREMENT_NODES;
                        /* realloc */
                       // *nodes = (TreeNode *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode));
                         *nodes = (pll_unode_t *) realloc (*nodes, *numNodes  * (long) sizeof(pll_unode_t));
                        if (!(*nodes))
                        {
                            fprintf (stderr, "Could not reallocate nodes (%lu bytes)\n", *numNodes  * (long) sizeof(TreeNode));
                            exit (-1);
                        }
                        
                    }
                }
                else  {//origin of oldest pop reached
                    nodeIdAncesterMRCA=*nextAvailable-1;//for the last population, nodeIdAncesterMRCA is the MRCA instead of ancester of MRCA
                    r = *nodes + *nextAvailable-1;//popI->idsActiveGametes[0]
                     u= (TreeNode *)(r->data);
                    u->indexOldClone = u->indexCurrentClone = index;
                    u->orderCurrentClone = order;
                    MRCA= r;
                    
                    
                }
            }
        }
        if (programOptions->noisy > 3)
        {
            fprintf (stderr, "\nActive nodes (%d):",  this->numActiveGametes);
            for (i = 0; i < this->numActiveGametes; i++)
                fprintf (stderr, " %d", idsActiveGametes[i]);
            fprintf (stderr, "\t|\tNext node available = %d", *nextAvailable);
        }
    }
    if (programOptions->noisy > 1)
        fprintf (stderr, "\n\nEvolutionary history of clone %d is completed \n", index);
    
}
void population::ChooseRandomIndividual(int *firstInd,   int numClones,   int *secondInd, long *seed, int choosePairIndividuals)
{
    double random;
    int k, w;
    double *cumPopulPart = (double *) malloc((numActiveGametes + 1)* (long) sizeof(double));
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
    //    for (k = 0; k < numClones; k++)
    //    {
    //        pop = *(populations + k);
    //        w = w + (pop->numActiveGametes);
    //        //fprintf (stderr, "\nClone %d with %d gametes. w=%d ", k, pop->numActiveGametes, w);
    //    }
    //    if (w != numActiveGametes)
    //    {
    //        fprintf (stderr, "\nError. The sum of partial active gametes is different to the total number of gametes, w %d != numActiveGametes %d. In Coalescence.", w, numActiveGametes);
    //        exit (-1);
    //    }
    random = RandomUniform(seed);
    //fprintf (stderr, "\nran = %lf ", ran);
    *firstInd = bbinClones(random, cumPopulPart, numActiveGametes)-1;
    w = 0;
    
    if (*firstInd >= numActiveGametes || *firstInd < 0 ) /* checking */
    {
        fprintf (stderr, "\n\nERROR: firstInd out of range!\n");
        exit (-1);
    }
    
    if (choosePairIndividuals== YES && numActiveGametes > 1) {
        
        do//choose randomly another individual to coalesce
        {
            random = RandomUniform(seed);
            *secondInd = bbinClones(random, cumPopulPart, numActiveGametes)-1;
            
        } while (*firstInd == *secondInd  );
    }
    free (cumPopulPart);
    cumPopulPart=NULL;
}

void population::MakeCoalescenceEvent(pll_unode_t **nodes, int numClones, long int* seed, int noisy,  int *numActiveGametes,   int* nextAvailable,
                                      int*labelNodes, double *currentTime, int *numNodes)
{
    int k, w;
    double rand, ran;
    // TreeNode  *p, *q, *r;
    pll_unode_t  *p, *q, *r, *r1, *r2, *r3 ;
    TreeNode *u1, *u2, *u3;
    int firstInd, i, j, secondInd=0, newInd=0;
    int choosePairIndividuals = YES;
    
    ChooseRandomIndividual(&firstInd, numClones,   &secondInd, seed, choosePairIndividuals);
    
    newInd = *nextAvailable;
    if (noisy > 1)
        fprintf (stderr, "Coalescence involving %d and %d to create node %d (in clone %d)", idsActiveGametes[firstInd], idsActiveGametes[secondInd], newInd, index);
    /*  set pointers between nodes */
    p = *nodes + idsActiveGametes[firstInd];
    q = *nodes + idsActiveGametes[secondInd];
    
 
    r1 = *nodes + newInd;   /* new nodelet ancester */
    r2 = *nodes + (newInd + 1);   /* new  nodelet ancester */
    r3 = *nodes + (newInd + 2);   /* new nodelet ancester */
    
    r1->node_index = newInd;
    r2->node_index = newInd +1;
    r3->node_index = newInd +2;
    
    r1->data =  malloc(sizeof(TreeNode));
    r2->data =  malloc(sizeof(TreeNode));
    r3->data =  malloc(sizeof(TreeNode));
    
    u1= (TreeNode *)(r1->data);
    u2= (TreeNode *)(r2->data);
    u3= (TreeNode *)(r3->data);
    
    u1->index =  u2->index=  u3->index= *nextAvailable;
    u1->label = u2->label =u3->label =*labelNodes;
    *labelNodes=*labelNodes+1;
    u1->indexOldClone = u1->indexCurrentClone = index;//here the clone number is updated
     u2->indexOldClone = u2->indexCurrentClone = index;//here the clone number is updated
     u3->indexOldClone = u3->indexCurrentClone = index;//here the clone number is updated
    // r->indexCurrentClone = p->indexCurrentClone;
    // r->orderCurrentClone = p->orderCurrentClone;
    u1->orderCurrentClone = u2->orderCurrentClone =u3->orderCurrentClone =order;
    u1->effectPopSize=u2->effectPopSize=u3->effectPopSize=effectPopSize;
    u1->nodeClass =u2->nodeClass =u3->nodeClass = 4;
    // link the nodes
    r1->back = p;
    p->back= r1;
    r1->next=r2;
    r2->next=r3;
    r3->next=r1;
    r2->back=q;
    r3->back=NULL;
    
//    r->left = p;
//    r->right = q;
//    p->anc1 = r;
//    q->anc1 = r;
    u1->time = u2->time = u3->time =*currentTime;
    u1->timePUnits = u2->timePUnits = u3->timePUnits =*currentTime * (effectPopSize);
    
    //fprintf (stderr, "\n r->index = %d, r->time = %lf, ClonePopSizeMeffectBegin[ThisCloneNumber] = %lf, ThisCloneNumber = %d\n", r->index, r->time, ClonePopSizeMeffectBegin[ThisCloneNumber], ThisCloneNumber);
    if (noisy > 1)
        fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", u1->timePUnits);
    /* readjust active nodes */
    
    // idsActiveGametes[firstInd] = newInd;
    idsActiveGametes[firstInd] = r3->node_index;// the r3 node
    idsActiveGametes[secondInd] = idsActiveGametes[this->numActiveGametes - 1];
    *numActiveGametes = *numActiveGametes - 1; /* less 1 active node */
    
    //update list ids nodes
    idsGametes[numGametes] = newInd;
    numGametes = numGametes +1;
    
   // *nextAvailable=*nextAvailable+1; /* 1 node more is available */
    *nextAvailable=*nextAvailable+3; /* 1 node more is available */
    
    CoalescentEventTimes[ numCompletedCoalescences]=  u1->time;
    numActiveGametes = numActiveGametes - 1; /* now this clone
                                                          has 1 less node */
    
    numCompletedCoalescences= numCompletedCoalescences+1;
    /* memory for number of nodes */
    if (*nextAvailable >= *numNodes)  /* if there aren't enough nodes it go into and it addition more */
    {
        /* ReallocNodes(&numNodes, activeGametes); */
        if (noisy == 4)
            fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
        *numNodes += INCREMENT_NODES;
        /* realloc */
        //*nodes = (TreeNode *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode) );
        *nodes = (pll_unode_t *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode) );
        if (!(*nodes))
        {
            fprintf (stderr, "Could not reallocate nodes (%lu bytes)\n", *numNodes  * (long) sizeof(TreeNode));
            exit (-1);
        }
        idsActiveGametes = (int *) realloc (idsActiveGametes, *numNodes * (long) sizeof(int));
        if (!(idsActiveGametes))
        {
            fprintf (stderr, "Could not reallocate idsActiveGametes for the current population(%lu bytes)\n", *numNodes * (long) sizeof(int));
            exit (-1);
        }
    }

}
