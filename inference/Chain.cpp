//
//  Chain.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 7/10/19.
//

#include "Chain.hpp"
using namespace std;
#include "libpll/pll_optimize.h"
#include "libpll/pll_tree.h"
#include "libpll/pllmod_algorithm.h"
#include "libpll/pll_msa.h"
#include "libpll/pll.h"
#include "data_utils.hpp"
#include "definitions.h"
//Parametrized Constructor

chain::chain( int chainNumber,
      int numClones,
      int gammaParam,
      int totalPopSize,
      double mutationRate,
      double seqErrorRate,
      double dropoutRate
             )

{
    if (chainNumber >=0)
        this->chainNumber=chainNumber;
   
    if (numClones >=0)
        this->numClones=numClones;
    
    if (gammaParam >=0)
        this->gammaParam=gammaParam;
    
    if (mutationRate >=0)
        this->mutationRate=mutationRate;
    
    if (seqErrorRate >=0)
        this->seqErrorRate=seqErrorRate;
    
    if (dropoutRate >=0)
        this->dropoutRate=dropoutRate;
    
    currentNumberIerations = 0;
    
    if (totalPopSize >=0)
        this->totalPopSize= totalPopSize;
    
     currentlogConditionalLikelihoodTree=0;
     currentlogConditionalLikelihoodSequences=0;
    
    
}
int chain::setIntitialTree(char * NewickString){
    
    this->initialTree = pll_utree_parse_newick_string_unroot(NewickString);
    if (this->initialTree ==NULL)
        return 1;
    else
        return 0;
    
}
void chain::MakeCoalescenceEvent(population *population, pll_unode_t **nodes, int numClones, long int* seed, int noisy,  int *numActiveGametes,   int* nextAvailable,
                                      int*labelNodes, double *currentTime, int *numNodes)
{
    int k, w;
    double rand, ran;
    // TreeNode  *p, *q, *r;
    pll_unode_t  *p, *q, *r, *r1, *r2, *r3 ;
    TreeNode *u1, *u2, *u3;
    int firstInd, i, j, secondInd=0, newInd=0;
    int choosePairIndividuals = YES;
    
    population->ChooseRandomIndividual(&firstInd, numClones,   &secondInd, seed, choosePairIndividuals);
    
    newInd = *nextAvailable;
    if (noisy > 1)
        fprintf (stderr, "Coalescence involving %d and %d to create node %d (in clone %d)", population->idsActiveGametes[firstInd], population->idsActiveGametes[secondInd], newInd, population->index);
    /*  set pointers between nodes */
    p = *nodes + population->idsActiveGametes[firstInd];
    q = *nodes + population->idsActiveGametes[secondInd];
    
    
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
    u1->indexOldClone = u1->indexCurrentClone = population->index;//here the clone number is updated
    u2->indexOldClone = u2->indexCurrentClone = population->index;//here the clone number is updated
    u3->indexOldClone = u3->indexCurrentClone = population->index;//here the clone number is updated
    // r->indexCurrentClone = p->indexCurrentClone;
    // r->orderCurrentClone = p->orderCurrentClone;
    u1->orderCurrentClone = u2->orderCurrentClone =u3->orderCurrentClone =population->order;
    u1->effectPopSize=u2->effectPopSize=u3->effectPopSize=population->effectPopSize;
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
    u1->timePUnits = u2->timePUnits = u3->timePUnits =*currentTime * (population->effectPopSize);
    
    //fprintf (stderr, "\n r->index = %d, r->time = %lf, ClonePopSizeMeffectBegin[ThisCloneNumber] = %lf, ThisCloneNumber = %d\n", r->index, r->time, ClonePopSizeMeffectBegin[ThisCloneNumber], ThisCloneNumber);
    if (noisy > 1)
        fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", u1->timePUnits);
    /* readjust active nodes */
    
    // idsActiveGametes[firstInd] = newInd;
    population->idsActiveGametes[firstInd] = r3->node_index;// the r3 node
    population->idsActiveGametes[secondInd] = population->idsActiveGametes[population->numActiveGametes - 1];
    *numActiveGametes = *numActiveGametes - 1; /* less 1 active node */
    
    //update list ids nodes
    population->idsGametes[population->numGametes] = newInd;
    population->numGametes = population->numGametes +1;
    
    // *nextAvailable=*nextAvailable+1; /* 1 node more is available */
    *nextAvailable=*nextAvailable+3; /* 1 node more is available */
    
    population->CoalescentEventTimes[ population->numCompletedCoalescences]=  u1->time;
    numActiveGametes = numActiveGametes - 1; /* now this clone
                                              has 1 less node */
    
    population->numCompletedCoalescences= population->numCompletedCoalescences+1;
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
        population->idsActiveGametes = (int *) realloc (population->idsActiveGametes, *numNodes * (long) sizeof(int));
        if (!(population->idsActiveGametes))
        {
            fprintf (stderr, "Could not reallocate idsActiveGametes for the current population(%lu bytes)\n", *numNodes * (long) sizeof(int));
            exit (-1);
        }
    }
    
}

void chain::SimulatePopulation(population *pop,ProgramOptions *programOptions,
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
    
    numParcialActiveGametes = pop->numActiveGametes;
    
    int numMigrations = pop->numIncomingMigrations; //taking into account also the    time of origin as a migration
    
    double timeNextMigration;
    int indexNextMigration = 0;
    population *incommingPop;
    
    if (programOptions->noisy > 1)
        fprintf (stderr, "\n\n>> Simulating evolutionary history of clone %d (number active gametes %d, original time to origin %lf)\n", pop->index, pop->numActiveGametes, pop->timeOriginInput);
    
    if (programOptions->noisy > 1)
        fprintf (stderr, "\n\n> Simulating evolutionary history of clone  or order  %d ..\n", pop->order);
    
    *currentTime=0;
    while (indexNextMigration < numMigrations) {
        timeNextMigration = (double)(pop->migrationTimes)[indexNextMigration];
        //fprintf (stderr, "\n\n> numParcialActiveGametes= %d \n", numParcialActiveGametes);
        if ( pop->numActiveGametes >= 2) {
            ThisRateCA = (double)  pop->numActiveGametes * ((double) pop->numActiveGametes - 1) / 2.0;
            ThisTimeCA_W = RandomExponential (ThisRateCA, seed) ;
            ThisTimeCA_V1 = population::FmodelTstandard (*currentTime, pop->timeOriginSTD, pop->delta);
            ThisTimeCA_V1 = ThisTimeCA_V1 + ThisTimeCA_W;
            // from standard time to model time, GstandardTmodel(V, T, delta)
            ThisTimeCA_V2 = GstandardTmodel(ThisTimeCA_V1, pop->timeOriginSTD, pop->delta);
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
            whichClone = pop->index;
            *currentTime = ThisTimeCA_V2; // update current time in model time
            eventTime = *currentTime;
            
            if (programOptions->noisy > 1)
            {
                fprintf (stderr, "\n\n*** Event %3d *** currentTime (model time) = %lf, currentTime (standard time) = %lf\n", *eventNum, ThisTimeCA_V2, ThisTimeCA_V1 );
                fprintf (stderr, "\n\n*** Event %3d *** currentTime (input units) = %lf\n", *eventNum, ThisTimeCA_V2);
            }
            if (programOptions->noisy == 4)
                fprintf (stderr, "* Coalescence *\n");
            MakeCoalescenceEvent( pop, nodes, numClones, seed, programOptions->noisy, numActiveGametes,  nextAvailable,
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
                r2 = *nodes + (newInd + 1);   /* new  nodelet ancester */
                r3 = *nodes + (newInd + 2);   /* new nodelet ancester */


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

                u1->indexCurrentClone = u2->indexCurrentClone= u3->indexCurrentClone= pop->index;
                u1->indexCurrentClone = u2->indexCurrentClone = u3->indexCurrentClone = pop->index;
                u1->orderCurrentClone =u2->orderCurrentClone= u3->orderCurrentClone= pop->order;
                u1->nodeClass = u2->nodeClass= u3->nodeClass= 4;
                
                
                //    p = *nodes + MatrixMigrationIDnodeMRCA[ThisCloneNumber][doAmigration]; // root of younger clone
                incommingPop = *((pop->immigrantsPopOrderedModelTime) + indexNextMigration );
                
                if (!incommingPop){
                    fprintf (stderr, "\nError. The incoming population to  poulation %d is empty \n", pop->index);
                    exit (-1);
                }
                //r->indexOldClone = incommingPop->index;
                p = *nodes + (incommingPop->nodeIdAncesterMRCA); // root of younger clone
                indexNextMigration = indexNextMigration + 1;
                
                v= (TreeNode *)p;
                
                v->indexCurrentClone = pop->index;
                v->indexOldClone = incommingPop->index;
                v->orderCurrentClone = pop->order;
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
                
              
                //p->time = *currentTime;
                // p->timePUnits = *currentTime * (popI->effectPopSize);
                
                
                u1->time =u2->time= u3->time= *currentTime;// this is equal to the time of the migration
                u1->timePUnits = u2->timePUnits= u3->timePUnits= *currentTime * (pop->effectPopSize);
                *nextAvailable=*nextAvailable+3; /* 3 nodelets more are available */
                
                k = v->indexCurrentClone;
                incommingPop->numActiveGametes = incommingPop->numActiveGametes - 1; /* now the other clone has 1 less node */
                // remove node from old clone in list of active gametes and add the new node of the current clone
                //popI->idsActiveGametes[popI->numActiveGametes]=r->index;//adding the superfluos node
                pop->idsActiveGametes[pop->numActiveGametes]=p->node_index;//adding the superfluos node
                pop->numActiveGametes = pop->numActiveGametes + 1; /* now this clone has 1 more node */
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
                    fprintf (stderr, "Clone origin %d at time (model units) = %lf\n", pop->index, *currentTime);
                if (pop->order < numClones - 1) // do not do it for the last clone
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
                    u1->indexOldClone =u1->indexCurrentClone = pop->index;
                    u2->indexOldClone =u2->indexCurrentClone = pop->index;
                    u3->indexOldClone =u3->indexCurrentClone = pop->index;
                    
                    u1->orderCurrentClone = u2->orderCurrentClone= u3->orderCurrentClone= pop->order;
                    u1->effectPopSize= u2->effectPopSize=  u3->effectPopSize= pop->effectPopSize;
                    pop->nodeIdAncesterMRCA=newInd;
                    
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
                    u1->timePUnits =  u2->timePUnits =  u3->timePUnits = *currentTime * pop->effectPopSize;
                    pop->MRCA = p;
                    
                   
                    
                    //fprintf (stderr, "\n r->index = %d, r->time = %lf\n", r->index, r->time);
                    /* readjust active nodes */
                    *nextAvailable=*nextAvailable+1; /* 1 node more is available */
                    pop->idsActiveGametes[0] = newInd;//always will be in the 0  position because there is only one left
                    
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
                    pop->nodeIdAncesterMRCA=*nextAvailable-1;//for the last population, nodeIdAncesterMRCA is the MRCA instead of ancester of MRCA
                    r = *nodes + *nextAvailable-1;//popI->idsActiveGametes[0]
                    u= (TreeNode *)(r->data);
                    u->indexOldClone = u->indexCurrentClone = pop->index;
                    u->orderCurrentClone = pop->order;
                    pop->MRCA= r;
                    
                    
                }
            }
        }
        if (programOptions->noisy > 3)
        {
            fprintf (stderr, "\nActive nodes (%d):",  pop->numActiveGametes);
            for (i = 0; i < pop->numActiveGametes; i++)
                fprintf (stderr, " %d", pop->idsActiveGametes[i]);
            fprintf (stderr, "\t|\tNext node available = %d", *nextAvailable);
        }
    }
    if (programOptions->noisy > 1)
        fprintf (stderr, "\n\nEvolutionary history of clone %d is completed \n", pop->index);
    
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

void chain::BuildTree(population *olderPopulation, long int *seed,
               ProgramOptions *programOptions,
               pll_unode_t    **nodes,
               pll_unode_t   **treeTips,
               pll_unode_t    **treeRootInit,
               int *nextAvailable,
               int *newInd,
               double *currentTime,
               int *labelNodes
               )
{
    int i, j, k;
    int indexCurrentTip;
    int  foundSuperflousNode;
    pll_unode_t *p, *q, *r;
    TreeNode *u, *anc, *healthyR, *u1, *u2, *u3;
    /********** BUILDING TREES ***********/
    if (programOptions->noisy > 1)
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
        for (i = 0; i < *nextAvailable; i++) // available all the nodes
        {
            p = *nodes + i;
            
            //fprintf (stderr, "\n\np->index = %d", p->index);
            
            //if (p->left == NULL && p->right == NULL && p->anc1 == NULL)
            if (p->next == NULL && p->back == NULL )
            {
                // nothing to do with this node because it is not connected to anything
                //fprintf (stderr, "\n * nothing to do with this node because it is not connected to anything");
            }
           // else if (p->left == NULL && p->right == NULL && p->anc1 != NULL)
            else if (p->back == NULL && p->next != NULL && p->next->back == NULL  )//first nodelet in post order of a leaf node
            {
                // (*treeTips)[indexCurrentTip]=p;
                if(indexCurrentTip <  programOptions->TotalNumSequences){
                    //treeTips[indexCurrentTip]=p;
                    treeTips[indexCurrentTip]=p->next->next;
                    indexCurrentTip++;
                }
               
                // do not do anything with this node because it is a tip
                //fprintf (stderr, "\n * do not do anything with this node because it is a tip");
            }
            else if (p->back == NULL && p->next != NULL && p->next->back != NULL)//second nodelet in post order of a leaf node
            {
                //do nothing because the tip is already added
                
            }
            //else if (p->left != NULL && p->right == NULL && p->anc1 != NULL)
            else if (p->back != NULL && p->next != NULL && p->next->back == NULL && p->next->next!=NULL)//first nodelet in post order of a superflous node
            {
                // this is a superflous node and can be removed(this superfluos nodes are the MRCA nodes of the demes
                foundSuperflousNode = YES;
                //q = p->left;
                //r = p->anc1;
                q = p->back;
                r = p->next->next->back;
                
                //if (p->anc1->left == p)  // p->anc up, p->left down, total: up and down for left
               // {
                    r->back = q;
                    q->back = r;
                p->next=NULL;
                p->back = NULL;
                p->next->next=NULL;
                p->next->back = NULL;
                p->next->next->next=NULL;
                p->next->next->back = NULL;
                   // r->left = q;
                    //q->anc1 = r;
                    //p->left = NULL;
                    //p->anc1 = NULL;
               
               // }
              //  else
               // {
//                    r->right = q;
//                    q->anc1 = r;
//                    p->left = NULL;
//                    p->anc1 = NULL;
//
              //  }
                
                //fprintf (stderr, "\n - this is a superflous node and can be removed (1)");
            }
           // else if (p->left == NULL && p->right != NULL && p->anc1 != NULL)
            else if (p->back != NULL && p->next != NULL && p->next->back != NULL && p->next->next!=NULL && p->next->next->back==NULL)
                //second nodelet in post order of a superflous node
            {
                // this is a superflous node and can be removed
                foundSuperflousNode = YES;
                q = p->back;
                r = p->next->back;
                
                r->back = q;
                q->back = r;
                p->next=NULL;
                p->back = NULL;
                p->next->next=NULL;
                p->next->back = NULL;
                p->next->next->next=NULL;
                p->next->next->back = NULL;
                
                //fprintf (stderr, "\n - this is a superflous node and can be removed (2)");
            }
           // else if (p->left != NULL && p->right != NULL && p->anc1 != NULL)
            else if (p->back != NULL && p->next!= NULL && p->next->back != NULL && p->next->next!= NULL && p->next->next->back!= NULL)
            {
               
                // this is an internal node formed by a coalescence event, do not touch
                //fprintf (stderr, "\n * this is an internal node formed by a coalescence event, do not touch");
            }
            //else if (p->left != NULL && p->right != NULL && p->anc1 == NULL)
            else if (p->back != NULL && p->next!= NULL && p->next->back != NULL && p->next->next!= NULL && p->next->next->back== NULL)// first nodelet of the root
            {
               
                // this is the last (coalescence event) in the tree, MRCA
                //fprintf (stderr, "\n * this is the last (coalescence event) in the tree, MRCA");
            }
            else if (p->back != NULL && p->next!= NULL && p->next->back == NULL && p->next->next!= NULL && p->next->next->back!= NULL)// second nodelet of the root
            {
                
                
            }
            else if (p->back == NULL && p->next!= NULL && p->next->back != NULL && p->next->next!= NULL && p->next->next->back!= NULL)// third nodelet of the root
            {
                
            }
            //else if (p->left != NULL && p->right == NULL && p->anc1 == NULL)
            else if (p->back != NULL && p->next!=NULL && p->next->back == NULL &&  p->next->next!= NULL && p->next->next->back == NULL)
            {
                // Seems to be the last coalescent event among sequences with non-ancestral material
                // it is not superfluous, we just remove it
                p->back->back = NULL;
                //p->left->anc1 = NULL;
                //fprintf (stderr, "\n - this is a superflous node and can be removed (3)");
            }
           // else if (p->left == NULL && p->right != NULL && p->anc1 == NULL)
            else if (p->back == NULL && p->next != NULL && p->next->back !=NULL && p->next->next!= NULL &&  p->next->next->back == NULL)
            {
                // not clear what this node could be doing, but we will remove it anyway
                fprintf (stderr, "strange\n");
                p->back=NULL;
                p->next->back =NULL;
                //p->left = NULL;
               // p->right->anc1 = NULL;
                //fprintf (stderr, "\n - this is a superflous node and can be removed (4)");
            }
            else
            {
                fprintf (stderr, "You should not be here, I think\n");
                //fprintf (stderr, "%d %d-- %d %d %d\n", Index(p), j, Index(p->left), Index(p->right), Index(p->anc1));
            }
            if (p->back != NULL)
            {//update length field
                
                //   p->length = p->anc1->time- p->time;
                //p->length = (p->anc1->timePUnits- p->timePUnits);
                 u = (TreeNode *)(p->data);
                 anc = (TreeNode *)(p->back->data);
                u->length =(anc->timePUnits- u->timePUnits);
                p->length = (anc->timePUnits- u->timePUnits);
                //*mutationRate;
                u->lengthModelUnits = (anc->time- u->time);
                //*mutationRate;
               
            }
        }
        //fprintf (stderr, "\n");
    }//while
    
    /* about the MRCA */
    *newInd=*nextAvailable-1;
    p = *nodes + *newInd; /* because the last one event is the last coalescence */
    //fprintf (stderr, "\n\n\n>> newInd = %d\n", newInd);
    
    if (programOptions->thereisOutgroup == NO)
    {
        p = *nodes + *newInd;
        u = (TreeNode *)(p->data);
        u->nodeClass = 5;
        *treeRootInit = p;
        //treeRootInit[0] = p;
        u->anc1 = NULL;
    }
    if (programOptions->thereisOutgroup == YES && programOptions->outgroupSelection > 0)  /*** Root and outgroup ***/
    {
        p = *nodes + *newInd; // MRCA
        u = (TreeNode *)(p->data);
        u->nodeClass = 4;
        
        if (programOptions->noisy > 1)
            fprintf (stderr, "\n\n>> Attaching outgroup .. ");
        
        //fprintf (stderr, "\n>> ThisCloneNumber = %d, ListMigrationTimesInitial[ThisCloneNumber] = %lf, ClonePopSizeMeffectBegin[ThisCloneNumber] = %lf \n", ThisCloneNumber, ListMigrationTimesInitial[ThisCloneNumber], ClonePopSizeMeffectBegin[ThisCloneNumber]);
        
        if (programOptions->outgroupSelection == 1)  /*** Root 2 times and outgroup ***/
            *currentTime = olderPopulation->timeOriginSTD; // origin of the clone; // currentTime + (outgroupBranchLength_Root1Root2 / mutationRate); // set time of the new root (from which the MRCA and outgroup nodes are derived)
        else if (programOptions->outgroupSelection == 2) { /*** Root 2 times and outgroup ***/
            *currentTime = olderPopulation->timeOriginSTD + (programOptions->outgroupBranchLength_Root1Root2 / olderPopulation->effectPopSize) ; // origin of the clone + time given by the user
            
        }
        else
        {
            fprintf (stderr, "\n\nError simulationg the outgroup. Check input settings\n");
            PrintUsage();
        }
        //TreeNode*       healthyRoot = *nodes + *nextAvailable;
        pll_unode_t*       healthyRoot = *nodes + *nextAvailable;
        u = (TreeNode *)(p->data);
        healthyR = (TreeNode *)(healthyRoot->data);
        healthyR->index = *nextAvailable;
        healthyR->label = *labelNodes;
        healthyR->effectPopSize= u->effectPopSize;
        *labelNodes=*labelNodes+1;
        //healthyRoot->left = p;//coalTreeMRCA;
        healthyRoot->back = p;//coalTreeMRCA;
        //        coalTreeMRCA->anc = healthyRoot;
        //p->anc1 = healthyRoot;
        p->back = healthyRoot;
        
       // healthyRoot->timePUnits = p->timePUnits * healthyRoot->effectPopSize;
        healthyR->timePUnits = u->timePUnits * healthyR->effectPopSize;
        healthyR->nodeClass = 5;
        //        coalTreeMRCA->length = transformingBranchLength/mutationRate;
        p->length = 0;
        u->length = 0;
        //        coalTreeMRCA->branchLength = transformingBranchLength;
        u->lengthModelUnits = 0;
        
        //        healthyRoot->time = currentTime +  transformingBranchLength/mutationRate;
        healthyR->time = *currentTime  ;
        
        int transformingBranchLength=1.001;
        // healthyRoot->time = p->time * transformingBranchLength ;
        healthyR->timePUnits = *currentTime * healthyR->effectPopSize;
        //p->length = (p->anc1->timePUnits- p->timePUnits);
         anc = (TreeNode *)(p->back->data);
        p->length = (anc->timePUnits- u->timePUnits);
        u->length = (anc->timePUnits- u->timePUnits);
        //*mutationRate;
        u->lengthModelUnits = (anc->time- u->time);
        //*mutationRate;
        
        healthyRoot->length = 0;
        healthyR->length = 0;
        //        healthyRoot->length = 0;
        
        //        if (noisy > 2)
        //            fprintf (stderr, "DONE");
        //
        (*nextAvailable)++;
        //
        //        /* connect the healthy ancestral cell with the tip healthy cell*/
        //        if (noisy > 2)
        //            fprintf (stderr, "\n>> Adding healthy tip ... ");
        //TreeNode* healthyTip = *nodes + *nextAvailable;
        pll_unode_t* healthyTip1 = *nodes + *nextAvailable;
        pll_unode_t* healthyTip2 = *nodes + (*nextAvailable
        +1);
        pll_unode_t* healthyTip3 = *nodes + (*nextAvailable + 2);
        
        healthyTip1->back = healthyTip2->back = NULL;
        healthyTip1->next = healthyTip2 ;
        healthyTip2->next = healthyTip3;
        healthyTip3->next =healthyTip1;
       
        healthyTip1->data =  malloc(sizeof(TreeNode));
        healthyTip2->data =  malloc(sizeof(TreeNode));
        healthyTip3->data =  malloc(sizeof(TreeNode));
        
        u1= (TreeNode *)(healthyTip1->data);
        u2= (TreeNode *)(healthyTip2->data);
        u3= (TreeNode *)(healthyTip3->data);
        
        
        u1->effectPopSize= u2->effectPopSize=u3->effectPopSize=healthyR->effectPopSize;
        
     
        healthyTip3->back =healthyRoot;
        healthyRoot->back=healthyTip3;
        //healthyTip->anc1 = healthyRoot;
       // healthyRoot->right = healthyTip;
        u1->time = u2->time =u3->time =0;
        u1->timePUnits =  u2->timePUnits = u3->timePUnits =0;
        double  healthyTipBranchLengthRatio =1;
        
        healthyTip1->length = healthyTip2->length=healthyTip3->length=(healthyR->timePUnits- u1->timePUnits);
        u1->length = u2->length=u3->length=(healthyR->timePUnits- u1->timePUnits);
        
        u1->lengthModelUnits =  u2->lengthModelUnits = u3->lengthModelUnits =(healthyR->time- u1->time);
        
        u1->isOutgroup=u2->isOutgroup=u3->isOutgroup= YES;
        
       
        
        
        *treeRootInit=healthyRoot;
        
    }
    
    int intLabel = 0;
    if (programOptions->noisy > 1)
        fprintf (stderr, "\n\n>> Relabeling nodes on tree... \n\n");
    if (programOptions->thereisOutgroup == YES)
        intLabel = programOptions->TotalNumSequences + 1;
    else
        intLabel = programOptions->TotalNumSequences;
    
   // RelabelNodes(*treeRootInit, treeRootInit, &intLabel );
    
}
void chain::MakeCoalescenceTree (long int *seed,
                           int *numNodes,
                           int numClones,
                           ProgramOptions *programOptions,
                           double      cumNumCA,
                           double meanNumCA,
                           double cumNumMIG,
                           double meanNumMIG,
                           int  *numMIG,
                           int  *numCA,
                           double *numEventsTot,
//                           pll_unode_t** nodes,
//                           pll_unode_t** treeTips,
//                           pll_unode_t    **treeRootInit,
                           char* ObservedCellNames[],
                           int *sampleSizes
                           ) {
    
    int      c, d, i, j, w, k, m, cumIndivid, *activeGametes = NULL, isCoalescence, whichInd,
    firstInd, secondInd, newInd, eventNum, numActiveGametes, foundSuperflousNode,
    isMigration, whichClone, currentNumberAliveClones;
    int     labelNodes;
    double    currentTime, eventTime;
    pll_unode_t  *p, *q, *r;
    double    ran;
    double    *cumPopulPart;
    int     *numParcialActiveGametes, *CumSamNodes;
    int         *fromClone;
    int         ThisCloneNumber, ThisOriginCloneNumber;
    double      minValue;
    double      ThisRateCA;
    double      ThisTimeCA_W;
    double      ThisTimeCA_V1;
    double      ThisTimeCA_V2;
    
    int         doAmigration;
    int         ThisCloneNumberMigrations, ThisM;
    int nextAvailable;
    population *pop;
    
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
    
    //*numNodes = 2 * TotalNumSequences * numClones+ 1; // (2 * TotalNumSequences) + numClones(superfluos) - 1, but let's allocate some more..
    currentNumberAliveClones = numClones;
    
    //InitListPossibleMigrations(populations, numClones);
    for (i = 0; i < numClones; i++){
        pop = populations[i];
        pop->resetMigrationsList();
    }
    //resetMigrationsList( populations,  numClones);
    
    //allocate memory for the nodes
    nodes = (pll_unode_t *) malloc (3*(programOptions->numNodes + 1)* (sizeof(pll_unode_t) + 2*MAX_NAME * sizeof(char) + sizeof(TreeNode) )); /* nodes */
    if (!(nodes))
    {
        fprintf (stderr, "Could not allocate nodes (%lu bytes)\n", 3*(programOptions->numNodes+ 1)  * (long) sizeof(pll_unode_t));
        exit (1);
    }
    
    for (i=0; i< programOptions->TotalNumSequences; i++){
        treeTips[i]=NULL;
    }

    for (i = 0; i < *numNodes; i++)
    {
        //p = (*nodes + i);
        p = &(nodes[i]);
        p->next = NULL;
        p->back = NULL;
        p->data = NULL;
        p->length = 0;
        p->node_index = 0;
        //p->label = 0;

    }
  
    
    AssignSequencesToPopulations( programOptions, *numNodes, programOptions->noisy, programOptions->TotalNumSequences, &numActiveGametes,  &nextAvailable,
                                       &labelNodes, ObservedCellNames, programOptions->doUseObservedCellNames, sampleSizes);
    population *currentPop;
    population *fatherPop;
    i=0;
    currentTime=0.0;
    while (i < numClones) {
        //currentPop = *(populations + i);
        currentPop = *(populations + i);
        SimulatePopulation( currentPop, programOptions, seed,
                           &(programOptions->numNodes),
                           numClones,
                           cumNumCA,
                           meanNumCA,
                           cumNumMIG,
                           meanNumMIG,
                           numMIG,
                           numCA,
                           numEventsTot,
                           &nodes,
                           &nextAvailable ,
                           &numActiveGametes,
                           &labelNodes,
                           &currentTime,
                           &eventNum);
        if (i< numClones-1)   //if it is not the last one
        {
            //choose the father population from which the population i came
            fatherPop= ChooseFatherPopulation( numClones, currentPop, seed,  programOptions->noisy);
            currentPop->FatherPop = fatherPop;
            //update list of migrant times
            population::UpdateListMigrants( numClones, currentPop, fatherPop);
        }
        i = i + 1;
    }
  
     BuildTree(currentPop,seed,
                   programOptions,
                   &nodes,
                   treeTips,
                   &root,
                   &nextAvailable,
                   &newInd,
                   &currentTime,
                   &labelNodes
                   );

    if (programOptions->noisy > 1)
        fprintf (stderr, "\n\n>> Relabeling nodes on tree... \n\n");
   
    
    // RelabelNodes(treeRootInit[0], treeRootInit, &intLabel );
    }
population* chain::ChooseFatherPopulation( int numClones, population  *PopChild,  long int *seed, int noisy) {
    
    population *pOrigin, *pTarget, *p;
    double pij, ran;
    double *ptr;
    int i, j, k;
    double cumProb[numClones  - (int)(PopChild->order)];
    cumProb[0] = 0.0;
    for (j = PopChild->order + 1; j < numClones; j++)
    {
        cumProb[j - PopChild->order] = 0.0;
        p = *(populations + j);
        pij= PopChild->ProbabilityComeFromPopulation ( p, populations,  numClones);
        //pij = ProbabilityCloneiFromClonej2(PopChild, p, populations, numClones);
        cumProb[j - PopChild->order] = cumProb[j - 1 - PopChild->order] + pij;
        
    }
    // now selecting the ancestral clone
    ran = RandomUniform(seed);
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
    population *result =  *(populations + PopChild->order + w);
    if (noisy > 3)
        fprintf (stderr, "\nClone %d derived from clone %d\n", PopChild->index, result->index); // clone ThisCloneNumber (i) is originated from clone ThisOriginCloneNumber (j)
    /* Update list of migation times considering that clone i comes from clone j */
    if (noisy > 1)
        fprintf (stderr, "\n*** Updating list of migration times (considering that clone %d comes from clone %d) ..\n", PopChild->index, result->index);
    return (result); //the father population has  order  (PopChild->order) + w
}
void chain::AssignSequencesToPopulations( ProgramOptions* programOptions,
                                         int numNodes, int noisy,  int TotalNumSequences, int *numActiveGametes, int* nextAvailable,
                                        int *labelNodes, char* ObservedCellNames[], int doUseObservedCellNames, int *sampleSizes)
{
    population *pop;
    pll_unode_t *p;
    TreeNode *u;
    int i, j;
    double currentSampleSize;
    *numActiveGametes = 0;
    int indexFirstObservedCellName;
    
    int*  CumSamNodes = (int *) malloc ((numClones + 1)* (long) sizeof(int)); /* cumulative of samples in clones */
    if (!CumSamNodes)
    {
        fprintf (stderr, "Could not allocate CumSamNodes (%lu bytes)\n", (numNodes + 1) * (long) sizeof(int));
        exit (1);
    }
    CumSamNodes[0] = 0;
    
    for (i = 1; i <= numClones; i++)
    {
        pop = *(populations + i-1);
        CumSamNodes[i] = 0;
        currentSampleSize = pop->sampleSize;
        //pop->idsActiveGametes= (int*) malloc( (numNodes) * sizeof( int));
        pop->numCompletedCoalescences=0;
        pop->nodeIdAncesterMRCA =0;
        //pop->idsActiveGametes= (int*) malloc( (TotalNumSequences) * sizeof( int));
        pop->idsActiveGametes=NULL;
        pop->idsActiveGametes= (int*) malloc( (pop->sampleSize + pop->numPossibleMigrations) * sizeof( int));
        
        if (!( pop->idsActiveGametes))
        {
            fprintf (stderr, "Could not allocate  pop->idsActiveGametes (%lu bytes)\n", (pop->sampleSize + pop->numPossibleMigrations ) * (long) sizeof(int));
            exit (1);
        }
        pop->idsGametes=NULL;
        pop->idsGametes= (int*) malloc( (2* pop->sampleSize + pop->numPossibleMigrations) * sizeof( int));
        if (!( pop->idsGametes))
        {
            fprintf (stderr, "Could not allocate  pop->idsGametes (%lu bytes)\n", (2* pop->sampleSize + pop->numPossibleMigrations) * (long) sizeof(int));
            exit (1);
        }
        pop ->numActiveGametes=0;
        pop ->numGametes=0;
        CumSamNodes[i] = CumSamNodes[i - 1] + currentSampleSize;
    }
    if (noisy > 1)
        fprintf (stderr, "\n Initial relation nodes-clones:");
    int  cumIndivid = 0;
    int currentPopIndex;
    *numActiveGametes = 0;
    
    for (j = 0; j < TotalNumSequences; j++)
    {
        cumIndivid++;
        p = nodes + j;
        //activeGametes[*numActiveGametes] = j;
        p->node_index = j;
        //p->label = j;
        
        *labelNodes = j;
        
         u= (TreeNode *)(p->data);
         u->nodeClass = 1;
        // fprintf (stderr,"\nIn cumIndivid=%d j=%d", cumIndivid, j);
        for (i = 1; i <= numClones; i++)
        {
            pop = *(populations + i - 1);
            currentPopIndex = pop->index;
            indexFirstObservedCellName= pop->indexFirstObservedCellName;
            // Identify to which clone belongs this node (sample)
            if (cumIndivid <= CumSamNodes[i] && cumIndivid > CumSamNodes[i - 1])
            {
                //fprintf (stderr,"\ncumIndivid=%d <= CumSamNodes[i]=%d, in clone %d\n", cumIndivid, CumSamNodes[i], pop->index);
                pop->idsActiveGametes[pop->numActiveGametes]=j;
                
                
                pop->idsGametes[pop->numGametes]=j;
                pop->numGametes=pop->numGametes+1;
                
                u->indexOldClone = currentPopIndex;
                u->indexCurrentClone = currentPopIndex;
                u->effectPopSize= pop->effectPopSize;
                u->orderCurrentClone = pop->order;
            
                
                if (programOptions->doUseObservedCellNames)
                    strcpy( u->observedCellName,ObservedCellNames[indexFirstObservedCellName + pop->numActiveGametes ]);
                pop->numActiveGametes=pop->numActiveGametes+1;
                
                break;
            }
        }
        
        //        if(doUseObservedCellNames == YES)
        //            strcpy( p->observedCellName,ObservedCellNames[j]);
        
        
        if (noisy > 1)
            fprintf (stderr,"\n > The node %d(%d) belongs to clone %d", u->index, cumIndivid, u->indexOldClone);
        *numActiveGametes = *numActiveGametes + 1;
    }
    //AssignObservedCellNamestoTips(nodes, populations, sampleSizes,ObservedCellNames,  programOptions);
    // AssignObservedCellNamestoTips2(nodes, populations, sampleSizes,ObservedCellNames,  programOptions);
    free(CumSamNodes);
    CumSamNodes=NULL;
    *nextAvailable = *numActiveGametes;
    *labelNodes = *labelNodes + 1;
}
void chain::SetPopulationsBirthRate( double lambda){
    int i;
    population *popI;
    for( i = 0 ; i < numClones; i++)
    {
        popI=*(populations + i);
        if (popI != NULL)
            popI->birthRate = lambda;
    }
}
void chain::GenerateEffectPopSizesFromPriors2( int noisy,   long int *seed,  int doGenerateProportionsVector){
    int i, j, k, l;
    population *popI, *popJ, *popL, *popIPlus1;
    double W, V,P, aboveTerm, belowTerm, temp;
    double sum;
    double rand;
    
    if (noisy > 1){
        fprintf (stderr, "Estimation of pop sizes of clones ..\n");
        fflush (stdout);
    }
    double *outputVector;
    if (doGenerateProportionsVector == YES)
    {
        for (j = 0; j < numClones; j++)
        {
            oldproportionsVector[j]=proportionsVector[j] ;
        }
        
        RandomDirichlet(gammaParam,  numClones, &(proportionsVector), seed );
        //chain->proportionsVector = outputVector;
        
    }
   
    if (proportionsVector == NULL)
        fprintf (stderr, "ERROR: the proportions vector is null..\n");
    for (j = 0; j < numClones; j++)
    {
        popJ = *(populations + j);
        //temp=*(outputVector + i);
        popJ->oldeffectPopSize = popJ->effectPopSize;
        popJ->effectPopSize = proportionsVector[j] * totalPopSize;
    
    }
}
void chain::FillChainPopulationsFromPriors( ProgramOptions *programOptions,  MCMCoptions *mcmcOptions, int *sampleSizes, long int *seed)
{
    population* popI;
    int i;
    if ( populations !=NULL)
    {
        int i;
        population *popI;
        double randomDelta;
        mutationRate = RandomLogUniform(mcmcOptions->MutRatefrom, mcmcOptions->MutRateto, seed);
        
        totalPopSize= RandomLogUniform(mcmcOptions->totalEffectPopSizefrom, mcmcOptions->totalEffectPopSizeto, seed);
        
        double lambda = 1;
        
        SetPopulationsBirthRate(  lambda);
     
        
        GenerateEffectPopSizesFromPriors2( programOptions->noisy,   seed, YES);
        if (programOptions->populationSampleSizesKnown == NO)
        {
            InitPopulationSampleSizes(populations, programOptions->numCells, programOptions->numClones,proportionsVector, seed);
        }
        //else fill the sample  sizes
        else{
            setChainPopulationSampleSizes( sampleSizes, programOptions);
        }
        //generate the population  sizes
        
       for( i = 0 ; i < numClones; i++)
       {
           popI=*(populations + i);
           popI->InitCoalescentEvents(numClones);
       }
        if (programOptions->populationSampleSizesKnown ==YES)
        {
            for( i = 0 ; i < numClones; i++)
            {
                popI=*(populations + i);
                do {
                    randomDelta = RandomLogUniform(mcmcOptions->Deltafrom, mcmcOptions->Deltato, seed);
                    //popI->delta = chain->proportionsVector[i] * randomDelta;
                    popI->delta =  randomDelta;
                    popI->growthRate =popI->delta  / popI->effectPopSize;
                    popI->popSize=popI->effectPopSize * popI->birthRate;
                    popI->deathRate= popI->birthRate - popI->growthRate;
                }
                while(popI->popSize < popI->sampleSize);
            }
        }
        ListClonesAccordingTimeToOrigin2();
        GenerateTimesFromPriorsOriginal(programOptions->noisy,  seed);
        ListClonesAccordingTimeToOrigin2();
    
        for( i = 0 ; i < numClones; i++)
        {
            popI=*(populations + i);
            popI->InitListPossibleMigrations();
        }
      
    }
}
void chain::setChainPopulationSampleSizes(int *sampleSizes,  ProgramOptions *programOptions)
{
    int i;
    population *popI;
    int cumSampleSize=0;
    
    if (sampleSizes== NULL)
        fprintf (stderr, "ERROR: the sample sizes  vector is null..\n");
    if (numClones >0 && populations!=NULL)
    {
        for( i = 0 ; i < numClones; i++)
        {
            popI=*(populations + i);
            popI->sampleSize=sampleSizes[i];
            popI->indexFirstObservedCellName =cumSampleSize;
            cumSampleSize +=sampleSizes[i];
        }
    }
}
void chain::ListClonesAccordingTimeToOrigin2() {
    
    
    qsort(populations, numClones, sizeof(Population* ), comparePopulationsByTimeOrigin);
}
void chain::GenerateTimesFromPriorsOriginal(int noisy,  long int *seed) {
    double    *Uvector;
    int i, j, z, m,  l;
    double    TotalProbability, LocalProbability, aboveTerm, belowTerm, ranHere;
    int  AttemptsAcceptation = 0;
    TotalProbability = 0.0;
    LocalProbability = 0.0;
    aboveTerm = 0.0;
    belowTerm = 0.0;
    ranHere = 0.0;
    if (noisy > 1)
        printf("Estimation of times of origin of clones ..\n");
    Uvector = (double *) malloc((numClones )* (long) sizeof(double));
    if (!Uvector)
    {
        fprintf (stderr, "Could not allocate Xvector (%lu bytes)\n", (numClones ) * (long) sizeof(double));
        exit (-1);
    }
    int doAcceptTimes = NO;
    int doEstimateTimesOriginClones=YES;
    population *popJ, *popI, *popL;
    double rand;
    if (doEstimateTimesOriginClones == YES)
    {
        while (doAcceptTimes == NO)
        {
            AttemptsAcceptation++;
            // Calculate t and T
            if (noisy > 2)
                printf("\nProposed times of origin of clones ..");
            for (j = 0; j < numClones; j++)
            {
                popJ=*(populations + j);
                Uvector[j] = 0;
                Uvector[j] = RandomUniform (seed);
                rand = RandomUniform(seed);
                popJ->timeOriginSTD= (1/ popJ->delta)*log(1-(popJ->delta / log(rand)));
                popJ->timeMigrationSTDCurrentPop= popJ->timeOriginSTD ;
                // popJ->timeOriginInput=proportionsVector[j]*popJ->timeOriginSTD;
                popJ->timeOriginInput=popJ->effectPopSize*popJ->timeOriginSTD;
            }
            
            if (numClones == 1) {
                doAcceptTimes = YES;
                break;
            }
            ListClonesAccordingTimeToOrigin2();
            // Calculate probabilities P
            LocalProbability = 0.0;
            TotalProbability = 0.0;
            // printf("Number of clones = %d\n", numClones);
            for (i = 0; i < numClones - 1; i++)
            {
                popI=*(populations + i);
                printf("\ni  = %d ", i);
                LocalProbability = 0.0;
                aboveTerm = 0.0;
                belowTerm = 0.0;
                m = popI->order;
                // printf("\n\nCalculating P for clone = %d ", m);
                for (l = i + 1; l < numClones; l++)
                {
                    popL=*(populations + l);
                    j = popL->order ;
                    aboveTerm = aboveTerm + (popL->popSize * population::CalculateH(popI->timeOriginSTD * popI->effectPopSize / popL->effectPopSize, popL->timeOriginSTD, popL->delta));
                    belowTerm = belowTerm + popL->popSize;
                    // printf("\nClones %d and %d ", m, j);
                }
                LocalProbability = aboveTerm / belowTerm;
                //printf("\naboveTerm = %lf    belowTerm = %lf", aboveTerm, belowTerm);
                if (i == 1)
                    TotalProbability = 1.0 * LocalProbability;
                else
                    TotalProbability = TotalProbability * LocalProbability;
                // printf("\nLocalProbability = %lf    TotalProbability = %lf", LocalProbability, TotalProbability);
            }
            // printf("\nTotalProbability = %lf \n", TotalProbability);
            // Random number
            ranHere = RandomUniform (seed);
            // printf("\n\nranHere = %lf \n", ranHere);
            // free memory
            //free (CloneNameBeginOrderByTimes);
            // free (CloneNameBeginOrderByTimesOnlyControl);
            // check accept or reject attempt
            if (ranHere <= TotalProbability)
                doAcceptTimes = YES;
            if (noisy > 2)
                printf("\nProbability = %lf (random number = %lf) [#attempt = %d]\n", TotalProbability, ranHere, AttemptsAcceptation);
        }
    }
    if (noisy > 2)
    {
        printf("\nTimes accepted .. \n");
        printf("\n\nDONE! ranHere = %lf / TotalProbability = %lf [total attempts = %d]\n", ranHere, TotalProbability, AttemptsAcceptation);
    }
}
void chain::InitChainPopulations( int noisy,  int TotalNumSequences  ) {
    int z;
    //struct Population* pops = malloc(numClones * (sizeof(struct Population)+TotalNumSequences * sizeof( int) + numClones * sizeof(double) ));
    population* pops =(population*) malloc(numClones * (sizeof(struct population) ));
    for (z = 0; z <= (numClones - 1); z++)
    {
        pops[z].index =z;
        pops[z].order=0;
        pops[z].birthRate =0.0;
        pops[z].deathRate = 0;
        pops[z].growthRate = 0;
        pops[z].sampleSize = 0;
        pops[z].popSize =0;
        pops[z].effectPopSize =0;
        pops[z].delta =0 ;
        pops[z].numActiveGametes = 0;
        pops[z].isAlive = 1;
        pops[z].nodeIdAncesterMRCA = 0;
        pops[z].numCompletedCoalescences = 0;
        pops[z].nextAvailableIdInmigrant = 0;
        pops[z].numIncomingMigrations = 0;
        pops[z].numPossibleMigrations = 0;
        pops[z].doEstimateTimeOrigin = NO;
        pops[z].timeOriginInput = 0;
        pops[z].timeOriginSTD =0;
        pops[z].timeMigrationSTDCurrentPop = 0;
        pops[z].FatherPop =NULL;
        
        *(populations + z) = &pops[z];
        
    }
}
