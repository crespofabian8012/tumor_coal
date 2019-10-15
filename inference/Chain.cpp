//
//  Chain.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 7/10/19.
//

#include "chain.hpp"
using namespace std;

#include "data_utils.hpp"
#include "definitions.hpp"
#include "tree_node.hpp"
#include "random.h"
#include "constants.hpp"

#include <libpll/pll_msa.h> //for pllmod_msa_empirical_frequencies
#include "eigen.hpp"
//#include "libpll/pllmod_common.h"



//Parametrized Constructor

Chain::Chain( int chainNumber,
             int numClones,
             int gammaParam,
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
    

    
    currentlogConditionalLikelihoodTree=0;
    currentlogConditionalLikelihoodSequences=0;
    
    
}
int Chain::setIntitialTree(char * NewickString){
    
    this->initialTree = pll_utree_parse_newick_string_unroot(NewickString);
    if (this->initialTree ==NULL)
        return 1;
    else
        return 0;
    
}
//void Chain::MakeCoalescenceEvent(Population *population, pll_unode_t **nodes, int numClones, long int* seed, int noisy,  int *numActiveGametes,   int* nextAvailable,
//                                 int*labelNodes, double *currentTime, int *numNodes)
void Chain::MakeCoalescenceEvent( Population *population, vector<pll_unode_t *> &nodes, int numClones, long int* seed, int noisy,   int &numActiveGametes, int &nextAvailable,
                          int &labelNodes, double &currentTime, int &numNodes)
{
    int k, w;
    double rand, ran;
    // TreeNode  *p, *q, *r;
    pll_unode_t  *p, *q, *r, *r1, *r2, *r3 ;
    TreeNode *u1, *u2, *u3;
    int firstInd, i, j, secondInd=0, newInd=0;
    int choosePairIndividuals = YES;
    
    population->ChooseRandomIndividual(&firstInd, numClones,   &secondInd, seed, choosePairIndividuals);
    
    newInd = nextAvailable;
    if (noisy > 1)
        fprintf (stderr, "Coalescence involving %d and %d to create node %d (in clone %d)", population->idsActiveGametes[firstInd], population->idsActiveGametes[secondInd], newInd, population->index);
    /*  set pointers between nodes */
    p = nodes[population->idsActiveGametes[firstInd]];
    q = nodes[population->idsActiveGametes[secondInd]];
    
    
    r1 = nodes[newInd];   /* new nodelet ancester */
    r2 = nodes[newInd + 1];   /* new  nodelet ancester */
    r3 = nodes[newInd + 2];   /* new nodelet ancester */
    
    r1->node_index = newInd;
    r2->node_index = newInd +1;
    r3->node_index = newInd +2;
    
//    r1->data =  malloc(sizeof(TreeNode));
//    r2->data =  malloc(sizeof(TreeNode));
//    r3->data =  malloc(sizeof(TreeNode));
    
    r1->data =  new TreeNode();;
    r2->data =  new TreeNode();
    r3->data = new TreeNode();;
    
    u1= (TreeNode *)(r1->data);
    u2= (TreeNode *)(r2->data);
    u3= (TreeNode *)(r3->data);
    
    u1->index =  u2->index=  u3->index= nextAvailable;
    u1->label = u2->label =u3->label =labelNodes;
    labelNodes=labelNodes+1;
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
    u1->time = u2->time = u3->time =currentTime;
    u1->timePUnits = u2->timePUnits = u3->timePUnits =currentTime * (population->effectPopSize);
    
    //fprintf (stderr, "\n r->index = %d, r->time = %lf, ClonePopSizeMeffectBegin[ThisCloneNumber] = %lf, ThisCloneNumber = %d\n", r->index, r->time, ClonePopSizeMeffectBegin[ThisCloneNumber], ThisCloneNumber);
    if (noisy > 1)
        fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", u1->timePUnits);
    /* readjust active nodes */
    
    // idsActiveGametes[firstInd] = newInd;
    population->idsActiveGametes[firstInd] = r3->node_index;// the r3 node
    population->idsActiveGametes[secondInd] = population->idsActiveGametes[population->numActiveGametes - 1];
    numActiveGametes = numActiveGametes - 1; /* less 1 active node */
    
    //update list ids nodes
    population->idsGametes[population->numGametes] = newInd;
    population->numGametes = population->numGametes +1;
    
    // *nextAvailable=*nextAvailable+1; /* 1 node more is available */
    nextAvailable=nextAvailable+3; /* 1 node more is available */
    
    population->CoalescentEventTimes[ population->numCompletedCoalescences]=  u1->time;
    numActiveGametes = numActiveGametes - 1; /* now this clone
                                              has 1 less node */
    
    population->numCompletedCoalescences= population->numCompletedCoalescences+1;
    /* memory for number of nodes */
    if (nextAvailable >= numNodes)  /* if there aren't enough nodes it go into and it addition more */
    {
        /* ReallocNodes(&numNodes, activeGametes); */
        if (noisy == 4)
            fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
        numNodes += INCREMENT_NODES;
        /* realloc */
        //*nodes = (TreeNode *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode) );
//        *nodes = (pll_unode_t *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode) );
//        if (!(*nodes))
//        {
//            fprintf (stderr, "Could not reallocate nodes (%lu bytes)\n", *numNodes  * (long) sizeof(TreeNode));
//            exit (-1);
//        }
//        population->idsActiveGametes = (int *) realloc (population->idsActiveGametes, *numNodes * (long) sizeof(int));
//        if (!(population->idsActiveGametes))
//        {
//            fprintf (stderr, "Could not reallocate idsActiveGametes for the current population(%lu bytes)\n", *numNodes * (long) sizeof(int));
//            exit (-1);
//        }
    }
    
}

//void Chain::SimulatePopulation(Population *pop,ProgramOptions *programOptions,
//                               long int *seed,
//                               int *numNodes,
//                               int numClones,
//                               double      cumNumCA,
//                               double meanNumCA,
//                               double cumNumMIG,
//                               double meanNumMIG,
//                               int  *numMIG,
//                               int  *numCA,
//                               double *numEventsTot,
//                               pll_unode_t    **nodes,
//                               int *nextAvailable,
//                               int*  numActiveGametes,
//                               int* labelNodes,
//                               double *currentTime,
//                               int* eventNum
//                               )
void Chain::SimulatePopulation( Population *popI, vector<Population*> &populations,
                        ProgramOptions &programOptions,
                        long int *seed,
                        int &numNodes,
                        int numClones,
                        double      cumNumCA,
                        double meanNumCA,
                        double cumNumMIG,
                        double meanNumMIG,
                        int  &numMIG,
                        int  &numCA,
                        double &numEventsTot,
                        vector<pll_unode_t *> &nodes,
                        int &nextAvailable,
                        int &numActiveGametes,
                        int &labelNodes,
                        double &currentTime,
                        int &eventNum)
{
    int  c, d, i, j, w, k, m, cumIndivid, isCoalescence, whichInd,
    firstInd, secondInd, newInd,  foundSuperflousNode,
    isMigration, whichClone, currentNumberAliveClones;
    double     eventTime;
    pll_unode_t  *p, *q, *r, *r1, *r2, *r3 ;
    TreeNode  *u,*v, *u1, *u2, *u3;
    double    ran;
    double    *cumPopulPart;
    eventNum = 0;
    int numSimClones = 0;
    double ThisRateCA = 0.0;
    double ThisTimeCA_W = 0.0;
    double ThisTimeCA_V1 = 0.0;
    double ThisTimeCA_V2 = 0.0;
    int numParcialActiveGametes;
    int choosePairIndividuals;
    
    numParcialActiveGametes = popI->numActiveGametes;
    
    int numMigrations = popI->numIncomingMigrations; //taking into account also the    time of origin as a migration
    
    double timeNextMigration;
    int indexNextMigration = 0;
    Population *incomingPop;
    
    if (programOptions.noisy > 1)
        fprintf (stderr, "\n\n>> Simulating evolutionary history of clone %d (number active gametes %d, original time to origin %lf)\n", popI->index, popI->numActiveGametes, popI->timeOriginInput);
    
    if (programOptions.noisy > 1)
        fprintf (stderr, "\n\n> Simulating evolutionary history of clone  or order  %d ..\n", popI->order);
    
    currentTime=0;
    while (indexNextMigration < numMigrations) {
        timeNextMigration = popI->immigrantsPopOrderedByModelTime[indexNextMigration].first;
        //fprintf (stderr, "\n\n> numParcialActiveGametes= %d \n", numParcialActiveGametes);
        if ( popI->numActiveGametes >= 2) {
            ThisRateCA = (double)  popI->numActiveGametes * ((double) popI->numActiveGametes - 1) / 2.0;
            ThisTimeCA_W = RandomExponential (ThisRateCA, seed) ;
            ThisTimeCA_V1 = Population::FmodelTstandard (currentTime, popI->timeOriginSTD, popI->delta);
            ThisTimeCA_V1 = ThisTimeCA_V1 + ThisTimeCA_W;
            // from standard time to model time, GstandardTmodel(V, T, delta)
            ThisTimeCA_V2 = Population::GstandardTmodel(ThisTimeCA_V1, popI->timeOriginSTD, popI->delta);
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
            numCA = numCA + 1;
            eventNum= eventNum +1;
            whichClone = popI->index;
            currentTime = ThisTimeCA_V2; // update current time in model time
            eventTime = currentTime;
            
            if (programOptions.noisy > 1)
            {
                fprintf (stderr, "\n\n*** Event %3d *** currentTime (model time) = %lf, currentTime (standard time) = %lf\n", eventNum, ThisTimeCA_V2, ThisTimeCA_V1 );
                fprintf (stderr, "\n\n*** Event %3d *** currentTime (input units) = %lf\n", eventNum, ThisTimeCA_V2);
            }
            if (programOptions.noisy == 4)
                fprintf (stderr, "* Coalescence *\n");
             MakeCoalescenceEvent( popI, nodes, numClones, seed, programOptions.noisy, numActiveGametes, nextAvailable, labelNodes, currentTime,  programOptions.numNodes);
//            MakeCoalescenceEvent( pop, nodes, numClones, seed, programOptions->noisy, numActiveGametes,  nextAvailable,
//                                 labelNodes, currentTime,  &(programOptions->numNodes));
            
            
        }
        else
        {
            if (indexNextMigration < numMigrations - 1) //indexNextMigration corresponds to one of the true migrations
            {
                isCoalescence = NO;
                isMigration = YES;
                numMIG = numMIG + 1;
                eventNum= eventNum +1;
                currentTime = timeNextMigration;
                eventTime = currentTime;
                //update migration times in model time
                if (programOptions.noisy > 1)
                {
                    fprintf (stderr, "\n\n*** Event %3d *** *currentTime (model units) = %lf\n", eventNum, ThisTimeCA_V2);
                }
                if (programOptions.noisy == 4)
                {fprintf (stderr, "* Migration *\n");}
//                newInd = nextAvailable;
//
//
//                r1 = nodes[newInd];   /* new nodelet ancester */
//                r2 = nodes[newInd + 1];   /* new  nodelet ancester */
//                r3 = nodes[newInd + 2];   /* new nodelet ancester */
//
//
//                r1->node_index = newInd;
//                r2->node_index = newInd+1;
//                r3->node_index = newInd+2;
//
//                r1->data =  new TreeNode();
//                r2->data =  new TreeNode();
//                r3->data =  new TreeNode();
//
//                u1= (TreeNode *)(r1->data);
//                u2= (TreeNode *)(r2->data);
//                u3= (TreeNode *)(r3->data);
//
//                u1->label = u2->label= u3->label= *labelNodes;
//                *labelNodes=*labelNodes+1;
//
//                u1->indexCurrentClone = u2->indexCurrentClone= u3->indexCurrentClone= pop->index;
//                u1->indexCurrentClone = u2->indexCurrentClone = u3->indexCurrentClone = pop->index;
//                u1->orderCurrentClone =u2->orderCurrentClone= u3->orderCurrentClone= pop->order;
//                u1->nodeClass = u2->nodeClass= u3->nodeClass= 4;
//
                
                //    p = *nodes + MatrixMigrationIDnodeMRCA[ThisCloneNumber][doAmigration]; // root of younger clone
                incomingPop = popI->immigrantsPopOrderedByModelTime[indexNextMigration].second;
                
               
                //r->indexOldClone = incommingPop->index;
                p = nodes[incomingPop->nodeIdAncestorMRCA]; // root of younger clone
                
                
                indexNextMigration = indexNextMigration + 1;
                
                v= (TreeNode *)p;
                
                printf( "\n The incoming population %d to  population %d with node %d and nodelet %d time %lf", incomingPop->order, popI->order, v->index, p->node_index, v->timePUnits );
                
                v->indexCurrentClone = popI->index;
                v->indexOldClone = incomingPop->index;
                v->orderCurrentClone = popI->order;
                // link the nodes
                //r->left = p;
//                r1->back = p;
//                r1->next = r2;
//                r2->next = r3;
//                r3->next = r1;
//
//                //r->right = NULL;
//                r2->back = NULL;
//                r3->back = NULL;
                //choosePairIndividuals = NO;
                
                //ChooseRandomIndividual(&firstInd, numClones, popI,  &secondInd, seed, choosePairIndividuals);
                //q=*nodes + firstInd;
                //r->right = q;//choose another random living individual of the population
                
               // p->back = r1;
                
                //p->anc1 = r;
                //q->anc1 = r;
                
                
                //p->time = *currentTime;
                // p->timePUnits = *currentTime * (popI->effectPopSize);
                
                
              //  u1->time =u2->time= u3->time= currentTime;// this is equal to the time of the migration
              //  u1->timePUnits = u2->timePUnits= u3->timePUnits= currentTime * (pop->effectPopSize);
                //nextAvailable= nextAvailable+3; /* 3 nodelets more are available */
                
                k = v->indexCurrentClone;
                incomingPop->numActiveGametes = incomingPop->numActiveGametes - 1; /* now the other clone has 1 less node */
                // remove node from old clone in list of active gametes and add the new node of the current clone
                //popI->idsActiveGametes[popI->numActiveGametes]=r->index;//adding the superfluos node
                popI->idsActiveGametes[popI->numActiveGametes]=p->node_index;//adding the superfluos node
                popI->numActiveGametes = popI->numActiveGametes + 1; /* now this clone has 1 more node */
                //                if (noisy > 1)
                //                    fprintf (stderr, "Migration, creating node %d (clone %d) derived from node %d (clone %d)", newInd, popI->index, incommingPop->nodeIdAncesterMRCA, k);
                //fprintf (stderr, "Migration, creating node %d (clone %d) derived from node %d (clone %d)", newInd, ThisCloneNumber, MatrixMigrationIDnodeMRCA[ThisCloneNumber][doAmigration], k);
                
                for(int i=0; i < popI->idsActiveGametes.size();i++)
                    fprintf (stderr, "\n pop of order  %d Active gamete id %d", popI->order, popI->idsActiveGametes[i]);
                
                if (programOptions.noisy > 1)
                    fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", v->timePUnits);
                // fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", r->timePUnits);
                /* memory for number of nodes */
                if (nextAvailable >= numNodes)  /* if there aren't enough nodes it go into and it addition more */
                {
                    /* ReallocNodes(&numNodes, activeGametes); */
                    if (programOptions.noisy == 4)
                        fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
                    numNodes += INCREMENT_NODES;
                    /* realloc */
                    // *nodes = (TreeNode *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode));
//                    nodes = (pll_unode_t *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode));
//                    if (!(*nodes))
//                    {
//                        fprintf (stderr, "Could not reallocate nodes (%lu bytes)\n", *numNodes  * (long) sizeof(pll_unode_t));
//                        exit (-1);
//                    }
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
                currentTime = timeNextMigration;
                indexNextMigration = indexNextMigration + 1;
                isCoalescence = NO;
                isMigration = NO;
                eventTime = currentTime;
                if (programOptions.noisy > 1)
                    fprintf (stderr, "\n\n*** Clone origin ***\n");
                if (programOptions.noisy == 4)
                    fprintf (stderr, "Clone origin %d at time (model units) = %lf\n", popI->index, currentTime);
                if (popI->order < numClones - 1) // do not do it for the last clone
                {
                    newInd = nextAvailable;
                    
                    r1 = nodes[newInd];   /* new nodelet ancester */
                    r1 = nodes[newInd + 1];   /* new  nodelet ancester */
                    r1 = nodes[newInd + 2];   /* new nodelet ancester */
                    
                    r1->node_index = newInd;
                    r2->node_index = newInd +1;
                    r3->node_index = newInd +2;
                    
                    r1->data =  new TreeNode();
                    r2->data =   new TreeNode();
                    r3->data =   new TreeNode();
                    
                    u1= (TreeNode *)(r1->data);
                    u2= (TreeNode *)(r2->data);
                    u3= (TreeNode *)(r3->data);
                    
                    u1->index =u2->index= u3->index= nextAvailable;
                    u1->label =u2->label= u3->label= labelNodes;
                    labelNodes=labelNodes+1;
                    u1->indexOldClone =u1->indexCurrentClone = popI->index;
                    u2->indexOldClone =u2->indexCurrentClone = popI->index;
                    u3->indexOldClone =u3->indexCurrentClone = popI->index;
                    
                    u1->orderCurrentClone = u2->orderCurrentClone= u3->orderCurrentClone= popI->order;
                    u1->effectPopSize= u2->effectPopSize=  u3->effectPopSize= popI->effectPopSize;
                    popI->nodeIdAncestorMRCA=newInd;
                    
                    u1->nodeClass = u2->nodeClass =u3->nodeClass =4;
                    
                    //firstInd = *nextAvailable - 1;
                    firstInd = nextAvailable - 1;
                    //p = nodes + activeGametes[firstInd]; // descendant node (previously generated node, nextAvailable - 1)
                    p = nodes[firstInd]; // descendant node (previously generated node, nextAvailable - 1)
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
                    
                    u1->time = u2->time= u3->time= currentTime;
                    u1->timePUnits =  u2->timePUnits =  u3->timePUnits = currentTime * popI->effectPopSize;
                    popI->nodeletMRCA = p;
                    
                    
                    
                    //fprintf (stderr, "\n r->index = %d, r->time = %lf\n", r->index, r->time);
                    /* readjust active nodes */
                    nextAvailable=nextAvailable+1; /* 1 node more is available */
                    popI->idsActiveGametes[0] = newInd;//always will be in the 0  position because there is only one left
                    
                    //popI->idsGametes[popI->numGametes] = newInd; r is a superflous node and it will be removed so no need to add it
                    //popI->numGametes = popI->numGametes +1;
                    
                    if (programOptions.noisy > 1)
                        fprintf (stderr, "Creating origin node, it creates node %d derived from node %d", newInd, firstInd);
                    if (programOptions.noisy > 1)
                        fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", u1->timePUnits);
                    /* memory for number of nodes */
                    if (nextAvailable >= numNodes)  /* if there aren't enough nodes it go into and it addition more */
                    {
                        /* ReallocNodes(&numNodes); */
                        if (programOptions.noisy == 4)
                            fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
                        numNodes += INCREMENT_NODES;
                        /* realloc */
                        // *nodes = (TreeNode *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode));
//                        *nodes = (pll_unode_t *) realloc (*nodes, *numNodes  * (long) sizeof(pll_unode_t));
//                        if (!(*nodes))
//                        {
//                            fprintf (stderr, "Could not reallocate nodes (%lu bytes)\n", *numNodes  * (long) sizeof(TreeNode));
//                            exit (-1);
//                        }
                        
                    }
                }
                else  {//origin of oldest pop reached
                    popI->nodeIdAncestorMRCA=nextAvailable-1;//for the last population, nodeIdAncesterMRCA is the MRCA instead of ancester of MRCA
                    r = nodes[nextAvailable-1];//popI->idsActiveGametes[0]
                    u= (TreeNode *)(r->data);
                    u->indexOldClone = u->indexCurrentClone = popI->index;
                    u->orderCurrentClone = popI->order;
                    popI->nodeletMRCA= r;
                    
                    
                }
            }
        }
        if (programOptions.noisy > 3)
        {
            fprintf (stderr, "\nActive nodes (%d):",  popI->numActiveGametes);
            for (i = 0; i < popI->numActiveGametes; i++)
                fprintf (stderr, " %d", popI->idsActiveGametes[i]);
            fprintf (stderr, "\t|\tNext node available = %d", nextAvailable);
        }
    }
    if (programOptions.noisy > 1)
        fprintf (stderr, "\n\nEvolutionary history of clone %d is completed \n", popI->index);
    
}

pll_unode_t* Chain::BuildTree(vector<Population* > &populations,
                       Population *CurrentPop,
                       long int *seed,
                       ProgramOptions &programOptions,
                       vector<pll_unode_t *> &nodes,
                       vector<pll_unode_t *> &treeTips,
                       pll_unode_t *tumour_mrca,
                       int &nextAvailable,
                       int &newInd,
                       double &currentTime,
                       int &labelNodes)
{
    int i, j, k;
    int indexCurrentTip;
    int  foundSuperflousNode;
    pll_unode_t *p, *q, *r;
    TreeNode *u, *anc, *healthyR, *u1, *u2, *u3;
    pll_unode_t *treeRootInit;
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
            p = nodes[i];
            
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
                if(indexCurrentTip <  programOptions.TotalNumSequences){
                    //treeTips[indexCurrentTip]=p;
                    treeTips.push_back(p->next->next);
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
                p->next= NULL;
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
    newInd=nextAvailable-1;
    p = nodes[newInd]; /* because the last one event is the last coalescence */
    //fprintf (stderr, "\n\n\n>> newInd = %d\n", newInd);
    
    if (programOptions.thereisOutgroup == NO)
    {
        p = nodes[newInd];
        u = (TreeNode *)(p->data);
        u->nodeClass = 5;
        tumour_mrca = p;
        //treeRootInit[0] = p;
        u->anc1 = NULL;
        treeRootInit =tumour_mrca;
    }
    if (programOptions.thereisOutgroup == YES && programOptions.outgroupSelection > 0)  /*** Root and outgroup ***/
    {
        p = nodes[newInd]; // MRCA
        u = (TreeNode *)(p->data);
        u->nodeClass = 4;
        
        if (programOptions.noisy > 1)
            fprintf (stderr, "\n\n>> Attaching outgroup .. ");
        
        //fprintf (stderr, "\n>> ThisCloneNumber = %d, ListMigrationTimesInitial[ThisCloneNumber] = %lf, ClonePopSizeMeffectBegin[ThisCloneNumber] = %lf \n", ThisCloneNumber, ListMigrationTimesInitial[ThisCloneNumber], ClonePopSizeMeffectBegin[ThisCloneNumber]);
        
        if (programOptions.outgroupSelection == 1)  /*** Root 2 times and outgroup ***/
            currentTime = CurrentPop->timeOriginSTD; // origin of the clone; // currentTime + (outgroupBranchLength_Root1Root2 / mutationRate); // set time of the new root (from which the MRCA and outgroup nodes are derived)
        else if (programOptions.outgroupSelection == 2) { /*** Root 2 times and outgroup ***/
            currentTime = CurrentPop->timeOriginSTD + (programOptions.outgroupBranchLength_Root1Root2 / CurrentPop->effectPopSize) ; // origin of the clone + time given by the user
            
        }
        else
        {
            fprintf (stderr, "\n\nError simulationg the outgroup. Check input settings\n");
            PrintUsage();
        }
        //TreeNode*       healthyRoot = *nodes + *nextAvailable;
        pll_unode_t*       healthyRoot = nodes[nextAvailable];
        u = (TreeNode *)(p->data);
        healthyR = (TreeNode *)(healthyRoot->data);
        healthyR->index = nextAvailable;
        healthyR->label = labelNodes;
        healthyR->effectPopSize= u->effectPopSize;
        labelNodes= labelNodes+1;
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
        healthyR->time = currentTime  ;
        
        int transformingBranchLength=1.001;
        // healthyRoot->time = p->time * transformingBranchLength ;
        healthyR->timePUnits = currentTime * healthyR->effectPopSize;
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
        nextAvailable++;
        //
        //        /* connect the healthy ancestral cell with the tip healthy cell*/
        //        if (noisy > 2)
        //            fprintf (stderr, "\n>> Adding healthy tip ... ");
        //TreeNode* healthyTip = *nodes + *nextAvailable;
        pll_unode_t* healthyTip1 = nodes[nextAvailable];
        pll_unode_t* healthyTip2 = nodes[nextAvailable
                                             +1];
        pll_unode_t* healthyTip3 = nodes[nextAvailable + 2];
        
        healthyTip1->back = healthyTip2->back = NULL;
        healthyTip1->next = healthyTip2 ;
        healthyTip2->next = healthyTip3;
        healthyTip3->next =healthyTip1;
        
        healthyTip1->data =  new TreeNode();
        healthyTip2->data =  new TreeNode();
        healthyTip3->data =   new TreeNode();
        
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
        
    
        
        treeRootInit=healthyRoot;
        
    }
    
    int intLabel = 0;
    if (programOptions.noisy > 1)
        fprintf (stderr, "\n\n>> Relabeling nodes on tree... \n\n");
    if (programOptions.thereisOutgroup == YES)
        intLabel = programOptions.TotalNumSequences + 1;
    else
        intLabel = programOptions.TotalNumSequences;
    
     RelabelNodes( treeRootInit, intLabel );
        
    return treeRootInit;
}
pll_unode_t * Chain::MakeCoalescenceTree (long int *seed,
                                 int &numNodes,
                                 int numClones,
                                 ProgramOptions &programOptions,
                                 double      cumNumCA,
                                 double meanNumCA,
                                 double cumNumMIG,
                                 double meanNumMIG,
                                 int  &numMIG,
                                 int  &numCA,
                                 double &numEventsTot,
                                 //                           pll_unode_t** nodes,
                                 //                           pll_unode_t** treeTips,
                                 //                           pll_unode_t    **treeRootInit,
                                 char* ObservedCellNames[],
                                 vector<int> &sampleSizes
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
    Population *pop;
    
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
        pop =populations[i];
        pop->resetMigrationsList();
    }
    //resetMigrationsList( populations,  numClones);
    
    //allocate memory for the nodes
//    nodes = (pll_unode_t *) malloc (3*(programOptions->numNodes + 1)* (sizeof(pll_unode_t) + 2*MAX_NAME * sizeof(char) + sizeof(TreeNode) )); /* nodes */
//    if (!(nodes))
//    {
//        fprintf (stderr, "Could not allocate nodes (%lu bytes)\n", 3*(programOptions->numNodes+ 1)  * (long) sizeof(pll_unode_t));
//        exit (1);
//    }
    
    for (i=0; i< programOptions.TotalNumSequences; i++){
        treeTips[i]=NULL;
    }
    
    for (i = 0; i < numNodes; i++)
    {
        //p = (*nodes + i);
        p = new pll_unode_t();
        nodes.push_back(p);
        
    }
   
    
    AssignSequencesToPopulations( programOptions, numNodes, programOptions.noisy, programOptions.TotalNumSequences, numActiveGametes,  nextAvailable,
                                 labelNodes, ObservedCellNames, programOptions.doUseObservedCellNames, sampleSizes);
    Population *currentPop;
    Population *fatherPop;
    i=0;
    currentTime=0.0;
    while (i < numClones) {
        //currentPop = *(populations + i);
        currentPop = populations[i];
        SimulatePopulation( currentPop, populations, programOptions, seed,
                          programOptions.numNodes,
                           numClones,
                           cumNumCA,
                           meanNumCA,
                           cumNumMIG,
                           meanNumMIG,
                           numMIG,
                           numCA,
                           numEventsTot,
                           nodes,
                           nextAvailable ,
                           numActiveGametes,
                           labelNodes,
                           currentTime,
                           eventNum);
        if (i< numClones-1)   //if it is not the last one
        {
            //choose the father population from which the population i came
            fatherPop= ChooseFatherPopulation( numClones, currentPop, seed,  programOptions.noisy);
            currentPop->FatherPop = fatherPop;
            //update list of migrant times
            Population::UpdateListMigrants( numClones, currentPop, fatherPop);
        }
        i = i + 1;
    }
 
    pll_unode_t *root =  BuildTree(populations, currentPop,seed,
              programOptions,
              nodes,
              treeTips,
              root,
              nextAvailable,
              newInd,
              currentTime,
              labelNodes
              );
    
    if (programOptions.noisy > 1)
        fprintf (stderr, "\n\n>> Relabeling nodes on tree... \n\n");
    
    
    // RelabelNodes(treeRootInit[0], treeRootInit, &intLabel );
    return root;
}
Population* Chain::ChooseFatherPopulation( int numClones, Population  *PopChild,  long int *seed, int noisy) {
    
    Population *pOrigin, *pTarget, *p;
    double pij, ran;
    double *ptr;
    int i, j, k;
    double cumProb[numClones  - (int)(PopChild->order)];
    cumProb[0] = 0.0;
    for (j = PopChild->order + 1; j < numClones; j++)
    {
        cumProb[j - PopChild->order] = 0.0;
        p = populations[ j];
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
    Population *result =  populations[ PopChild->order + w];
    if (noisy > 3)
        fprintf (stderr, "\nClone %d derived from clone %d\n", PopChild->index, result->index); // clone ThisCloneNumber (i) is originated from clone ThisOriginCloneNumber (j)
    /* Update list of migation times considering that clone i comes from clone j */
    if (noisy > 1)
        fprintf (stderr, "\n*** Updating list of migration times (considering that clone %d comes from clone %d) ..\n", PopChild->index, result->index);
    return (result); //the father population has  order  (PopChild->order) + w
}
void Chain::AssignSequencesToPopulations(
                                        ProgramOptions &programOptions,
                                        int numNodes, int noisy,  int TotalNumSequences,
                                        int &numActiveGametes, int &nextAvailable,
                                        int &labelNodes, char* ObservedCellNames[], int doUseObservedCellNames, vector<int> &sampleSizes)

{
    Population *pop;
    pll_unode_t *p;
    TreeNode *u;
    int i, j;
    double currentSampleSize;
    numActiveGametes = 0;
    int indexFirstObservedCellName;
    
  vector<int> CumSumNodes(numClones+1);
    CumSumNodes[0] = 0;
    
    for (i = 1; i <= numClones; i++)
    {
        pop = populations[i-1];
        pop->resetActiveGametes();
        currentSampleSize = pop->sampleSize;
        CumSumNodes[i] = CumSumNodes[i - 1] + currentSampleSize;
    }
    if (noisy > 1)
        fprintf (stderr, "\n Initial relation nodes-clones:");
    int  cumIndivid = 0;
    int currentPopIndex;
    numActiveGametes = 0;
    
    for (j = 0; j < TotalNumSequences; j++)
    {
        cumIndivid++;
          p = nodes[j];
        //activeGametes[*numActiveGametes] = j;
        p->node_index = j;
        //p->label = j;
        
        labelNodes = j;
        
        u= (TreeNode *)(p->data);
        u->nodeClass = 1;
        // fprintf (stderr,"\nIn cumIndivid=%d j=%d", cumIndivid, j);
        for (i = 1; i <= numClones; i++)
        {
            pop = populations[ i - 1];
            currentPopIndex = pop->index;
            indexFirstObservedCellName= pop->indexFirstObservedCellName;
            // Identify to which clone belongs this node (sample)
            if (cumIndivid <= CumSumNodes[i] && cumIndivid > CumSumNodes[i - 1])
            {
                //fprintf (stderr,"\ncumIndivid=%d <= CumSamNodes[i]=%d, in clone %d\n", cumIndivid, CumSamNodes[i], pop->index);
                pop->idsActiveGametes.push_back(j);
                
                
                pop->idsGametes.push_back(j);
                pop->numGametes=pop->numGametes+1;
                
                u->indexOldClone = currentPopIndex;
                u->indexCurrentClone = currentPopIndex;
                u->effectPopSize= pop->effectPopSize;
                u->orderCurrentClone = pop->order;
                
                
                if (programOptions.doUseObservedCellNames)
                    strcpy( u->observedCellName,ObservedCellNames[indexFirstObservedCellName + pop->numActiveGametes ]);
                pop->numActiveGametes=pop->numActiveGametes+1;
                
                break;
            }
        }
        
        //        if(doUseObservedCellNames == YES)
        //            strcpy( p->observedCellName,ObservedCellNames[j]);
        
        
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
void Chain::SetPopulationsBirthRate( double lambda){
    int i;
    Population *popI;
    for( i = 0 ; i < numClones; i++)
    {
        popI=populations[ i];
        if (popI != NULL)
            popI->birthRate = lambda;
    }
}
void Chain::GenerateEffectPopSizesFromPriors2( int noisy,   long int *seed,  int doGenerateProportionsVector){
    int i, j, k, l;
    Population *popI, *popJ, *popL, *popIPlus1;
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
        
        RandomDirichlet(gammaParam,  numClones, proportionsVector, seed );
        //chain->proportionsVector = outputVector;
        
    }
    
//    if (proportionsVector == NULL)
//        fprintf (stderr, "ERROR: the proportions vector is null..\n");
    for (j = 0; j < numClones; j++)
    {
        popJ = populations [j];
        //temp=*(outputVector + i);
        popJ->oldeffectPopSize = popJ->effectPopSize;
        popJ->effectPopSize = proportionsVector[j] * totalPopSize;
        
    }
}
void Chain::FillChainPopulationsFromPriors( ProgramOptions &programOptions,  MCMCoptions &mcmcOptions, vector<int> &sampleSizes, long int *seed)
{
    Population* popI;
    int i;
  
  
        double randomDelta;
        mutationRate = RandomLogUniform(mcmcOptions.MutRatefrom, mcmcOptions.MutRateto, seed);
        
        totalPopSize= RandomLogUniform(mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto, seed);
        
        double lambda = 1;
        
        SetPopulationsBirthRate(  lambda);
        
        
        GenerateEffectPopSizesFromPriors2( programOptions.noisy,   seed, YES);
        if (programOptions.populationSampleSizesKnown == NO)
        {
            InitPopulationSampleSizes(populations, programOptions.numCells, programOptions.numClones,proportionsVector, seed);
        }
        //else fill the sample  sizes
        else{
            setChainPopulationSampleSizes( sampleSizes, programOptions);
        }
        //generate the population  sizes
        
        for( i = 0 ; i < numClones; i++)
        {
            popI=populations[i];
            popI->InitCoalescentEvents(numClones);
        }
        if (programOptions.populationSampleSizesKnown ==YES)
        {
            for( i = 0 ; i < numClones; i++)
            {
                popI=populations[i];
                do {
                    randomDelta = RandomLogUniform(mcmcOptions.Deltafrom, mcmcOptions.Deltato, seed);
                    //popI->delta = chain->proportionsVector[i] * randomDelta;
                    popI->delta =  randomDelta;
                    popI->growthRate =popI->delta  / popI->effectPopSize;
                    popI->popSize=popI->effectPopSize * popI->birthRate;
                    popI->deathRate= popI->birthRate - popI->growthRate;
                }
                while(popI->popSize < popI->sampleSize);
            }
        }
        ListClonesAccordingTimeToOrigin(populations);
        GenerateTimesFromPriorsOriginal(programOptions.noisy,  seed);
        ListClonesAccordingTimeToOrigin(populations);
        
        for( i = 0 ; i < numClones; i++)
        {
            popI=populations[  i];
            popI->InitListPossibleMigrations(i);
        }
        
    
}
void Chain::setChainPopulationSampleSizes(vector<int > &sampleSizes,  ProgramOptions &programOptions)
{
    int i;
    Population *popI;
    int cumSampleSize=0;
    
    if (sampleSizes.size() == numClones)
        fprintf (stderr, "ERROR: the sample sizes  vector is null..\n");
   
        for( i = 0 ; i < numClones; i++)
        {
            popI=populations[i];
            popI->sampleSize = sampleSizes[i];
            popI->indexFirstObservedCellName =cumSampleSize;
            cumSampleSize +=sampleSizes[i];
        }
   
}
void Chain::ListClonesAccordingTimeToOrigin(vector<Population *> &populations)
{
    sort(populations.begin(), populations.end(), comparePopulationsByTimeOrigin);
    //    qsort(populations, numClones, sizeof(Population), comparePopulationsByTimeOrigin);
}

void Chain::GenerateTimesFromPriorsOriginal(int noisy,  long int *seed) {
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
    Population *popJ, *popI, *popL;
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
                popJ = populations[j];
                Uvector[j] = 0;
                Uvector[j] = RandomUniform (seed);
                rand = RandomUniform(seed);
                popJ->timeOriginSTD= (1/ popJ->delta)*log(1-(popJ->delta / log(rand)));
                //popJ->timeMigrationSTDCurrentPop= popJ->timeOriginSTD ;
                // popJ->timeOriginInput=proportionsVector[j]*popJ->timeOriginSTD;
                popJ->timeOriginInput=popJ->effectPopSize*popJ->timeOriginSTD;
            }
            
            if (numClones == 1) {
                doAcceptTimes = YES;
                break;
            }
            ListClonesAccordingTimeToOrigin(populations);
            // Calculate probabilities P
            LocalProbability = 0.0;
            TotalProbability = 0.0;
            // printf("Number of clones = %d\n", numClones);
            for (i = 0; i < numClones - 1; i++)
            {
                popI=populations[i];
                printf("\ni  = %d ", i);
                LocalProbability = 0.0;
                aboveTerm = 0.0;
                belowTerm = 0.0;
                m = popI->order;
                // printf("\n\nCalculating P for clone = %d ", m);
                for (l = i + 1; l < numClones; l++)
                {
                    popL=populations[ l];
                    j = popL->order ;
                    aboveTerm = aboveTerm + (popL->popSize * Population::CalculateH(popI->timeOriginSTD * popI->effectPopSize / popL->effectPopSize, popL->timeOriginSTD, popL->delta));
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
void Chain::InitChainPopulations( int noisy,  int TotalNumSequences  ) {
    int z;
    //struct Population* pops = malloc(numClones * (sizeof(struct Population)+TotalNumSequences * sizeof( int) + numClones * sizeof(double) ));
//    population* pops =(population*) malloc(numClones * (sizeof( population) ));
    for (z = 0; z <= (numClones - 1); z++)
    {
        int ind = z;
        int ord = 0;
        double timeOriginInput = 0.0;
        int sampleSize = 0;
        int popSize =0;
        double birthRate = 1.0;
        double deathRate = 0.99;
        
        auto pop = new Population(ind, ord, timeOriginInput, sampleSize, popSize, birthRate, deathRate, NO);
        populations.push_back(pop);
//        pops[z].index =z;
//        pops[z].order=0;
//        pops[z].birthRate =0.0;
//        pops[z].deathRate = 0;
//        pops[z].growthRate = 0;
//        pops[z].sampleSize = 0;
//        pops[z].popSize =0;
//        pops[z].effectPopSize =0;
//        pops[z].delta =0 ;
//        pops[z].numActiveGametes = 0;
//        pops[z].isAlive = 1;
//        pops[z].nodeIdAncesterMRCA = 0;
//        pops[z].numCompletedCoalescences = 0;
//        pops[z].nextAvailableIdInmigrant = 0;
//        pops[z].numIncomingMigrations = 0;
//        pops[z].numPossibleMigrations = 0;
//        pops[z].doEstimateTimeOrigin = NO;
//        pops[z].timeOriginInput = 0;
//        pops[z].timeOriginSTD =0;
//        pops[z].timeMigrationSTDCurrentPop = 0;
//        pops[z].FatherPop =NULL;
//
//        *(populations + z) = &pops[z];
//
    }
}
/**************** RelabelNodes **************/

void Chain::RelabelNodes(pll_unode_t *p, int &intLabel)
{
    if (p != NULL)
    {
      
        TreeNode *u= (TreeNode *)(p->data);
        TreeNode *u1= (TreeNode *)(p->next->data);
        TreeNode *u2= (TreeNode *)(p->next->next->data);
        
        RelabelNodes (p->next->back, intLabel);
        RelabelNodes (p->next->next->back, intLabel);
        /*RelabelNodes (p->outgroup);*/
       // if (p->left == NULL && p->right == NULL) /* is tip */
        if (p->back != NULL && p->next->back == NULL) /* is tip */
        {
            // p->label = intLabel++;
            //  p->label = (*intLabel);
            // *intLabel=*intLabel+1;
            //p->label = p->index ;
            
            u->label =u1->label= u2->label= u->index ;
        }
        else                  /* all ancester */
        {
            //p->label = intLabel++;
           // p->label = intLabel;
            u->label =u1->label= u2->label= intLabel ;
            intLabel=intLabel+1;
        }
    }
}
void Chain::InitPopulationSampleSizes(vector<Population*> populations, int TotalSampleSize, int numClones, vector<double> &proportionsVector, long int *seed)
{
    int i,j;
    Population *popI, *popJ;
    double rand;
    double    *cumSum = (double *) calloc((numClones +1), (long) sizeof(double));
    if (!cumSum)
    {
        fprintf (stderr, "Could not allocate Xvector (%lu bytes)\n", (numClones +1 ) * (long) sizeof(double));
        exit (-1);
    }
    cumSum[0]=0.0;
    for (j = 1; j <= numClones; j++)
    {
        cumSum[j]=0;
        cumSum[j]=cumSum[j-1]+proportionsVector[j-1];
        popJ = populations[j - 1];
        popJ->sampleSize=0;
    }
    for (i = 0; i < TotalSampleSize; i++)
    {
        rand=RandomUniform(seed);
        for (j = 1; j <= numClones;j ++)
        {
            popJ = populations[ j-1 ];
            if (rand <= cumSum[j] && rand > cumSum[j - 1])
            {
                popJ->sampleSize=popJ->sampleSize +1;
                break;
            }
        }
    }
    free(cumSum);
    cumSum=NULL;
}
double Chain::SumBranches(pll_unode_t *p, double mutationRate){
    
    static double sum;
    
    if (p != NULL)
    {
        if (p->back  == NULL && p->next !=NULL && p->next->back !=NULL)
            sum = 0;
        else{
            //sum += (p->anc1->time- p->time)* mutationRate;//p->lengthModelUnits;
            sum += p->length;
            //            sum += p->length;
        }
        //            sum += p->lengthModelUnits;//length;
        SumBranches (p->next->back,  mutationRate);
        SumBranches (p->next->next->back,   mutationRate);
    }
    
    return sum;
    
    
}

char * Chain::toNewickString ( pll_unode_t *p, double mutationRate,     int doUseObservedCellNames)
{
    char buffer[1024];
    //char *newickString = malloc( size);
    char *newickString =NULL;
    char *left=NULL;
    char *right=NULL;
    char *outgroup =NULL;

    if (p != NULL)
    {
        TreeNode *u= (TreeNode *)(p->data);
        if (u->isOutgroup == YES)     /* Outgroup */
        {
            strcpy( u->cellName,"healthycell");
            strcpy( p->label,"healthycell");
            //strcat(newickString, "healthycell:%10.9lf", (p->anc1->timePUnits - p->timePUnits) * mutationRate);
            if (asprintf(&newickString,  "healthycell:%10.9lf",  (u->anc1->timePUnits - u->timePUnits) * mutationRate)<0)
                return NULL;
            // snprintf(newickString,  size,  "healthycell:%10.9lf",  (p->anc1->timePUnits - p->timePUnits) * mutationRate);
            return newickString;
        }
       // else if (p->left == NULL && p->right == NULL)   /* tip of the tree */
        else if (p->next != NULL && p->next->back == NULL && p->next->next != NULL && p->next->next->back == NULL)   /* tip of the tree */
        {
            snprintf(buffer, sizeof(buffer), "tip_i%05d_C%d_%d", u->index,u->indexOldClone,u->indexCurrentClone);
            strcpy( u->cellName,buffer);
            strcpy( p->label,buffer);
            //            if (p->isOutgroup == YES)     /* Outgroup */
            //            {
            //                //strcat(newickString, "healthycell:%10.9lf", (p->anc1->timePUnits - p->timePUnits) * mutationRate);
            //                snprintf(newickString,  size,  "healthycell:%10.9lf",  (p->anc1->timePUnits - p->timePUnits) * mutationRate);
            //                return newickString;
            //            }
            //else{
            if  (doUseObservedCellNames == YES)
            {
                if (asprintf(&newickString,   "%s:%10.9lf",  u->observedCellName, (u->anc1->timePUnits - u->timePUnits)*mutationRate)<0)
                    return NULL;
                //snprintf(newickString,  size,  "%s:%10.9lf",  p->observedCellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
                return newickString;
            }
            else{
                if (asprintf(&newickString,   "%s:%10.9lf",  u->cellName, (u->anc1->timePUnits - u->timePUnits)*mutationRate)<0)
                    return NULL;
                //snprintf(newickString,  size,  "%s:%10.9lf",  p->observedCellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
                
                return newickString;
            }
            //  }
        }
        else
        {
            // fprintf (fpTrees2, "(");
           // if ( p->left != NULL  )
            if ( p->next != NULL  &&   p->next->back!=NULL )
            {
                left = toNewickString (p->next->back, mutationRate,   doUseObservedCellNames);
                //right = toNewickString2 (p->right, mutationRate,   cellNames);
                //snprintf(newickString, size, "(%s,%s):%10.9lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate );
                //                free(left);
                //                left = NULL;
                //                free(right);
                //                right=NULL;
                //                return newickString;
            }
            if ( p->next->next != NULL  &&   p->next->next->back!=NULL )
            {
                right = toNewickString( p->next->next->back, mutationRate,   doUseObservedCellNames);
                //                snprintf(newickString, size, "(%s,%s):%10.9lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate );
                //                free(left);
                //                left = NULL;
                //                free(right);
                //                right=NULL;
                //                return newickString;
            }
            //outgroup =toNewickString (u->outgroup, mutationRate,   doUseObservedCellNames);
            if(left!= NULL && right!= NULL && p->back != NULL)
            {
                snprintf(buffer, sizeof(buffer), "int_i%05d_C%d_%d",  u->index,u->indexOldClone,u->indexCurrentClone);
                strcpy( u->cellName,buffer);
                
                if (asprintf(&newickString, "(%s,%s):%10.9lf", left, right,  (u->anc1->timePUnits - u->timePUnits)*mutationRate )<0)
                    return NULL;
                //snprintf(newickString, size, "(%s,%s):%10.9lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate );
                free(left);
                left = NULL;
                free(right);
                right=NULL;
                free(outgroup);
                outgroup=NULL;
            }
            else if (left != NULL &&  right!= NULL  && p->back == NULL)
            {
                snprintf(buffer, sizeof(buffer), "root_i%05d_C%d_%d",  u->index,u->indexOldClone,u->indexCurrentClone);
                strcpy( u->cellName,buffer);
                // left = toNewickString2 (p->left, mutationRate,   cellNames);
                //right = toNewickString2 (p->right, mutationRate,   cellNames);
                if (asprintf(&newickString,  "(%s,%s);", left, right)<0)
                    return NULL;
                //snprintf(newickString, size, "(%s,%s);", left, outgroup);
                free(left);
                left = NULL;
                free(right);
                right=NULL;
                free(outgroup);
                outgroup=NULL;
                // return newickString;
            }
            return newickString;
        }
    }
    return newickString;
}
/************************ LogConditionalLikelihoodTree ***********************/
/*  LogConditionalLikelihoodTree */
double Chain::LogConditionalLikelihoodTree(pll_unode_t  *tree, ProgramOptions &programOptions  )
{
    Population* popI;
    Population* popJ;
    Population* fatherPop;
    double product=0;
    int i, j;
    double temp;
    for ( i = 0; i < numClones; i++)
    {
        popI=populations[i];
        product = product + log( popI->DensityTime( popI->timeOriginSTD));
    }
    for ( j = 0; j < numClones - 1; j++)
    {
        popJ = populations[j ];
        product = product + log( popJ->popSize);
        fatherPop = popJ -> FatherPop;
        temp=popJ->timeOriginSTD * popJ->effectPopSize / fatherPop->effectPopSize;
        temp=Population::CalculateH(popJ->timeOriginSTD * popJ->effectPopSize / fatherPop->effectPopSize, fatherPop->timeOriginSTD, fatherPop->delta);
        product = product  + log( temp);
    }
    //for ( i = 0; i < numClones; i++)
    //  {
    //     popI=*(populations + i );
    product = product + LogDensityCoalescentTimesForPopulation(tree);
    //}
    return product;
}
/************************ LogDensityCoalescentTimesForPopulation ***********************/
/*  LogDensityCoalescentTimesForPopulation */
double Chain::LogDensityCoalescentTimesForPopulation(pll_unode_t  *tree)
{
    double result =0;
    int i, k;
    Population *popI;
    int numberLeftCoalescences;
    int numberLeftMigrations;
    int numberAliveCells;
    int currentCoalescentEvent=0;
    int currentMigrationEvent=0;
    double temp;
    for ( i = 0; i < numClones; i++){
        popI = populations[i];
        currentCoalescentEvent=0;
        currentMigrationEvent=0;
        numberAliveCells= popI->sampleSize;
        // numberLeftCoalescences =  popI->numCompletedCoalescences; //in the numCompletedCoalescences we are considering also the migrations
        numberLeftCoalescences =  popI->numCompletedCoalescences - (popI->numIncomingMigrations-1);
        numberLeftMigrations = popI->numIncomingMigrations-1;
        //we are not counting time of origin as a migration
        if (numberLeftCoalescences ==0)
            return  result;
        while(numberLeftMigrations > 0)
        {
            while(popI->CoalescentEventTimes[currentCoalescentEvent] < popI->immigrantsPopOrderedByModelTime[currentMigrationEvent].first
                 && numberAliveCells > 1)
            {
                temp=log(numberAliveCells * (numberAliveCells-1)/2);
                result= result + temp;
                temp = log(1/Population::CalculateH(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta));
                result= result + temp;
                temp =  (numberAliveCells/2)* (numberAliveCells-1)*(Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta)-Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent+1], popI->timeOriginSTD, popI->delta));
                result= result -temp;
                currentCoalescentEvent++;
                numberLeftCoalescences--;
                numberAliveCells--;
            }
            //if (numberLeftMigrations > 0 && numberAliveCells > 1)// if there are migrations
            if (numberLeftMigrations > 0 )// if there are migrations
            {  temp= popI->LogProbNoCoalescentEventBetweenTimes(popI->CoalescentEventTimes[currentCoalescentEvent],popI->immigrantsPopOrderedByModelTime[currentMigrationEvent].first, numberAliveCells );
                result= result+ temp;
                numberLeftMigrations--;
                currentMigrationEvent++;
            }
        }
        //here there are only coalescents events left(al least one event)
        while(numberLeftCoalescences > 0 && numberAliveCells > 1)
        {   temp = log(numberAliveCells * (numberAliveCells-1)/2);
            result= result + temp;
            temp = log(1/Population::CalculateH(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta));
            result= result + temp;
            temp=( numberAliveCells/2)* (numberAliveCells-1)*(Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta)-Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent+1], popI->timeOriginSTD, popI->delta));
            result= result -  temp;
            currentCoalescentEvent++;
            numberLeftCoalescences--;
            numberAliveCells--;
        }
    }
    
    return result;
}
//double  Chain::LogConditionalLikelihoodSequences(pll_msa_t * msa, char* NewickString, ProgramOptions &programOptions, double seqError,double dropoutError)
//{
//
//    FILE    *outputShell;
//    char script[80];
//    char buf[1000];
//
//    //pll_state_t pll_map_gt10_2[256];
//    if (NewickString == NULL)
//    {
//        fprintf (stderr, "\nERROR: The newick representation of the tree cannot be empty\n\n");
//        PrintUsage();
//        return 0;
//    }
//
//    unsigned int i;
//    unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
//    pll_partition_t * partition;
//
//    pll_utree_t *unrootedTree = pll_utree_parse_newick_string_unroot(NewickString);
//
//    /* compute node count information */
//    tip_nodes_count = unrootedTree->tip_count;
//    inner_nodes_count = unrootedTree->inner_count;
//    nodes_count = inner_nodes_count + tip_nodes_count;
//    branch_count = unrootedTree->edge_count;
//
////    pllmod_treeinfo_t * treeinfo = pllmod_treeinfo_create((pll_unode_s*)unrootedTree->vroot,
////                                                          tip_nodes_count,
////                                                          1,
////                                                          PLLMOD_COMMON_BRLEN_LINKED);
//    pllmod_treeinfo_t * treeinfo = pllmod_treeinfo_create(unrootedTree->vroot,
//                                                          tip_nodes_count,
//                                                          1, PLLMOD_COMMON_BRLEN_LINKED);
//
//    pll_unode_t ** tipnodes = unrootedTree->nodes;
//
//    /* create a libc hash table of size tip_nodes_count */
//    hcreate(tip_nodes_count);
//
//    /* populate a libc hash table with tree tip labels */
//    unsigned int * data = (unsigned int *)malloc(tip_nodes_count *
//                                                 sizeof(unsigned int));
//    char *label;
//    for (i = 0; i < tip_nodes_count; ++i)
//    {
//        data[i] = tipnodes[i]->clv_index;
//        ENTRY entry;
//        label= tipnodes[i]->label;
//        entry.key = tipnodes[i]->label;
//        entry.data = (void *)(data+i);
//        hsearch(entry, ENTER);
//    }
//
//
//    if (msa !=NULL)
//        printf("Original sequence (alignment) length : %d\n", msa->length);
//    else{
//        fprintf (stderr, "\nERROR: The multiple sequence alignment is empty\n\n");
//        PrintUsage();
//        return 0;
//    }
//
//
//
//    pllmod_subst_model_t * model = pllmod_util_model_info_genotype(GT_MODEL);
//
//    /* create the PLL partition instance
//
//     tip_nodes_count : the number of tip sequences we want to have
//     inner_nodes_count : the number of CLV buffers to be allocated for inner nodes
//     model->states : the number of states that our data have
//     1 : number of different substitution models (or eigen decomposition)
//     to use concurrently (i.e. 4 for LG4)
//     branch_count: number of probability matrices to be allocated
//     RATE_CATS : number of rate categories we will use
//     inner_nodes_count : how many scale buffers to use
//     PLL_ATTRIB_ARCH_AVX : list of flags for hardware acceleration
//     */
//    partition = pll_partition_create(tip_nodes_count,
//                                     inner_nodes_count,
//                                     model->states,
//                                     (unsigned int)(msa->length),
//                                     1,
//                                     branch_count,
//                                     RATE_CATS,
//                                     inner_nodes_count,
//                                     PLL_ATTRIB_ARCH_AVX);
//
//    set_partition_tips_costum( partition, msa, programOptions,  seqError, dropoutError);
//
//    treeinfo = pllmod_treeinfo_create(unrootedTree->vroot,
//                                      tip_nodes_count,
//                                      1,
//                                      PLLMOD_COMMON_BRLEN_LINKED);
//    int params_to_optimize = 0;
//    unsigned int params_indices[RATE_CATS] = {0};
//
//    int retval = pllmod_treeinfo_init_partition(treeinfo,
//                                                0,
//                                                partition,
//                                                params_to_optimize,
//                                                PLL_GAMMA_RATES_MEAN,
//                                                1.0, /* alpha*/
//                                                params_indices, /* param_indices */
//                                                model->rate_sym /* subst matrix symmetries*/
//                                                );
//
//    double * empirical_frequencies;
//     empirical_frequencies = pllmod_msa_empirical_frequencies(partition);
//    
//     double * empirical_subst_rates = pllmod_msa_empirical_subst_rates( partition);
//
//    unsigned int * weight = pll_compress_site_patterns(msa->sequence,
//                                                       pll_map_gt10,
//                                                       tip_nodes_count,
//                                                       &(msa->length));
//    printf("Number of unique site patterns: %d\n\n", msa->length);
//
//    /* initialize the array of base frequencies  AA CC GG TT AC/CA AG/GA AT/TA CG/GC CT/TC GT/TG  */
////    double user_freqs[10] = { 0.254504, 0.157238, 0.073667, 0.287452, 0.027336,
////        0.058670, 0.041628, 0.024572, 0.064522, 0.010412 };
//
//
//    /* substitution rates: for GTR4 model those are 6 "regular" DNA susbt. rates + 1 rate
//     * for "unlikely" double substitutions (eg A/A -> C/T) */
//    double unique_subst_rates[7] = { 0.001000, 0.101223, 0.001000, 0.001000, 1.000000,
//        0.001000, 0.447050 };
//
//    /* get full above-diagonal half-matrix */
//    double * user_subst_rates = expand_uniq_rates(model->states, unique_subst_rates,
//                                                  model->rate_sym);
//
//    double rate_cats[RATE_CATS] = {0};
//
//    /* compute the discretized category rates from a gamma distribution
//     with alpha shape 1 and store them in rate_cats  */
//    pll_compute_gamma_cats(1, RATE_CATS, rate_cats, PLL_GAMMA_RATES_MEAN);
//
//    /* set frequencies at model with index 0 (we currently have only one model) */
//    //pll_set_frequencies(partition, 0, model->freqs ? model->freqs : user_freqs);
//    pll_set_frequencies(partition, 0, model->freqs ? model->freqs : empirical_frequencies);
//
//    /* set substitution parameters at model with index 0 */
//   // pll_set_subst_params(partition, 0, model->rates ? model->rates : user_subst_rates);
//     pll_set_subst_params(partition, 0, model->rates ? model->rates : empirical_subst_rates);
//    free(user_subst_rates);
//
//    /* set rate categories */
//    pll_set_category_rates(partition, rate_cats);
//
//    /* set pattern weights and free the weights array */
//    pll_set_pattern_weights(partition, weight);
//    free(weight);
//
//    //set_partition_tips(chain, partition, msa, programOptions);
//
//    // pll_msa_destroy(msa);
//
//    /* destroy hash table */
//    hdestroy();
//    /* we no longer need these two arrays (keys and values of hash table... */
//    free(data);
//
//    if (!retval)
//        fprintf(stderr, "Error initializing partition!");
//    /* Compute initial LH of the starting tree */
//    double loglh = pllmod_treeinfo_compute_loglh(treeinfo, 1);
//    pllmod_treeinfo_destroy(treeinfo);
//    double * clv ;
//    int s, clv_index, j, k;
//    int scaler_index = PLL_SCALE_BUFFER_NONE;
//    unsigned int * scaler = (scaler_index == PLL_SCALE_BUFFER_NONE) ?
//    NULL : partition->scale_buffer[scaler_index];
//    unsigned int states = partition->states;
//    unsigned int states_padded = partition->states_padded;
//    unsigned int rates = partition->rate_cats;
//    double prob;
//    unsigned int *site_id = 0;
//    int index;
//    unsigned int float_precision;
//
//    pll_partition_destroy(partition);
//    /* we will no longer need the tree structure */
//    //pll_utree_destroy(unrootedTree, NULL);
//    destroyTree(unrootedTree, NULL);
//    pllmod_util_model_destroy(model);
//
//    return loglh;
//}
void Chain::set_partition_tips_costum( pll_partition_t * partition, pll_msa_t * msa, ProgramOptions &programOptions, double seqError, double dropoutError)
{
    //pll_state_t pll_map_gt10_2[256];
    int states =10;
    int i, currentState;
    int from, to;
    
    
    unsigned int state;
    //double * _freqs;
    
    /* find sequences in hash table and link them with the corresponding taxa */
    for (i = 0; i < msa->count; ++i)
    {
        ENTRY query;
        query.key = msa->label[i];
        ENTRY * found = NULL;
        
        found = hsearch(query,FIND);
        
        if (!found)
            fprintf(stderr,"Sequence with header %s does not appear in the tree", msa->label[i]);
        
        unsigned int tip_clv_index = *((unsigned int *)(found->data));
        
        if (programOptions.doUseGenotypes == NO)
        {
            pll_set_tip_states(partition, tip_clv_index, pll_map_gt10, msa->sequence[i]);
        }
        else
        {
            
            //            for ( currentState = 0; currentState < states; currentState++)
            //            {
            //                from=1;
            //                to=1;
            this->set_tipclv_custom_error_model( partition,
                        tip_clv_index,
                        pll_map_gt10,
                        msa->sequence[i],seqError,dropoutError);
            //                 compute_state_probs( msa->sequence[i],  &(partition->clv[tip_clv_index]),  states, from, to,
            //                                chain->seqErrorRate,
            //                                 chain->dropoutRate);
            //      }
            
            //pll_set_tip_clv(partition, tip_clv_index, tipCLV, PLL_FALSE);
            
            
        }
    }
    
}


int Chain::set_tipclv_custom_error_model(pll_partition_t * partition,
                unsigned int tip_index,
                const pll_state_t * map,
                const char * sequence,
                double _seq_error_rate,
                double _dropout_rate
                )
{

    pll_state_t c;
    unsigned int i,j;
    double * tipclv = partition->clv[tip_index];
    unsigned int state_id ;
    
    pll_repeats_t * repeats = partition->repeats;
    unsigned int use_repeats = pll_repeats_enabled(partition);
    unsigned int ids = use_repeats ?
    repeats->pernode_ids[tip_index] : partition->sites;
    
    static const double one_3 = 1. / 3.;
    static const double one_6 = 1. / 6.;
    static const double one_8 = 1. / 8.;
    static const double three_8 = 3. / 8.;
    static const double one_12 = 1. / 12.;
    
    // TODO: move it out of here
    unsigned int undef_state = (unsigned int) (pow(2, partition->states)) - 1;
    
    double sum_lh = 0.;
    
    /* iterate through sites */
    for (i = 0; i < ids; ++i)
    {
        unsigned int index = use_repeats ?
        repeats->pernode_id_site[tip_index][i] : i;
        if ((c = map[(int)sequence[index]]) == 0)
        {
            pll_errno = PLL_ERROR_TIPDATA_ILLEGALSTATE;
            snprintf(pll_errmsg, 200, "Illegal state code in tip \"%c\"", sequence[index]);
            return PLL_FAILURE;
        }
        
        /* decompose basecall into the encoded residues and set the appropriate
         positions in the tip vector */
        state_id = __builtin_ctz(c);
        for (j = 0; j < partition->states; ++j)
        {
            if (c == undef_state)
                tipclv[j] = 1.;
            else
            {
                if (j == state_id)
                {
                    /* 0 letters away */
                    if (HOMO(state_id))
                        tipclv[j] = 1. - _seq_error_rate + 0.5 * _seq_error_rate * _dropout_rate;
                    else
                        tipclv[j] =  (1. - _dropout_rate ) * (1. - _seq_error_rate) + one_12 * _seq_error_rate * _dropout_rate;
                }
                else if (mut_dist[state_id][j] == 1)
                {
                    /* 1 letter away */
                    if (HOMO(j))
                    {
                        tipclv[j] = one_12 * _seq_error_rate * _dropout_rate +
                        one_3  * (1. - _dropout_rate) * _seq_error_rate;
                    }
                    else
                    {
                        if (HOMO(state_id))
                        {
                            tipclv[j] = 0.5 * _dropout_rate + one_6 * _seq_error_rate -
                            three_8 * _seq_error_rate * _dropout_rate;
                        }
                        else
                        {
                            tipclv[j]= one_6 * _seq_error_rate -
                            one_8 * _seq_error_rate * _dropout_rate;
                        }
                    }
                }
                else
                {
                    /* 2 letters away */
                    if (HOMO(state_id))
                        tipclv[j] = one_12 * _seq_error_rate * _dropout_rate;
                    else
                        tipclv[j] = 0.;
                }
                sum_lh += tipclv[j];
            }
            // tipclv[j] = c & 1;
            //  c >>= 1;
        }
        
        /* fill in the entries for the other gamma values */
        tipclv += partition->states_padded;
        for (j = 0; j < partition->rate_cats - 1; ++j)
        {
            memcpy(tipclv, tipclv - partition->states_padded,
                   partition->states * sizeof(double));
            tipclv += partition->states_padded;
        }
    }
    
    /* if asc_bias is set, we initialize the additional positions */
    if (partition->asc_bias_alloc)
    {
        for (i = 0; i < partition->states; ++i)
        {
            for (j = 0; j < partition->states; ++j)
            {
                tipclv[j] = j==i;
            }
            
            /* fill in the entries for the other gamma values */
            tipclv += partition->states_padded;
            for (j = 0; j < partition->rate_cats - 1; ++j)
            {
                memcpy(tipclv, tipclv - partition->states_padded,
                       partition->states * sizeof(double));
                tipclv += partition->states_padded;
            }
        }
    }
    
    return PLL_SUCCESS;
}
/************************ destroyTree ***********************/
/*  destroyTree */
void  Chain::destroyTree(pll_utree_t * tree, void (*cb_destroy)(void *))
{
    
    unsigned int i;
    
    /* deallocate tip nodes */
    for (i = 0; i < tree->tip_count; ++i)
    {
        dealloc_data_costum(tree->nodes[i], cb_destroy);
        
        free(tree->nodes[i]);
    }
    /* deallocate inner nodes */
    for (i = tree->tip_count; i < tree->tip_count + tree->inner_count; ++i)
    {
        pll_unode_t * first = tree->nodes[i];
        assert(first);
        if (first->label)
            free(first->label);
        
        pll_unode_t * node = first;
        do
        {
            pll_unode_t * next = node->next;
            dealloc_data_costum(node, cb_destroy);
            free(node);
            node = next;
        }
        while(node && node != first);
    }
    
    /* deallocate tree structure */
    free(tree->nodes);
    free(tree);
    
}
void Chain::dealloc_data_costum(pll_unode_t * node, void (*cb_destroy)(void *))
{
    if (node->data)
    {
        if (cb_destroy)
            cb_destroy(node->data);
    }
}

void Chain::InitializeChains(vector<Chain*> &chains,   ProgramOptions &programOptions,  MCMCoptions &mcmcOptions, vector<int> &sampleSizes, long int *seed, char* ObservedCellNames[], pll_msa_t *msa, pll_utree_t * initialTree)
{
    int chainNumber;
    Chain *currentChain=NULL;
    double totalTreeLength;
    int        numCA, numMIG;
    double      cumNumCA, meanNumCA, cumNumMIG, meanNumMIG, numEventsTot;
    char *newickString2;
    for(chainNumber=0; chainNumber< mcmcOptions.numChains;chainNumber++)
    {
        
        //chains[chainNumber] = (Chain *)malloc(sizeof(Chain));
        auto chain= new Chain(  chainNumber,
                              programOptions.numClones,
                              1,
                              programOptions.mutationRate,
                              programOptions.seqErrorRate,
                              programOptions.dropoutRate
                              );
        
        if (programOptions.numberClonesKnown)
        {
            chain->numNodes = programOptions.numNodes;
            
            
        }
        
        chain->InitChainPopulations( programOptions.noisy, programOptions.TotalNumSequences );
        
        
        chain->FillChainPopulationsFromPriors( programOptions,mcmcOptions, sampleSizes, seed );
        
        if (programOptions.doUseFixedTree == NO)
        {
            chain->root = chain->MakeCoalescenceTree (seed,
                                                      chain->numNodes,
                                                      chain->numClones,
                                                      programOptions,
                                                      cumNumCA,
                                                      meanNumCA,
                                                      cumNumMIG,
                                                      meanNumMIG,
                                                      numMIG,
                                                      numCA,
                                                      numEventsTot,
                                                      ObservedCellNames, sampleSizes
                                                      ) ;
            
        }
        else {
            chain->initialTree = initialTree;
            
        }
        totalTreeLength = chain->SumBranches(chain->root, chain->mutationRate);
        //        cumNumMUperTree=0;
        
        newickString2=NULL;
        newickString2 = chain->toNewickString ( chain->root, chain->mutationRate,     programOptions.doUseObservedCellNames);
        printf("\n newick = %s  \n", newickString2);
        
        
        chain->currentlogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(chain->root, programOptions);
        printf ( "Initial likelihood of the tree of chain %d is:  %lf \n",chainNumber, chain->currentlogConditionalLikelihoodTree );
        
       //chain->currentlogConditionalLikelihoodSequences= chain->LogConditionalLikelihoodSequences( msa,  newickString2, programOptions, chain->seqErrorRate, chain->dropoutRate);
        
        
        fprintf (stderr, "Initial likelihood of the sequences of chain %d  is = %lf  \n", chainNumber,chain->currentlogConditionalLikelihoodSequences );
        
        free(newickString2);
        newickString2=NULL;
        
        chains.push_back(chain);
        
    }
}


