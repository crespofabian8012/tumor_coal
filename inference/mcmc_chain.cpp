//
//  mcmc_chain.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 7/10/19.
//
//#include "mcmc_chain.hpp"
#include <map>
#include <string>
#include <queue>
#include <unordered_set>
#include <numeric>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

extern "C"
{
#include <libpll/pll_msa.h> //for pllmod_msa_empirical_frequencies
#include <libpll/pll_tree.h>
#include <libpll/pllmod_util.h>
#include <libpll/pllmod_common.h>
#include <libpll/pllmod_algorithm.h>
}

#include <boost/algorithm/string.hpp>

#include "mcmc_chain.hpp"
#include "data_utils.hpp"
#include "definitions.hpp"
#include "eigen.hpp"
#include "tree_node.hpp"
#include "random.h"
#include "constants.hpp"
#include "mcmc_move.hpp"

using namespace std;

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
int Chain::setInitialTreeFromNewick(char * NewickString){
    
    
    this->initialUnrootedTree = pll_utree_parse_newick_string_unroot(NewickString);
    this->numNodes = initialUnrootedTree->tip_count +  initialUnrootedTree->inner_count;
   
    if (this->initialUnrootedTree ==NULL)
        return 1;
    else
        return 0;
    
}
int Chain::setInitialTreeUnrootedTree(pll_utree_t *unrootedTree){
    if (unrootedTree !=NULL)
    {
        initialUnrootedTree = unrootedTree;
        
        numNodes = initialUnrootedTree->tip_count +  initialUnrootedTree->inner_count;
        return 0;
    }
    else
        return 1;
}

//void Chain::MakeCoalescenceEvent(Population *population, pll_unode_t **nodes, int numClones, long int* seed, int noisy,  int *numActiveGametes,   int* nextAvailable,
//                                 int*labelNodes, double *currentTime, int *numNodes)
void Chain::MakeCoalescenceEvent( Population *population, vector<pll_unode_t *> &nodes, int numClones, long int* seed, int noisy,   int &numActiveGametes, int &nextAvailable,
                          int &labelNodes, double &currentTime, int &numNodes)
{
 
    // TreeNode  *p, *q, *r;
    pll_unode_t  *p, *q,  *r1, *r2, *r3 ;
    TreeNode *u1, *u2, *u3;
    int firstInd,  secondInd=0, newInd=0;
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
    int   i, j, w, k, isCoalescence, whichInd,
    firstInd,  newInd,
    isMigration, whichClone;
    double     eventTime;
    pll_unode_t  *p, *r, *r1, *r2, *r3 ;
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
    
    int      c,  i,  m, cumIndivid, *activeGametes = NULL, isCoalescence, whichInd,
    newInd, eventNum, numActiveGametes,
    isMigration, whichClone, currentNumberAliveClones;
    int     labelNodes;
    double    currentTime, eventTime;
    pll_unode_t  *p;
   
 
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
    
    Population  *p;
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
void Chain::GenerateEffectPopSizesFromPriors2( int noisy, int doGenerateProportionsVector){
    int i, j;
    Population  *popJ;
 
    if (doGenerateProportionsVector == YES)
    {
        initProportionsVector();
        double alpha[numClones];
        if (proportionsVector.size()==0)
           std::fill_n(alpha, numClones, 1.0);
        
        generateProportionsVectorFromDirichlet(alpha);
    }
    updateEffectPopSizesCurrentProportionsVector();
}
void Chain::updateEffectPopSizesCurrentProportionsVector(){
    int j=0;
    Population  *popJ;
    for (j = 0; j < numClones; j++)
    {
        popJ = populations [j];
        //temp=*(outputVector + i);
        popJ->oldeffectPopSize = popJ->effectPopSize;
        popJ->effectPopSize = proportionsVector[j] * totalEffectPopSize;
    }
}
void Chain::initProportionsVector(){
    for (size_t i = 0; i < numClones; i++) {
        proportionsVector.push_back(0);
        oldproportionsVector.push_back(0);
    }
 
}
void Chain::generateProportionsVectorFromDirichlet(double alpha[]){
    
    double theta[numClones];
//    for (unsigned int i = 0; i < numClones; ++i){
//        theta[i]=0;
//    }
    randomDirichletFromGsl(numClones, alpha, theta);
    for (int i=0; i< numClones; i++)
    {
        proportionsVector.at(i)=theta[i];
        oldproportionsVector.at(i)=theta[i];
    }
    //std::copy(proportionsVector.begin(), proportionsVector.end(), theta);
}
void Chain::initTotalEffectivePopulationSize(MCMCoptions &mcmcOptions, long int *seed)
{
    totalEffectPopSize= RandomLogUniform(mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
}
void Chain::FillChainPopulationsFromPriors( ProgramOptions &programOptions,  MCMCoptions &mcmcOptions, vector<int> &sampleSizes, long int *seed)
{
    Population* popI;
    int i;
    double randomDelta;
    mutationRate = RandomLogUniform(mcmcOptions.MutRatefrom, mcmcOptions.MutRateto);
    initTotalEffectivePopulationSize(mcmcOptions, seed);
    initProportionsVector();
    //double lambda = 1;
    SetPopulationsBirthRate(mcmcOptions.fixedLambda);
    
       // GenerateEffectPopSizesFromPriors2( programOptions.noisy,    YES);
    if (programOptions.doUseFixedTree ==NO)
    {
        double alpha[numClones];
        std::fill_n(alpha, numClones, 1.0);
        generateProportionsVectorFromDirichlet(alpha);
      if (programOptions.populationSampleSizesKnown == NO)
         {
            InitPopulationSampleSizes(populations, programOptions.numCells, programOptions.numClones,proportionsVector, seed);
         }
        //else fill the sample  sizes
      else{
            setChainPopulationSampleSizes( sampleSizes, programOptions);
          }
        //generate the population  sizes
        initEffectPopulationSizesFromProportionsVector();
        initPopulationsCoalTimes();
        
         if (programOptions.populationSampleSizesKnown ==YES)
         {
            for( i = 0 ; i < numClones; i++)
            {
                popI=populations[i];
                do {
                    randomDelta = RandomLogUniform(mcmcOptions.Deltafrom, mcmcOptions.Deltato);
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
        
        initPopulationMigration();
    }
    
}
void Chain::initEffectPopulationSizesFromProportionsVector(){
    int i;
    Population *popI;
    for( i = 0 ; i < numClones; i++)
    {
        popI=populations[i];
        popI->effectPopSize = totalEffectPopSize * proportionsVector[i];
    }
}
void Chain::initTimeOriginSTD(){
    int i;
    Population *popI;
    for( i = 0 ; i < numClones; i++)
    {
        popI=populations[i];
        popI->timeOriginSTD = popI->timeOriginInput /  popI->effectPopSize;
    }
}
void Chain::initPopulationsCoalTimes(){
    int i;
    Population *popI;
    for( i = 0 ; i < numClones; i++)
    {
        popI=populations[i];
        popI->InitCoalescentEvents(numClones);
    }
}
void Chain::initPopulationMigration(){
    int i=0;
    Population *popI;
    InitListPossibleMigrations(populations,numClones);
//    for( i = 0 ; i < numClones; i++)
//    {
//        popI=populations[i];
//        popI->InitListPossibleMigrations(i);
//    }
}
void Chain::setChainPopulationSampleSizes(vector<int > &sampleSizes,  ProgramOptions &programOptions)
{
    int i;
    Population *popI;
    int cumSampleSize=0;
    
    if (sampleSizes.size() == numClones)
       
   
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
       // if (p->back  == NULL && p->next !=NULL && p->next->back !=NULL)
        if (p->next  == NULL )
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
double Chain::SumBranches2(pll_rnode_t *p, double mutationRate){
    
    static double sum;
    
    if (p != NULL)
    {
        // if (p->back  == NULL && p->next !=NULL && p->next->back !=NULL)
        if (p->left  == NULL )
            sum = 0;
        else{
            //sum += (p->anc1->time- p->time)* mutationRate;//p->lengthModelUnits;
            sum += p->length;
            //            sum += p->length;
        }
        //            sum += p->lengthModelUnits;//length;
        SumBranches2 (p->left,  mutationRate);
        SumBranches2 (p->right,   mutationRate);
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
double Chain::LogConditionalLikelihoodTree( ProgramOptions &programOptions  )
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
    product = product + LogDensityCoalescentTimesForPopulation();
    //}
    return product;
}
/************************ LogDensityCoalescentTimesForPopulation ***********************/
/*  LogDensityCoalescentTimesForPopulation */
double Chain::LogDensityCoalescentTimesForPopulation()
{
    double result =0;
    int i;
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
double  Chain::LogConditionalLikelihoodSequences(pll_msa_t * msa, char* NewickString, ProgramOptions &programOptions, double seqError,double dropoutError)
{

  

    //pll_state_t pll_map_gt10_2[256];
    if (NewickString == NULL)
    {
        fprintf (stderr, "\nERROR: The newick representation of the tree cannot be empty\n\n");
        PrintUsage();
        return 0;
    }

    unsigned int i;
    unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
    pll_partition_t * partition;

    pll_utree_t *unrootedTree = pll_utree_parse_newick_string_unroot(NewickString);

    /* compute node count information */
    tip_nodes_count = unrootedTree->tip_count;
    inner_nodes_count = unrootedTree->inner_count;
    nodes_count = inner_nodes_count + tip_nodes_count;
    branch_count = unrootedTree->edge_count;

//    pllmod_treeinfo_t * treeinfo = pllmod_treeinfo_create((pll_unode_s*)unrootedTree->vroot,
//                                                          tip_nodes_count,
//                                                          1,
//                                                          PLLMOD_COMMON_BRLEN_LINKED);
    pllmod_treeinfo_t * treeinfo = pllmod_treeinfo_create(unrootedTree->vroot,
                                                          tip_nodes_count,
                                                          1, PLLMOD_COMMON_BRLEN_LINKED);

    pll_unode_t ** tipnodes = unrootedTree->nodes;

    /* create a libc hash table of size tip_nodes_count */
    hcreate(tip_nodes_count);

    /* populate a libc hash table with tree tip labels */
    unsigned int * data = (unsigned int *)malloc(tip_nodes_count *
                                                 sizeof(unsigned int));
    char *label;
    for (i = 0; i < tip_nodes_count; ++i)
    {
        data[i] = tipnodes[i]->clv_index;
        ENTRY entry;
        label= tipnodes[i]->label;
        entry.key = tipnodes[i]->label;
        entry.data = (void *)(data+i);
        hsearch(entry, ENTER);
    }


    if (msa !=NULL)
        printf("Original sequence (alignment) length : %d\n", msa->length);
    else{
//        fprintf (stderr, "\nERROR: The multiple sequence alignment is empty\n\n");
//        PrintUsage();
//        return 0;
    }

    pllmod_subst_model_t * model = pllmod_util_model_info_genotype(JC_MODEL);

    /* create the PLL partition instance

     tip_nodes_count : the number of tip sequences we want to have
     inner_nodes_count : the number of CLV buffers to be allocated for inner nodes
     model->states : the number of states that our data have
     1 : number of different substitution models (or eigen decomposition)
     to use concurrently (i.e. 4 for LG4)
     branch_count: number of probability matrices to be allocated
     RATE_CATS : number of rate categories we will use
     inner_nodes_count : how many scale buffers to use
     PLL_ATTRIB_ARCH_AVX : list of flags for hardware acceleration
     */
    partition = pll_partition_create(tip_nodes_count,
                                     inner_nodes_count,
                                     model->states,
                                     (unsigned int)(msa->length),
                                     1,
                                     branch_count,
                                     RATE_CATS,
                                     inner_nodes_count,
                                     PLL_ATTRIB_ARCH_AVX);

    set_partition_tips_costum( partition, msa, programOptions,  seqError, dropoutError);

    treeinfo = pllmod_treeinfo_create(unrootedTree->vroot,
                                      tip_nodes_count,
                                      1,
                                      PLLMOD_COMMON_BRLEN_LINKED);
    int params_to_optimize = 0;
    unsigned int params_indices[RATE_CATS] = {0};

    int retval = pllmod_treeinfo_init_partition(treeinfo,
                                                0,
                                                partition,
                                                params_to_optimize,
                                                PLL_GAMMA_RATES_MEAN,
                                                1.0, /* alpha*/
                                                params_indices, /* param_indices */
                                                model->rate_sym /* subst matrix symmetries*/
                                                );

    double * empirical_frequencies;
     empirical_frequencies = pllmod_msa_empirical_frequencies(partition);
    
     double * empirical_subst_rates = pllmod_msa_empirical_subst_rates( partition);

    unsigned int * weight = pll_compress_site_patterns(msa->sequence,
                                                       pll_map_gt10,
                                                       tip_nodes_count,
                                                       &(msa->length));
    printf("Number of unique site patterns: %d\n\n", msa->length);

    /* initialize the array of base frequencies  AA CC GG TT AC/CA AG/GA AT/TA CG/GC CT/TC GT/TG  */
//    double user_freqs[10] = { 0.254504, 0.157238, 0.073667, 0.287452, 0.027336,
//        0.058670, 0.041628, 0.024572, 0.064522, 0.010412 };


    /* substitution rates: for GTR4 model those are 6 "regular" DNA susbt. rates + 1 rate
     * for "unlikely" double substitutions (eg A/A -> C/T) */
    double unique_subst_rates[7] = { 0.001000, 0.101223, 0.001000, 0.001000, 1.000000,
        0.001000, 0.447050 };

    /* get full above-diagonal half-matrix */
    double * user_subst_rates = Chain::expand_uniq_rates(model->states, unique_subst_rates,
                                                  model->rate_sym);

    double rate_cats[RATE_CATS] = {0};

    /* compute the discretized category rates from a gamma distribution
     with alpha shape 1 and store them in rate_cats  */
    pll_compute_gamma_cats(1, RATE_CATS, rate_cats, PLL_GAMMA_RATES_MEAN);

    /* set frequencies at model with index 0 (we currently have only one model) */
    //pll_set_frequencies(partition, 0, model->freqs ? model->freqs : user_freqs);
    
    
    //pll_set_frequencies(partition, 0, model->freqs ? model->freqs : empirical_frequencies);
    if (model->freqs)
        pll_set_frequencies(partition, 0,  model->freqs);
    else
        pll_set_frequencies(partition, 0,  empirical_frequencies);
    /* set substitution parameters at model with index 0 */
   // pll_set_subst_params(partition, 0, model->rates ? model->rates : user_subst_rates);
     //pll_set_subst_params(partition, 0, model->rates ? model->rates : empirical_subst_rates);
    if (model->rates)
        pll_set_subst_params(partition, 0, model->rates);
    else
        pll_set_subst_params(partition, 0, empirical_subst_rates);
    free(user_subst_rates);

    /* set rate categories */
    pll_set_category_rates(partition, rate_cats);

    /* set pattern weights and free the weights array */
    pll_set_pattern_weights(partition, weight);
    free(weight);

    //set_partition_tips(chain, partition, msa, programOptions);

    // pll_msa_destroy(msa);

    /* destroy hash table */
    hdestroy();
    /* we no longer need these two arrays (keys and values of hash table... */
    free(data);

    if (!retval)
        fprintf(stderr, "Error initializing partition!");
    /* Compute initial LH of the starting tree */
    double loglh = pllmod_treeinfo_compute_loglh(treeinfo, 1);
    pllmod_treeinfo_destroy(treeinfo);
    double * clv ;
    int s, clv_index, j, k;
    int scaler_index = PLL_SCALE_BUFFER_NONE;
    unsigned int * scaler = (scaler_index == PLL_SCALE_BUFFER_NONE) ?
    NULL : partition->scale_buffer[scaler_index];
    unsigned int states = partition->states;
    unsigned int states_padded = partition->states_padded;
    unsigned int rates = partition->rate_cats;
    double prob;
    unsigned int *site_id = 0;
    int index;
    unsigned int float_precision;

    pll_partition_destroy(partition);
    /* we will no longer need the tree structure */
    //pll_utree_destroy(unrootedTree, NULL);
    destroyTree(unrootedTree, NULL);
    pllmod_util_model_destroy(model);

    return loglh;
}
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

Chain *Chain::initializeChains(vector<Chain*> &chains,   ProgramOptions &programOptions,  MCMCoptions &mcmcOptions, vector<int> &sampleSizes, long int *seed, char* ObservedCellNames[], pll_msa_t *msa, pll_utree_t * initialTree, pll_rtree_t * initialRootedTree, string& healthyTipLabel)
{
    double totalTreeLength;
    int        numCA, numMIG;
    double      cumNumCA, meanNumCA, cumNumMIG, meanNumMIG, numEventsTot;
    char *newickString2;
    char * rootedNewick2;
    Chain *chain=new Chain(0, programOptions.numClones, 1, programOptions.mutationRate, programOptions.seqErrorRate, programOptions.dropoutRate);
    
    if (programOptions.numberClonesKnown)
    {
        chain->numNodes = programOptions.numNodes;
    }
    chain->InitChainPopulations( programOptions.noisy, programOptions.TotalNumSequences );
    if (programOptions.doUsefixedMutationRate)
        chain->mutationRate = programOptions.mutationRate;
    else
       chain->mutationRate = RandomLogUniform(mcmcOptions.MutRatefrom, mcmcOptions.MutRateto);
    chain->initTotalEffectivePopulationSize(mcmcOptions, seed);
    chain->SetPopulationsBirthRate(1);
    //chain->FillChainPopulationsFromPriors( programOptions,mcmcOptions, sampleSizes, seed );
    
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
    else
    {
        chain->initialUnrootedTree = initialTree;
        chain->initialRootedTree = initialRootedTree;
        chain->rescaleRootedTreeBranchLengths(chain->mutationRate);
        //chain->root = initialRootedTree->root;
        chain->rootRootedTree = initialRootedTree->root;
        chain->root =  chain->initialUnrootedTree->nodes[chain->initialUnrootedTree->tip_count + chain->initialUnrootedTree->inner_count - 1];
        chain->initPopulationsTipsFromTree(initialTree, NO);
        chain->initPopulationsTipsFromRootedTree(initialRootedTree, NO);
        chain->initNodeDataFromTree();
        chain->initNodeDataFromRootedTree();
        chain->initNumberTipsSubTree(chain->rootRootedTree);
       // pll_utree_traverse_apply(chain->initialUnrootedTree->vroot, NULL, Chain::computeNumberTipsSubTree, chain->root->data);
        //std::map<pll_unode_t*, Population*> mrcaOfPopulation=chain->chooseTimeOfOriginsOnTree( seed);
        bool existsZeroSampleSizePop=false;
        double alpha[chain->numClones];
        int numberPoints = chain->numClones -1;
        //std::map<pll_rnode_t*, Population*>  rmrcaOfPopulation;
        do {
            chain->rMRCAPopulation.clear();
            chain->rMRCAPopulation = chain->initTimeOfOriginsOnRootedTree(numberPoints,  healthyTipLabel);
            
            chain->initPopulationsSampleSizes( chain->rMRCAPopulation);
           
            for (unsigned int i = 0; i < chain->numClones; ++i){
                auto pop =  chain->populations[i];
                alpha[i]= pop->sampleSize;
                if (pop->sampleSize == 0)
                    existsZeroSampleSizePop=true;
            }
            existsZeroSampleSizePop=false;
        }
        while(existsZeroSampleSizePop);
            
         int totalSampleSize=chain->initialRootedTree->tip_count-1;//not the healthytip
        std::transform(alpha, alpha + chain->numClones , alpha,std::bind2nd(std::divides<double>(),totalSampleSize));
        chain->initProportionsVector();
        chain->generateProportionsVectorFromDirichlet(alpha);
        chain->initEffectPopulationSizesFromProportionsVector();
        chain->initTimeOriginSTD();
        chain->initPopulationMigration();//after setting the timeSTD
        chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, healthyTipLabel);
        chain->filterSortPopulationsCoalescentEvents();
        chain->SetPopulationsBirthRate(mcmcOptions.fixedLambda);
        chain->samplePopulationDeltaFromPriors(mcmcOptions, seed );
        //chain->initPopulationsCoalescentAndMigrationEvents(mrcaOfPopulation);
        
        //chain->initializeCoalescentEventTimes(chain->initialUnrootedTree, sampleSizes);
    }
    totalTreeLength = chain->SumBranches2(chain->rootRootedTree, chain->mutationRate);
    //        cumNumMUperTree=0;

    rootedNewick2 =  pll_rtree_export_newick(initialRootedTree->root,NULL);
    //newickString2 = chain->toNewickString ( chain->root, chain->mutationRate,     programOptions.doUseObservedCellNames);
    printf("\n newick = %s  \n", rootedNewick2);
    
    chain->currentlogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree( programOptions);
    printf ( "Initial likelihood of the tree of chain %d is:  %lf \n",0, chain->currentlogConditionalLikelihoodTree );
    
    chain->currentlogConditionalLikelihoodSequences= chain->LogConditionalLikelihoodSequences( msa,  rootedNewick2, programOptions, chain->seqErrorRate, chain->dropoutRate);
    fprintf (stderr, "Initial likelihood of the sequences of chain %d  is = %lf  \n", 0,chain->currentlogConditionalLikelihoodSequences );
    free(rootedNewick2);
    rootedNewick2=NULL;
    return chain;
//    for(chainNumber=0; chainNumber< mcmcOptions.numChains;chainNumber++)
//    {
//        auto chain= new Chain(chainNumber,
//                              programOptions.numClones,
//                              1,
//                              programOptions.mutationRate,
//                              programOptions.seqErrorRate,
//                              programOptions.dropoutRate
//                              );
    
//        if (programOptions.numberClonesKnown)
//        {
//            chain->numNodes = programOptions.numNodes;
//
//
//        }
//
//        chain->InitChainPopulations( programOptions.noisy, programOptions.TotalNumSequences );
//
//
//        chain->FillChainPopulationsFromPriors( programOptions,mcmcOptions, sampleSizes, seed );
//
//        if (programOptions.doUseFixedTree == NO)
//        {
//            chain->root = chain->MakeCoalescenceTree (seed,
//                                                      chain->numNodes,
//                                                      chain->numClones,
//                                                      programOptions,
//                                                      cumNumCA,
//                                                      meanNumCA,
//                                                      cumNumMIG,
//                                                      meanNumMIG,
//                                                      numMIG,
//                                                      numCA,
//                                                      numEventsTot,
//                                                      ObservedCellNames, sampleSizes
//                                                      ) ;
//
//        }
//        else
//        {
//            chain->initialUnrootedTree = initialTree;
//            chain->initialRootedTree = initialRootedTree;
//            //chain->root = initialRootedTree->root;
//            chain->rootRootedTree = initialRootedTree->root;
//            chain->root =  chain->initialUnrootedTree->nodes[chain->initialUnrootedTree->tip_count + chain->initialUnrootedTree->inner_count - 1];
//
//            chain->initializeCoalescentEventTimes(chain->initialUnrootedTree, sampleSizes);
//        }
//        totalTreeLength = chain->SumBranches(chain->root, chain->mutationRate);
//        //        cumNumMUperTree=0;
//
//        newickString2=NULL;
//        newickString2 = chain->toNewickString ( chain->root, chain->mutationRate,     programOptions.doUseObservedCellNames);
//        printf("\n newick = %s  \n", newickString2);
//
//
//        chain->currentlogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(chain->root, programOptions);
//        printf ( "Initial likelihood of the tree of chain %d is:  %lf \n",chainNumber, chain->currentlogConditionalLikelihoodTree );
//
//       //chain->currentlogConditionalLikelihoodSequences= chain->LogConditionalLikelihoodSequences( msa,  newickString2, programOptions, chain->seqErrorRate, chain->dropoutRate);
//
//
//        fprintf (stderr, "Initial likelihood of the sequences of chain %d  is = %lf  \n", chainNumber,chain->currentlogConditionalLikelihoodSequences );
//
//        free(newickString2);
//        newickString2=NULL;
//
//        chains.push_back(chain);
//
//    }
}
int* Chain::computeNumberTipsSubTree(pll_unode_t *node, void *data)
{
    TreeNode * treeNode= (TreeNode *) data;
    if (node->next ==NULL)//tip
        treeNode->numberOfTipsSubTree =0;
    else {
         TreeNode * treeNodeNext= (TreeNode *) (node->next->back->data);
         TreeNode * treeNodeNextNext= (TreeNode *) (node->next->next->back->data);
        treeNode->numberOfTipsSubTree = treeNodeNext->numberOfTipsSubTree + treeNodeNextNext->numberOfTipsSubTree;
        }
}
void  Chain::initNumberTipsSubTree(pll_rnode_t *node)
{
      TreeNode* treeNode= (TreeNode *)(node->data);
     if (node->left ==NULL)//tip
        treeNode->numberOfTipsSubTree =1;
      else {
          initNumberTipsSubTree(node->left);
          initNumberTipsSubTree(node->right);
          TreeNode * treeNodeLeft= (TreeNode *) (node->left->data);
          TreeNode * treeNodeRight= (TreeNode *) (node->right->data);
          treeNode->numberOfTipsSubTree = treeNodeLeft->numberOfTipsSubTree + treeNodeRight->numberOfTipsSubTree;
      }
}
void Chain::runChain(   MCMCoptions &mcmcOptions,  long int *seed,  FilePaths &filePaths, Files &files,  ProgramOptions &programOptions,
              char* ObservedCellNames[], pll_msa_t * msa, vector<int> &sampleSizes
              ){
    
    vector<double>  varTimeGMRCA;
    //Population** populations;
    //TreeNode** nodes; TreeNode** treeTips;
    if (numClones <= 0)
    {
        fprintf (stderr, "\nERROR: The number of clones cannot be negative.");
        PrintUsage();
        
    }
    Population *popI;
    double      cumNumCA, meanNumCA, cumNumMIG, meanNumMIG, numEventsTot, countTMRCA;
    int        numCA, numMIG;
    int i, j, k;
    double    *proportionsVector;
    long numCells;
    int currentIteration;
    double randomDelta;
    double meanNumSNVs, meanNumMU, meanNumDEL, meanNumCNLOH;
    double   cumNumSNVs, cumNumMU, cumNumDEL, cumNumCNLOH, cumCountMLgenotypeErrors;
    double cumNumMUSq, cumNumSNVsSq, cumNumDELSq, cumNumCNLOHSq;
    double varNumMU, varNumSNVs, varNumDEL, varNumCNLOH;
    double expNumMU, expVarNumMU;
    double logConditionalLikelihoodTree;
    double totalTreeLength;
    cumNumCA=0;
    meanNumCA=0;
    cumNumMIG=0;
    meanNumMIG=0;
   
    
    
    if (programOptions.doPrintSeparateReplicates == YES)
        PrepareSeparateFiles(chainNumber, currentIteration, 1,  filePaths, programOptions, files);
    
    numCA = numMIG = 0;
    numEventsTot=0;
    countTMRCA = 0.0;
    
    // order of proposals:
    fprintf (stderr, "\n>> Current log conditional Likelihood tree of the chain %d is = %lf  \n", chainNumber,currentlogConditionalLikelihoodTree );
    
    // 1- Update the topology using Wilson Balding move
//    if (programOptions->doUseFixedTree == NO)
//        WilsonBaldingMove( seed,  chain->currentlogConditionalLikelihoodSequences, programOptions,  ObservedCellNames,  msa);
    
    //2- Update M_T
    NewTotalEffectPopSizeMove *newTotalEffectPopSizeMove= new NewTotalEffectPopSizeMove(this, "new Total Effect Pop Size Move");
    newTotalEffectPopSizeMove->move(programOptions, mcmcOptions);
//    newTotalEffectivePopulationSizeMove(programOptions, ObservedCellNames,  msa,  opt, sampleSizes);
    
    //3-Update the vector of proportions \theta(this will change M_i)
    
    NewProportionsVectorMove *newProportionsVector= new NewProportionsVectorMove(this, "new Proportions Vector Move");
    newProportionsVector->move(programOptions, mcmcOptions);
  //newProportionsVectorMove(  programOptions, ObservedCellNames, msa, opt, sampleSizes);
    
    //        4- Update the Delta_i
    NewGrowthRateMoveForPopulation *newGrowthRateMoveForPopulation;
    for( i = 0 ; i < numClones; i++)
    {
        popI=populations[i];
        newGrowthRateMoveForPopulation= new NewGrowthRateMoveForPopulation(this, "new Growth Rate Move for population", popI);
        newGrowthRateMoveForPopulation->move(programOptions, mcmcOptions);
        //newScaledGrowthRateMoveforPopulation( popI, seed,  programOptions,ObservedCellNames, msa, opt, sampleSizes);
    }
    //        5A-Update  T_i
    //this move will update the time of origin of a population of order i. For a fixed tree, the move will try to choose another available edge
    NewTimeOriginOnTreeforPopulationMove *newTimeOriginOnTreeforPopulationMove ;
    for( i = 0 ; i < numClones; i++)
    {
        popI=populations[i];
        newTimeOriginOnTreeforPopulationMove= new NewTimeOriginOnTreeforPopulationMove(this, "new Time of Origin Move for population", popI);
        newTimeOriginOnTreeforPopulationMove->move(programOptions, mcmcOptions);
        //newScaledGrowthRateMoveforPopulation( popI, seed,  programOptions,ObservedCellNames, msa, opt, sampleSizes);
    }
    //        5B-Update  T_i
    //this move will update the time of origin of a population of order i. For a fixed tree, the move will try to change the time of origin within the same edge
    
    
    
    //        6-Update the sub tree_i for sub population i by changing the topology of the sub tree of the time of some internal nodes.
    
    //slidingWindowMove
    
    
    //proposalChangeCoalTimeInternalNodePopulation(chain , programOptions, root, nodes, Population *pop, int numNodes, long int *seed, double mutationRate, pll_msa_t * msa);
    
    
    //   }
    
}
//void Chain::WilsonBaldingMove( long int *seed, double currentLogLikelihoodSequences, ProgramOptions &programOptions, char *ObservedCellNames[], pll_msa_t * msa)
//{
//    if (numClones > 1)
//    {
//        int i,j;
//        double currentProb, proposalProb;
//        pll_utree_rb_t * rb = (pll_utree_rb_t *)  malloc( (long)sizeof(pll_utree_rb_t));;
//        int min =0;
//        int max=numClones -2;// except the oldest population
//        int orderPopulationToDisconnect = rand() % (max + 1 - min) + min;// generate random integer between 0 and max -min = max-0=max= numClones -2
//        int orderPopulationToReAttach;
//        int indexPopulationToReAttach;
//        int indexNodeToReAttach =0;
//        TreeNode *tempNode;
//        Population *popI=*(chain->populations + orderPopulationToDisconnect);
//        Population *proposalPop;
//        Population *currentFatherPop;
//        TreeNode * MRCA=  popI->MRCA;
//
//        currentFatherPop= popI->FatherPop;
//
//        double timeOrigin = popI->timeOriginSTD;
//        TreeNode *nodeReAttach;
//        Population *tempPopulation;
//        double timeOriginScaledOtherUnits;
//        double distanceModelTimeProposalEdge;
//        double acceptanceProbability=0;
//        int numCandidateNodes=0;
//        double sumNumerators=0;
//        double sumDenominators=0;
//
//        TreeNode *listCandidateNodes[chain->numNodes];
//        //finding branches candidates
//        for( i = orderPopulationToDisconnect+1 ; i < chain->numClones; i++)
//        {
//            tempPopulation = *(chain->populations + i);
//            for (j = tempPopulation->sampleSize  ; j < tempPopulation->numGametes ; j++)
//            {
//                // tempNode= *(nodes) + tempPopulation->idsGametes[j];
//                tempNode= chain->nodes + tempPopulation->idsGametes[j];
//                timeOriginScaledOtherUnits = (timeOrigin * popI->effectPopSize ) /(tempPopulation->effectPopSize);
//                if ( tempNode->index != MRCA->anc1->index && tempNode->index != MRCA->index   && timeOriginScaledOtherUnits > tempNode->time )
//                    // if tempNode is older in time backwards and different from MRCA
//                {
//                    if (( tempNode->anc1 != NULL ) && (tempNode->anc1->orderCurrentClone == tempNode->orderCurrentClone)  && (tempNode->anc1->index != MRCA->anc1->index) && ( timeOriginScaledOtherUnits <= tempNode->anc1->time ))
//                        //if tempNode has father node and the father node is older than  the rescaled time of origin
//                    {
//                        listCandidateNodes[numCandidateNodes]= tempNode;
//                        numCandidateNodes ++;
//                    }
//                }
//            }
//        }
//        if (numCandidateNodes ==0)
//        {
//            fprintf (stderr, "\n: No candidate edges  to Wilson Balding move.\n");
//            return;
//        }
//
//        indexNodeToReAttach = rand() % (numCandidateNodes  - min) + min;
//
//        nodeReAttach = listCandidateNodes[indexNodeToReAttach];
//        indexPopulationToReAttach = nodeReAttach->indexCurrentClone;
//        orderPopulationToReAttach = nodeReAttach->orderCurrentClone;
//
//        if (orderPopulationToReAttach <= (chain->numClones- 1))
//            proposalPop =  *(chain->populations + orderPopulationToReAttach);
//        else{
//            fprintf (stderr, "\n: Index of population outside range.\n");
//            return;
//        }
//        // proposalPop = *(populations + orderPopulationToReAttach);
//
//        pll_tree_edge_t  * edgeReAttach = nodeReAttach->edgeBack;
//
//        pll_tree_rollback_t * rollback_info= (pll_tree_rollback_t *)  malloc( (long)sizeof(pll_tree_rollback_t));
//
//        if (currentFatherPop !=NULL && proposalPop !=NULL )
//        {
//            if (currentFatherPop->order !=proposalPop->order)
//            {
//                proposalProb=ProbabilityCloneiFromClonej2(popI, proposalPop, chain->populations, chain->numClones);
//                currentProb=ProbabilityCloneiFromClonej2(popI, currentFatherPop , chain->populations, chain->numClones);
//                sumNumerators= sumNumerators + log(proposalProb);
//                sumDenominators= sumDenominators + log(currentProb);
//
//                popI->oldFatherPop =currentFatherPop;
//                popI->FatherPop = proposalPop;
//
//            }
//        }
//        distanceModelTimeProposalEdge=(nodeReAttach->anc1->time - (timeOrigin * popI->effectPopSize ) /(proposalPop->effectPopSize));
//        double time_from_r_to_p = RandomUniform(seed)*(distanceModelTimeProposalEdge *proposalPop->effectPopSize)  ;
//        double effect_pop_size =proposalPop->effectPopSize ;
//
//        if (distanceModelTimeProposalEdge >0)
//            sumNumerators= sumNumerators + log(distanceModelTimeProposalEdge);
//
//        double distanceModelTimeCurrentEdge =MRCA->anc1->time- (timeOrigin * popI->effectPopSize ) /(currentFatherPop->effectPopSize);
//        if (distanceModelTimeCurrentEdge >0)
//            sumDenominators= sumDenominators +  log(distanceModelTimeCurrentEdge);
//
//        char *oldNewickString=NULL;
//        oldNewickString = toNewickString2 ( chain->root, chain->mutationRate, programOptions->doUseObservedCellNames);
//        printf("\n before newick = %s  \n", oldNewickString);
//        oldNewickString=NULL;
//        //oldNewickString = toNewickString3 ( treeRootInit[0], treeRootInit[0]->nodeBack, mutationRate, ObservedCellNames);
//        //  oldNewickString=NULL;
//        // printf("\n  before newick = %s  \n", oldNewickString);
//        TreeNode * container_of_u;
//        TreeNode * container_of_v;
//
//        if (! pllmod_utree_spr1(MRCA->nodeBack->back, nodeReAttach->nodeBack->back, rollback_info,rb,  &(chain->nodes), time_from_r_to_p, currentFatherPop->effectPopSize,
//                                proposalPop->effectPopSize, MRCA->anc1, nodeReAttach->anc1,
//                                MRCA,//mrca
//                                nodeReAttach,//the other end point of the edge to reconnect
//                                container_of_u,
//                                container_of_v
//                                ))//first argument is the back nodelet of the ancester of the MRCA and
//        {
//            fprintf (stderr, "\nERROR: Wilson Balding move cannot be applied.\n");
//
//            popI->FatherPop = popI->oldFatherPop;
//
//            free(oldNewickString);
//            oldNewickString= NULL;
//            return;
//            //exit (1);
//        }
//        char *newNewickString=NULL;
//        //   newNewickString=NULL;
//        //   newNewickString = toNewickString3 ( treeRootInit[0], treeRootInit[0]->nodeBack, mutationRate, ObservedCellNames);
//        // printf("\n after newick 1 = %s  \n", newNewickString);
//        newNewickString = toNewickString2 ( root, chain->mutationRate, programOptions->doUseObservedCellNames);
//        printf("\n after newick 2 = %s  \n", newNewickString);
//
//        //compute the new likelihood
//        double newlogConditionalLikelihoodSequences;
//        if (newNewickString !=NULL)
//            newlogConditionalLikelihoodSequences = LogConditionalLikelihoodSequences(chain, msa,  newNewickString, programOptions);
//
//
//        free(newNewickString);
//        newNewickString=NULL;
//        //        logl = pllmod_utree_compute_lk(partition,
//        //                                       tree,
//        //                                       params_indices,
//        //                                       1,
//        //                                       1);
//
//        /* compute marginal likelihoods */
//        //        printf("\nMarginal likelihoods:\n");
//        //        logl = pll_compute_root_loglikelihood (partition,
//        //                                               tree->clv_index,
//        //                                               tree->scaler_index,
//        //                                               params_indices,
//        //                                               NULL);
//        //        printf ("  Log-L Partial at %s: %f\n", tree->label, logl);
//        //        logl = pll_compute_root_loglikelihood (partition,
//        //                                               tree->back->clv_index,
//        //                                               tree->back->scaler_index,
//        //                                               params_indices,
//        //                                               NULL);
//        //        printf ("  Log-L Partial at %s: %f\n", tree->back->label, logl);
//        //
//        //        /* compute global likelihood */
//
//        sumNumerators = sumNumerators + newlogConditionalLikelihoodSequences;
//        sumDenominators = sumDenominators + currentLogLikelihoodSequences;
//
//
//        acceptanceProbability =sumNumerators- sumDenominators;
//
//        double randomNumber= RandomUniform (seed);
//
//        double LogAcceptanceRate = (sumNumerators - sumDenominators) >0? (sumNumerators - sumDenominators) :0;
//
//        if (log(randomNumber) < LogAcceptanceRate )
//        {
//            // accept the WB move
//            printf("\n Accepted Wilson Balding move\n");
//            //int result= pll_utree_tbr(MRCA->nodeBack,edgeReAttach,rollback_info);
//            chain->currentlogConditionalLikelihoodSequences =newlogConditionalLikelihoodSequences;
//            double  currentLikelihoodTree= chain->currentlogConditionalLikelihoodTree;
//            //update the coalescentEventTimes
//            double newLogConditionalLikelihoodTree= LogConditionalLikelihoodTree(chain->root, chain->nodes, chain->populations,  chain->numClones);
//            chain->currentlogConditionalLikelihoodTree =newLogConditionalLikelihoodTree;
//
//        }
//        else {
//            //rollback the move
//            if (rollback_info !=NULL){
//                printf("\n Rejected Wilson Balding move\n");
//
//                pllmod_utree_rollback_spr1(rollback_info,rb, &(chain->nodes),
//                                           currentFatherPop->effectPopSize,
//                                           proposalPop->effectPopSize,
//                                           MRCA->anc1,
//                                           nodeReAttach->anc1,
//                                           MRCA,//mrca
//                                           nodeReAttach,
//                                           container_of_u,
//                                           container_of_v
//                                           );
//                free(rollback_info);
//                free(rb);
//                rollback_info = NULL;
//                rb=NULL;
//            }
//
//        }
//
//}
void Chain::newScaledGrowthRateMoveforPopulation( Population *popI, long int *seed,  ProgramOptions &programOptions, char *ObservedCellNames[], pll_msa_t * msa, MCMCoptions & mcmcOptions, vector<int> &sampleSizes)
{
    ListClonesAccordingTimeToOrigin(populations);

    
    //popI=*(populations + index);
    double randomDelta = RandomLogUniform(mcmcOptions.Deltafrom, mcmcOptions.Deltato);
    //popI->delta = chain->proportionsVector[i] * randomDelta;
    popI->olddelta= popI->delta;
    popI->delta =  randomDelta;
    popI->growthRate =popI->delta  / popI->effectPopSize;
    popI->deathRate= popI->birthRate - popI->growthRate;
    
    ListClonesAccordingTimeToOrigin(populations);
    
    double newLogConditionalLikelihoodTree= LogConditionalLikelihoodTree(programOptions);
    
    fprintf (stderr, "\n>> New log conditional Likelihood tree after the new total growth rate for population %d,  of the chain %d is = %lf  \n", popI->index, chainNumber, newLogConditionalLikelihoodTree );
    
//    char *newickString2;
//    newickString2=NULL;
//    newickString2 = toNewickString2 ( oldroot, programOptions.mutationRate,     programOptions.doUseObservedCellNames);
//    printf("\n newick after move= %s  \n", newickString2);
    
//    newickString2 = toNewickString4 ( oldroot, programOptions.mutationRate,     programOptions.doUseObservedCellNames);
//    printf("\n newick after move= %s  \n", newickString2);
    
//    free(newickString2);
//    newickString2=NULL;
    
    
    double priorDensityNewScaledGrowthRate =  LogUniformDensity(randomDelta, mcmcOptions.Deltafrom, mcmcOptions.Deltato);
    
    double priorDensitScaledGrowthRate =  LogUniformDensity(popI->olddelta, mcmcOptions.Deltafrom, mcmcOptions.Deltato);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewScaledGrowthRate;
    
    double sumLogDenominators=currentlogConditionalLikelihoodTree +priorDensitScaledGrowthRate;
    
    double randomNumber= RandomUniform (seed);
    
    double LogAcceptanceRate = (sumLogNumerators - sumLogDenominators) >0? (sumLogNumerators - sumLogDenominators) :0;
    
    if (log(randomNumber) < LogAcceptanceRate )
    {//accept the move
        printf("\n Accepted new growth rate move\n");
        currentlogConditionalLikelihoodTree =newLogConditionalLikelihoodTree;
    }
    else
    {
        //reject the move
        printf("\n Rejected new growth rate move\n");
        popI->olddelta= popI->delta;
    }
}

int Chain::getPopulationIndex(char * label)
{
    string nodeLabelString = std::string(label);
    vector<string> strs;
    boost::split(strs,nodeLabelString,boost::is_any_of("_"));
    if (strs.size() == 4) {
        int indexPopulation = stoi(strs[2].substr(1, strs[2].size()));
        return indexPopulation;
    }
    
    return -1;
}


bool is_leaf(pll_unode_t *node)
{
    // not complete -- should handle all cases
   // if (node->back == NULL && node->next !=NULL && node->next->back == NULL && node->next->next!=NULL
    //    && node->next->next->back!=NULL)
    if (node->next == NULL)
    {
        return true;
    }
    return false;
}

void updateCoalTimes(Population *pop)
{
    vector<pll_unode_t *> tips = pop->tips;
    deque<pll_unode_t *> q(tips.begin(), tips.end());
    unordered_map<pll_unode_t *, double> nodes_visited;
    for (auto it = tips.begin(); it != tips.end(); ++it) {
        nodes_visited[*it] = 0;
    }
    while (q.size() != 0) {
        pll_unode_t *node = q.front();
        q.pop_front();
        double time = nodes_visited[node];
        double length = node->length;
        double coal_time = time + length;
        
        pll_unode_t *parent = node->back;
        if (nodes_visited.count(parent) == 0) {
            nodes_visited[parent] = coal_time;
            q.push_back(parent);
        }
        
        pop->CoalescentEventTimes.push_back(coal_time);
    }
}

void Chain::initializeCoalescentEventTimesFormSampleSizes(pll_utree_t *utree, vector<int > &sampleSizes )
{
    int totalSampleSize=0;
    for (unsigned int i = 0; i < sampleSizes.size(); ++i)
    {
        totalSampleSize += sampleSizes[i];
    }
    if ( utree->tip_count - 1 != totalSampleSize )
    {
        
        printf("\n The number of tips have to be equal to the total sample size\n");
        exit(1);
    }
    
    initPopulationsTipsFromTree(utree, YES);
    
    for (unsigned int i = 0; i < numClones; ++i)
    {
        auto pop = populations[i];
        updateCoalTimes(pop);
    }
}
void Chain::initPopulationsTipsFromTree(pll_utree_t *utree, bool assignationKnown )
{
    if (utree != NULL)
    {
    for (unsigned int i = 0; i < utree->inner_count + utree->tip_count; ++i)
    {
        auto node = utree->nodes[i];
        
       
        
        if (node->label != 0) // only the tips have labels
        {
            cout << node->label << ": " << node->node_index << ", back: " << node->back->node_index << ", next: " << node->next << endl;
            treeTips.push_back(node);
            if (assignationKnown)
            {
             int popIdx = getPopulationIndex(node->label);
             populations[popIdx]->tips.push_back(node);
            }
        }
        else{
            
            cout << "interior: " << node->node_index << ", " << node->clv_index << ", back: " << node->back->node_index << ", next:" << node->next->node_index << ", nextnext:" << node->next->next->node_index<< endl;
            
        }
      }
    }
}
void Chain::initPopulationsTipsFromRootedTree(pll_rtree_t *rtree, bool assignationKnown )
{
    if (rtree != NULL)
    {
        for (unsigned int i = 0; i < rtree->inner_count + rtree->tip_count; ++i)
        {
            auto node = rtree->nodes[i];
            
            
            
            if (node->label != 0) // only the tips have labels
            {
                cout << node->label << ": " << node->node_index << ", parent: " <<node->parent->node_index   << endl;
                rtreeTips.push_back(node);
                if (assignationKnown)
                {
                    int popIdx = getPopulationIndex(node->label);
                    populations[popIdx]->rtips.push_back(node);
                }
            }
            else{
                if (node->parent!=NULL)
                    cout << "interior: " << node->node_index << ", " << node->clv_index << ", parent: " << node->parent->node_index << ", left:" << node->left->node_index << ", right:" << node->right->node_index<< endl;
                else
                    cout << "interior: " << node->node_index << ", " << node->clv_index << ", parent: " << " " << ", left:" << node->left->node_index << ", right:" << node->right->node_index<< endl;
                
            }
        }
    }
}
void Chain::initializeMapPopulationAssignFromTree()
{
    if (initialUnrootedTree)
    {
        pll_unode_t * node;
        char *nodeLabel;
        std::string nodeLabelString;
        int indexFirstMatch=0;
        int indexPopulation=0;
        Population *pop;
      for (unsigned int i = 0; i < initialUnrootedTree->inner_count + initialUnrootedTree->tip_count; ++i)
     {
        node = initialUnrootedTree->nodes[i];
      //  if (node->back ==NULL && node->next!=NULL && node->next->back ==NULL && node->next->next!=NULL
      //      && node->next->next->back!=NULL)
        if (node->next ==NULL)
        { // first nodelet of a tip!
            nodeLabel =node->label;
            nodeLabelString = std::string(nodeLabel);
            if (nodeLabelString.find("C_"))
            {
                indexFirstMatch= nodeLabelString.find_first_of("C_");
                indexPopulation = nodeLabelString.at(indexFirstMatch +1)-'0';
                pop= getPopulationbyIndex(indexPopulation);
                tipsAssign[node] = pop;
                labelsAssign[nodeLabelString] =node;
//                tipsAssign.insert({ node, pop});
//                labelsAssign.insert({ nodeLabelString, node});
                
            }
         // std::map<pll_unode_t, Population> tipsAssign;
       
        }
      }
    }
    
}
void Chain::initNodeDataFromTree()
{
    if (treeTips.size()>0)
    {
        TreeNode* u1, *u2, *u3;
         deque<pll_unode_t *> q(treeTips.begin(), treeTips.end());
          unordered_map<pll_unode_t *, double> nodes_visited;
   
         for (auto it = treeTips.begin(); it != treeTips.end(); ++it)
          {
            nodes_visited[*it] = 0;
          }
        
        while (q.size() != 0) {
            pll_unode_t *node = q.front();
            q.pop_front();
            double time = nodes_visited[node];
            double length = node->length;
            double coal_time = time + length;
            if (node->data == NULL)
            {
                node->data =  new TreeNode();
                 u1 = (TreeNode *)(node->data);
                 u1->initNumberTipsVector(numClones);
                u1->timePUnits=time;
                printf("updated data for node %d \n", node->node_index);
            }
            if (node->next != NULL &&  node->next->data ==NULL)
            {
                node->next->data =  new TreeNode();
                u2 = (TreeNode *)(node->next->data);
                 u2->initNumberTipsVector(numClones);
                 u2->timePUnits = time;
                printf("updated data for node %d \n", node->next->node_index);
            }
            if (node->next != NULL &&  node->next->next!=NULL && node->next->next->data ==NULL)
            {
                node->next->next->data = new TreeNode();
                u3 = (TreeNode *)(node->next->next->data);
                 u3->initNumberTipsVector(numClones);
                 u3->timePUnits = time;
                printf("updated data for node %d \n", node->next->next->node_index);
            }
           // printf("timePUnits : %lf \n", ((TreeNode *) node->data)->timePUnits);
            
//            if (node->label != 0)
//            {
//                 u->cellName = node->label ;
//                 u->observedCellName = node->label ;
//            }
            pll_unode_t *parent;
         
            if (node->next == NULL)//tip
            {
                if (nodes_visited.count(node->back) == 0)
                {
                    parent = node->back;
                    nodes_visited[parent] = coal_time;
                    q.push_back(parent);
                }
            }
            else {//internal nodelet

                if ( nodes_visited.count(node->next->back) == 0)
                {  parent= node->next->back;
                   nodes_visited[parent] = coal_time;
                    q.push_back(parent);
                }
                if ( nodes_visited.count(node->next->next->back) == 0)
                {
                    parent= node->next->next->back;
                    nodes_visited[parent] = coal_time;
                    q.push_back(parent);
                }
            }
        }
    }
}
void Chain::initNodeDataFromRootedTree()
{
    if (rtreeTips.size()>0)
    {
        TreeNode* u1;
        deque<pll_rnode_t *> q(rtreeTips.begin(), rtreeTips.end());
        unordered_map<pll_rnode_t *, double> nodes_visited;
        
        for (auto it = rtreeTips.begin(); it != rtreeTips.end(); ++it)
        {
            nodes_visited[*it] = 0;
        }
        while (q.size() != 0) {
            pll_rnode_t *node = q.front();
            q.pop_front();
            double time = nodes_visited[node];
            double length = node->length;
            double coal_time = time + length;
            if (node->data == NULL)
            {
                node->data =  new TreeNode();
                u1 = (TreeNode *)(node->data);
                u1->timePUnits=time;
                u1->initNumberTipsVector(numClones);
                printf("updated data for node %d \n", node->node_index);
            }
         
           
            pll_rnode_t *parent = node->parent;
            if (parent!=NULL && nodes_visited.count(parent) == 0) {
                nodes_visited[parent] = coal_time;
                q.push_back(parent);
            }
            
       
        }
    }
}
std::map<pll_unode_t*, Population*>  Chain::chooseTimeOfOriginsOnTree( long int *seed)
{
   
    int numberPoints= numClones -1;
    std::unordered_set<pll_unode_t *>  ancestorsOfTimeOfOrigins(numberPoints);
    TreeNode *u, *v;
    vector<double> cumulativeBranchLengths;
    vector<double> branchLengths;
    pll_tree_edge * edge;
    vector<pll_tree_edge_t *> edges;
    set<pll_tree_edge_t *> visitedEdges;
    set<pll_unode_t *> visitedNodes;
    cumulativeBranchLengths.push_back(0);
    double cumBranchLength=0;
    double random;
    double totalBranchLength1=0;
    std::map<pll_unode_t*, Population*> mrcaOfPopulation;
    
    
    for (unsigned int i = 0; i < initialUnrootedTree->tip_count+ initialUnrootedTree->inner_count; ++i)
    {
        auto node =initialUnrootedTree->nodes[i];
        

        if (node->back && visitedNodes.count(node) == 0  )
        {
            edge = new pll_tree_edge_t();
            u= (TreeNode *)(node->data);
            v= (TreeNode *)(node->back->data);
            edge->edge.utree.parent = u->timePUnits > v->timePUnits ? node :  node->back;
            edge->edge.utree.child= u->timePUnits < v->timePUnits ? node :  node->back;;
            edge->length = node->length;
            visitedNodes.insert(node->back);
            totalBranchLength1 += edge->length;
            
           // if(visitedEdges.count(edge) ==0)
           // {
                branchLengths.push_back(edge->length);
                edges.push_back(edge);
                cumBranchLength=edge->length + cumulativeBranchLengths.at(cumulativeBranchLengths.size()-1);
                cumulativeBranchLengths.push_back(cumBranchLength);
            //}
        }
         visitedNodes.insert(node);
    }
    
    if (edges.size()>0)
    {
        //first normalize
        double maxCumBranchLength = *std::max_element(cumulativeBranchLengths.begin(), cumulativeBranchLengths.end());
        double maxBranchLength = *std::max_element(branchLengths.begin(), branchLengths.end());
//        std::vector<int>::iterator it = std::find(branchLengths.begin(), branchLengths.end(), maxBranchLength);
        branchLengths.erase(std::remove(branchLengths.begin(), branchLengths.end(), maxBranchLength), branchLengths.end());
        
        vector<double>  cumBranchLengths( branchLengths.size());
        std::partial_sum(branchLengths.begin(), branchLengths.end(), cumBranchLengths.begin(), plus<double>());
        
        cumBranchLengths.insert(cumBranchLengths.begin(), 0);
        
        maxCumBranchLength = *std::max_element(cumBranchLengths.begin(), cumBranchLengths.end());
        
        std::transform(cumBranchLengths.begin(), cumBranchLengths.end(), cumBranchLengths.begin(),std::bind2nd(std::divides<double>(),maxCumBranchLength));
        
        double cumBranchLengthsArray[cumBranchLengths.size()];
        std::copy(cumBranchLengths.begin() ,cumBranchLengths.end()
                  ,cumBranchLengthsArray );
        int  nextEvent;
        do{
            //random = RandomUniform(seed);
            random =randomUniformFromGsl();
            nextEvent = bbinClones(random, cumBranchLengthsArray, cumBranchLengths.size());
            pll_unode_t * ancestorMRCA;
            if (edges.at(nextEvent)->length != branchLengths.at(nextEvent -1))
                printf( "something wrong");
                
            ancestorMRCA =edges.at(nextEvent)->edge.utree.parent;
            if (ancestorsOfTimeOfOrigins.count(ancestorMRCA)==0)
            {
                ancestorsOfTimeOfOrigins.insert(ancestorMRCA);
            }
        }
        while(ancestorsOfTimeOfOrigins.size() < numberPoints);
    
//        vector<double> cumProbPop(numClones+1);
//        cumProbPop.push_back(0.0);
//        for (unsigned int i = 1; i <= numClones  ; ++i)
//            cumProbPop.push_back(i* 1.0/numClones);
//
//        double cumProbPop2[numClones+1];
//        std::copy(cumProbPop.begin() ,cumProbPop.end()
//                  ,cumProbPop2 );
//        int popIdx;
        int k=0;
        pll_unode_t * ancestorMRCA;
        for ( auto it = ancestorsOfTimeOfOrigins.begin(); it != ancestorsOfTimeOfOrigins.end(); ++it )
        {
             ancestorMRCA = *it;
            auto pop=getPopulationbyIndex(k);
            pop->nodeletMRCA = ancestorMRCA->back;
            
            mrcaOfPopulation[ancestorMRCA->back]=pop;
            u= (TreeNode *)(ancestorMRCA->data);
            v= (TreeNode *)(ancestorMRCA->back->data);
            double proposedTime= v->timePUnits+ (u->timePUnits- v->timePUnits)*randomUniformFromGsl();
            if (proposedTime<=0)
                printf( "time of origin is positive");
            pop->timeOriginInput =proposedTime;
            k++;
            
        }
        Population* pop=getPopulationbyIndex(k);
        pop->nodeletMRCA =initialUnrootedTree->nodes[initialUnrootedTree->tip_count+ initialUnrootedTree->inner_count-1];//the last population has MRCA the root of the tree
        mrcaOfPopulation[initialUnrootedTree->nodes[initialUnrootedTree->tip_count+ initialUnrootedTree->inner_count-1]]=pop;
    }
   return(mrcaOfPopulation);
}
void Chain::initBranches(string& healthyCellLabel,vector<double> &branchLengths, vector<pll_tree_edge_t *> &edges)
{
    TreeNode *u, *v;
    vector<double> cumulativeBranchLengths;
    pll_tree_edge * edge;
    //vector<pll_tree_edge_t *> edges;
    set<pll_tree_edge_t *> visitedEdges;
    set<pll_rnode_t *> visitedNodes;
    cumulativeBranchLengths.push_back(0);
    double cumBranchLength=0;
    double totalBranchLength1=0;
    int totalNumberNodes  =initialRootedTree->tip_count+ initialRootedTree->inner_count;
    for (unsigned int i = 0; i < totalNumberNodes; ++i)
    {
        auto node =initialRootedTree->nodes[i];
        if (node->parent !=NULL && visitedNodes.count(node) == 0  )
        {
            edge = new pll_tree_edge_t();
            u= (TreeNode *)(node->data);
            v= (TreeNode *)(node->parent->data);
            edge->edge.rtree.parent = u->timePUnits > v->timePUnits ? node :  node->parent;
            edge->edge.rtree.child= u->timePUnits < v->timePUnits ? node :  node->parent;;
            edge->length = node->length;
            //visitedNodes.insert(node->parent);
            totalBranchLength1 += edge->length;
            // if(visitedEdges.count(edge) ==0)
            // {
            if (node->label && std::string(node->label).compare(healthyCellLabel)!=0)//not  the branch from the healthyTip
            {
                branchLengths.push_back(edge->length);
            }
            else if(!node->label && node->parent->parent !=NULL)//not the branch to the root
            {
                branchLengths.push_back(edge->length);
            }
            edges.push_back(edge);
            cumBranchLength = edge->length + cumulativeBranchLengths.at(cumulativeBranchLengths.size()-1);
            cumulativeBranchLengths.push_back(cumBranchLength);
            //}
        }
        visitedNodes.insert(node);
    }
}
std::map<pll_rnode_t*, Population*>  Chain::initTimeOfOriginsOnRootedTree( int numberPoints, string &healthyCellLabel)
{
    std::unordered_set<pll_rnode_t *>  ancestorsOfTimeOfOrigins(numberPoints);
    std::unordered_set<pll_rnode_t *>  MRCAs(numberPoints);
    TreeNode *u, *v;
    //vector<double> cumulativeBranchLengths;
    vector<double> branchLengths;
    //vector<pll_tree_edge_t *> edges;
    set<pll_tree_edge_t *> visitedEdges;
    set<pll_rnode_t *> visitedNodes;
    //cumulativeBranchLengths.push_back(0);
    double random;
    std::map<pll_rnode_t*, Population*> mrcaOfPopulation;
    
     initBranches(healthyCellLabel, branchLengths, edges);
    
    if (branchLengths.size()>0)
    {
      //double maxCumBranchLength = *std::max_element(cumulativeBranchLengths.begin(), cumulativeBranchLengths.end());
         //first erase the branch from the healthyTip
       // double maxBranchLength = *std::max_element(branchLengths.begin(), branchLengths.end());
        //        std::vector<int>::iterator it = std::find(branchLengths.begin(), branchLengths.end(), maxBranchLength);
     //   branchLengths.erase(std::remove(branchLengths.begin(), branchLengths.end(), maxBranchLength), branchLengths.end());
        double maxCumBranchLength ;
        vector<double>  cumBranchLengths( branchLengths.size());
        std::partial_sum(branchLengths.begin(), branchLengths.end(), cumBranchLengths.begin(), plus<double>());
        
        cumBranchLengths.insert(cumBranchLengths.begin(), 0);
        
        maxCumBranchLength = *std::max_element(cumBranchLengths.begin(), cumBranchLengths.end());
        
        std::transform(cumBranchLengths.begin(), cumBranchLengths.end(), cumBranchLengths.begin(),std::bind2nd(std::divides<double>(),maxCumBranchLength));
        
        double cumBranchLengthsArray[cumBranchLengths.size()];
        std::copy(cumBranchLengths.begin() ,cumBranchLengths.end()
                  ,cumBranchLengthsArray );
        int  nextEvent;
        vector<int> eventIds;
        pll_rnode_t * ancestorMRCA;
        pll_rnode_t * MRCA;
        do{
            //random = RandomUniform(seed);
            do{
            random =randomUniformFromGsl();
            nextEvent = bbinClones(random, cumBranchLengthsArray, cumBranchLengths.size());
            }
            while(std::find(eventIds.begin(), eventIds.end(), nextEvent) != eventIds.end());
            
            if (edges.at(nextEvent)->length != branchLengths.at(nextEvent -1))
                printf( "the lengths of the branches %d does not match: %lf, %lf ",nextEvent, edges.at(nextEvent)->length, branchLengths.at(nextEvent -1));
            
            ancestorMRCA =edges.at(nextEvent)->edge.rtree.parent;
            MRCA =edges.at(nextEvent)->edge.rtree.child;
           // if (ancestorsOfTimeOfOrigins.count(ancestorMRCA)==0)
             if (ancestorsOfTimeOfOrigins.count(ancestorMRCA)==0)//the parents of the 2 events are differents
            {
                eventIds.push_back(nextEvent);
                ancestorsOfTimeOfOrigins.insert(ancestorMRCA);
                MRCAs.insert(MRCA);
            }
            else if(ancestorsOfTimeOfOrigins.count(ancestorMRCA)>0 && ancestorMRCA->parent->parent!= NULL)//the parents of the 2 events are igual, but this same parent is not the left child of the root
            {
                eventIds.push_back(nextEvent);
                ancestorsOfTimeOfOrigins.insert(ancestorMRCA);
                 MRCAs.insert(MRCA);
            }
        }
        while(ancestorsOfTimeOfOrigins.size() < numberPoints);
        
        int k=0;
        
        for ( auto it = MRCAs.begin(); it != MRCAs.end(); ++it )
        {
            MRCA = *it;
            auto pop=getPopulationbyIndex(k);
            pop->rMRCA = MRCA;
            mrcaOfPopulation[MRCA]=pop;
            u= (TreeNode *)(MRCA->data);
            v= (TreeNode *)(MRCA->parent->data);
            double proposedTime= u->timePUnits+ (v->timePUnits- u->timePUnits)*randomUniformFromGsl();
            if (proposedTime<=0)
                printf( "time of origin is positive");
            pop->timeOriginInput =proposedTime;
            printf( "\n MRCA node id %d with time %lf anf parent with time %lf was assigned to pop %d and time of origin %lf \n", MRCA->node_index, u->timePUnits, v->timePUnits, pop->index, pop->timeOriginInput  );
            k++;
        }
        Population* pop=getPopulationbyIndex(k);
        pop->rMRCA =initialRootedTree->root;//the last population has MRCA node the root  of the tree
        u= (TreeNode *)(initialRootedTree->root->data);
        pop->timeOriginInput=u->timePUnits;
        mrcaOfPopulation[pop->rMRCA]=pop;
    }
    return(mrcaOfPopulation);
}

void Chain::initPopulationCoalescentAndMigrationEventsFromRootNodeOnTree(pll_unode_t *p, Population *currentPopulation, std::map<pll_unode_t*, Population*> mrcaOfPopulation )
{
    if (p!=NULL)
    {
        if (!p->next)
        {
            currentPopulation->tips.push_back(p);
        }
        else if(p->next->back !=NULL && mrcaOfPopulation.count(p->next->back) != 0 )
        { // in this branch we have  migration
            TreeNode * u= (TreeNode *)p->data;
            currentPopulation->CoalescentEventTimes.push_back(u->timePUnits );
            currentPopulation->UpdateListMigrants(numClones, currentPopulation, mrcaOfPopulation[p->next->back]);
            
            initPopulationCoalescentAndMigrationEventsFromRootNodeOnTree(p->next->back, mrcaOfPopulation[p->next->back], mrcaOfPopulation);
        }
        else if(p->next->next->back !=NULL && mrcaOfPopulation.count(p->next->next->back) != 0 )
        { // in this branch we have  migration
            currentPopulation->UpdateListMigrants(numClones, currentPopulation, mrcaOfPopulation[p->next->next->back]);
            
            initPopulationCoalescentAndMigrationEventsFromRootNodeOnTree(p->next->next->back, mrcaOfPopulation[p->next->next->back], mrcaOfPopulation);
        }
        else{
            initPopulationCoalescentAndMigrationEventsFromRootNodeOnTree(p->next->back, currentPopulation, mrcaOfPopulation);
            
            initPopulationCoalescentAndMigrationEventsFromRootNodeOnTree(p->next->next->back, currentPopulation, mrcaOfPopulation);
        }
    }
}
void Chain::initPopulationsSampleSizes( std::map<pll_rnode_t*, Population*>  rmrcaOfPopulation){
    
      if (rmrcaOfPopulation.count(initialRootedTree->root)==0)
          printf("Error, the root node is not the MRCA of any population");
    
    initPopulationSampleSizesFromNodeOnRootedTree(initialRootedTree->root, rmrcaOfPopulation[initialRootedTree->root], rmrcaOfPopulation);
   
    TreeNode *treeNode;
     for (unsigned int i = 0; i < numClones; ++i)
        {
            auto pop=populations[i];
            treeNode = (TreeNode *)(pop->rMRCA->data);
            if (pop->rMRCA->parent !=NULL)//not the oldest population
            { pop->sampleSize =treeNode->numberTipsByPopulation[i];
            }
            else{
                if (treeNode->numberTipsByPopulation[i] -1 >0)
                  pop->sampleSize =treeNode->numberTipsByPopulation[i] -1;
                else
                   printf("Error, the last population must have at least one tip");
            }
        }

}

void Chain::initPopulationsCoalescentAndMigrationEvents(std::map<pll_unode_t*, Population*> mrcaOfPopulation ){
    
    for (unsigned int i = 0; i < numClones; ++i)
    {
        auto pop = getPopulationbyIndex(i);
        pop->resetMigrationsList();//this requires  that you have set first timeOriginSTD in the population and effective population size
        initPopulationCoalescentAndMigrationEventsFromRootNodeOnTree(pop->nodeletMRCA, pop, mrcaOfPopulation );
    }
}
void Chain::initPopulationsCoalescentAndMigrationEventsFromRootedTree(std::map<pll_rnode_t*, Population*> rmrcaOfPopulation, string& healthyTipLabel){
    
    for (unsigned int i = 0; i < numClones; ++i)
    {
        auto pop = getPopulationbyIndex(i);
        pop->resetMigrationsList();//this requires  that you have set first timeOriginSTD in the population and effective population size
        //pop->InitCoalescentEvents(numClones);
    }
    initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(rootRootedTree, rmrcaOfPopulation[rootRootedTree], rmrcaOfPopulation , healthyTipLabel);
}
void Chain::initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(pll_rnode_t *p, Population *currentPopulation, std::map<pll_rnode_t*, Population*> rmrcaOfPopulation, string& healthyTipLabel ){
    
    if (p!=NULL)
    {
        if (p->left == NULL && std::string(p->label).compare(healthyTipLabel)!=0 )//tip different from healthy tip
        {
            currentPopulation->rtips.push_back(p);
        }
        else
        {
            TreeNode * u= (TreeNode *)p->data;
            currentPopulation->CoalescentEventTimes.push_back(u->timePUnits / currentPopulation->effectPopSize  );
            if(p->left !=NULL && rmrcaOfPopulation.count(p->left) != 0 )
             { // in this branch we have  migration
             
               currentPopulation->UpdateListMigrants(numClones, rmrcaOfPopulation[p->left], currentPopulation);
                 rmrcaOfPopulation[p->left]->FatherPop=currentPopulation;
                 initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(p->left, rmrcaOfPopulation[p->left], rmrcaOfPopulation, healthyTipLabel);
             }
            else{
                initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(p->left, currentPopulation, rmrcaOfPopulation, healthyTipLabel);
              }
            if(p->right !=NULL && rmrcaOfPopulation.count(p->right) != 0 )
             { // in this branch we have  migration
              currentPopulation->UpdateListMigrants(numClones,  rmrcaOfPopulation[p->right], currentPopulation);
           
                rmrcaOfPopulation[p->right]->FatherPop=currentPopulation;
                 initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(p->right, rmrcaOfPopulation[p->right], rmrcaOfPopulation, healthyTipLabel);
              }
            else{
           initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(p->right, currentPopulation, rmrcaOfPopulation, healthyTipLabel);
               }
        }
    }
}
//void Chain::initPopulationSampleSizesFromRootNodeOnTree(pll_unode_t *p, Population *population )
//{
//
//    if (p!=NULL)
//    {
//        if (p->next == NULL)
//        {
//            population->sampleSize=  population->sampleSize + 1;
//        }
//
//        else{
//            TreeNode *pData = (TreeNode *)p->data;
//            TreeNode *datapNextBack = (TreeNode *)p->next->back->data;
//            TreeNode *datapNextNextBack = (TreeNode *)p->next->next->back->data;
//            if (datapNextBack->timePUnits <= pData->timePUnits)
//                initPopulationSampleSizesFromRootNodeOnTree(p->next->back, population);
//            if (datapNextNextBack->timePUnits <= pData->timePUnits)
//                initPopulationSampleSizesFromRootNodeOnTree(p->next->next->back, population);
//        }
//
//    }
//}
void Chain::initPopulationSampleSizesFromNodeOnRootedTree(pll_rnode_t *p, Population *currentPopulation, std::map<pll_rnode_t*, Population*>  rmrcaOfPopulation)
{
    if (p!=NULL)
    {
        TreeNode * treeNode=(TreeNode *)(p->data);
        if (p->left == NULL)//tip
        {
            
            int idx= currentPopulation->index;
            if (idx <= numClones -1)
            {
                treeNode->numberTipsByPopulation.at(currentPopulation->index) = 1;
                return;
            }
            else{
                printf("error the index is out of bounds");
            }
        }
        else
        {
            if (rmrcaOfPopulation.count(p->left)>0 )
           {
            std::map<pll_rnode_t*, Population*>::iterator it;
            it = rmrcaOfPopulation.find(p->left);
            if (rmrcaOfPopulation.count(p->left) != 1)
                printf("error the MRCA(left) node appear more than once");
            initPopulationSampleSizesFromNodeOnRootedTree(p->left, rmrcaOfPopulation[p->left], rmrcaOfPopulation );
           }
            else{
               initPopulationSampleSizesFromNodeOnRootedTree(p->left, currentPopulation, rmrcaOfPopulation );
               }
        
           if (rmrcaOfPopulation.count(p->right)>0 )
           {
            std::map<pll_rnode_t*, Population*>::iterator it;
            it = rmrcaOfPopulation.find(p->right);
            if (rmrcaOfPopulation.count(p->right) != 1)
                printf("error the MRCA(right) node appear more than once");
            initPopulationSampleSizesFromNodeOnRootedTree(p->right, rmrcaOfPopulation[p->right], rmrcaOfPopulation );
            }
           else
           {
            initPopulationSampleSizesFromNodeOnRootedTree(p->right, currentPopulation, rmrcaOfPopulation);
            }
             TreeNode * treeNodeLeft=(TreeNode *)(p->left->data);
             TreeNode * treeNodeRight=(TreeNode *)(p->right->data);
            std::transform (treeNodeLeft->numberTipsByPopulation.begin(), treeNodeLeft->numberTipsByPopulation.end(), treeNodeRight->numberTipsByPopulation.begin(), treeNode->numberTipsByPopulation.begin(), std::plus<int>());
            
            printf("node %d has %d tips below \n", p->node_index,treeNode->numberOfTipsSubTree );
            printf("node %d has %d accumulated tips below \n", p->node_index,accumulate(treeNode->numberTipsByPopulation.begin(),treeNode->numberTipsByPopulation.end(),0));
        }
    }
}
Population * Chain::getPopulationbyIndex(int indexPopulation)
{
    Population *pop;
    for (unsigned int i = 0; i < numClones; ++i){
        pop = populations[i];
        if (pop->index ==indexPopulation)
            return pop;
    }
    return NULL;
}
void Chain::filterSortPopulationsCoalescentEvents()
{
     Population *pop;
    for (unsigned int i = 0; i < numClones; ++i){
        pop = populations[i];
        pop->filterAndSortCoalescentEvents();
    }
}
void Chain::samplePopulationDeltaFromPriors(MCMCoptions &mcmcOptions, long int *seed )
{
    int i;
    Population *popI;
    double randomDelta;
    for( i = 0 ; i < numClones; i++)
    {
            popI=populations[i];
            randomDelta = RandomLogUniform(mcmcOptions.Deltafrom, mcmcOptions.Deltato);
            //popI->delta = chain->proportionsVector[i] * randomDelta;
            popI->delta =  randomDelta;
            popI->growthRate =popI->delta  / popI->effectPopSize;
            popI->popSize=popI->effectPopSize * popI->birthRate;
            popI->deathRate= popI->birthRate - popI->growthRate;
    }
}
void Chain::rescaleRootedTreeBranchLengths(double mutationRate)
{
    if (mutationRate <=0)
        printf("Error, the mutation rate must be positive \n");
    if (initialRootedTree !=NULL )
    {
        for (unsigned int i = 0; i < initialRootedTree->inner_count + initialRootedTree->tip_count; ++i)
        {
            auto node = initialRootedTree->nodes[i];
            if (node->parent!=NULL)
            {
                node->length = node->length / mutationRate;
            }
        }
    }
}
void Chain::newTotalEffectivePopulationSizeMove( ProgramOptions &programOptions, char *ObservedCellNames[], pll_msa_t * msa, MCMCoptions &mcmcOptions, vector<int> &sampleSizes)
{
    int i,j;
    Population *popJ;
    double totalTreeLength;
    int        numCA, numMIG;
    double      cumNumCA, meanNumCA, cumNumMIG, meanNumMIG, numEventsTot;
 
    Population *popI;
    // save the current values of
    oldTotalEffectPopSize = totalEffectPopSize;
    double newTotalEffectPopulationSize= proposalSlidingWindow(oldTotalEffectPopSize,  mcmcOptions.slidingWindowSizeTotalEffectPopSize);
    
    totalEffectPopSize = (int)newTotalEffectPopulationSize;
    // RandomLogUniform(mcmcOptions->totalEffectPopSizefrom,mcmcOptions->totalEffectPopSizeto);
    // chain->totalPopSize= RandomLogUniform(7,13, seed);
   // GenerateEffectPopSizesFromPriors2(programOptions.noisy, NO);
    
    updateEffectPopSizesCurrentProportionsVector();
  
    //setChainPopulationSampleSizes(chain, sampleSizes, programOptions);
    //InitPopulationsCoalescentEvents( numClones,  populations) ;
    
    if (programOptions.populationSampleSizesKnown ==YES)
    {
        for( i = 0 ; i < numClones; i++)
        {
            popI=populations[i];
            do
            {
                popI->growthRate =popI->delta  / popI->effectPopSize;
                popI->popSize=popI->effectPopSize * popI->birthRate;
                popI->deathRate= popI->birthRate - popI->growthRate;
            }
            while(popI->popSize < popI->sampleSize);
        }
    }
    ListClonesAccordingTimeToOrigin(populations);
    
    double newLogConditionalLikelihoodTree= LogConditionalLikelihoodTree(programOptions);
    
    fprintf (stderr, "\n>> New log conditional Likelihood tree after the new total effective population size,  of the chain %d is = %lf  \n", chainNumber,newLogConditionalLikelihoodTree );
//    char *newickString2;
//    newickString2=NULL;
//    newickString2 = toNewickString2 ( oldroot, programOptions.mutationRate,     programOptions.doUseObservedCellNames);
//    printf("\n newick after move= %s  \n", newickString2);
//    newickString2 = toNewickString4 ( oldroot, programOptions.mutationRate,     programOptions.doUseObservedCellNames);
//    printf("\n newick after move= %s  \n", newickString2);
//    free(newickString2);
//    newickString2=NULL;
    double priorDensityNewTotalEffectivePopulationSize= LogUniformDensity(newTotalEffectPopulationSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    double priorDensityCurrentTotalEffectivePopulationSize= LogUniformDensity(oldTotalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewTotalEffectivePopulationSize;
    
    double sumLogDenominators=currentlogConditionalLikelihoodTree +priorDensityCurrentTotalEffectivePopulationSize;
    
    double randomNumber= randomUniformFromGsl();
    
    double LogAcceptanceRate = (sumLogNumerators - sumLogDenominators) >0? (sumLogNumerators - sumLogDenominators) :0;
    if (log(randomNumber) < LogAcceptanceRate )
    {//accept the move
        printf("\n Accepted new total effective population size move\n");
    }
    else
    {
        //reject the move
        printf("\n Rejected new total effective population size move\n");
        totalEffectPopSize = oldTotalEffectPopSize;
        for (j = 0; j < numClones; j++)
        {
            popJ = populations[j];
            popJ->effectPopSize = popJ->oldeffectPopSize;
        }
    }
}
/********************* proposalSlidingWindow **********************/
/* proposalSlidingWindow*/
double  Chain::proposalSlidingWindow( double oldvalue,  double windowSize)
{
    double newvalue =0;
    newvalue = oldvalue + (randomUniformFromGsl()-0.5) * windowSize ;
    if (newvalue <0)
        newvalue = -  newvalue;
    return newvalue;
}
void Chain::newProportionsVectorMove(ProgramOptions &programOptions, char *ObservedCellNames[], pll_msa_t * msa, MCMCoptions &mcmcOptions, vector<int> &sampleSizes)
{
    int i,j,k;
    Population *popJ;
    
    Population *popI;
    proposalProportionsVector(oldproportionsVector, mcmcOptions.tuningParameter);
    
   // GenerateEffectPopSizesFromPriors2(chain, programOptions->noisy,  chain->numClones, chain->populations,  seed,  chain->gammaParam, chain->totalPopSize, YES);
    if (programOptions.populationSampleSizesKnown == NO)
    {
        //InitPopulationSampleSizes(populations, programOptions.numCells, programOptions.numClones, proportionsVector, seed);
    }
    //else fill the sample  sizes
    else{
        setChainPopulationSampleSizes(sampleSizes, programOptions);
    }
    if (programOptions.populationSampleSizesKnown ==YES)
    {
        for( i = 0 ; i < numClones; i++)
        {
            popI=populations [i];
            do {
                popI->growthRate =popI->delta  / popI->effectPopSize;
                popI->popSize=popI->effectPopSize * popI->birthRate;
                popI->deathRate= popI->birthRate - popI->growthRate;
            }
            while(popI->popSize < popI->sampleSize);
        }
    }
    ListClonesAccordingTimeToOrigin(populations);
    double newLogConditionalLikelihoodTree= LogConditionalLikelihoodTree(programOptions);
    fprintf (stderr, "\n>> New log conditional Likelihood tree after the new total effective population size,  of the chain %d is = %lf  \n", chainNumber,newLogConditionalLikelihoodTree );
//    char *newickString2;
//    newickString2=NULL;
//    newickString2 = toNewickString2 ( chain->oldroot, programOptions->mutationRate,     programOptions->doUseObservedCellNames);
//    printf("\n newick after move= %s  \n", newickString2);
//    newickString2 = toNewickString4 ( chain->oldroot, programOptions->mutationRate,     programOptions->doUseObservedCellNames);
//    printf("\n newick after move= %s  \n", newickString2);
//    free(newickString2);
//    newickString2=NULL;
    double priorDensityNewProportionsVector =  DirichletDensity( proportionsVector, oldproportionsVector , numClones);
    
    double priorDensityCurrentProportionsVector = DirichletDensity(oldproportionsVector, proportionsVector, numClones);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewProportionsVector;
    
    double sumLogDenominators=currentlogConditionalLikelihoodTree +priorDensityCurrentProportionsVector;
    
    double randomNumber= randomUniformFromGsl();
    
    double LogAcceptanceRate = (sumLogNumerators - sumLogDenominators) >0? (sumLogNumerators - sumLogDenominators) :0;
    
    if (log(randomNumber) < LogAcceptanceRate )
    {//accept the move
        printf("\n Accepted new proportions vector move\n");
//        nodes = oldnodes;
//        root = oldroot;
        currentlogConditionalLikelihoodTree =newLogConditionalLikelihoodTree;
    }
    else
    {
        //reject the move
        printf("\n Rejected new proportions vector move\n");
       totalEffectPopSize = oldTotalEffectPopSize;
        for (j = 0; j <numClones; j++)
        {
            popJ =populations[j];
            popJ->effectPopSize = popJ->oldeffectPopSize;
        }
    }
}
void Chain::proposalProportionsVector(vector<double > &newProportionsvector, double tuningParameter )
{
    int i=0;
    double vectorForGamma[numClones];
    for( i = 0 ; i < numClones; i++)
    {
        vectorForGamma[i]= tuningParameter * proportionsVector[i];
    }
    double *proportionsVectorArray=&proportionsVector[0];
 randomDirichletFromGsl(numClones, proportionsVectorArray, &newProportionsvector[0]);

}
double Chain::DirichletDensity(vector<double> &proportionsVector,  vector<double> &concentrationVector, int sizeVector)
{
    int i;
    double sum=0;
    double logResult=0;
    for( i = 0 ; i < sizeVector; i++)
    {
        sum = sum +concentrationVector[i];
        logResult= logResult+(concentrationVector[i]-1)*log(proportionsVector[i]);
        logResult= logResult-lgamma(concentrationVector[i]);
    }
    logResult = logResult+lgamma(concentrationVector[i]);
    return logResult;
}
double * Chain::expand_uniq_rates(int states, const double * uniq_rates, const int * rate_sym)
{
    unsigned int i;
    
    unsigned int num_rates = states * (states-1) / 2;
    //double * subst_rates = calloc(num_rates, sizeof(double));
    double * subst_rates =(double *) calloc(num_rates, sizeof(double));
    for (i = 0; i < num_rates; ++i)
        subst_rates[i] = rate_sym ? uniq_rates[rate_sym[i]] : uniq_rates[i];
    
    return subst_rates;
}
int Chain::totalSampleSize()
{    int i;
    Population *popI;
    int totalSampleSize=0;
    for( i = 0 ; i < numClones; i++){
        popI = populations[i];
        totalSampleSize+=popI->sampleSize;
    }
    return totalSampleSize;
}
 std::map<pll_rnode_t*, Population*> Chain::chooseAvailableEdgeOnRootedTreeForPopulation(Population *pop, std::map<pll_rnode_t*, Population*> &mrcaOfPopulation, string &healthyCellLabel)
{
    std::unordered_set<pll_rnode_t *>  ancestorsOfTimeOfOrigins;
    std::unordered_set<pll_rnode_t *>  MRCAs;
    TreeNode *u, *v;
    //vector<double> cumulativeBranchLengths;
    vector<double> branchLengths;
    vector<pll_tree_edge_t *> availableEdges;
    set<pll_tree_edge_t *> visitedEdges;
    set<pll_rnode_t *> visitedNodes;
    //cumulativeBranchLengths.push_back(0);
    double random;
    pll_tree_edge_t *edge;
    pll_rnode_t *parent, *child;
    double cumBranchLength;
    
    for (unsigned int i = 0; i < edges.size(); ++i)
    {
        auto edge =edges.at(i);
        parent = edge->edge.rtree.parent ;
        child = edge->edge.rtree.child;
        if (parent !=NULL && edge->length >0 && mrcaOfPopulation.count(child) == 0  )
        {
            //totalBranchLength1 += edge->length;
            availableEdges.push_back(edge);
            branchLengths.push_back(edge->length);
//            cumBranchLength = edge->length + cumulativeBranchLengths.at(cumulativeBranchLengths.size()-1);
//            cumulativeBranchLengths.push_back(cumBranchLength);
        }
    }
    std::map<pll_rnode_t*, Population*> copyMRCAOfPopulation;
    ////////////////////////////////////
    if (availableEdges.size()>0)
    {
        double maxCumBranchLength ;
        vector<double>  cumBranchLengths( branchLengths.size());
        std::partial_sum(branchLengths.begin(), branchLengths.end(), cumBranchLengths.begin(), plus<double>());
        
        cumBranchLengths.insert(cumBranchLengths.begin(), 0);
        
        maxCumBranchLength = *std::max_element(cumBranchLengths.begin(), cumBranchLengths.end());
        
        std::transform(cumBranchLengths.begin(), cumBranchLengths.end(), cumBranchLengths.begin(),std::bind2nd(std::divides<double>(),maxCumBranchLength));
        
        double cumBranchLengthsArray[cumBranchLengths.size()];
        std::copy(cumBranchLengths.begin() ,cumBranchLengths.end()
                  ,cumBranchLengthsArray );
        int  nextEvent;
        vector<int> eventIds;
        pll_rnode_t * ancestorMRCA;
        pll_rnode_t * MRCA;
        bool existsZeroSampleSizePop=true;
        double alpha[numClones];
       // copyMRCAOfPopulation.insert(mrcaOfPopulation.begin(), mrcaOfPopulation.end());
        do
        {
            copyMRCAOfPopulation.clear();
            copyMRCAOfPopulation.insert(mrcaOfPopulation.begin(), mrcaOfPopulation.end());
            do
            {
                random =randomUniformFromGsl();
                nextEvent = bbinClones(random, cumBranchLengthsArray, cumBranchLengths.size());
            }
            while(std::find(eventIds.begin(), eventIds.end(), nextEvent) != eventIds.end());
//            if (edges.at(nextEvent)->length != branchLengths.at(nextEvent -1))
//                printf( "the lengths of the branches %d does not match: %lf, %lf ",nextEvent, edges.at(nextEvent)->length, branchLengths.at(nextEvent -1));
            ancestorMRCA =edges.at(nextEvent)->edge.rtree.parent;
            MRCA =edges.at(nextEvent)->edge.rtree.child;
            copyMRCAOfPopulation[MRCA]= pop;
            initPopulationsSampleSizes( copyMRCAOfPopulation);
            existsZeroSampleSizePop= (pop->sampleSize == 0)? true: false;
        }
        while(existsZeroSampleSizePop);
    //////////////////////////////////////
            pop->rMRCA = MRCA;
            mrcaOfPopulation[MRCA]=pop;
            u= (TreeNode *)(MRCA->data);
            v= (TreeNode *)(MRCA->parent->data);
            double proposedTime= u->timePUnits+ (v->timePUnits- u->timePUnits)*randomUniformFromGsl();
            if (proposedTime<=0)
                printf( "time of origin must be positive");
            pop->timeOriginInput =proposedTime;
            printf( "\n MRCA node id %d with time %lf anf parent with time %lf was assigned to pop %d and time of origin %lf \n", MRCA->node_index, u->timePUnits, v->timePUnits, pop->index, pop->timeOriginInput  );
    }
    return copyMRCAOfPopulation;
}
double Chain::sumAvailableBranchLengths(std::map<pll_rnode_t*, Population*> currentMRCAPopulation)
{
    double result=0.0;
    for (unsigned int i = 0; i < edges.size(); ++i)
    {
        auto edge =edges.at(i);
        auto parent = edge->edge.rtree.parent ;
        auto child = edge->edge.rtree.child;
        if (parent !=NULL && edge->length >0 && currentMRCAPopulation.count(child) != 0  )
        {
            result+= result;
        }
    }
    return result;
}
