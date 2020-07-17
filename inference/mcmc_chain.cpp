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
 * MCMC chain
 */
//#include "mcmc_chain.hpp"
#include <map>
#include <iterator>
#include <string>
#include <vector>
#include <queue>
#include <unordered_set>
#include <numeric>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <functional>
#include <random>
#include <chrono>
#include <set>
#include <boost/regex.h>
#include <boost/random.hpp>

#include <sys/stat.h>

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
#include "utils.hpp"
#include "constants.hpp"
#include "mcmc_move.hpp"

using namespace std;
using namespace std::placeholders;
//Parametrized Constructor

Chain::Chain( int chainNumber,
             int numClones,
             int gammaParam,
             long double  mutationRate,
             long double  seqErrorRate,
             long double  dropoutRate
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
    totalAccepted=0;
    totalRejected=0;
    
    theta = 0.0;
    oldtheta = 0.0;
    initProportionsVector();
   
    
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

void Chain::MakeCoalescenceEvent( Population *population, vector<pll_unode_t *> &nodes, int numClones, long int* seed, int noisy,   int &numActiveGametes, int &nextAvailable,
                                 int &labelNodes, long double  &currentTime, int &numNodes)
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
    
    r1->data =  new TreeNode(0);;
    r2->data =  new TreeNode(0);
    r3->data = new TreeNode(0);;
    
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
        fprintf (stderr, "\t|\tCurrentTime (input units) = %Lf", u1->timePUnits);
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
        
    }
    
}

void Chain::SimulatePopulation( Population *popI, vector<Population*> &populations,
                               ProgramOptions &programOptions,
                               long int *seed,
                               int &numNodes,
                               int numClones,
                               long double       cumNumCA,
                               long double  meanNumCA,
                               long double  cumNumMIG,
                               long double  meanNumMIG,
                               int  &numMIG,
                               int  &numCA,
                               long double  &numEventsTot,
                               vector<pll_unode_t *> &nodes,
                               int &nextAvailable,
                               int &numActiveGametes,
                               int &labelNodes,
                               long double  &currentTime,
                               int &eventNum)
{
    int   i,  k, isCoalescence,
    firstInd,  newInd,
    isMigration, whichClone;
    long double      eventTime;
    pll_unode_t  *p, *r, *r1, *r2, *r3 ;
    TreeNode  *u,*v, *u1, *u2, *u3;
    
    eventNum = 0;
    
    long double  ThisRateCA = 0.0;
    long double  ThisTimeCA_W = 0.0;
    long double  ThisTimeCA_V1 = 0.0;
    long double  ThisTimeCA_V2 = 0.0;
    int numParcialActiveGametes;
    
    
    numParcialActiveGametes = popI->numActiveGametes;
    
    int numMigrations = popI->numIncomingMigrations; //taking into account also the    time of origin as a migration
    
    long double  timeNextMigration;
    int indexNextMigration = 0;
    Population *incomingPop;
    
    if (programOptions.noisy > 1)
        fprintf (stderr, "\n\n>> Simulating evolutionary history of clone %d (number active gametes %d, original time to origin %Lf)\n", popI->index, popI->numActiveGametes, popI->timeOriginInput);
    
    if (programOptions.noisy > 1)
        fprintf (stderr, "\n\n> Simulating evolutionary history of clone  or order  %d ..\n", popI->order);
    
    currentTime=0;
    while (indexNextMigration < numMigrations) {
        timeNextMigration = popI->immigrantsPopOrderedByModelTime[indexNextMigration].first;
        //fprintf (stderr, "\n\n> numParcialActiveGametes= %d \n", numParcialActiveGametes);
        if ( popI->numActiveGametes >= 2) {
            ThisRateCA = (double)  popI->numActiveGametes * ((double) popI->numActiveGametes - 1) / 2.0;
            ThisTimeCA_W = Random::RandomExponential (ThisRateCA, seed, true) ;
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
                fprintf (stderr, "\n\n*** Event %3d *** currentTime (model time) = %Lf, currentTime (standard time) = %Lf\n", eventNum, ThisTimeCA_V2, ThisTimeCA_V1 );
                fprintf (stderr, "\n\n*** Event %3d *** currentTime (input units) = %Lf\n", eventNum, ThisTimeCA_V2);
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
                    fprintf (stderr, "\n\n*** Event %3d *** *currentTime (model units) = %Lf\n", eventNum, ThisTimeCA_V2);
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
                
                printf( "\n The incoming population %d to  population %d with node %d and nodelet %d time %Lf", incomingPop->order, popI->order, v->index, p->node_index, v->timePUnits );
                
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
                    fprintf (stderr, "\t|\tCurrentTime (input units) = %Lf", v->timePUnits);
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
                    fprintf (stderr, "Clone origin %d at time (model units) = %Lf\n", popI->index, currentTime);
                if (popI->order < numClones - 1) // do not do it for the last clone
                {
                    newInd = nextAvailable;
                    
                    r1 = nodes[newInd];   /* new nodelet ancester */
                    r2 = nodes[newInd + 1];   /* new  nodelet ancester */
                    r3 = nodes[newInd + 2];   /* new nodelet ancester */
                    
                    r1->node_index = newInd;
                    r2->node_index = newInd +1;
                    r3->node_index = newInd +2;
                    
                    r1->data =  new TreeNode(0);
                    r2->data =   new TreeNode(0);
                    r3->data =   new TreeNode(0);
                    
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
                        fprintf (stderr, "\t|\tCurrentTime (input units) = %Lf", u1->timePUnits);
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
                              long double   &currentTime,
                              int &labelNodes)
{
    int i, j;
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
        
        // int transformingBranchLength=1.001;
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
        
        healthyTip1->data =  new TreeNode(0);
        healthyTip2->data =  new TreeNode(0);
        healthyTip3->data =   new TreeNode(0);
        
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
        //long double   healthyTipBranchLengthRatio =1;
        
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
                                          gsl_rng * randomGenerator,
                                          pll_msa_t *msa,
                                          int &numNodes,
                                          int numClones,
                                          ProgramOptions &programOptions,
                                          long double       cumNumCA,
                                          long double  meanNumCA,
                                          long double  cumNumMIG,
                                          long double  meanNumMIG,
                                          int  &numMIG,
                                          int  &numCA,
                                          long double  &numEventsTot,
                                          //                           pll_unode_t** nodes,
                                          //                           pll_unode_t** treeTips,
                                          //                           pll_unode_t    **treeRootInit,
                                          std::vector<std::vector<int> > ObservedData,
                                            char* ObservedCellNames[],
                                          vector<int> &sampleSizes
                                          ) {
    
    int      c,  i,  m, cumIndivid, isCoalescence, whichInd,
    newInd, eventNum, numActiveGametes,
    isMigration, whichClone, currentNumberAliveClones;
    int     labelNodes;
    long double     currentTime, eventTime;
    pll_unode_t  *p;
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
        pop->sampleSize= sampleSizes[i];
        pop->resetMigrationsList();
    }
    //resetMigrationsList( populations,  numClones);

    
    for (i=0; i< msa->count; i++){
        p=new pll_unode_t();
        p->label = msa->label[i];
        p->next =NULL;
        p->back=NULL;
        t = new TreeNode(msa->length);
       
        t->genotypeSequence = ObservedData[i];
        
        p->data = new TreeNode(msa->length);
        treeTips.push_back(p);
        nodes.push_back(p);
    }
    
    for (i = 0; i < (numNodes - programOptions.TotalNumSequences); i++)
    {
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
    pll_unode_t *root=NULL;
    root =  BuildTree(populations, currentPop,seed,
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
    long double  pij, ran;
    int  j, k;
    long double  cumProb[numClones  - (int)(PopChild->order)];
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
    ran =  Random::randomUniformFromGsl();
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
    long double  currentSampleSize;
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
void Chain::SetPopulationsBirthRate( long double  lambda){
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
    
    
    if (doGenerateProportionsVector == YES)
    {
        initProportionsVector();
        long double  alpha[numClones];
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
        popJ->effectPopSize = popJ->x * totalEffectPopSize;
    }
}
void Chain::initProportionsVector(){
    for (size_t i = 0; i < numClones; i++) {
        proportionsVector.push_back(0);
        oldproportionsVector.push_back(0);
    }
    
}
void Chain::initProportionsVectorFromSampleSizes(vector<int> sampleSizes)
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
void Chain::generateProportionsVectorFromDirichlet(long double  alpha[]){
    
    long double  theta[numClones];
    vector<long double> outputVector;
    for (unsigned int i = 0; i < numClones; ++i){
        theta[i]=0;
    }
    vector<long double> alphaVector(alpha, alpha + numClones);
    Random::randomDirichletFromVector(alphaVector, outputVector);
    //randomDirichletFromGsl(numClones, alpha, theta);
    for (int i=0; i< numClones; i++)
    {
        proportionsVector.at(i)=outputVector.at(i);
        //oldproportionsVector.at(i)=outputVector.at(i);
    }
    //std::copy(proportionsVector.begin(), proportionsVector.end(), theta);
}
void Chain::copyProportionsVector(long double  alpha[]){
    
    
    for (int i=0; i< numClones; i++)
    {
        proportionsVector.at(i)=alpha[i];
    }
    
}
void Chain::initTotalEffectivePopulationSize(MCMCoptions &mcmcOptions, gsl_rng* randomGenerator)
{
    
    int currentTotalSampleSize = mcmcOptions.numberTumorCells;
    
    if(mcmcOptions.priorsType==0)
        totalEffectPopSize = Random::RandomLogUniform(mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto, mcmcOptions.useGSLRandomGenerator, randomGenerator);
    else if(mcmcOptions.priorsType==1)
        totalEffectPopSize = Random::RandomExponentialStartingFrom(mcmcOptions.lambdaExponentialPriorTotalEffectivePopSize, currentTotalSampleSize/mcmcOptions.fixedLambda, mcmcOptions.useGSLRandomGenerator, randomGenerator);
    else
        totalEffectPopSize = Random::RandomPowerLawDistribution(mcmcOptions.parameterPowerLawDistributionTotalEffectPopSize, currentTotalSampleSize/mcmcOptions.fixedLambda, mcmcOptions.useGSLRandomGenerator);
    //printf("\n initial total effect pop size %d \n", totalEffectPopSize);
}
void Chain::FillChainPopulationsFromPriors( ProgramOptions &programOptions,  MCMCoptions &mcmcOptions, vector<int> &sampleSizes, gsl_rng* randomGenerator)
{
    Population* popI;
    int i;
    long double  randomDelta;
    mutationRate = Random::RandomLogUniform(mcmcOptions.MutRatefrom, mcmcOptions.MutRateto, mcmcOptions.useGSLRandomGenerator, NULL);
    initTotalEffectivePopulationSize(mcmcOptions,  randomGenerator);
    initProportionsVector();
    //long double  lambda = 1;
    SetPopulationsBirthRate(mcmcOptions.fixedLambda);
    
    // GenerateEffectPopSizesFromPriors2( programOptions.noisy,    YES);
    if (programOptions.doUseFixedTree ==NO)
    {
        long double  alpha[numClones];
        std::fill_n(alpha, numClones, 1.0);
        generateProportionsVectorFromDirichlet(alpha);
        if (programOptions.populationSampleSizesKnown == NO)
        {
            InitPopulationSampleSizes(populations, programOptions.numCells, programOptions.numClones,proportionsVector,  randomGenerator);
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
                    randomDelta = Random::RandomLogUniform(mcmcOptions.Deltafrom, mcmcOptions.Deltato, mcmcOptions.useGSLRandomGenerator, NULL);
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
        GenerateTimesFromPriorsOriginal(programOptions.noisy,  randomGenerator);
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
        popI->effectPopSize = totalEffectPopSize * popI->x;
        popI->popSize = totalEffectPopSize * popI->x;
        //popI->effectPopSize = totalEffectPopSize * proportionsVector[i];
        //printf("\n initial effect pop size  for pop index %d  is %lf \n", popI->index,popI->effectPopSize);
    }
}
void Chain::initTimeOriginSTDYoungerPopulations(MCMCoptions &mcmcOptions){
    int i;
    Population *popI;
    Population* oldestPop=getPopulationbyIndex(numClones -1);
    for( i = 0 ; i < numClones; i++)
    {
        popI=populations[i];
        if (popI!= oldestPop && popI->x !=0)
            popI->timeOriginSTD = popI->scaledtimeOriginInput /  popI->x;
        if (popI->x ==0 && mcmcOptions.verbose>=2)
            printf("\n error effective  pop size proportion  size is 0 \n");
        if (mcmcOptions.verbose>=2)
            printf("\n population %d with order %d  has  time origin std of: %Lf and input: %Lf\n", popI->index, popI->order, popI->timeOriginSTD, popI->timeOriginInput );
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
    
    InitListPossibleMigrations(populations,numClones);
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

void Chain::GenerateTimesFromPriorsOriginal(int noisy,  gsl_rng* randomGenerator) {
    long double     *Uvector;
    int i, j, m,  l;
    long double     TotalProbability, LocalProbability, aboveTerm, belowTerm, ranHere;
    int  AttemptsAcceptation = 0;
    TotalProbability = 0.0;
    LocalProbability = 0.0;
    aboveTerm = 0.0;
    belowTerm = 0.0;
    ranHere = 0.0;
    if (noisy > 1)
        printf("Estimation of times of origin of clones ..\n");
    Uvector = (long double  *) malloc((numClones )* (long) sizeof(double));
    if (!Uvector)
    {
        fprintf (stderr, "Could not allocate Xvector (%lu bytes)\n", (numClones ) * (long) sizeof(double));
        exit (-1);
    }
    int doAcceptTimes = NO;
    int doEstimateTimesOriginClones=YES;
    Population *popJ, *popI, *popL;
    long double  rand;
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
                Uvector[j] =  Random::randomUniformFromGsl();
                rand = Random::randomUniformFromGsl();
                popJ->timeOriginSTD= (1/ popJ->delta)*log(1-(popJ->delta / log(rand)));
                //popJ->timeMigrationSTDCurrentPop= popJ->timeOriginSTD ;
                // popJ->timeOriginInput=proportionsVector[j]*popJ->timeOriginSTD;
                popJ->timeOriginInput = popJ->effectPopSize*popJ->timeOriginSTD;
            }
            
            if (numClones == 1) {
                doAcceptTimes = YES;
                break;
            }
            ListClonesAccordingTimeToOrigin(populations);
            // Calculate probabilities P
            LocalProbability = 0.0;
            TotalProbability = 0.0;
            
            for (i = 0; i < numClones - 1; i++)
            {
                popI=populations[i];
                printf("\ni  = %d ", i);
                LocalProbability = 0.0;
                aboveTerm = 0.0;
                belowTerm = 0.0;
                m = popI->order;
                
                for (l = i + 1; l < numClones; l++)
                {
                    popL=populations[ l];
                    j = popL->order ;
                    aboveTerm = aboveTerm + (popL->popSize * Population::CalculateH(popI->timeOriginSTD * popI->effectPopSize / popL->effectPopSize, popL->timeOriginSTD, popL->delta));
                    belowTerm = belowTerm + popL->popSize;
                }
                LocalProbability = aboveTerm / belowTerm;
                if (i == 1)
                    TotalProbability = 1.0 * LocalProbability;
                else
                    TotalProbability = TotalProbability * LocalProbability;
                
            }
            ranHere =  Random::randomUniformFromGsl2(randomGenerator);
            if (ranHere <= TotalProbability)
                doAcceptTimes = YES;
            if (noisy > 2)
                printf("\nProbability = %Lf (random number = %Lf) [#attempt = %d]\n", TotalProbability, ranHere, AttemptsAcceptation);
        }
    }
    if (noisy > 2)
    {
        printf("\nTimes accepted .. \n");
        printf("\n\nDONE! ranHere = %Lf / TotalProbability = %Lf [total attempts = %d]\n", ranHere, TotalProbability, AttemptsAcceptation);
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
        long double  timeOriginInput = 0.0;
        int sampleSize = 0;
        int popSize =0;
        long double  birthRate = 1.0;
        long double  deathRate = 0.99;
        
        
        auto pop = new Population(ind, ord, timeOriginInput, sampleSize, popSize, birthRate, deathRate, NO);
        populations.push_back(pop);
        
    }
}
/**************** RelabelNodes **************/

void Chain::RelabelNodes(pll_unode_t *p, int &intLabel)
{
    if (p != NULL)
    {
        TreeNode *u=NULL;
        TreeNode *u1=NULL;
        TreeNode *u2=NULL;
        if (p->data !=NULL)
          u= (TreeNode *)(p->data);
        if (p->next  !=NULL && p->next->data!=NULL )
           u1= (TreeNode *)(p->next->data);
        if (p->next !=NULL && p->next->next!=NULL && p->next->next->data!=NULL )
           u2= (TreeNode *)(p->next->next->data);
        
        if (u1!=NULL)
           RelabelNodes (p->next->back, intLabel);
        if (u2!=NULL)
           RelabelNodes (p->next->next->back, intLabel);
        /*RelabelNodes (p->outgroup);*/
        // if (p->left == NULL && p->right == NULL) /* is tip */
        if (p->back != NULL && p->next->back == NULL) /* is tip */
        {
            // p->label = intLabel++;
            //  p->label = (*intLabel);
            // *intLabel=*intLabel+1;
            //p->label = p->index ;
            if (u!=NULL)
            {
                u->label =u->index;
                if (u1!=NULL)
                    u1->label=u->index ;
                if (u2!=NULL)
                  u2->label= u->index ;
            }
        }
        else                  /* all ancester */
        {
            //p->label = intLabel++;
            // p->label = intLabel;
            if (u!=NULL)
            {
                u->label =u->index;
                if (u1!=NULL)
                    u1->label=u->index ;
                if (u2!=NULL)
                    u2->label= u->index ;
            }
            intLabel=intLabel+1;
        }
    }
}
void Chain::InitPopulationSampleSizes(vector<Population*> populations, int TotalSampleSize, int numClones, vector<long double> &proportionsVector, gsl_rng* randomGenerator)
{
    int i,j;
    Population  *popJ;
    double rand;
    long double     *cumSum = (long double  *) calloc((numClones +1), (long) sizeof(double));
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
        rand= Random::randomUniformFromGsl2(randomGenerator);
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
long double  Chain::SumBranches(pll_unode_t *p, long double  mutationRate){
    
    static long double  sum;
    
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
long double  Chain::SumBranches2(pll_rnode_t *p, long double  mutationRate){
    
    static long double  sum;
    
    if (p != NULL)
    {
        // if (p->back  == NULL && p->next !=NULL && p->next->back !=NULL)
        if (p->parent  == NULL )
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
char * Chain::toNewickString ( pll_unode_t *p, long double  mutationRate,     int doUseObservedCellNames)
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
            if (asprintf(&newickString,  "healthycell:%10.9Lf",  (u->anc1->timePUnits - u->timePUnits) * mutationRate)<0)
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
                if (asprintf(&newickString,   "%s:%10.9Lf",  u->observedCellName, (u->anc1->timePUnits - u->timePUnits)*mutationRate)<0)
                    return NULL;
                //snprintf(newickString,  size,  "%s:%10.9lf",  p->observedCellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
                return newickString;
            }
            else{
                if (asprintf(&newickString,   "%s:%10.9Lf",  u->cellName, (u->anc1->timePUnits - u->timePUnits)*mutationRate)<0)
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
                
                if (asprintf(&newickString, "(%s,%s):%10.9Lf", left, right,  (u->anc1->timePUnits - u->timePUnits)*mutationRate )<0)
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
long double Chain::SumLogDensitiesTimeOriginSTDPopulations() {
    Population* popI;
    
    long double temp;
    long double result = 0;
    unsigned int i;
    
    for ( i = 0; i < numClones; i++)
    {
        popI=populations[i];
        temp=popI->LogDensityTimeSTDFrom( popI->timeOriginSTD, 0);
        if (isnan(temp) || isinf(temp))
            fprintf (stderr, "\n isNan temp LogDensityTime, Delta: %Lf, T:%Lf \n",popI->delta ,popI->timeOriginSTD);
        result = result + temp;
        // fprintf (stderr, "\n Product log Density Time   = %lf after pop order %d,  popI->LogDensityTime( popI->timeOriginSTD) %lf,  popI->timeOriginSTD: %lf, popI->delta: %lf\n", product, popI->order, popI->LogDensityTime( popI->timeOriginSTD), popI->timeOriginSTD,popI->delta );
    }
    return result;
}

long double Chain::SumLogProbFatherPopulations() {
    
    Population *fatherPop;
    Population *popJ;
    long double  result=0;
    long double temp;
    for ( unsigned int j = 0; j < numClones-1 ; j++)
    {
        popJ = populations[j ];
        if (popJ->FatherPop !=NULL)
        {
            fatherPop = popJ->FatherPop;
            result = result + log( fatherPop->popSize);
            temp=popJ->timeOriginSTD * popJ->x / fatherPop->x;
            temp=Population::LogCalculateH(temp, fatherPop->timeOriginSTD, fatherPop->delta);
            result = result  + temp;
            // fprintf (stderr, "\n Product calculate H    = %lf after pop order %d \n", product, popJ->order);
        }
        if(popJ->order < numClones-1 && popJ->FatherPop == NULL)
            fprintf (stderr, "\n the pop with order %d has father pop null\n",popJ->order );
    }
    return(result);
}

/************************ LogConditionalLikelihoodTree ***********************/
/*  LogConditionalLikelihoodTree */
long double  Chain::LogConditionalLikelihoodTree( ProgramOptions &programOptions, MCMCoptions &mcmcOptions  )
{

    long double result;
    result =SumLogDensitiesTimeOriginSTDPopulations();
    result= result + SumLogProbFatherPopulations();
        // fprintf (stderr, "\n Product after  = %lf \n", result);
    result = result + LogDensityCoalescentTimesForPopulation2();

    return result;
}
/************************ LogDensityCoalescentTimesForPopulation ***********************/
/*  LogDensityCoalescentTimesForPopulation */
long double  Chain::LogDensityCoalescentTimesForPopulation()
{
    long double  result =0;
    int i;
    Population *popI;
    int numberLeftCoalescences;
    int numberLeftMigrations;
    int numberAliveCells;
    int currentCoalescentEvent=0;
    int currentMigrationEvent=0;
    long double  temp;
    long double  termOnlyAfterFirstCoalEvent=0.0;
    long double  lastEventTimeBeforeMigration=0;
    for ( i = 0; i < numClones; i++){
        popI = populations[i];
        fprintf (stderr, "\n Likelihood population index %d order %d and sample size %d \n", popI->index, popI->order, popI->sampleSize);
        currentCoalescentEvent=0;
        currentMigrationEvent=0;
        lastEventTimeBeforeMigration=0;
        numberAliveCells= popI->sampleSize;
        if (popI->sampleSize <=1)
            continue;
        // numberLeftCoalescences =  popI->numCompletedCoalescences; //in the numCompletedCoalescences we are considering also the migrations
        //numberLeftCoalescences =  popI->numCompletedCoalescences - (popI->numIncomingMigrations-1);
        numberLeftCoalescences = popI->CoalescentEventTimes.size() - (popI->immigrantsPopOrderedByModelTime.size()-1);
        //numberLeftMigrations = popI->numIncomingMigrations-1;
        fprintf (stderr, "\n The number of completed coalescences %d, %lu, migrations %d,%lu  for population order %d and sample size %d \n",popI->numCompletedCoalescences, popI->CoalescentEventTimes.size(), popI->numIncomingMigrations -1, popI->immigrantsPopOrderedByModelTime.size() -1, popI->order, popI->sampleSize );
        if (popI->numCompletedCoalescences != popI->CoalescentEventTimes.size()){
            fprintf (stderr, "\n wrong \n");
        }
        numberLeftMigrations = popI->immigrantsPopOrderedByModelTime.size()-1;
        //we are not counting time of origin as a migration
        if (numberLeftCoalescences <=0)
            continue;
        while(numberLeftMigrations > 0)
        {
            while(popI->CoalescentEventTimes[currentCoalescentEvent] < popI->immigrantsPopOrderedByModelTime[currentMigrationEvent].first
                  && numberAliveCells > 1 &&  (currentCoalescentEvent) <= (popI->CoalescentEventTimes.size()-1))
            {
                fprintf (stderr, "\n Coalescesce between migrations %d time %Lf\n",currentCoalescentEvent,popI->CoalescentEventTimes[currentCoalescentEvent] );
                temp=log(numberAliveCells * (numberAliveCells-1.0)/2.0);
                if (isnan(temp) || isinf(temp) )
                    fprintf (stderr, "\n isNan product\n");
                result= result + temp;
                
                temp = -1.0 * Population::LogCalculateH(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta);
                if (isnan(temp) || isinf(temp) )
                    fprintf (stderr, "\n isNan product\n");
                result= result + temp;
                termOnlyAfterFirstCoalEvent =(currentCoalescentEvent == 0)?0:Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent-1], popI->timeOriginSTD, popI->delta);
                temp =  (numberAliveCells/ 2.0)* (numberAliveCells-1)*(Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta)-termOnlyAfterFirstCoalEvent);
                if (isnan(temp) || isinf(temp) )
                    fprintf (stderr, "\n isNan product\n");
                result= result -temp;
                if (isnan(result) || isinf(result))
                    fprintf (stderr, "\n isNan product\n");
                lastEventTimeBeforeMigration=popI->CoalescentEventTimes[currentCoalescentEvent];
                currentCoalescentEvent++;
                numberLeftCoalescences--;
                numberAliveCells--;
            }
            //if (numberLeftMigrations > 0 && numberAliveCells > 1)// if there are migrations
            if (numberLeftMigrations > 0 &&  currentMigrationEvent <= (popI->immigrantsPopOrderedByModelTime.size()-1))// if there are migrations
            {
                fprintf (stderr, "\n Migration  %d time %Lf from pop order %d \n",currentMigrationEvent,popI->immigrantsPopOrderedByModelTime[currentMigrationEvent].first, popI->immigrantsPopOrderedByModelTime[currentMigrationEvent].second->order);
                temp= popI->LogProbNoCoalescentEventBetweenTimes(lastEventTimeBeforeMigration,popI->immigrantsPopOrderedByModelTime[currentMigrationEvent].first, numberAliveCells );
                if (isnan(temp) || isinf(temp) )
                    fprintf (stderr, "\n isNan product\n");
                lastEventTimeBeforeMigration=popI->immigrantsPopOrderedByModelTime[currentMigrationEvent].first;
                result= result+ temp;
                numberAliveCells++;
                numberLeftMigrations--;
                currentMigrationEvent++;
            }
        }
        fprintf (stderr, "\n More coalescences \n" );
        //here there are only coalescents events left(al least one event)
        while(numberLeftCoalescences > 0 && numberAliveCells > 1 && currentCoalescentEvent <= (popI->CoalescentEventTimes.size()-1) )
        {
            fprintf (stderr, "\n Coalescesce after all true migrations  %d time %Lf\n",currentCoalescentEvent,popI->CoalescentEventTimes[currentCoalescentEvent] );
            
            temp = log(numberAliveCells * (numberAliveCells-1.0)/2.0);
            if (isnan(temp) || isinf(temp) )
                fprintf (stderr, "\n isNan temp\n");
            result= result + temp;
            temp = -1.0 * Population::LogCalculateH(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta);
            
            if (isnan(temp) || isinf(temp) )
                fprintf (stderr, "\n isNan temp\n");
            
            temp = -1.0 * Population::LogCalculateH(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta);
            result= result + temp;
            
            termOnlyAfterFirstCoalEvent =(currentCoalescentEvent == 0)?0:Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent-1], popI->timeOriginSTD, popI->delta);
            temp=( numberAliveCells/ 2.0)* (numberAliveCells-1.0)*(Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta)-termOnlyAfterFirstCoalEvent);
            if (isnan(temp) || isinf(temp) )
                fprintf (stderr, "\n isNan temp\n");
            temp=( numberAliveCells/ 2.0)* (numberAliveCells-1.0)*(Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta)-termOnlyAfterFirstCoalEvent);
            result= result -  temp;
            
            currentCoalescentEvent++;
            numberLeftCoalescences--;
            numberAliveCells--;
        }
        //fprintf (stderr, "\n Result = %lf after population order %d \n", result, popI->order);
    }
    //fprintf (stderr, "\n Result = %lf \n", result);
    
    if (isnan(result) || isinf(result))
        fprintf (stderr, "\n isNan product\n");
    
    return result;
}
long double  Chain::LogDensityCoalescentTimesForPopulation2(){
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
                    temp = -1.0 * Population::LogCalculateH(timeCurrentEvent,popI->timeOriginSTD, popI->delta);
                    result= result + temp;
                    termOnlyAfterFirstCoalEvent =(currentCoalescentEventInThisEpoch == 0)?Population::FmodelTstandard(lastEventTimeBeforeMigration, popI->timeOriginSTD, popI->delta):Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEventInThisEpoch-1], popI->timeOriginSTD, popI->delta);//if no coalescent
                    temp =  (numberAliveCells / 2.0)* (numberAliveCells - 1.0)*(Population::FmodelTstandard(timeCurrentEvent,popI->timeOriginSTD, popI->delta)-termOnlyAfterFirstCoalEvent);
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
                    temp= popI->LogProbNoCoalescentEventBetweenTimes(lastEventTimeBeforeMigration,timeCurrentEvent, numberAliveCells );
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
long double  Chain::LogConditionalLikelihoodSequences(pll_msa_t * msa, char* NewickString, ProgramOptions &programOptions, long double seqError,long double dropoutError)
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
    {
        
        // printf("Original sequence (alignment) length : %d\n", msa->length);
        
        
    }
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
    
    double  * empirical_frequencies;
    empirical_frequencies = pllmod_msa_empirical_frequencies(partition);
    
    double  * empirical_subst_rates = pllmod_msa_empirical_subst_rates( partition);
    
    unsigned int * weight = pll_compress_site_patterns(msa->sequence,
                                                       pll_map_gt10,
                                                       tip_nodes_count,
                                                       &(msa->length));
    
    //  printf("Number of unique site patterns: %d\n\n", msa->length);
    
    /* initialize the array of base frequencies  AA CC GG TT AC/CA AG/GA AT/TA CG/GC CT/TC GT/TG  */
    //    double user_freqs[10] = { 0.254504, 0.157238, 0.073667, 0.287452, 0.027336,
    //        0.058670, 0.041628, 0.024572, 0.064522, 0.010412 };
    
    
    /* substitution rates: for GTR4 model those are 6 "regular" DNA susbt. rates + 1 rate
     * for "unlikely" double substitutions (eg A/A -> C/T) */
    double  unique_subst_rates[7] = { 0.001000, 0.101223, 0.001000, 0.001000, 1.000000,
        0.001000, 0.447050 };
    
    /* get full above-diagonal half-matrix */
    double  * user_subst_rates = Chain::expand_uniq_rates(model->states, unique_subst_rates,
                                                          model->rate_sym);
    
    double  rate_cats[RATE_CATS] = {0};
    
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
    long double  loglh = pllmod_treeinfo_compute_loglh(treeinfo, 1);
    pllmod_treeinfo_destroy(treeinfo);
    //long double  * clv ;
    
    //int scaler_index = PLL_SCALE_BUFFER_NONE;
    //unsigned int * scaler = (scaler_index == PLL_SCALE_BUFFER_NONE) ?NULL : partition->scale_buffer[scaler_index];
    //unsigned int states = partition->states;
    //unsigned int states_padded = partition->states_padded;
    //unsigned int rates = partition->rate_cats;
    //unsigned int *site_id = 0;
    
    pll_partition_destroy(partition);
    /* we will no longer need the tree structure */
    //pll_utree_destroy(unrootedTree, NULL);
    destroyTree(unrootedTree, NULL);
    pllmod_util_model_destroy(model);
    
    return loglh;
}
void Chain::set_partition_tips_costum( pll_partition_t * partition, pll_msa_t * msa, ProgramOptions &programOptions, long double  seqError, long double  dropoutError)
{
    //pll_state_t pll_map_gt10_2[256];
    //int states =10;
    int i;
    //int currentState;
    //int from, to;
    
    
    //unsigned int state;
    //long double  * _freqs;
    
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
                                         long double  _seq_error_rate,
                                         long double  _dropout_rate
                                         )
{
    
    pll_state_t c;
    unsigned int i,j;
    double  * tipclv = partition->clv[tip_index];
    unsigned int state_id ;
    
    pll_repeats_t * repeats = partition->repeats;
    unsigned int use_repeats = pll_repeats_enabled(partition);
    unsigned int ids = use_repeats ?
    repeats->pernode_ids[tip_index] : partition->sites;
    
    static const long double  one_3 = 1. / 3.;
    static const long double  one_6 = 1. / 6.;
    static const long double  one_8 = 1. / 8.;
    static const long double  three_8 = 3. / 8.;
    static const long double  one_12 = 1. / 12.;
    
    // TODO: move it out of here
    unsigned int undef_state = (unsigned int) (pow(2, partition->states)) - 1;
    
    long double  sum_lh = 0.;
    
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

void Chain::initMutationRate( MCMCoptions &mcmcOptions, ProgramOptions &programOptions, gsl_rng * randomGenerator) {
    if (programOptions.doUsefixedMutationRate)
        theta = programOptions.mutationRate;
    else
    {
        if(mcmcOptions.priorsType==0)
        {
            theta = Random::RandomLogUniform(mcmcOptions.MutRatefrom, mcmcOptions.MutRateto, mcmcOptions.useGSLRandomGenerator, randomGenerator);
            
        }
        else if(mcmcOptions.priorsType==1)
        {
            theta = Random::RandomExponentialStartingFrom( mcmcOptions.lambdaExponentialPriorMutationRate, 0, mcmcOptions.useGSLRandomGenerator, randomGenerator);
        }
        
        else{
            
            theta = Random::RandomPowerLawDistribution(mcmcOptions.parameterPowerLawDistributionMutationRate, 0, mcmcOptions.useGSLRandomGenerator);
            
        }
    }
}

void Chain::initLogLikelihoods(MCMCoptions &mcmcOptions, pll_msa_t *msa, ProgramOptions &programOptions) {
    
    char *rootedNewick2;
    rootedNewick2 =  pll_rtree_export_newick(initialRootedTree->root,NULL);
    
    //newickString2 = chain->toNewickString ( chain->root, chain->mutationRate,     programOptions.doUseObservedCellNames);
    //printf("\n newick = %s  \n", rootedNewick2);
    
    if (!mcmcOptions.noData)
    {
        long double result=0;
        currentlogSumDensitiesTimeOriginSTDPopulations =SumLogDensitiesTimeOriginSTDPopulations();
        result=currentlogSumDensitiesTimeOriginSTDPopulations;
    currentlogProbFatherPopulations=SumLogProbFatherPopulations();
        result= result + currentlogProbFatherPopulations;
        currentlogDensityCoalescentTimesForPopulation= LogDensityCoalescentTimesForPopulation2();
        result = result + currentlogDensityCoalescentTimesForPopulation;
        currentlogConditionalLikelihoodTree= result;
        // if (mcmcOptions.verbose>=2)
        printf ( "Initial log likelihood of the tree of chain %d is:  %Lf \n",0, currentlogConditionalLikelihoodTree );
        if (mcmcOptions.useSequencesLikelihood ==1)
        {
            currentlogConditionalLikelihoodSequences= LogConditionalLikelihoodSequences( msa,  rootedNewick2, programOptions, seqErrorRate, dropoutRate);
            free(rootedNewick2);
            rootedNewick2=NULL;
        }
        else
        {
            currentlogConditionalLikelihoodSequences=0;
        }
        if (mcmcOptions.verbose>=2)
            fprintf (stderr, "Initial log likelihood of the sequences of chain %d  is = %Lf  \n", 0,currentlogConditionalLikelihoodSequences );
    }
}

void Chain::initPopulationsThetaDelta() {
    for (unsigned int i = 0; i < numClones; ++i)
    {
        auto popI =  populations[i];
        popI->x = proportionsVector[i];
        popI->theta = popI->x * theta;
        popI->delta = popI->deltaT * popI->x;
    }
}

void Chain::initChainTree( std::string &healthyTipLabel, MCMCoptions &mcmcOptions, pll_msa_t *msa, ProgramOptions &programOptions ) {
    
    if (initialRootedTree->inner_count > 1)
        initialUnrootedTree = pll_rtree_unroot(initialRootedTree);
    //chain->initialUnrootedTree = pll_utree_parse_newick(treefileName);
    
    rootRootedTree = initialRootedTree->root;
    root =initialUnrootedTree->nodes[initialUnrootedTree->tip_count + initialUnrootedTree->inner_count - 1];
    
    mcmcOptions.numberTumorCells =initialUnrootedTree->tip_count-1;
  
    pll_utree_reset_template_indices(root, initialUnrootedTree->tip_count);
    
    initPopulationsTipsFromTree(initialUnrootedTree, NO, programOptions.healthyTipLabel);
    initPopulationsTipsFromRootedTree(initialRootedTree, NO, programOptions.healthyTipLabel);
    initNodeDataFromTree();
    rescaleNodeDataFromRootedTree(theta );
    initNumberTipsSubTree(rootRootedTree);
    
    // pll_utree_traverse_apply(chain->initialUnrootedTree->vroot, NULL, Chain::computeNumberTipsSubTree, chain->root->data);
    
    
    initBranches(healthyTipLabel, edgeLengths, edges);
   
    //std::transform(alpha, alpha + chain->numClones , alpha,std::bind2nd(std::divides<double>(),totalSampleSize));
    // std::transform(alpha, alpha + chain->numClones , alpha,[totalSampleSize](long double  a) {return a /totalSampleSize; } );
}
vector<long double> Chain::initVectorSampleSizes(std::string &healthyTipLabel, MCMCoptions &mcmcOptions, ProgramOptions &programOptions){
    
    bool existsZeroSampleSizePop=false;
    vector<long double> alpha(numClones);
    
    vector<double> branchLengths;
    existsZeroSampleSizePop=false;
    
    
    existsZeroSampleSizePop= initPopulationsSampleSizes( rMRCAPopulation, healthyTipLabel);
    
    //int totalSampleSize=chain->initialRootedTree->tip_count-1;//not including the healthytip
    
    for (unsigned int i = 0; i < numClones; ++i)
    {
        auto popI =  populations[i];
        if (existsZeroSampleSizePop)
            alpha.at(i)= popI->sampleSize+1;//this is for avoid to pass  a sample size 0 to the Dirichlet and get a effective population size of 0
        else
            alpha.at(i)= popI->sampleSize;
        if (popI->sampleSize == 0)
            printf("\n population %d and order %d with sample size 0:  %d \n", popI->index,popI->order,popI->sampleSize);
    }
    return alpha;
}
void Chain::initMRCAOldestPopulation(string& healthyTipLabel)
{
    Population* oldestPop=getPopulationbyIndex(numClones -1);
    TreeNode *u;
    if (std::string(initialRootedTree->root->right->label).compare(healthyTipLabel)==0)
    {
        oldestPop->rMRCA =initialRootedTree->root->left;
    }
    else  if (std::string(initialRootedTree->root->left->label).compare(healthyTipLabel)==0)
    {
        oldestPop->rMRCA =initialRootedTree->root->right;
    }
    else{
        //this should not happen
        oldestPop->rMRCA =initialRootedTree->root;
        u= (TreeNode *)(initialRootedTree->root->data);
        printf( "\n Error, the healthy cell is not a direct child of the tree root \n"  );
    }
}

Chain *Chain::initializeChain(   ProgramOptions &programOptions,  MCMCoptions &mcmcOptions, vector<int> &sampleSizes, gsl_rng * randomGenerator, std::vector<std::vector<int> > ObservedData,char* ObservedCellNames[], pll_msa_t *msa, pll_rtree_t * initialRootedTree, string& healthyTipLabel)
{
  
    int        numCA, numMIG;
    long double       cumNumCA = 0.0, meanNumCA = 0.0, cumNumMIG = 0.0, meanNumMIG = 0.0, numEventsTot;
    long int *seed;
    char * rootedNewick2;
    Chain *chain=new Chain(0, programOptions.numClones, 1, programOptions.mutationRate, programOptions.seqErrorRate, programOptions.dropoutRate);
    
    if (programOptions.numberClonesKnown)
    {
        chain->numNodes = programOptions.numNodes;
    }
    chain->InitChainPopulations( programOptions.noisy, programOptions.TotalNumSequences );
    
    chain->initMutationRate( mcmcOptions, programOptions, randomGenerator);
    //if (mcmcOptions.verbose>=2)
        fprintf(stderr, "\n initial theta %.10Lf \n",chain->theta );
    chain->samplePopulationGrowthRateFromPriors(mcmcOptions, randomGenerator );
    
    chain->SetPopulationsBirthRate(mcmcOptions.fixedLambda);
    
    if (programOptions.doUseFixedTree == NO)
    {
        Population* oldestPop=chain->getPopulationbyIndex(chain->numClones -1);
        oldestPop->timeOriginSTD  =   Random::RandomDensityModelTimeOrigin (oldestPop->delta, mcmcOptions.useGSLRandomGenerator, 0, randomGenerator );
        chain->initProportionsVectorFromSampleSizes(sampleSizes);
        chain->initPopulationsThetaDelta();
        chain->root = chain->MakeCoalescenceTree (seed,
                                                 randomGenerator,
                                                  msa,
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
                                                  ObservedData,
                                                  ObservedCellNames,
                                                  sampleSizes
                                                  ) ;
        
        //generate the newick representation of the simulated tree
        //
        rootedNewick2 = chain->toNewickString ( chain->root, chain->mutationRate,     programOptions.doUseObservedCellNames);
        printf("\n newick = %s  \n", rootedNewick2);
        chain->initialRootedTree = pll_rtree_parse_newick(rootedNewick2);
        chain->initChainTree( healthyTipLabel, mcmcOptions, msa, programOptions);
        free(rootedNewick2);
        rootedNewick2=NULL;
    }
    else
    {
        chain->initialRootedTree = initialRootedTree;
        if (chain->initialRootedTree ==NULL)
        {
            fprintf(stderr, "\n There was an error reading the tree file! \n" );
            exit(-1);
        }
         Population* oldestPop=chain->getPopulationbyIndex(chain->numClones -1);
       
        chain->initChainTree( healthyTipLabel, mcmcOptions, msa, programOptions);
        vector<long double> alpha;
        
        chain->rMRCAPopulation = chain->initTimeOfOriginsOnRootedTree( chain->edgeLengths, programOptions.numClones -1,  healthyTipLabel, mcmcOptions, randomGenerator);
        alpha=chain->initVectorSampleSizes(healthyTipLabel, mcmcOptions, programOptions);
        Random::randomDirichletFromVector (alpha, chain->proportionsVector);
        chain->initPopulationsThetaDelta();
        
        TreeNode *u= (TreeNode *)(oldestPop->rMRCA->data);
        oldestPop->setLowerBoundTimeOriginInput(u->timeInputTreeUnits);
        
        oldestPop->timeOriginSTD =   Random::RandomDensityModelTimeOrigin (oldestPop->delta, mcmcOptions.useGSLRandomGenerator, oldestPop->lowerBoundTimeOriginInput/ (chain->theta * oldestPop->x), randomGenerator);
        chain->initOriginTimeOldestPopulation(healthyTipLabel, mcmcOptions);

        chain->ListClonesAccordingTimeToOrigin(chain->populations);
        chain->initTimeOriginSTDYoungerPopulations(mcmcOptions);
        chain->initPopulationMigration();//after setting the timeSTD
        chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, healthyTipLabel);
        chain->filterSortPopulationsCoalescentEvents();
        
        // totalTreeLength = SumBranches2(rootRootedTree, mutationRate);
        
        chain->initLogLikelihoods(mcmcOptions, msa, programOptions);
        
    }//else
    return chain;
}
void Chain::computeNumberTipsSubTree(pll_unode_t *node, void *data)
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
    if (node->left ==NULL && treeNode!=NULL)//tip
        treeNode->numberOfTipsSubTree =1;
    else if(node->left ==NULL && treeNode==NULL){
        node->data =  new TreeNode(0);
        ((TreeNode *)(node->data))->numberOfTipsSubTree =1;
    }
    else {
        initNumberTipsSubTree(node->left);
        initNumberTipsSubTree(node->right);
        TreeNode * treeNodeLeft= (TreeNode *) (node->left->data);
        TreeNode * treeNodeRight= (TreeNode *) (node->right->data);
        treeNode->numberOfTipsSubTree = treeNodeLeft->numberOfTipsSubTree + treeNodeRight->numberOfTipsSubTree;
    }
}
void  Chain::initListMoves()
{
    moves.clear();
    Population *popI;
    if (numClones > 1){
        NewProportionsVectorMove *newProportionsVector= new NewProportionsVectorMove(this, "new Proportions Vector Move");
        moves.push_back(newProportionsVector );
    }
    NewGlobalScaledMutationRateMove *newGlobalScaledMutationRateMove;
    NewGlobalScaledGrowthRateForPopulationMove *newGlobalScaledGrowthRateMove;
    NewTimeOriginOnTreeforPopulationMove *newTimeOriginOnTreeforPopulationMove ;
    
    //Pair of theta and r move
    newGlobalScaledMutationRateMove= new NewGlobalScaledMutationRateMove(this, "new Pair Global scaled mutation rate(theta) and scaled growth rate(deltaT) Move", NULL);
    moves.push_back(newGlobalScaledMutationRateMove);
    for( int i = 0 ; i < numClones; i++)
    {
        popI=populations[i];
        
        newGlobalScaledGrowthRateMove =  new NewGlobalScaledGrowthRateForPopulationMove(this, "new global scaled growth rate(deltaT)", popI);
        moves.push_back(newGlobalScaledGrowthRateMove);
        
        newTimeOriginOnTreeforPopulationMove= new NewTimeOriginOnTreeforPopulationMove(this, "new Time of Origin Move on Tree for population", popI);
        moves.push_back(newTimeOriginOnTreeforPopulationMove);
    }
}

void Chain::runChain(   MCMCoptions &mcmcOptions,  gsl_rng *randomGenerator,  FilePaths &filePaths, Files &files,  ProgramOptions &programOptions,
                     char* ObservedCellNames[], pll_msa_t * msa, vector<int> &sampleSizes, int currentIteration
                     ){
    
    std::random_device rng;
    std::mt19937 urng(rng());
    std::shuffle ( moves.begin(), moves.end(), urng );
    vector<MCMCmove*>::iterator it;
    MCMCmove* currentMove;
    for ( it=moves.begin(); it!=moves.end(); ++it)
    {
        currentMove=  *it;
        //fprintf (stderr, "\n>> started  move %s ", currentMove->name().c_str());
        currentMove->move(programOptions, mcmcOptions, randomGenerator);
        // fprintf (stderr, "\n>> finished  move %s \n", currentMove->name().c_str());
    }
    //moves.clear();
    //   }
}

//void Chain::WilsonBaldingMove( long int *seed, long double  currentLogLikelihoodSequences, ProgramOptions &programOptions, char *ObservedCellNames[], pll_msa_t * msa)
//{
//    if (numClones > 1)
//    {
//        int i,j;
//        long double  currentProb, proposalProb;
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
//        long double  timeOrigin = popI->timeOriginSTD;
//        TreeNode *nodeReAttach;
//        Population *tempPopulation;
//        long double  timeOriginScaledOtherUnits;
//        long double  distanceModelTimeProposalEdge;
//        long double  acceptanceProbability=0;
//        int numCandidateNodes=0;
//        long double  sumNumerators=0;
//        long double  sumDenominators=0;
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
//        long double  time_from_r_to_p = RandomUniform(seed)*(distanceModelTimeProposalEdge *proposalPop->effectPopSize)  ;
//        long double  effect_pop_size =proposalPop->effectPopSize ;
//
//        if (distanceModelTimeProposalEdge >0)
//            sumNumerators= sumNumerators + log(distanceModelTimeProposalEdge);
//
//        long double  distanceModelTimeCurrentEdge =MRCA->anc1->time- (timeOrigin * popI->effectPopSize ) /(currentFatherPop->effectPopSize);
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
//        long double  newlogConditionalLikelihoodSequences;
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
//        long double  randomNumber= RandomUniform (seed);
//
//        long double  LogAcceptanceRate = (sumNumerators - sumDenominators) >0? (sumNumerators - sumDenominators) :0;
//
//        if (log(randomNumber) < LogAcceptanceRate )
//        {
//            // accept the WB move
//            printf("\n Accepted Wilson Balding move\n");
//            //int result= pll_utree_tbr(MRCA->nodeBack,edgeReAttach,rollback_info);
//            chain->currentlogConditionalLikelihoodSequences =newlogConditionalLikelihoodSequences;
//            long double   currentLikelihoodTree= chain->currentlogConditionalLikelihoodTree;
//            //update the coalescentEventTimes
//            long double  newLogConditionalLikelihoodTree= LogConditionalLikelihoodTree(chain->root, chain->nodes, chain->populations,  chain->numClones);
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
        long double  time = nodes_visited[node];
        long double  length = node->length;
        long double  coal_time = time + length;
        
        pll_unode_t *parent = node->back;
        if (nodes_visited.count(parent) == 0) {
            nodes_visited[parent] = coal_time;
            q.push_back(parent);
        }
        
        pop->CoalescentEventTimes.push_back(coal_time);
    }
}

void Chain::initializeCoalescentEventTimesFormSampleSizes(pll_utree_t *utree, vector<int > &sampleSizes, string &healthyCellLabel )
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
    
    initPopulationsTipsFromTree(utree, YES, healthyCellLabel);
    
    for (unsigned int i = 0; i < numClones; ++i)
    {
        auto pop = populations[i];
        updateCoalTimes(pop);
    }
}
void Chain::initPopulationsTipsFromTree(pll_utree_t *utree, bool assignationKnown, string &healthyCellLabel )
{
    if (utree != NULL)
    {
        for (unsigned int i = 0; i < utree->inner_count + utree->tip_count; ++i)
        {
            auto node = utree->nodes[i];
            
            
            
            if (node->label != 0 && node->label && std::string(node->label).compare(healthyCellLabel)!=0) // only the tips have labels
            {
                //cout << node->label << ": " << node->node_index << ", back: " << node->back->node_index << ", next: " << node->next << endl;
                treeTips.push_back(node);
                if (assignationKnown)
                {
                    int popIdx = getPopulationIndex(node->label);
                    populations[popIdx]->tips.push_back(node);
                }
            }
            else{
                
                //  cout << "interior: " << node->node_index << ", " << node->clv_index << ", back: " << node->back->node_index << ", next:" << node->next->node_index << ", nextnext:" << node->next->next->node_index<< endl;
                
            }
        }
    }
}
void Chain::initPopulationsTipsFromRootedTree(pll_rtree_t *rtree, bool assignationKnown, string &healthyCellLabel )
{
    if (rtree != NULL)
    {
        for (unsigned int i = 0; i < rtree->inner_count + rtree->tip_count; ++i)
        {
            auto node = rtree->nodes[i];
            
            
            
            if (node->label != 0 && node->label && std::string(node->label).compare(healthyCellLabel)!=0) // only the tips have labels and distinct from healthy tip
            {
                //cout << node->label << ": " << node->node_index << ", parent: " <<node->parent->node_index   << endl;
                rtreeTips.push_back(node);
                if (assignationKnown)
                {
                    int popIdx = getPopulationIndex(node->label);
                    populations[popIdx]->rtips.push_back(node);
                }
            }
            else{
                //                if (node->parent!=NULL)
                //                    cout << "interior: " << node->node_index << ", " << node->clv_index << ", parent: " << node->parent->node_index << ", left:" << node->left->node_index << ", right:" << node->right->node_index<< endl;
                //                else
                //                    cout << "interior: " << node->node_index << ", " << node->clv_index << ", parent: " << " " << ", left:" << node->left->node_index << ", right:" << node->right->node_index<< endl;
                
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
            long double  time = nodes_visited[node];
            long double  length = node->length;
            long double  coal_time = time + length;
            if (node->data == NULL)
            {
                node->data =  new TreeNode(0);
                u1 = (TreeNode *)(node->data);
                u1->initNumberTipsVector(numClones);
                u1->timePUnits = time;
                u1->timeInputTreeUnits = time;
                //printf("updated data for node %d \n", node->node_index);
            }
            if (node->next != NULL &&  node->next->data ==NULL)
            {
                node->next->data =  new TreeNode(0);
                u2 = (TreeNode *)(node->next->data);
                u2->initNumberTipsVector(numClones);
                u2->timePUnits = time;
                u2->timeInputTreeUnits=time;
                //printf("updated data for node %d \n", node->next->node_index);
            }
            if (node->next != NULL &&  node->next->next!=NULL && node->next->next->data ==NULL)
            {
                node->next->next->data = new TreeNode(0);
                u3 = (TreeNode *)(node->next->next->data);
                u3->initNumberTipsVector(numClones);
                u3->timePUnits = time;
                u3->timeInputTreeUnits=time;
                //printf("updated data for node %d \n", node->next->next->node_index);
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
            long double  time = nodes_visited[node];
            long double  length = node->length;
            long double  coal_time = time + length;
            if (node->data == NULL)
            {
                node->data =  new TreeNode(0);
                u1 = (TreeNode *)(node->data);
                u1->timePUnits = time;
                u1->initNumberTipsVector(numClones);
                //printf("updated data for node %d \n", node->node_index);
            }
            
            
            pll_rnode_t *parent = node->parent;
            if (parent!=NULL && nodes_visited.count(parent) == 0) {
                nodes_visited[parent] = coal_time;
                q.push_back(parent);
            }
        }
    }
}
void Chain::rescaleNodeDataFromRootedTree(long  double  scale)
{
    long double time;
    long double time_input_tree_units;
    long double length;
    long double coal_time;
    long double coal_time_input_tree_units;
    if (rtreeTips.size()>0)
    {
        TreeNode* u1;
        deque<pll_rnode_t *> q(rtreeTips.begin(), rtreeTips.end());
        unordered_map<pll_rnode_t *, long double> nodes_visited;
        unordered_map<pll_rnode_t *, long double> nodes_visited_input_tree_units;
        for (auto it = rtreeTips.begin(); it != rtreeTips.end(); ++it)
        {
            nodes_visited[*it] = 0;
            nodes_visited_input_tree_units[*it] = 0;
        }
        while (q.size() != 0) {
            pll_rnode_t *node = q.front();
            q.pop_front();
            time = nodes_visited[node];
            time_input_tree_units = nodes_visited_input_tree_units[node];
            length = node->length ;
            
            coal_time = time + (length / scale);
            coal_time_input_tree_units = time_input_tree_units + length ;
            if (node->data == NULL)
            {
                node->data =  new TreeNode(0);
                u1 = (TreeNode *)(node->data);
                u1->timePUnits= time ;
                u1->timeInputTreeUnits = time_input_tree_units;
                u1->initNumberTipsVector(numClones);
                // printf("updated data for node %d and time input %Lf and scaled time %Lf by theta %Lf \n", node->node_index, u1->timeInputTreeUnits, u1->timePUnits, scale );
            }
            pll_rnode_t *parent = node->parent;
            if (parent!=NULL && nodes_visited.count(parent) == 0) {
                nodes_visited[parent] = coal_time;
                nodes_visited_input_tree_units[parent] = coal_time_input_tree_units;
                q.push_back(parent);
            }
        }
    }
}
void Chain::updateNodeScaledTimeForRootedTree(long double  newTheta)
{
    for (unsigned int i = 0; i < initialRootedTree->tip_count+ initialRootedTree->inner_count; ++i)
    {
        auto node =initialRootedTree->nodes[i];
        ((TreeNode *) node->data)->timePUnits =  ((TreeNode *) node->data)->timeInputTreeUnits / newTheta;
        
    }
}



void Chain::safeTreeNodeCurrentTimePUnits(){
    
    TreeNode *u;
    for (unsigned int i = 0; i < initialRootedTree->tip_count+ initialRootedTree->inner_count; ++i)
    {
        auto node =initialRootedTree->nodes[i];
        u = (TreeNode *) node->data;
        u->oldtimePUnits = u->timePUnits;
        
    }
}

void Chain::rollbackTreeNodeCurrentTimePUnits(){
    
    TreeNode *u;
    for (unsigned int i = 0; i < initialRootedTree->tip_count+ initialRootedTree->inner_count; ++i)
    {
        auto node =initialRootedTree->nodes[i];
        u = (TreeNode *) node->data;
        u->timePUnits = u->oldtimePUnits;
    }
}
std::map<pll_unode_t*, Population*>  Chain::chooseTimeOfOriginsOnTree( gsl_rng * randomGenerator)
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
    long double  cumBranchLength=0;
    long double  random;
    long double  totalBranchLength1=0;
    std::map<pll_unode_t*, Population*> mrcaOfPopulation;
    
    
    for (unsigned int i = 0; i < initialUnrootedTree->tip_count+ initialUnrootedTree->inner_count; ++i)
    {
        auto node =initialUnrootedTree->nodes[i];
        
        
        if (node->back && visitedNodes.count(node) == 0  )
        {
            edge = new pll_tree_edge_t();
            u= (TreeNode *)(node->data);
            v= (TreeNode *)(node->back->data);
            edge->edge.utree.parent = u->timeInputTreeUnits > v->timeInputTreeUnits ? node :  node->back;
            edge->edge.utree.child= u->timeInputTreeUnits < v->timeInputTreeUnits ? node :  node->back;;
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
        long double  maxCumBranchLength = *std::max_element(cumulativeBranchLengths.begin(), cumulativeBranchLengths.end());
        long double  maxBranchLength = *std::max_element(branchLengths.begin(), branchLengths.end());
        //        std::vector<int>::iterator it = std::find(branchLengths.begin(), branchLengths.end(), maxBranchLength);
        branchLengths.erase(std::remove(branchLengths.begin(), branchLengths.end(), maxBranchLength), branchLengths.end());
        
        vector<long double>  cumBranchLengths( branchLengths.size());
        std::partial_sum(branchLengths.begin(), branchLengths.end(), cumBranchLengths.begin(), plus<double>());
        
        cumBranchLengths.insert(cumBranchLengths.begin(), 0);
        
        maxCumBranchLength = *std::max_element(cumBranchLengths.begin(), cumBranchLengths.end());
        
        
        std::transform(cumBranchLengths.begin(), cumBranchLengths.end(), cumBranchLengths.begin(),[maxCumBranchLength](long double  a) {return a /maxCumBranchLength; } );
        
        long double  cumBranchLengthsArray[cumBranchLengths.size()];
        std::copy(cumBranchLengths.begin() ,cumBranchLengths.end()
                  ,cumBranchLengthsArray );
        int  nextEvent;
        do{
            
            random =Random::randomUniformFromGsl2(randomGenerator);
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
        //        long double  cumProbPop2[numClones+1];
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
            long double  proposedTime= v->timeInputTreeUnits+ (u->timeInputTreeUnits- v->timeInputTreeUnits)*Random::randomUniformFromGsl2(randomGenerator);
            
            if (proposedTime<=0)
                printf( "time of origin is positive");
            pop->timeOriginInput = proposedTime;
            pop->scaledtimeOriginInput = pop->timeOriginInput / theta;
            k++;
            
        }
        Population* pop=getPopulationbyIndex(k);
        pop->nodeletMRCA =initialUnrootedTree->nodes[initialUnrootedTree->tip_count+ initialUnrootedTree->inner_count-1];//the last population has MRCA the root of the tree
        mrcaOfPopulation[initialUnrootedTree->nodes[initialUnrootedTree->tip_count+ initialUnrootedTree->inner_count-1]]=pop;
    }
    return(mrcaOfPopulation);
}
void Chain::initBranches(string& healthyCellLabel,vector<pair<double, pll_tree_edge_t *> > &edgeLengths, vector<pll_tree_edge_t *> &edges)
{
    TreeNode *u, *v;
    vector<double> cumulativeBranchLengths;
    pll_tree_edge * edge;
    //vector<pll_tree_edge_t *> edges;
    set<pll_tree_edge_t *> visitedEdges;
    set<pll_rnode_t *> visitedNodes;
    cumulativeBranchLengths.push_back(0);
    long double  cumBranchLength=0;
    long double  totalBranchLength1=0;
    int totalNumberNodes  =initialRootedTree->tip_count+ initialRootedTree->inner_count;
    edgeLengths.clear();
    for (unsigned int i = 0; i < totalNumberNodes; ++i)
    {
        auto node =initialRootedTree->nodes[i];
        if (node->parent !=NULL && visitedNodes.count(node) == 0  )
        {
            edge = new pll_tree_edge_t();
            u= (TreeNode *)(node->data);
            v= (TreeNode *)(node->parent->data);
            edge->edge.rtree.parent = u->timeInputTreeUnits >= v->timeInputTreeUnits ? node :  node->parent;
            edge->edge.rtree.child= u->timeInputTreeUnits < v->timeInputTreeUnits ? node :  node->parent;
            if (abs(abs(u->timeInputTreeUnits - v->timeInputTreeUnits ) - node->length)> 0.000001)
            {
                // printf( "\n the lengths are different! %s \n", node->label );
            }
            edge->length = node->length;
            //visitedNodes.insert(node->parent);
            totalBranchLength1 += edge->length;
            // if(visitedEdges.count(edge) ==0)
            // {
            if (node->label && std::string(node->label).compare(healthyCellLabel)!=0)//not  the branch from the healthyTip
            {
                //branchLengths.push_back(edge->length);
                //                if (edgesCumulativeLengths.size()==0)
                //                    edgesCumulativeLengths.push_back(make_pair(0, edge));
                //                else
                //                {
                edgeLengths.push_back(make_pair(edge->length, edge));
                //}
                //edges.push_back(edge);
            }
            else if(!node->label
                    //&& node->parent->parent !=NULL
                    )//not the branch to the root
            {
                //branchLengths.push_back(edge->length);
                edgeLengths.push_back(make_pair(edge->length, edge));
                // edges.push_back(edge);
            }
            edges.push_back(edge);
            cumBranchLength = edge->length + cumulativeBranchLengths.at(cumulativeBranchLengths.size()-1);
            cumulativeBranchLengths.push_back(cumBranchLength);
            //}
        }
        visitedNodes.insert(node);
    }
}
std::map<pll_rnode_t*, vector<Population*> >  Chain::initTimeOfOriginsOnRootedTree( vector<pair<double, pll_tree_edge_t *> > edgeLengths, int numberPoints, string &healthyCellLabel, MCMCoptions &mcmcOptions, gsl_rng * randomGenerator)
{
    std::unordered_set<pll_rnode_t *>  ancestorsOfTimeOfOrigins(numberPoints);
    std::unordered_set<pll_rnode_t *>  MRCAs(numberPoints);
    TreeNode *u, *v;
    long double  proposedTime;
    //vector<double> cumulativeBranchLengths;
    
    //vector<pll_tree_edge_t *> edges;
    set<pll_tree_edge_t *> visitedEdges;
    set<pll_rnode_t *> visitedNodes;
    //cumulativeBranchLengths.push_back(0);
    long double  random;
    std::map<pll_rnode_t*, vector<Population*>> mrcaOfPopulation;
    
    vector<double> branchLengths;
    
    if (numberPoints >0  )
    {
        std::transform(std::begin(edgeLengths), std::end(edgeLengths),
                       std::back_inserter(branchLengths), [](auto const& pair){return pair.first;}
                       );
        
        long double  maxCumBranchLength ;
        vector<double>  cumBranchLengths( branchLengths.size(),0);
        std::partial_sum(branchLengths.begin(), branchLengths.end(), cumBranchLengths.begin(),  plus<double>());
        
        cumBranchLengths.insert(cumBranchLengths.begin(), 0);
        
        maxCumBranchLength = cumBranchLengths.back(); //*std::max_element(cumBranchLengths.begin(), cumBranchLengths.end());
        
        std::transform(cumBranchLengths.begin(), cumBranchLengths.end(), cumBranchLengths.begin(),[maxCumBranchLength](long double  a) {return a /maxCumBranchLength; } );
        
        long double  cumBranchLengthsArray[cumBranchLengths.size()];
        std::copy(cumBranchLengths.begin() ,cumBranchLengths.end()
                  ,cumBranchLengthsArray );
        int  nextEvent;
        vector<int> eventIds;
        pll_rnode_t * ancestorMRCA;
        pll_rnode_t * MRCA;
        long double  proportionInsideChosenEdge=0.0;
        int assignedPop=0;
        while(assignedPop < numberPoints)
        {
            //do{
            random =Random::randomUniformFromGsl2(randomGenerator);
            nextEvent = bbinClones(random, cumBranchLengthsArray, cumBranchLengths.size());
            //}
            // while(std::find(eventIds.begin(), eventIds.end(), nextEvent) != eventIds.end());
            
            printf( "\n attempted edge id: %d out of  %lu and length %lf \n",nextEvent-1, cumBranchLengths.size()-1, branchLengths.at(nextEvent-1) );
            
            proportionInsideChosenEdge = (random - cumBranchLengths[nextEvent-1]) /(cumBranchLengths[nextEvent]-cumBranchLengths[nextEvent-1]);
            
            if (proportionInsideChosenEdge < 0 || proportionInsideChosenEdge > 1 )
                printf( "\n bad proportion inside proposed edge  \n" );
            
            if (assignedPop ==0)
                nextEvent =15;
            else if (assignedPop ==1)
                nextEvent =15;
            
            ancestorMRCA =edgeLengths.at(nextEvent-1).second->edge.rtree.parent;
            MRCA =edgeLengths.at(nextEvent-1).second->edge.rtree.child;
            
            eventIds.push_back(nextEvent);
            ancestorsOfTimeOfOrigins.insert(ancestorMRCA);
            
            auto popI=getPopulationbyIndex(MRCAs.size());
            MRCAs.insert(MRCA);
            popI->rMRCA = MRCA;
            
            u= (TreeNode *)(MRCA->data);
            v= (TreeNode *)(MRCA->parent->data);
            // proposedTime= u->timeInputTreeUnits+ (v->timeInputTreeUnits- u->timeInputTreeUnits)*proportionInsideChosenEdge;
            proposedTime= u->timeInputTreeUnits+ (v->timeInputTreeUnits- u->timeInputTreeUnits)*Random::randomUniformFromGsl2(randomGenerator);
            
            if (proposedTime<=0)
                printf( "time of origin is positive");
            popI->timeOriginInput = proposedTime  ;
            popI->scaledtimeOriginInput = proposedTime / theta;
            
            if (mcmcOptions.verbose>=1)
                printf( "\n MRCA node id %d with input time %Lf and scaled %Lf and parent with input time %Lf and scaled %Lf was assigned to the oldest pop %d and time of origin %Lf, scaled %Lf \n", MRCA->node_index,u->timeInputTreeUnits, u->timePUnits, v->timeInputTreeUnits,  v->timePUnits, popI->index, popI->timeOriginInput, popI->scaledtimeOriginInput  );
            
            mrcaOfPopulation[MRCA].push_back(popI);
            assignedPop++;
            
            if (mrcaOfPopulation[MRCA].size()>1)
                sort(mrcaOfPopulation[MRCA].begin(), mrcaOfPopulation[MRCA].end(), comparePopulationsByTimeOrigin);
        }//while
    }
    addOldestPopulation(mrcaOfPopulation,  healthyCellLabel, mcmcOptions);
    return(mrcaOfPopulation);
}
void Chain::addOldestPopulation(std::map<pll_rnode_t*, vector<Population*> > &mrcaOfPopulation, string &healthyCellLabel, MCMCoptions &mcmcOptions)
{
    TreeNode *u, *v;
    //long double  proposedTime;
    Population* oldestPop=getPopulationbyIndex(numClones -1);//for the oldest population
    if (std::string(initialRootedTree->root->right->label).compare(healthyCellLabel)==0)
    {
        oldestPop->rMRCA =initialRootedTree->root->left;
        
    }
    else  if (std::string(initialRootedTree->root->left->label).compare(healthyCellLabel)==0)
    {
        oldestPop->rMRCA =initialRootedTree->root->right;
        
    }
    else{
        //this should not happen
        oldestPop->rMRCA =initialRootedTree->root;
        u= (TreeNode *)(initialRootedTree->root->data);
        printf( "\n Error, the healthy cell is not a direct child of the tree root \n"  );
    }
    u= (TreeNode *)(oldestPop->rMRCA->data);
    v= (TreeNode *)(oldestPop->rMRCA->parent->data);
    mrcaOfPopulation[oldestPop->rMRCA].push_back(oldestPop);
    if (mrcaOfPopulation[oldestPop->rMRCA].size()>1)
        sort(mrcaOfPopulation[oldestPop->rMRCA].begin(), mrcaOfPopulation[oldestPop->rMRCA].end(), comparePopulationsByTimeOrigin);
    
}
void Chain::initOriginTimeOldestPopulation( string &healthyCellLabel, MCMCoptions &mcmcOptions){
    TreeNode *u, *v;
 
    Population* oldestPop=getPopulationbyIndex(numClones -1);//for the oldest population
  
     oldestPop->scaledtimeOriginInput = oldestPop->timeOriginSTD * oldestPop->x;
    oldestPop->timeOriginInput =oldestPop->scaledtimeOriginInput *theta;
    u= (TreeNode *)(oldestPop->rMRCA->data);
    v= (TreeNode *)(oldestPop->rMRCA->parent->data);
    if (mcmcOptions.verbose>=1)
       printf( "\n MRCA node id %d with input time %.10Lf and scaled %.10Lf and parent with input time %.10Lf and scaled %.10Lf was assigned to the oldest pop %d and time of origin %.10Lf, scaled %.10Lf std %.10Lf \n", oldestPop->rMRCA->node_index,u->timeInputTreeUnits, u->timePUnits, v->timeInputTreeUnits,  v->timePUnits, oldestPop->index, oldestPop->timeOriginInput, oldestPop->scaledtimeOriginInput,  oldestPop->timeOriginSTD  );
    //update the father node of the MRCA if needed
    if (v->timeInputTreeUnits < oldestPop->timeOriginInput){
        printf( "\n updated node from %.10Lf to %.10Lf", v->timeInputTreeUnits, oldestPop->timeOriginSTD);
        v->timeInputTreeUnits = oldestPop->timeOriginInput;
        v->timePUnits =  v->timeInputTreeUnits / theta;
        oldestPop->rMRCA->parent->length =  v->timeInputTreeUnits - u->timeInputTreeUnits ;
        oldestPop->rMRCA->length= v->timeInputTreeUnits - u->timeInputTreeUnits;
    }
    
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
            currentPopulation->numCompletedCoalescences++;
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
bool Chain::initPopulationsSampleSizes( std::map<pll_rnode_t*, vector<Population*>>  rmrcaOfPopulation, string &healthyTipLabel){
    bool existsZeroSampleSizePop=false;
    //      if (rmrcaOfPopulation.count(initialRootedTree->root)==0)
    //          printf("Error, the root node is not the MRCA of any population");
    initPopulationSampleSizesFromNodeOnRootedTree(initialRootedTree->root, NULL, rmrcaOfPopulation,healthyTipLabel );
    TreeNode *treeNode;
    for (unsigned int i = 0; i < numClones; ++i)
    {
        auto pop=populations[i];
        treeNode = (TreeNode *)(pop->rMRCA->data);
        pop->sampleSize =treeNode->numberTipsByPopulation[pop->index];
        existsZeroSampleSizePop =   pop->sampleSize <=0;
    }
    return existsZeroSampleSizePop;
}

void Chain::initPopulationsCoalescentAndMigrationEvents(std::map<pll_unode_t*, Population*> mrcaOfPopulation ){
    
    for (unsigned int i = 0; i < numClones; ++i)
    {
        auto pop = getPopulationbyIndex(i);
        pop->resetMigrationsList();//this requires  that you have set first timeOriginSTD in the population and effective population size
        initPopulationCoalescentAndMigrationEventsFromRootNodeOnTree(pop->nodeletMRCA, pop, mrcaOfPopulation );
    }
}
void Chain::initPopulationsCoalescentAndMigrationEventsFromRootedTree(std::map<pll_rnode_t*, vector<Population*>> rmrcaOfPopulation, string& healthyTipLabel){
    
    for (unsigned int i = 0; i < numClones; ++i)
    {
        auto pop = getPopulationbyIndex(i);
        pop->rtips.clear();
        pop->numCompletedCoalescences =0;
        // pop->resetMigrationsList();//this requires  that you have set first timeOriginSTD in the population and effective population size
        //pop->InitCoalescentEvents(numClones);
    }
    initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(rootRootedTree, NULL, rmrcaOfPopulation , healthyTipLabel);
}
void Chain::initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(pll_rnode_t *p, Population *currentPopulation, std::map<pll_rnode_t*, vector<Population*> > &rmrcaOfPopulation, string& healthyTipLabel ){
    
    if (p!=NULL)
    {
        if (p->left == NULL && std::string(p->label).compare(healthyTipLabel)!=0 )//tip different from healthy tip
        {
            currentPopulation->rtips.push_back(p);
        }
        else if (p->left == NULL && std::string(p->label).compare(healthyTipLabel)==0)
        {
            return;
            //nothing to do with the healthy cell
        }
        else
        {
            //initialize currentPopulation
            if (p->parent == NULL && (std::string(p->right->label).compare(healthyTipLabel)==0) && rmrcaOfPopulation.count(p->left)>0)
            {
                updateFatherPopOnSameEdge(rmrcaOfPopulation, p->left,NULL);
                
                currentPopulation = rmrcaOfPopulation[p->left].front();
                initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(p->left, currentPopulation, rmrcaOfPopulation , healthyTipLabel);
                return;
            }
            else if (p->parent == NULL && (std::string(p->left->label).compare(healthyTipLabel)==0) && rmrcaOfPopulation.count(p->right)>0)
            {
                updateFatherPopOnSameEdge(rmrcaOfPopulation, p->right,NULL);
                currentPopulation = rmrcaOfPopulation[p->right].front();
                initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(p->right, currentPopulation, rmrcaOfPopulation , healthyTipLabel);
                return;
            }
            //if is not the root ot the tree
            TreeNode * u= (TreeNode *)p->data;
            if (p->parent != NULL &&  currentPopulation != NULL )
            {
                //if is not the root of the tree(the parent of the MRCA of the oldest population
                currentPopulation->CoalescentEventTimes.push_back(u->timePUnits / (currentPopulation->x) );
                //printf("\n coalescent event node index %d for pop %d, time input tree %lf, scaled time %lf, model time  %lf\n ", p->clv_index,currentPopulation->index, u->timeInputTreeUnits, u->timePUnits, u->timePUnits / (currentPopulation->effectPopSize));
                currentPopulation->numCompletedCoalescences= currentPopulation->numCompletedCoalescences +1;
            }
            if(p->left !=NULL && rmrcaOfPopulation.count(p->left) != 0 )
            { // in this branch we have  migration
                
                updateFatherPopOnSameEdge(rmrcaOfPopulation, p->left,currentPopulation);
                
                initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(p->left, rmrcaOfPopulation[p->left].front(), rmrcaOfPopulation, healthyTipLabel);
            }
            else{
                initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(p->left, currentPopulation, rmrcaOfPopulation, healthyTipLabel);
            }
            if(p->right !=NULL && rmrcaOfPopulation.count(p->right) != 0 )
            { // in this branch we have  migration
                updateFatherPopOnSameEdge(rmrcaOfPopulation, p->right, currentPopulation);
                
                initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(p->right, rmrcaOfPopulation[p->right].front(), rmrcaOfPopulation, healthyTipLabel);
            }
            else
            {
                initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(p->right, currentPopulation, rmrcaOfPopulation, healthyTipLabel);
            }
        }
    }
}
void Chain::updateFatherPopOnSameEdge(std::map<pll_rnode_t*, vector<Population*> > &rmrcaOfPopulation, pll_rnode_t *p, Population * populationOfCurrentNode)
{
    int numberPop =((vector<Population*>)rmrcaOfPopulation[p]).size();
    Population *currentFatherPopulation;
    Population *currentPopulation;
    if (numberPop >0)
    {
        for (unsigned i=(numberPop -1); i >=1 ; i--)
        {
            currentPopulation=rmrcaOfPopulation[p].at(i-1);
            currentFatherPopulation =rmrcaOfPopulation[p].at(i);
            currentPopulation->FatherPop =currentFatherPopulation;
            Population::UpdateListMigrants(numClones, currentPopulation, currentFatherPopulation);
        }
        currentPopulation=rmrcaOfPopulation[p].back();
        currentFatherPopulation = populationOfCurrentNode;
        currentPopulation->FatherPop =currentFatherPopulation;
        if (currentFatherPopulation !=NULL)
        {     Population::UpdateListMigrants(numClones, currentPopulation, currentFatherPopulation);
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
void Chain::initPopulationSampleSizesFromNodeOnRootedTree(pll_rnode_t *p, Population *currentPopulation, std::map<pll_rnode_t*, vector<Population*>>  rmrcaOfPopulation, string &healthyTipLabel)
{
    if (p!=NULL)
    {
        TreeNode * treeNode=(TreeNode *)(p->data);
        //initialize currentPopulation
        if (p->parent == NULL && (std::string(p->right->label).compare(healthyTipLabel)==0) && rmrcaOfPopulation.count(p->left) >0)
        {
            
            currentPopulation = getYoungestPopulationOnEdge(p->left, rmrcaOfPopulation);
            initPopulationSampleSizesFromNodeOnRootedTree(p->left, currentPopulation, rmrcaOfPopulation, healthyTipLabel);
            return;
            
        }
        else if (p->parent == NULL && (std::string(p->left->label).compare(healthyTipLabel)==0) && rmrcaOfPopulation.count(p->right) >0)
        {
            currentPopulation = getYoungestPopulationOnEdge(p->right, rmrcaOfPopulation);
            initPopulationSampleSizesFromNodeOnRootedTree(p->right, currentPopulation, rmrcaOfPopulation, healthyTipLabel);
            return;
        }
        
        if (p->left == NULL)//tip
        {
            int idx= currentPopulation->index;
            if (idx <= numClones -1)
            {
                treeNode->resetNumberTipsVector(numClones);
                treeNode->numberTipsByPopulation.at(currentPopulation->index) = 1;
                return;
            }
            else
            {
                printf("error the index is out of bounds");
            }
        }
        else
        {
            if (rmrcaOfPopulation.count(p->left)>0 )
            {
                
                initPopulationSampleSizesFromNodeOnRootedTree(p->left, ((vector<Population*>)rmrcaOfPopulation[p->left]).front(), rmrcaOfPopulation, healthyTipLabel );
            }
            else{
                initPopulationSampleSizesFromNodeOnRootedTree(p->left, currentPopulation, rmrcaOfPopulation, healthyTipLabel );
            }
            if (rmrcaOfPopulation.count(p->right)>0 )
            {
                
                initPopulationSampleSizesFromNodeOnRootedTree(p->right, ((vector<Population*>)rmrcaOfPopulation[p->right]).front(), rmrcaOfPopulation, healthyTipLabel );
            }
            else
            {
                initPopulationSampleSizesFromNodeOnRootedTree(p->right, currentPopulation, rmrcaOfPopulation, healthyTipLabel);
            }
            TreeNode * treeNodeLeft=(TreeNode *)(p->left->data);
            TreeNode * treeNodeRight=(TreeNode *)(p->right->data);
            std::transform (treeNodeLeft->numberTipsByPopulation.begin(), treeNodeLeft->numberTipsByPopulation.end(), treeNodeRight->numberTipsByPopulation.begin(), treeNode->numberTipsByPopulation.begin(), std::plus<int>());
            // printf("\n node %d has %d tips below \n", p->node_index,treeNode->numberOfTipsSubTree );
            // printf("node %d has %d accumulated tips below \n", p->node_index,accumulate(treeNode->numberTipsByPopulation.begin(),treeNode->numberTipsByPopulation.end(),0));
        }
    }
}
Population* Chain::getYoungestPopulationOnEdge(pll_rnode_t* p, std::map<pll_rnode_t*, vector<Population*> >  rmrcaOfPopulation)
{
    Population* youngestPop=NULL;
    vector<Population*> vec;
    int numberOccurNode = rmrcaOfPopulation.count(p);
    int numberPoints;
    if (p!=NULL &&  numberOccurNode > 0)
    {
        vec = rmrcaOfPopulation[p];
        numberPoints =  vec.size();
        youngestPop = vec.front();
    }
    return youngestPop;
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
    long double  randomDelta;
    for( i = 0 ; i < numClones; i++)
    {
        popI=populations[i];
        randomDelta = Random::RandomLogUniform(mcmcOptions.Deltafrom, mcmcOptions.Deltato, mcmcOptions.useGSLRandomGenerator, NULL);
        //popI->delta = chain->proportionsVector[i] * randomDelta;
        popI->delta =  randomDelta;
        popI->growthRate =popI->delta  / popI->effectPopSize;
        popI->popSize=popI->effectPopSize * popI->birthRate;
        popI->deathRate= popI->birthRate - popI->growthRate;
    }
}
void Chain::samplePopulationGrowthRateFromPriors(MCMCoptions &mcmcOptions, gsl_rng * randomGenerator )
{
    int i;
    Population *popI;
    long double  randomGrowthRate;
    for( i = 0 ; i < numClones; i++)
    {
        popI=populations[i];
        if (mcmcOptions.priorsType ==0)
        {
            randomGrowthRate = Random::RandomLogUniform(mcmcOptions.GrowthRatefrom, mcmcOptions.GrowthRateto, mcmcOptions.useGSLRandomGenerator, randomGenerator);
        }
        if (mcmcOptions.priorsType ==2)
        {
            randomGrowthRate =Random::RandomPowerLawDistribution(mcmcOptions.parameterPowerLawDistributionGrowthRate, 0, mcmcOptions.useGSLRandomGenerator);
            
        }
        if (mcmcOptions.priorsType ==1)
        {
            randomGrowthRate =Random::RandomExponentialStartingFrom(mcmcOptions.lambdaExponentialPriorGrowthRate,0, mcmcOptions.useGSLRandomGenerator, randomGenerator);
        }
        
        popI->deltaT =  randomGrowthRate;
        
        if (mcmcOptions.verbose>=2)
            fprintf(stderr, "\n initial  r %.10Lf \n",popI->deltaT );
        
        
        if (mcmcOptions.verbose>=2)
            fprintf(stderr, "\n initial  delta  %.10Lf \n",popI->delta );
    }
}
void Chain::rescaleRootedTreeBranchLengths(long  double  mutationRate)
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
    
    
    
    Population *popI;
    // save the current values of
    oldTotalEffectPopSize = totalEffectPopSize;
    long double  newTotalEffectPopulationSize= proposalSlidingWindow(oldTotalEffectPopSize,  mcmcOptions.slidingWindowSizeTotalEffectPopSize);
    
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
    
    long double  newLogConditionalLikelihoodTree= LogConditionalLikelihoodTree(programOptions,mcmcOptions);
    
    //fprintf (stderr, "\n>> New log conditional Likelihood tree after the new total effective population size,  of the chain %d is = %lf  \n", chainNumber,newLogConditionalLikelihoodTree );
    //    char *newickString2;
    //    newickString2=NULL;
    //    newickString2 = toNewickString2 ( oldroot, programOptions.mutationRate,     programOptions.doUseObservedCellNames);
    //    printf("\n newick after move= %s  \n", newickString2);
    //    newickString2 = toNewickString4 ( oldroot, programOptions.mutationRate,     programOptions.doUseObservedCellNames);
    //    printf("\n newick after move= %s  \n", newickString2);
    //    free(newickString2);
    //    newickString2=NULL;
    long double  priorDensityNewTotalEffectivePopulationSize= Distributions::LogUniformDensity(newTotalEffectPopulationSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    long double  priorDensityCurrentTotalEffectivePopulationSize= Distributions::LogUniformDensity(oldTotalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    
    long double  sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewTotalEffectivePopulationSize;
    
    long double  sumLogDenominators=currentlogConditionalLikelihoodTree +priorDensityCurrentTotalEffectivePopulationSize;
    
    long double  randomNumber= Random::randomUniformFromGsl();
    
    long double  LogAcceptanceRate = (sumLogNumerators - sumLogDenominators) >0? (sumLogNumerators - sumLogDenominators) :0;
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
long double   Chain::proposalSlidingWindow( long double  oldvalue,  long double  windowSize)
{
    long double  newvalue =0;
    newvalue = oldvalue + (Random::randomUniformFromGsl()-0.5) * 0.5 * windowSize ;//windowSize is the distance  from the minimum to the  maximum value
    if (newvalue <0 )
        newvalue = -  newvalue;
    return newvalue;
}
void Chain::newProportionsVectorMove(ProgramOptions &programOptions, char *ObservedCellNames[], pll_msa_t * msa, MCMCoptions &mcmcOptions, vector<int> &sampleSizes)
{
    int i,j;
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
    long double  newLogConditionalLikelihoodTree= LogConditionalLikelihoodTree(programOptions,mcmcOptions);
    //fprintf (stderr, "\n>> New log conditional Likelihood tree after the new total effective population size,  of the chain %d is = %lf  \n", chainNumber,newLogConditionalLikelihoodTree );
    //    char *newickString2;
    //    newickString2=NULL;
    //    newickString2 = toNewickString2 ( chain->oldroot, programOptions->mutationRate,     programOptions->doUseObservedCellNames);
    //    printf("\n newick after move= %s  \n", newickString2);
    //    newickString2 = toNewickString4 ( chain->oldroot, programOptions->mutationRate,     programOptions->doUseObservedCellNames);
    //    printf("\n newick after move= %s  \n", newickString2);
    //    free(newickString2);
    //    newickString2=NULL;
    long double  priorDensityNewProportionsVector =  DirichletDensity( proportionsVector, oldproportionsVector , numClones);
    
    long double  priorDensityCurrentProportionsVector = DirichletDensity(oldproportionsVector, proportionsVector, numClones);
    
    long double  sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewProportionsVector;
    
    long double  sumLogDenominators=currentlogConditionalLikelihoodTree +priorDensityCurrentProportionsVector;
    
    long double  randomNumber= Random::randomUniformFromGsl();
    
    long double  LogAcceptanceRate = (sumLogNumerators - sumLogDenominators) >0? (sumLogNumerators - sumLogDenominators) :0;
    
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
void Chain::proposalProportionsVector(vector<long double  > &newProportionsvector, long double  tuningParameter )
{
    int i=0;
    long double  vectorForGamma[numClones];
    for( i = 0 ; i < numClones; i++)
    {
        vectorForGamma[i]= tuningParameter * proportionsVector[i];
    }
    double  *proportionsVectorArray=(double  *)&proportionsVector[0];
    Random::randomDirichletFromGsl(numClones, proportionsVectorArray, &newProportionsvector[0]);
    
}
long double  Chain::DirichletDensity(vector<long double> &proportionsVector,  vector<long double> &concentrationVector, int sizeVector)
{
    int i;
    long double  sum=0;
    long double  logResult=0;
    for( i = 0 ; i < sizeVector; i++)
    {
        sum = sum +concentrationVector[i];
        logResult= logResult+(concentrationVector[i]-1)*log(proportionsVector[i]);
        logResult= logResult-lgamma(concentrationVector[i]);
    }
    logResult = logResult+lgamma(sum);
    return logResult;
}
double  * Chain::expand_uniq_rates(int states, const  double  * uniq_rates, const int * rate_sym)
{
    unsigned int i;
    
    unsigned int num_rates = states * (states-1) / 2;
    //long double  * subst_rates = calloc(num_rates, sizeof(double));
    double  * subst_rates =(double  *) calloc(num_rates, sizeof(double));
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
void Chain::computeAvailableEdges( vector<pair<double, pll_tree_edge_t *> > &availableEdges, std::map<pll_rnode_t *, Population *> &currentMrcaOfPopulation, std::string &healthyCellLabel)
{
    pll_rnode_t *child;
    pll_rnode_t *parent;
    for (unsigned int i = 0; i < edgeLengths.size(); ++i)
    {
        auto edge = edgeLengths.at(i).second;
        parent = edge->edge.rtree.parent ;
        child = edge->edge.rtree.child;
        if (parent !=NULL && edge->length >0 && currentMrcaOfPopulation.count(child) == 0  )
        {
            //totalBranchLength1 += edge->length;
            if(child->left == NULL    && std::string(child->label).compare(healthyCellLabel)!=0 )//leaf different from healthy cell
            {
                availableEdges.push_back(make_pair(edge->length, edge));
            }
            else if(child->left != NULL)
            {//not a leaf
                availableEdges.push_back(make_pair(edge->length, edge));
            }
            //            cumBranchLength = edge->length + cumulativeBranchLengths.at(cumulativeBranchLengths.size()-1);
            //            cumulativeBranchLengths.push_back(cumBranchLength);
        }
    }
}

long double  Chain::computeAdjacentEdges( vector<pair<double, pll_tree_edge_t *> > &availableEdges, std::map<pll_rnode_t *, vector<Population *>> &currentMrcaOfPopulation, std::string &healthyCellLabel, Population *pop, pll_rnode_t * mrca)
{
    
    long double  sum =0.0;
    if (currentMrcaOfPopulation.count(mrca) ==0 )//if the dictionary doesnt have  the MRCA node we want to test then nothing to do (this should be an error
    {
        return sum;
    }
    else{// the MRCA of the population is in the dictionary
        availableEdges.clear();
        pll_tree_edge * edge;
        //upper branch
        if (mrca->parent != NULL &&  mrca->parent->parent != NULL)
            //not the edge from the tumor MRCA to the root
        {
            edge = new pll_tree_edge_t();
            edge->edge.rtree.child  =  mrca->parent;
            edge->edge.rtree.parent = mrca->parent->parent;
            edge->length = mrca->parent->length;
            sum += mrca->parent->length;
            availableEdges.push_back(make_pair(mrca->parent->length, edge));
        }
        else if(mrca->parent != NULL &&  mrca->parent->parent == NULL &&  currentMrcaOfPopulation[mrca].size() >=1 && pop->timeOriginInput == ((TreeNode *)mrca->parent->data)->timePUnits )// if it is the edge from the tumor MRCA to the root and  we have  at least 1 event and is the oldest population(the tiem of origin is the root of the tree)
        {
            //            edge = new pll_tree_edge_t();
            //            edge->edge.rtree.parent   =  mrca->parent;
            //
            //            edge->edge.rtree.parent = mrca->parent->parent;
            //            edge->length = mrca->parent->length;
            //            sum += mrca->parent->length;
            //            availableEdges.push_back(make_pair(mrca->parent->length, edge));
        }
        //left branch
        if (mrca->left != NULL && mrca->parent != NULL &&  mrca->parent->parent != NULL )//not a leaf and not the edge from the tumor MRCA to the root
        {
            edge = new pll_tree_edge_t();
            edge->edge.rtree.child  =  mrca->left;
            edge->edge.rtree.parent = mrca;
            edge->length = mrca->left->length;
            sum += mrca->left->length;
            availableEdges.push_back(make_pair(mrca->left->length, edge));
        }
        else if(mrca->left != NULL && mrca->parent != NULL &&  mrca->parent->parent == NULL  &&  currentMrcaOfPopulation[mrca].size() >=2 )// if it is the edge from the tumor MRCA to the root and  we have more than one event on that edge
        {
            edge = new pll_tree_edge_t();
            edge->edge.rtree.child  =  mrca->left;
            edge->edge.rtree.parent = mrca;
            edge->length = mrca->left->length;
            sum += mrca->left->length;
            availableEdges.push_back(make_pair(mrca->left->length, edge));
            
            
        }
        //right branch
        if (mrca->right != NULL && mrca->parent != NULL &&  mrca->parent->parent != NULL )//not a leaf and not the edge from the tumor MRCA to the root
        {
            edge = new pll_tree_edge_t();
            edge->edge.rtree.child  =  mrca->right;
            edge->edge.rtree.parent = mrca;
            edge->length = mrca->right->length;
            sum += mrca->right->length;
            availableEdges.push_back(make_pair(mrca->right->length, edge));
        }
        else if(mrca->right != NULL && mrca->parent != NULL &&  mrca->parent->parent == NULL &&  currentMrcaOfPopulation[mrca].size() >=2 && std::string(mrca->label).compare(healthyCellLabel)!=0)// if it is the edge from the tumor MRCA to the root and  we have more than one event on that edge
        {
            edge = new pll_tree_edge_t();
            edge->edge.rtree.child  =  mrca->right;
            edge->edge.rtree.parent = mrca;
            edge->length = mrca->right->length;
            sum += mrca->right->length;
            availableEdges.push_back(make_pair(mrca->right->length, edge));
        }
    }
    return sum;
}

std::map<pll_rnode_t*, vector<Population*>> Chain::chooseAvailableEdgeOnRootedTreeForPopulation(Population *pop, std::map<pll_rnode_t*, vector<Population*>> &currentMrcaOfPopulation, string &healthyCellLabel)
{
    std::unordered_set<pll_rnode_t *>  ancestorsOfTimeOfOrigins;
    std::unordered_set<pll_rnode_t *>  MRCAs;
    TreeNode *u, *v;
    //vector<double> cumulativeBranchLengths;
    set<pll_tree_edge_t *> visitedEdges;
    set<pll_rnode_t *> visitedNodes;
    //cumulativeBranchLengths.push_back(0);
    long double  random;
    
    vector<double> branchLengths;
    vector<pair<double, pll_tree_edge_t *> > availableEdges;
    
    computeAdjacentEdges( availableEdges, currentMrcaOfPopulation, healthyCellLabel, pop, pop->rMRCA);
    std::map<pll_rnode_t*, vector<Population*>> copyMRCAOfPopulation;
    branchLengths.clear();
    std::transform(std::begin(availableEdges), std::end(availableEdges),
                   std::back_inserter(branchLengths), [](auto const& pair){return pair.first;}
                   );
    ////////////////////////////////////
    if (availableEdges.size()>0)
    {
        long double  maxCumBranchLength ;
        vector<double>  cumBranchLengths( branchLengths.size());
        std::partial_sum(branchLengths.begin(), branchLengths.end(), cumBranchLengths.begin(), plus<double>());
        
        cumBranchLengths.insert(cumBranchLengths.begin(), 0);
        
        maxCumBranchLength = *std::max_element(cumBranchLengths.begin(), cumBranchLengths.end());
        
        std::transform(cumBranchLengths.begin(), cumBranchLengths.end(), cumBranchLengths.begin(),[maxCumBranchLength](long double  a) {return a /maxCumBranchLength; } );
        
        long double  cumBranchLengthsArray[cumBranchLengths.size()];
        std::copy(cumBranchLengths.begin() ,cumBranchLengths.end()
                  ,cumBranchLengthsArray );
        int  nextEvent;
        vector<int> eventIds;
        pll_rnode_t * ancestorMRCA;
        pll_rnode_t * MRCA;
        
        
        map<pll_rnode_t*,Population *>::iterator it;
        eventIds.clear();
        
        copyMRCAOfPopulation.clear();
        copyMRCAOfPopulation.insert(currentMrcaOfPopulation.begin(), currentMrcaOfPopulation.end());
        
        if (copyMRCAOfPopulation[pop->oldrMRCA].size() > 1)
        {
            copyMRCAOfPopulation[pop->oldrMRCA].erase(std::remove(copyMRCAOfPopulation[pop->oldrMRCA].begin(), copyMRCAOfPopulation[pop->oldrMRCA].end(), pop), copyMRCAOfPopulation[pop->oldrMRCA].end());
            
        }
        else if (copyMRCAOfPopulation[pop->oldrMRCA].size() == 1)
        {
            copyMRCAOfPopulation.erase(pop->oldrMRCA);
        }
        
        random =Random::randomUniformFromGsl();
        nextEvent = bbinClones(random, cumBranchLengthsArray, cumBranchLengths.size());
        MRCA =availableEdges.at(nextEvent-1).second->edge.rtree.child;
        
        ancestorMRCA =availableEdges.at(nextEvent-1).second->edge.rtree.parent;
        MRCA =availableEdges.at(nextEvent-1).second->edge.rtree.child;
        eventIds.push_back(MRCA->node_index);
        
        pop->rMRCA=MRCA;
        
        //after selecting the  edge we can use the same random value to
        // compute the new origin on the selected edge
        
        u= (TreeNode *)(MRCA->data);
        v= (TreeNode *)(MRCA->parent->data);
        
        long double  proposedTime;
        
        long double  proportionInsideChosenEdge = (random - cumBranchLengths[nextEvent-1]) /(cumBranchLengths[nextEvent]-cumBranchLengths[nextEvent-1]);
        
        if (proportionInsideChosenEdge < 0 || proportionInsideChosenEdge > 1 )
            printf( "\n bad proportion isside proposed edge  \n" );
        if (v->timeInputTreeUnits- u->timeInputTreeUnits >=0)
            proposedTime= u->timeInputTreeUnits+ (v->timeInputTreeUnits- u->timeInputTreeUnits)*proportionInsideChosenEdge;
        else
            proposedTime= v->timeInputTreeUnits+ (u->timeInputTreeUnits- v->timeInputTreeUnits)*(1-proportionInsideChosenEdge);
        
        /////////////////////////////////////
        
        pop->timeOriginInput =proposedTime;
        pop->scaledtimeOriginInput = pop->timeOriginInput / theta;
        
        copyMRCAOfPopulation[MRCA].push_back(pop);
        if (copyMRCAOfPopulation[pop->rMRCA].size()>1)
            sort(copyMRCAOfPopulation[pop->rMRCA].begin(), copyMRCAOfPopulation[pop->rMRCA].end(), comparePopulationsByTimeOrigin);
        printf( "\n MRCA node id %d with time %Lf and parent with time %Lf was assigned to pop %d with order %d and time of origin %Lf, scaled %Lf \n", MRCA->node_index, u->timePUnits, v->timePUnits, pop->index, pop->order,  pop->timeOriginInput, pop->scaledtimeOriginInput );
    }
    return copyMRCAOfPopulation;
}
bool Chain::isOldestPopulation(Population *pop, std::map<pll_rnode_t*, vector<Population*>> &currentMrcaOfPopulation)
{
    bool result = false ;
    if  (pop!= NULL && pop->rMRCA->parent!=NULL
         && pop->rMRCA->parent->parent == NULL )
    {
        if(currentMrcaOfPopulation[pop->rMRCA].size()>=1  && currentMrcaOfPopulation[pop->rMRCA].back() == pop )
        {
            result=true;
            
        }
    }
    return result;
}
long double  Chain::sumAvailableBranchLengths(std::map<pll_rnode_t*, vector<Population*>> currentMRCAPopulation)
{
    long double  result=0.0;
    pll_rnode_t* parent, *child;
    pll_tree_edge* edge;
    std::map<pll_rnode_t*, vector<Population*>>::iterator it;
    for (unsigned int i = 0; i < edges.size(); ++i)
    {
        edge =edges.at(i);
        parent = edge->edge.rtree.parent ;
        child = edge->edge.rtree.child;
        it= currentMRCAPopulation.find(child);
        if (parent !=NULL && edge->length >0 && it == currentMRCAPopulation.end() )
        {
            result+= result;
        }
    }
    return result;
}
long double  Chain::sumAdjacentEdges(std::map<pll_rnode_t*, vector<Population*>> currentMRCAPopulation, pll_rnode_t* MRCA, Population *pop)
{
    long double  result=0.0;
    if (currentMRCAPopulation.count(pop->rMRCA) ==0)//if the dictionary doesnt have  the MRCA of the population then nothing to do (this should be an error
    {
        return result;
    }
    else{// the MRCA of the population is in the dictionary
        pll_rnode_t * mrca= pop->rMRCA;
        result=mrca->length;
        
    }
    return result;
}
std::map<pll_rnode_t*, vector<Population*>> Chain::chooseNewTimeofOriginOnEdge(Population *pop, MCMCoptions &mcmcOptions, gsl_rng* randomGenerator)
{
    TreeNode * u, *v;
    double m;
    long double randomNumber;
    //long double currentTransformedTime, newTransformedTime;
    std::map<pll_rnode_t*, vector<Population*>> copyMRCAOfPopulation;
    copyMRCAOfPopulation.clear();
    copyMRCAOfPopulation.insert(rMRCAPopulation.begin(), rMRCAPopulation.end());
    
    
    u= (TreeNode *)( pop->rMRCA->data);
    v= (TreeNode *)( pop->rMRCA->parent->data);
    
    // pop->setLowerBoundTimeOriginInput(u->timeInputTreeUnits);
 
    long double proposedTime=0.0;
    
    if(mcmcOptions.kernelType ==1)
    {
        proposedTime= Random::randomNormalGreaterThan(pop->timeOriginInput, mcmcOptions.sigmaNormalKernelTimeofOrigin, u->timeInputTreeUnits);
    }
    else if(mcmcOptions.kernelType ==0)
    {
        if(mcmcOptions.useGSLRandomGenerator)
            randomNumber= Random::randomUniformFromGsl2(randomGenerator);
        else
            randomNumber= Random::randomUniformBoost();
        
        m = exp(2 * log(mcmcOptions.paramMultiplierTimeOriginOldestPop) *(randomNumber -0.5));
        
      
        proposedTime = m*  pop->timeOriginSTD;
        
    }
    //printf( "\n proposed time of origin input tree units %Lf,  old %Lf \n", proposedTime,pop->timeOriginInput );
   
    
    pop->timeOriginSTD =proposedTime;
    pop->scaledtimeOriginInput =pop->timeOriginSTD *pop->x;
    pop->timeOriginInput=pop->scaledtimeOriginInput * theta;
    
    //update the node if needed
    if (!mcmcOptions.noData && v->timeInputTreeUnits < pop->timeOriginInput)
    {
        updateNodeInfoOldestPopulation(pop, pop->timeOriginInput);
    }
    //copyMRCAOfPopulation[pop->oldrMRCA].push_back(pop);
    
    if (copyMRCAOfPopulation[pop->rMRCA].size()>1)
        sort(copyMRCAOfPopulation[pop->rMRCA].begin(), copyMRCAOfPopulation[pop->rMRCA].end(), comparePopulationsByTimeOrigin);
    
    return copyMRCAOfPopulation;
}
void Chain::PrepareFiles(const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files, int chainNumber)
{
    char File[MAX_NAME];
    char dir[MAX_NAME];
    
    if (chainNumber == 0)
        mkdir("Results", S_IRWXU); /* Create "Results" folder (with type S_IRWXU (read, write and execute)) */
    //mkdir("Results",0);
    
#ifdef MAC
    strcpy (dir, ":Results:"); /* Copy the string in char variable dir = Results (char), is different mac vs windows */
#else
    strcpy (dir, "Results/");
#endif
    
    sprintf(File,"%s/Chain%d", filePaths.resultsDir,  chainNumber );
    mkdir(File,S_IRWXU);
    //strcpy (resultsDir, dir);
    if (programOptions.doPrintTrees == YES)
    {
        sprintf(File,"%s/Chain%d/%s", filePaths.resultsDir,  chainNumber ,filePaths.treeDir);
        mkdir(File,S_IRWXU);
        sprintf(File,"%s/Chain%d/%s/%s.trees", filePaths.resultsDir, chainNumber, filePaths.treeDir,filePaths.treeFile);
        if (openFile(&files.fpTrees, File) == -1)
        {
            fprintf (stderr, "Can't open \"%s\"\n", File);
            exit(-1);
        }
        sprintf(File,"%s/Chain%d/%s/%s_2.trees", filePaths.resultsDir,chainNumber,  filePaths.treeDir, filePaths.treeFile);
        if (openFile(&files.fpTrees2, File) == -1)
        {
            fprintf(stderr, "Can't open %s.\n", File);
            exit(-1);
        }
    }
    if (programOptions.doPrintTimes == YES)
    {
        sprintf(File,"%s/Chain%d/%s", filePaths.resultsDir,chainNumber, filePaths.timesDir );
        mkdir(File,S_IRWXU);
        sprintf(File,"%s/Chain%d/%s/%s.txt", filePaths.resultsDir, chainNumber, filePaths.timesDir ,filePaths.timesFile);
        if (openFile(&files.fpTimes, File) == -1)
        {
            fprintf (stderr, "Can't open \"%s\"\n", File);
            exit(-1);
        }
        sprintf(File,"%s/%s", filePaths.resultsDir, filePaths.timesDir);
        sprintf(File,"%s/Chain%d/%s/%s_2.txt", filePaths.resultsDir,chainNumber, filePaths.timesDir, filePaths.timesFile);
        if (openFile(&files.fpTimes2, File) == -1)
        {
            fprintf(stderr, "Can't open %s.\n", File);
            exit(-1);
        }
    }
    //log file
    sprintf(File,"%s/Chain%d/%s.log", filePaths.resultsDir, chainNumber,  filePaths.logFile);
    if (openFile(&files.fplog, File) == -1)
    {
        fprintf (stderr, "Can't open \"%s\"\n", File);
        exit(-1);
    }
}
void Chain::writeMCMCState( int  currentIteration, const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files , MCMCoptions &mcmcOptions)
{
    //    fprintf (files.fplog, "%d  %lf %lf %lf effect_pop_size(pop1) effect_pop_size(pop2) effect_pop_size(pop3) sample_size(pop1)  sample_size(pop2) sample_size(pop3) time_origin(pop1)  time_origin(pop2) time_origin(pop3) \n", pop);
    string paramName;
    Population *popI;
    fprintf (files.fplog, "%d\t", currentIteration);
    fprintf (files.fplog, "%.5Lf\t", currentlogConditionalLikelihoodTree);
    fprintf (files.fplog, "%.5Lf\t", currentlogConditionalLikelihoodSequences);
    
    fprintf (files.fplog, "%.15Lf\t", theta);
    for (unsigned int i = 0; i < populations.size(); ++i){
        popI =populations[i];
        //fprintf (files.fplog, "%.8Lf\t", pop->effectPopSize);
        fprintf (files.fplog, "%.40Lf\t", popI->delta);
        //fprintf (files.fplog, "%.40Lf\t", pop->growthRate);
        fprintf (files.fplog, "%.10Lf\t", popI->x);
        fprintf (files.fplog, "%.40Lf\t", popI->deltaT);
        fprintf (files.fplog, "%.15Lf\t", popI->theta);
        // fprintf (files.fplog, "%.40Lf\t", mutationRate);
        fprintf (files.fplog, "%d\t", popI->sampleSize);
        fprintf (files.fplog, "%.15Lf\t", popI->timeOriginInput);//multiplied by mutation rate
        fprintf (files.fplog, "%.15Lf\t", popI->timeOriginSTD);//T
        fprintf (files.fplog, "%.15Lf\t", popI->scaledtimeOriginInput);//divided by theta_i
        fprintf (files.fplog, "%.15Lf\t", popI->timeOriginInput / mutationRate);//physical time(if lambda=1 is the time in generations)
        
    }
    
    if (currentIteration <= mcmcOptions.Niterations -1)
        fprintf (files.fplog, "\n");
}
void Chain::writeHeaderOutputChain(const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files )
{
    std::time_t timeNow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    fprintf (files.fplog, "#tumor_coal 1.0\n");
    fprintf (files.fplog, "# Generated %s\n", std::ctime(&timeNow));
    fprintf (files.fplog, "state\t");
    fprintf (files.fplog, "loglikTree\t");
    fprintf (files.fplog, "loglikSeq\t");
    string paramName;
    paramName =  "theta";
    fprintf (files.fplog, "%s\t", paramName.c_str());
    for (unsigned int i = 0; i < populations.size(); ++i){
        // paramName =  "effect_pop_size(pop" + std::to_string(i) + ")";
        // fprintf (files.fplog, "%s\t", paramName.c_str());
        paramName =  "delta(pop" + std::to_string(i) + ")";
        fprintf (files.fplog, "%s\t", paramName.c_str());
        // paramName =  "growth_rate(pop" + std::to_string(i) + ")";
        //fprintf (files.fplog, "%s\t", paramName.c_str());
        paramName =  "prop_effect_pop_size(pop" + std::to_string(i) + ")";
        fprintf (files.fplog, "%s\t", paramName.c_str());
        paramName =  "deltaT(pop" + std::to_string(i) + ")";
        fprintf (files.fplog, "%s\t", paramName.c_str());
        paramName =  "theta(pop" + std::to_string(i) + ")";
        fprintf (files.fplog, "%s\t", paramName.c_str());
        //paramName =  "mut_rate(pop" + std::to_string(i) + ")";
        //fprintf (files.fplog, "%s\t", paramName.c_str());
        paramName =  "sample_size(pop" + std::to_string(i) + ")";
        fprintf (files.fplog, "%s\t", paramName.c_str());
        paramName =  "time_origin_input(pop" + std::to_string(i) + ")";
        fprintf (files.fplog, "%s\t", paramName.c_str());
        paramName =  "time_origin_std(pop" + std::to_string(i) + ")";
        fprintf (files.fplog, "%s\t", paramName.c_str());
        paramName =  "scaled_time_origin(pop" + std::to_string(i) + ")";
        fprintf (files.fplog, "%s\t", paramName.c_str());
        paramName =  "physical_time(pop" + std::to_string(i) + ")";
        fprintf (files.fplog, "%s\t", paramName.c_str());
    }
    fprintf (files.fplog, "\n");
}
long double  Chain::autoCorrelation(int lag, vector<double> values)
{
    long double  result=0;
    long double  sumDenominator=0;
    long double  sumNumerator=0;
    if (values.size() >0)
    {
        long double  mean = accumulate( values.begin(), values.end(), 0.0)/ values.size();
        for (unsigned int i = 0; i <= values.size() -lag ; ++i)
        {
            sumNumerator = (values.at(i)- mean)*(values.at(i+lag)- mean);
            sumDenominator += pow(values.at(i)- mean,2);
        }
        for (unsigned int i = values.size() -lag +1; i <= values.size()  ; ++i)
        {
            sumDenominator += pow(values.at(i)- mean,2);
        }
        result = sumNumerator/ sumDenominator;
    }
    return result;
}
long double  Chain::ESS(int lag, vector<double> values)
{
    long double  result=0;
    long double  sumAutoCorrelations =0;
    for (unsigned int i = 0; i <= values.size()  ; ++i)
    {
        sumAutoCorrelations +=autoCorrelation(i, values);
    }
    result = values.size() / (1 + 2 * sumAutoCorrelations);
    return result;
}
bool Chain::checkMigrationsOrder()
{
    bool migrationsOrderCorrect = true;
    Population *popI;
    long double lastMigrationTime;
    for(unsigned  i = 0 ; i < numClones; i++)
    {
        popI = populations[i];
        lastMigrationTime= popI->immigrantsPopOrderedByModelTime.at(popI->immigrantsPopOrderedByModelTime.size()-1).first ;
        if (lastMigrationTime != popI->timeOriginSTD)//check the  last position
        {
            migrationsOrderCorrect=false;
            fprintf(stderr, "\n the last migration time  is not the time of origin \n");
            break;
        }
        for( unsigned i = 0 ; i < popI->immigrantsPopOrderedByModelTime.size()-1; i++)
        {
            if (popI->immigrantsPopOrderedByModelTime.at(i).first >= popI->immigrantsPopOrderedByModelTime.at(i+1).first )
            {
                migrationsOrderCorrect=false;
                fprintf(stderr,"\n bad order of migrations \n");
                break;
                
            }
        }
        if (!migrationsOrderCorrect)
            break;
    }
    
    return migrationsOrderCorrect;
}

void closeFiles(const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files, int chainNumber)
{
    fclose(files.fplog);
}
void Chain::updateNodeInfoOldestPopulation(Population * oldestPop, double newMRCAtimeInputTreeUnits){
    TreeNode *u, *v;
    u= (TreeNode *)( oldestPop->rMRCA->data);
    v= (TreeNode *)( oldestPop->rMRCA->parent->data);
     printf( "\n updated node from %.10Lf to %.10lf", v->timeInputTreeUnits, newMRCAtimeInputTreeUnits);
    v->timeInputTreeUnits = newMRCAtimeInputTreeUnits;
    
    v->timePUnits =  v->timeInputTreeUnits / theta;
    oldestPop->rMRCA->parent->length =  v->timeInputTreeUnits - u->timeInputTreeUnits ;
    oldestPop->rMRCA->length =  v->timeInputTreeUnits - u->timeInputTreeUnits ;
    
}
void Chain::drawModelTimeOriginFromConditionalDensity(Population * oldestPop, MCMCoptions &mcmcOptions,  gsl_rng * randomGenerator){
    
    long double modelTimeOriginOldestPop=Random::RandomDensityModelTimeOrigin(oldestPop->delta, mcmcOptions.useGSLRandomGenerator,  oldestPop->lowerBoundTimeOriginInput,  randomGenerator );
    
    oldestPop->timeOriginSTD = modelTimeOriginOldestPop;
    
    oldestPop->scaledtimeOriginInput=oldestPop->timeOriginSTD*(oldestPop->x);
    
    oldestPop->timeOriginInput =oldestPop->scaledtimeOriginInput * theta;
    
    
    updateNodeInfoOldestPopulation(oldestPop, oldestPop->timeOriginInput);
    
}
void Chain::printMovesSummary(){
    
    vector<MCMCmove*>::iterator it;
    MCMCmove* currentMove;
    string name;
    for ( it=moves.begin(); it!=moves.end(); ++it)
    {
        currentMove=  *it;
        name =currentMove->name();
        fprintf (stderr, "\n Number accepted moves of type %s: %d, number of rejected moves %d \n",name.c_str(), (int)(currentMove->numberAccepted()),(int)(currentMove->numberRejected()) );
    }
    
}
