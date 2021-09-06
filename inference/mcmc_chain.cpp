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
//

#include "mcmc_chain.hpp"

#include "data_utils.hpp"

#include "definitions.hpp"
#include "eigen.hpp"

#include "autocorrelation.hpp"
#include "tree_node.hpp"
#include "random.h"
#include "utils.hpp"


#include "treeLikelihood.hpp"


#include "constants.hpp"
#include "mcmc_move.hpp"
#include "output_functions.hpp"
#include "pll_utils.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/regex.h>
#include <boost/random.hpp>

#include <map>
#include <iterator>
#include <string>
#include <vector>
//#include <queue>
#include <unordered_set>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>


#include <functional>
#include <random>
#include <chrono>
#include <set>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <sys/stat.h>

extern "C"
    {
#include <libpll/pll_msa.h> //for pllmod_msa_empirical_frequencies
#include <libpll/pll_tree.h>
#include <libpll/pllmod_util.h>
#include <libpll/pllmod_common.h>
#include <libpll/pllmod_algorithm.h>
    
    //#include <gsl/gsl_sf_bessel.h>
    //#include <gsl/gsl_randist.h>
    //#include <gsl/gsl_cdf.h>
    
    }

//#include "data_utils.hpp"
//#include "definitions.hpp"
//#include "eigen.hpp"
//#include "tree_node.hpp"
//#include "random.h"
//#include "utils.hpp"
//#include "constants.hpp"
//#include "mcmc_move.hpp"
//#include "output_functions.hpp"



using namespace std;
using namespace std::placeholders;
//Parametrized Constructor

Chain::Chain( int chainNumber,
             int numClones,
             int gammaParam,
             long double  mutationRate,
             long double  seqErrorRate,
             long double  dropoutRate,
             MCMCoptions &mcmcOptions
             )

{
    converged = false;
    endReached =false;
    samplingFromPosteriorAfterConvergence=false;
    currentIdxMCMCParameters =0;
    long double zero =0.0;
    
    if (chainNumber >=0)
        this->chainNumber=chainNumber;
    
    if (numClones >=0)
        this->numClones=numClones;
    
    if (gammaParam >=0){
        this->gammaParam=gammaParam;
    }
    
    if (mutationRate >=0){
        this->mutationRate=mutationRate;
        mutationRatePar = new MCMCParameterWithKernel("ScaledMutRate",0.1,0.3,0.0);
        
        indexMCMCParameters.push_back(make_pair(mutationRatePar, currentIdxMCMCParameters));
        currentIdxMCMCParameters++;
        
    }
    
    if (seqErrorRate >0 && seqErrorRate <=1){
        this->seqErrorRate=seqErrorRate;
        seqErrorRatePar = new MCMCParameterWithKernel("SeqError",0.0, 0.3,0.0);
        
        indexMCMCParameters.push_back(make_pair(seqErrorRatePar, currentIdxMCMCParameters));
        currentIdxMCMCParameters++;
        
        long double alphaSeqError =1.2;
        long double betaSeqError =1;
        seqErrorRatePar->getParameter()->setPrior(&Distributions::LogBetaDensity);
        seqErrorRatePar->getParameter()->setThirdParLogPrior(betaSeqError);//beta parameter
        seqErrorRatePar->getParameter()->setSecondParLogPrior(alphaSeqError);//alpha parameter
        
        
    }
    
    if (dropoutRate >0 && dropoutRate <=1){
        this->dropoutRate=dropoutRate;
        dropoutRatePar = new MCMCParameterWithKernel("DropoutError",0.0, 0.3,0.0);
        
        indexMCMCParameters.push_back(make_pair(dropoutRatePar, currentIdxMCMCParameters));
        currentIdxMCMCParameters++;
        
        long double alphaDropout =1.2;
        long double betaDropout =1;
        dropoutRatePar->getParameter()->setPrior(&Distributions::LogBetaDensity);
        dropoutRatePar->getParameter()->setThirdParLogPrior(betaDropout);//beta parameter
        dropoutRatePar->getParameter()->setSecondParLogPrior(alphaDropout);//alpha parameter
        
    }
    
    
    currentNumberIerations = 0;
    
    currentlogConditionalLikelihoodTree=0;
    currentlogConditionalLikelihoodSequences=0;
    totalAccepted=0;
    totalRejected=0;
    
    theta = 0.0;
    oldtheta = 0.0;
    thetaPar = new MCMCParameterWithKernel("Theta",0.0, mcmcOptions.paramMultiplierMoveTheta,0.0);
    
    indexMCMCParameters.push_back(make_pair(thetaPar, currentIdxMCMCParameters));
    currentIdxMCMCParameters++;
    
    
    
    
    if(mcmcOptions.priorsType==0)
    {
        thetaPar->getParameter()->setPrior(&Distributions::LogUniformDensity);
        thetaPar->getParameter()->setSecondParLogPrior(mcmcOptions.MutRatefrom); //from parameter
        thetaPar->getParameter()->setThirdParLogPrior(mcmcOptions.MutRateto); //to parameter
    }
    if(mcmcOptions.priorsType==2)
    {
        thetaPar->getParameter()->setPrior(&Distributions::LogPowerLawDistibutionDensity);
        thetaPar->getParameter()->setSecondParLogPrior(mcmcOptions.parameterPowerLawDistributionMutationRate);
        thetaPar->getParameter()->setThirdParLogPrior(zero); //from parameter
    }
    
    if(mcmcOptions.priorsType==1) {
        thetaPar->getParameter()->setPrior(&Distributions::LogExponentialDensity);
        thetaPar->getParameter()->setThirdParLogPrior(mcmcOptions.lambdaExponentialPriorMutationRate);//lambda parameter
        thetaPar->getParameter()->setSecondParLogPrior(zero);//from parameter
    }
    
    
    currentNumberEdgesRootedTree=0;
    initProportionsVector();
    thinning=1;
    numberIndependentLongUpdates=0;
    initialRootedTree=NULL;
    
    files.fplog = new FilePath();
    files.fpTrees = new FilePath();
    files.fpTrees2 = new FilePath();
    files.fpTimes = new FilePath();
    files.fpTimes2 = new FilePath();
    
    Utils::init_to_empty_str(files.fplog->path);
    Utils::init_to_empty_str(files.fpTrees->path);
    Utils::init_to_empty_str(files.fpTrees2->path);
    Utils::init_to_empty_str(files.fpTimes->path);
    Utils::init_to_empty_str(files.fpTimes2->path);
}
int Chain::setInitialTreeFromNewick(char * NewickString){
    
    
    this->initialUnrootedTree = pll_utree_parse_newick_string_unroot(NewickString);
    this->numNodes = initialUnrootedTree->tip_count +  initialUnrootedTree->inner_count;
    
    if (this->initialUnrootedTree ==NULL)
        return 1;
    else
        return 0;
    
}
int Chain::setInitialRootedTreeFromNewick(char * NewickString){
    
    
    this->initialRootedTree = pll_rtree_parse_newick_string(NewickString);
    this->numNodes = initialRootedTree->tip_count +  initialRootedTree->inner_count;
    
    if (this->initialRootedTree ==NULL)
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
int Chain::setInitialTreeRootedTree(pll_rtree_t *rootedTree){
    if (initialRootedTree !=NULL)
    {
        initialRootedTree = rootedTree;
        
        numNodes = initialRootedTree->tip_count +  initialRootedTree->inner_count;
        rootRootedTree = initialRootedTree->root;
        
        pll_rtree_reset_template_indices(initialRootedTree->root, initialRootedTree->tip_count);
        return 0;
    }
    else
        return 1;
}
void Chain::MakeCoalescenceEvent( Population *population, vector<pll_unode_t *> &nodes, int numClones,const gsl_rng *randomGenerator, int noisy,   int &numActiveGametes, int &nextAvailable,
                                 int &labelNodes, long double  &currentTime, int &numNodes)
{
    
    // TreeNode  *p, *q, *r;
    pll_rnode_t  *p, *q,  *r ;
    TreeNode * tp, *tq, *tr ;
    int firstInd,  secondInd=0, newInd=0;
    int choosePairIndividuals = YES;
    
    population->ChooseRandomIndividual(&firstInd, numClones,   &secondInd, randomGenerator, choosePairIndividuals);
    
    newInd = nextAvailable;
    if (noisy > 1)
        fprintf (stderr, "Coalescence involving %d and %d to create node %d (in clone %d)", population->idsActiveGametes[firstInd], population->idsActiveGametes[secondInd], newInd, population->index);
    /*  set pointers between nodes */
    assert(firstInd< population->idsActiveGametes.size() );
    assert(secondInd< population->idsActiveGametes.size() );
    p = rnodes[population->idsActiveGametes[firstInd]];
    q = rnodes[population->idsActiveGametes[secondInd]];
    
    r = rnodes[newInd];
    
    r->node_index = nextAvailable;
    
    tr= (TreeNode *)(r->data);
    tp= (TreeNode *)(p->data);
    tq= (TreeNode *)(q->data);
    
    tr->index =   nextAvailable;
    tr->label = labelNodes;
    labelNodes=labelNodes+1;
    tr->indexOldClone = tr->indexCurrentClone = population->index;//here the clone number is updated
    
    // r->indexCurrentClone = p->indexCurrentClone;
    // r->orderCurrentClone = p->orderCurrentClone;
    tr->orderCurrentClone =population->order;
    //tr->effectPopSize=population->effectPopSize;
    tr->nodeClass  = 4;
    // link the nodes
    r->left = p;
    r->right = q;
    p->parent = r;
    q->parent = r;
    
    tr->time = currentTime;
    tr->scaledByThetaTimeInputTreeUnits =currentTime * theta;
    
    if (population->order==numClones-1)//oldest population
    {
        tr->timePUnits =currentTime * theta;
    }
    else{
        tr->timePUnits =currentTime * theta * population->x;
        
    }
    tr->timeInputTreeUnits = tr->timePUnits;
    
    p->length = tr->timePUnits -tp->timePUnits;
    q->length= tr->timePUnits -tq->timePUnits;
    
    //fprintf (stderr, "\n r->index = %d, r->time = %lf, ClonePopSizeMeffectBegin[ThisCloneNumber] = %lf, ThisCloneNumber = %d\n", r->index, r->time, ClonePopSizeMeffectBegin[ThisCloneNumber], ThisCloneNumber);
    if (noisy > 1)
        fprintf (stderr, "\t|\tCurrentTime (input units) = %Lf", tr->timePUnits);
    /* readjust active nodes */
    
    // idsActiveGametes[firstInd] = newInd;
    population->idsActiveGametes[firstInd] = r->node_index;// the r3 node
    population->idsActiveGametes[secondInd] = population->idsActiveGametes[population->numActiveGametes - 1];
    // for(int k=0; k < population ->idsActiveGametes.size();k++)
    //    fprintf (stderr, "\n pop of order  %d Active gamete id %d ", population->order, population->idsActiveGametes[k]);
    numActiveGametes = numActiveGametes - 1; /* less 1 active node */
    
    //update list ids nodes
    population->idsGametes[population->numGametes] = newInd;
    population->numGametes = population->numGametes +1;
    
    // *nextAvailable=*nextAvailable+1; /* 1 node more is available */
    nextAvailable=nextAvailable+1; /* 1 node more is available */
    
    population->CoalescentEventTimes[ population->numCompletedCoalescences]=  tr->time;
    population->numCompletedCoalescences= population->numCompletedCoalescences+1;
    
    population->numActiveGametes = population->numActiveGametes - 1; /* now this clone
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

void Chain::SimulatePopulation( Population *popI,
                               ProgramOptions &programOptions,
                               long int *seed,
                               const gsl_rng *randomGenerator,
                               int &numNodes,
                               int numClones,
                               int &nextAvailable,
                               int &numActiveGametes,
                               int &labelNodes,
                               long double  &currentTime,
                               int &eventNum, long double K)
{
    int   i,  k, isCoalescence,
    firstInd,  newInd,
    isMigration, whichClone;
    long double      eventTime;
    pll_rnode_t *  p, *r, *r1 ;
    TreeNode  *u,*v, *u1;
    
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
         std::cout<<  "\n\n>> Simulating evolutionary history of clone  " << popI->index << " (number active gametes " << popI->numActiveGametes << "original time to origin " << popI->timeOriginInput<< std::endl;
    
    if (programOptions.noisy > 1)
        std::cout<<  "\n\n> Simulating evolutionary history of clone  or order " << popI->order<< std::endl;
    
    currentTime=0;
    while (indexNextMigration < numMigrations) {
        timeNextMigration = popI->immigrantsPopOrderedByModelTime[indexNextMigration].first;
        //fprintf (stderr, "\n\n> numParcialActiveGametes= %d \n", numParcialActiveGametes);
        if ( popI->numActiveGametes >= 2) {
            ThisRateCA = (long double)  popI->numActiveGametes * ((long double) popI->numActiveGametes - 1) / 2.0;
            ThisTimeCA_W = Random::RandomExponentialStartingFrom (ThisRateCA,0,  true, randomGenerator, NULL) ;
            ThisTimeCA_V1 = Population::FmodelTstandard (currentTime, popI->timeOriginSTD, popI->delta, K);
            ThisTimeCA_V1 = ThisTimeCA_V1 + ThisTimeCA_W;
            // from standard time to model time, GstandardTmodel(V, T, delta)
            ThisTimeCA_V2 = Population::GstandardTmodel(ThisTimeCA_V1, popI->timeOriginSTD, popI->delta, K);
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
            whichClone = popI->index;
            currentTime = ThisTimeCA_V2; // update current time in model time
            eventTime = currentTime;
            
            if (programOptions.noisy > 1)
            {
                std::cout<<"\n\n*** Event "<< eventNum << " currentTime (model time) = "<< ThisTimeCA_V2 << " currentTime (standard time) = " << ThisTimeCA_V1 << std::endl;
                std::cout<< "\n\n*** Event " <<  eventNum << " currentTime (input units) =" << ThisTimeCA_V2<< std::endl;
            }
            if (programOptions.noisy == 4)
                std::cout << "* Coalescence *\n"<< std::endl;
            MakeCoalescenceEvent( popI, nodes, numClones, randomGenerator, programOptions.noisy, numActiveGametes, nextAvailable, labelNodes, currentTime,  programOptions.numNodes);
            
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
                     std::cout << "\n\n*** Event "<< eventNum << " currentTime (model units) = " << ThisTimeCA_V2<< std::endl;
                }
                if (programOptions.noisy == 4)
                {
                    std::cout <<"* Migration *\n"<< std::endl;
                    
                }
                
                incomingPop = popI->immigrantsPopOrderedByModelTime[indexNextMigration].second;
                
                //r->indexOldClone = incommingPop->index;
                p = rnodes[incomingPop->nodeIdAncestorMRCA]; // root of younger clone
                
                indexNextMigration = indexNextMigration + 1;
                
                v= (TreeNode *)(p->data);
                
                std::cout<< "\n The incoming population "<< incomingPop->order <<" to  population " << popI->order << " with node "<<v->index <<  " and nodelet "<< p->node_index <<" time " << v->timePUnits << std::endl;
                
                v->indexCurrentClone = popI->index;
                v->indexOldClone = incomingPop->index;
                v->orderCurrentClone = popI->order;
                
                
                k = v->indexCurrentClone;
                incomingPop->numActiveGametes = incomingPop->numActiveGametes - 1; /* now the other clone has 1 less node */
                // remove node from old clone in list of active gametes and add the new node of the current clone
                //popI->idsActiveGametes[popI->numActiveGametes]=r->index;//adding the superfluos node
                popI->idsActiveGametes[popI->numActiveGametes]=p->node_index;//adding the superfluos node
                
                popI->numActiveGametes = popI->numActiveGametes + 1; /* now this clone has 1 more node */
                //                if (noisy > 1)
                
                //for(int i=0; i < popI->idsActiveGametes.size();i++)
                //    fprintf (stderr, "\n pop of order  %d Active gamete id %d", popI->order, popI->idsActiveGametes[i]);
                
                if (programOptions.noisy > 1)
                    std::cout << "\t|\tCurrentTime (input units) =" <<  v->timePUnits << std::endl;
                // fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", r->timePUnits);
                /* memory for number of nodes */
                if (nextAvailable >= numNodes)  /* if there aren't enough nodes it go into and it addition more */
                {
                    /* ReallocNodes(&numNodes, activeGametes); */
                    if (programOptions.noisy == 4)
                       std::cout << "\n\n...Doing reallocation of nodes (Coalescence)\n"<< std::endl ;
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
                     std::cout <<  "\n\n*** Clone origin ***\n"<< std::endl;
                if (programOptions.noisy == 4)
                     std::cout <<  "Clone origin "<< popI->index << "at time " << currentTime << " in (model units) " << std::endl ;
                if (popI->order < numClones - 1) // do not do it for the last clone
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
                    u1->indexOldClone =u1->indexCurrentClone = popI->index;
                    
                    u1->orderCurrentClone = popI->order;
                    //u1->effectPopSize= popI->effectPopSize;
                    popI->nodeIdAncestorMRCA=newInd;
                    
                    u1->nodeClass  =4;
                    
                    //firstInd = *nextAvailable - 1;
                    firstInd = nextAvailable - 1;
                    //p = nodes + activeGametes[firstInd]; // descendant node (previously generated node, nextAvailable - 1)
                    p = rnodes[firstInd]; // descendant node (previously generated node, nextAvailable - 1)
                    // link the nodes
                    //r->left = p;
                    r1->left = p;
                    r1->right = NULL;
                    p->parent= r1;
                    
                    
                    //r->right = NULL;
                    //p->anc1 = r;
                    
                    u1->time = currentTime;
                    u1->scaledByThetaTimeInputTreeUnits = currentTime * theta;
                    
                    if (popI->order==numClones-1)//oldest population
                    {
                         u1->timePUnits =  currentTime * theta;
                         p->length = currentTime * theta ;
                    }
                    else{
                        
                        u1->timePUnits =  currentTime * theta * popI->x;
                        p->length = currentTime * theta * popI->x;
                    }
                        
                    
                    u1->timeInputTreeUnits = u1->timePUnits;
          
                    popI->rMRCA = p;
                    //insertMRCAMap(r,popI);
                    //rMRCAPopulation[p].push_back(popI);
                    
                    /* readjust active nodes */
                    nextAvailable=nextAvailable+1; /* 1 node more is available */
                    popI->idsActiveGametes[0] = newInd;//always will be in the 0  position because there is only one left
                    
                    // for(int k=0; k < popI->idsActiveGametes.size();k++)
                    //    fprintf (stderr, "\n pop of order  %d Active gamete id %d", popI->order, popI->idsActiveGametes[k]);
                    
                    if (programOptions.noisy > 1)
                        std::cout<< "Creating origin node, it creates node "<< newInd << " derived from node "<<  firstInd << std::endl;;
                    if (programOptions.noisy > 1)
                        std::cout<< "\t|\tCurrentTime (input units) =" <<  u1->timePUnits<< std::endl;
                    /* memory for number of nodes */
                    if (nextAvailable >= numNodes)  /* if there aren't enough nodes it go into and it addition more */
                    {
                        /* ReallocNodes(&numNodes); */
                        if (programOptions.noisy == 4)
                            std::cout << "\n\n...Doing reallocation of nodes (Coalescence)\n"<< std::endl;
                        numNodes += INCREMENT_NODES;
                    }
                }
                else  {//origin of oldest pop reached
                    popI->nodeIdAncestorMRCA=nextAvailable-1;//for the last population, nodeIdAncesterMRCA is the MRCA instead of ancester of MRCA
                    r = rnodes[nextAvailable-1];//popI->idsActiveGametes[0]
                    u= (TreeNode *)(r->data);
                    u->indexOldClone = u->indexCurrentClone = popI->index;
                    u->orderCurrentClone = popI->order;
                    popI->rMRCA= r;
                    //insertMRCAMap(r,popI);
                    //rMRCAPopulation[r].push_back(popI);
                }
            }
        }
        if (programOptions.noisy > 3)
        {
            std::cout << "\nActive nodes:" << popI->numActiveGametes << std::endl;
            for (i = 0; i < popI->numActiveGametes; i++)
                std::cout << popI->idsActiveGametes[i]<< std::endl;
            std::cout <<"\t|\tNext node available =" << nextAvailable<< std::endl;
        }
    }
    if (programOptions.noisy > 1)
        std::cout << "\n\nEvolutionary history of clone "<< popI->index<<   " is completed \n" << popI->index << std::endl;
    
}

pll_rnode_t* Chain::BuildTree(Population *CurrentPop,
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
        std::cout << "\n>> Building trees .."<< std::endl;
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
            std::cout << "\n\n>> Attaching outgroup .. "<<  std::endl;
        
        
        if (programOptions.outgroupSelection == 1)  /*** Root 2 times and outgroup ***/
            currentTime = CurrentPop->timeOriginSTD; // origin of the clone; // currentTime + (outgroupBranchLength_Root1Root2 / mutationRate); // set time of the new root (from which the MRCA and outgroup nodes are derived)
        else if (programOptions.outgroupSelection == 2) { /*** Root 2 times and outgroup ***/
            currentTime = CurrentPop->timeOriginSTD + (programOptions.outgroupBranchLength_Root1Root2 / (theta * CurrentPop->x)); // origin of the clone + time given by the user
            
        }
        else
        {
            std::cout << "\n\nError simulationg the outgroup. Check input settings" <<  std::endl;
            Output::PrintUsage();
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
        
        //healthyR->timePUnits = currentTime * theta *CurrentPop->x;
        healthyR->timePUnits = currentTime * theta;
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
        //        healthyRoot->length = 0;
        
        //        if (noisy > 2)
        //            fprintf (stderr, "DONE");
        //
        nextAvailable++;
        //
        //        /* connect the healthy ancestral cell with the tip healthy cell*/
        //        if (noisy > 2)
        //            fprintf (stderr, "\n>> Adding healthy tip ... ");
        
        pll_rnode_t * healthyTip1= rnodes[nextAvailable];
        if (healthyTip1->data==NULL)
            healthyTip1->data =  new TreeNode(0);
        //
        //        if (healthyTip ==NULL)
        //        {
        //            healthyTip1=rnodes[nextAvailable];
        //            if (healthyTip1->data==NULL)
        //               healthyTip1->data =  new TreeNode(0);
        //        }
        //       else
        //           healthyTip1= healthyTip;
        
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
            healthyTip1->label = new char[((string)(programOptions.healthyTipLabel)).length() + 1];
            healthyTip1->label= strcpy(healthyTip1->label, programOptions.healthyTipLabel.c_str());
        }
        
        healthyRoot->right = healthyTip1;
        healthyTip1->parent=healthyRoot;
        
        
        //u1->effectPopSize=healthyR->effectPopSize;
        
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
        std::cout <<  "\n\n>> Relabeling nodes on tree... " << std::endl ;
    if (programOptions.thereisOutgroup == YES)
        intLabel = programOptions.TotalTumorSequences + 1;
    else
        intLabel = programOptions.TotalTumorSequences;
    
    RelabelNodes2( treeRootInit, intLabel );
    
    return treeRootInit;
}
pll_rnode_t * Chain::MakeCoalescenceTree (long int *seed,
                                          const gsl_rng * randomGenerator,
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
                                          std::vector<std::vector<int> > &ObservedData,
                                          char* ObservedCellNames[],
                                          vector<int> &sampleSizes
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
    string nodeLabel;
    
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
        p=new pll_rnode_t();
        p->left =NULL;
        p->right=NULL;
        p->parent=NULL;
        p->label =new char[strlen(msa->label[i]) + 1]{};
        std::copy(msa->label[i], msa->label[i] + strlen(msa->label[i]), p->label);
        
        nodeLabel = p->label;
        
        t = new TreeNode(msa->length);
        t->initNumberTipsVector(numClones);
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
        t->initNumberTipsVector(numClones);
        p->data = t;
        rnodes.push_back(p);
    }
    
    AssignSequencesToPopulations( programOptions, numNodes, programOptions.noisy, programOptions.TotalTumorSequences, numActiveGametes,  nextAvailable,
                                 labelNodes, ObservedCellNames, programOptions.doUseObservedCellNames, sampleSizes);
    Population *currentPop;
    Population *fatherPop;
    i=0;
    currentTime=0.0;
    while (i < numClones) {
        //currentPop = *(populations + i);
        currentPop = populations[i];
        SimulatePopulation( currentPop,  programOptions, seed,randomGenerator,
                           programOptions.numNodes,
                           numClones,
                           nextAvailable ,
                           numActiveGametes,
                           labelNodes,
                           currentTime,
                           eventNum, programOptions.K);
        if (i< numClones-1)   //if it is not the last one
        {
            //choose the father population from which the population i came
            fatherPop= ChooseFatherPopulation( numClones, currentPop, randomGenerator,  programOptions.noisy, programOptions.K);
            currentPop->FatherPop = fatherPop;
            //update list of migrant times
            Population::UpdateListMigrants( numClones, currentPop, fatherPop);
        }
        i = i + 1;
    }
    pll_rnode_t *root=NULL;
    
    root =  BuildTree( currentPop,randomGenerator,
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
Population* Chain::ChooseFatherPopulation( int numClones, Population  *PopChild, const  gsl_rng *randomGenerator, int noisy, long double K) {
    
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
void Chain::AssignSequencesToPopulations(
                                         ProgramOptions &programOptions,
                                         int numNodes, int noisy,  int TotalTumorSequences,
                                         int &numActiveGametes, int &nextAvailable,
                                         int &labelNodes, char* ObservedCellNames[], int doUseObservedCellNames, vector<int> &sampleSizes)

{
    Population *pop;
    pll_rnode_t *p;
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
        //activeGametes[*numActiveGametes] = j;
        p->node_index = j;
        //p->label = j;
        
        labelNodes = j;
        
        if (p->data == NULL)
            p->data = new TreeNode(0);
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
                //pop->idsActiveGametes.push_back(j);
                pop->idsActiveGametes[pop->numActiveGametes]=j;
                //pop->rtips.push_back(p);
                pop->rtips[pop->numActiveGametes]=p;
                
                pop->numActiveGametes=pop->numActiveGametes+1;
                // for(int k=0; k < pop->idsActiveGametes.size();k++)
                //    fprintf (stderr, "\n pop of order  %d Active gamete id %d", pop->order, pop->idsActiveGametes[k]);
                // pop->idsGametes.push_back(j);
                pop->idsGametes[pop->numGametes]=j;
                pop->numGametes=pop->numGametes+1;
                
                
                
                u->indexOldClone = currentPopIndex;
                u->indexCurrentClone = currentPopIndex;
                //u->effectPopSize= pop->effectPopSize;
                u->orderCurrentClone = pop->order;
                
                //                if (programOptions.doUseObservedCellNames)
                //                    strcpy( u->observedCellName, ObservedCellNames[indexFirstObservedCellName + pop->numActiveGametes ]);
                
                
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
void Chain::GenerateEffectPopSizesFromPriors2( int noisy, int doGenerateProportionsVector, const gsl_rng* randomGenerator){
    
    
    if (doGenerateProportionsVector == YES)
    {
        initProportionsVector();
        long double  alpha[numClones];
        if (proportionsVector.size()==0)
            std::fill_n(alpha, numClones, 1.0);
        
        generateProportionsVectorFromDirichlet(alpha, randomGenerator);
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
        // popJ->effectPopSize = popJ->x * totalEffectPopSize;
    }
}
void Chain::initProportionsVector(){
    for (size_t i = 0; i < numClones; i++) {
        proportionsVector.push_back(0);
        oldproportionsVector.push_back(0);
    }
    
}
void Chain::initProportionsVectorFromSampleSizes(vector<int> &sampleSizes)
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
void Chain::generateProportionsVectorFromDirichlet(long double  alpha[], const gsl_rng* randomGenerator){
    
    long double  theta[numClones];
    vector<long double> outputVector;
    for (unsigned int i = 0; i < numClones; ++i){
        theta[i]=0;
    }
    vector<long double> alphaVector(alpha, alpha + numClones);
    Random::randomDirichletFromVector(alphaVector, outputVector,true, randomGenerator, NULL);
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
void Chain::FillChainPopulationsFromPriors( ProgramOptions &programOptions,  MCMCoptions &mcmcOptions, vector<int> &sampleSizes, const gsl_rng* randomGenerator, long double K)
{
    Population* popI;
    int i;
    long double  randomDelta;
    mutationRate = Random::RandomLogUniform(mcmcOptions.MutRatefrom, mcmcOptions.MutRateto, mcmcOptions.useGSLRandomGenerator, randomGenerator, NULL);
    
    initProportionsVector();
    //long double  lambda = 1;
    SetPopulationsBirthRate(mcmcOptions.fixedLambda);
    
    // GenerateEffectPopSizesFromPriors2( programOptions.noisy,    YES);
    if (programOptions.doUseFixedTree ==NO)
    {
        long double  alpha[numClones];
        std::fill_n(alpha, numClones, 1.0);
        generateProportionsVectorFromDirichlet(alpha, randomGenerator);
        if (programOptions.populationSampleSizesKnown == NO)
        {
            InitPopulationSampleSizes(populations, programOptions.numCells, programOptions.numClones,proportionsVector,  randomGenerator);
        }
        //else fill the sample  sizes
        else{
            setChainPopulationSampleSizes( sampleSizes, programOptions);
        }
        
        initPopulationsCoalTimes();
        
        if (programOptions.populationSampleSizesKnown ==YES)
        {
            for( i = 0 ; i < numClones; i++)
            {
                popI=populations[i];
                
                randomDelta = Random::RandomLogUniform(mcmcOptions.Deltafrom, mcmcOptions.Deltato, mcmcOptions.useGSLRandomGenerator, randomGenerator,NULL);
                //popI->delta = chain->proportionsVector[i] * randomDelta;
                popI->delta =  randomDelta;
                popI->growthRate =popI->delta  / popI->x;
                //popI->popSize=popI->x * popI->birthRate;
                popI->deathRate= popI->birthRate - popI->growthRate;
                
            }
        }
        ListClonesAccordingTimeToOrigin(populations);
        GenerateTimesFromPriorsOriginal(programOptions.noisy,  randomGenerator, K);
        ListClonesAccordingTimeToOrigin(populations);
        
        InitListPossibleMigrations();
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

void Chain::GenerateTimesFromPriorsOriginal(int noisy, const gsl_rng* randomGenerator, long double K) {
    vector<long double>    Uvector(numClones);
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
                Uvector[j] =  Random::randomUniformFromGsl2(randomGenerator);
                rand = Random::randomUniformFromGsl2(randomGenerator);
                popJ->timeOriginSTD= (1/ popJ->delta)*log(1-(popJ->delta / log(rand)));
                //popJ->timeMigrationSTDCurrentPop= popJ->timeOriginSTD ;
                // popJ->timeOriginInput=proportionsVector[j]*popJ->timeOriginSTD;
                popJ->timeOriginInput = popJ->x*popJ->timeOriginSTD;
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
                belowTerm = 1.0;
                m = popI->order;
                
                for (l = i + 1; l < numClones; l++)
                {
                    popL=populations[ l];
                    j = popL->order ;
                    //aboveTerm = aboveTerm + (popL->popSize * Population::CalculateH(popI->timeOriginSTD * popI->x / popL->x, popL->timeOriginSTD, popL->delta, K));
                    //belowTerm = belowTerm + popL->popSize;
                    aboveTerm = aboveTerm +  Population::CalculateH(popI->timeOriginSTD * popI->x / popL->x, popL->timeOriginSTD, popL->delta, K);
                    //belowTerm = belowTerm ;
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
void Chain::InitChainPopulations( int noisy,  int TotalNumSequences , MCMCoptions &mcmcOptions ) {
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
        
        auto pop = new Population(ind, ord, timeOriginInput, sampleSize, popSize, birthRate, deathRate, NO, mcmcOptions);
        populations.push_back(pop);
        
        indexMCMCParameters.push_back(std::make_pair(pop->DeltaT, currentIdxMCMCParameters));
        currentIdxMCMCParameters++;
        
        indexMCMCParameters.push_back(std::make_pair(pop->TimeOriginSTD, currentIdxMCMCParameters));
        currentIdxMCMCParameters++;
        
        
        //            Theta=new MCMCParameterWithKernel("Theta_pop"+std::to_string(index),0.0,mcmcOptions.paramMultiplierMoveTheta,0.0);
        //            long double doubleSampleSize= (long double) sampleSize;
        //
        //            X = new MCMCParameter<long double>("Proportion_pop"+std::to_string(index),doubleSampleSize, 0.0);
        
        resizeStoredMCMCparameters( mcmcOptions.numberWarmUpIterations,currentIdxMCMCParameters);
    }
}
void Chain::resizeStoredMCMCparametersAfterWarmUp(MCMCoptions &mcmcOptions){
    
    
    resizeStoredMCMCparameters( mcmcOptions.Niterations+mcmcOptions.numberWarmUpIterations,currentIdxMCMCParameters);
    
}
/**************** RelabelNodes2 **************/

void Chain::RelabelNodes2(pll_rnode_t *p, int &intLabel)
{
    if (p != NULL)
    {
        TreeNode *u;
        RelabelNodes2 (p->left, intLabel);
        RelabelNodes2 (p->right, intLabel);
        /*RelabelNodes (p->outgroup);*/
        if (p->left == NULL && p->right == NULL) /* is tip */
        {
            // p->label = intLabel++;
            //  p->label = (*intLabel);
            // *intLabel=*intLabel+1;
            u=(TreeNode*) p->data;
            u->label = u->index ;
        }
        else                  /* all ancester */
        {
            //p->label = intLabel++;
            u=(TreeNode*) p->data;
            u->label = intLabel;
            intLabel=intLabel+1;
        }
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
void Chain::InitPopulationSampleSizes(vector<Population*> &populations, int TotalSampleSize, int numClones, vector<long double> &proportionsVector, const gsl_rng* randomGenerator)
{
    int i,j;
    Population  *popJ;
    double rand;
    
    vector<long double> cumSum(numClones +1);
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
char * Chain::toNewickString2 ( pll_rnode_t *p,  string healthyTipLabel,    int doUseObservedCellNames)
{
    char buffer[1024];
    
    char *newickString =NULL;
    char *left=NULL;
    char *right=NULL;
    //char *outgroup =NULL;
    TreeNode *u, *anc;
    
    if (p != NULL)
    {
        u= (TreeNode*) (p->data);
        if (p->parent!=NULL)
            anc= (TreeNode*) (p->parent->data);
        
        if (u->isOutgroup == YES)     /* Outgroup */
        {
            strcpy( u->cellName, healthyTipLabel.c_str());
            
            
            if (asprintf(&newickString,  "%s:%.20Lf", healthyTipLabel.c_str(), (anc->timePUnits - u->timePUnits) )<0)
                return NULL;
            
            return newickString;
        }
        else if (p->left == NULL && p->right == NULL)   /* tip of the tree */
        {
            snprintf(buffer, sizeof(buffer), "tip_i%05d_C%d_%d", u->index,u->indexOldClone,u->indexCurrentClone);
            strcpy( u->cellName,buffer);
            
            if  (doUseObservedCellNames == YES)
            {
                if (asprintf(&newickString,   "%s:%.20Lf",  u->observedCellName, (anc->timePUnits - u->timePUnits))<0)
                    return NULL;
                
                return newickString;
            }
            else{
                if (asprintf(&newickString,   "%s:%.20Lf",  u->cellName, (anc->timePUnits - u->timePUnits))<0)
                    return NULL;
                
                
                return newickString;
            }
            //  }
        }
        else
        {
            // fprintf (fpTrees2, "(");
            if ( p->left != NULL  )
            {
                left = toNewickString2 (p->left, healthyTipLabel,  doUseObservedCellNames);
                
                if (left == NULL)
                    return NULL;
            }
            if ( p->right != NULL  )
            {
                right = toNewickString2 (p->right, healthyTipLabel,  doUseObservedCellNames);
                if (right == NULL)
                {
                    free(left);
                    return NULL;
                    
                }
            }
            // outgroup =toNewickString2 (u->outgroup, healthyTipLabel, doUseObservedCellNames);
            if(left!=NULL && right!=NULL && p->parent != NULL)
            {
                snprintf(buffer, sizeof(buffer), "int_i%05d_C%d_%d",  u->index,u->indexOldClone,u->indexCurrentClone);
                strcpy( u->cellName,buffer);
                
                if (asprintf(&newickString, "(%s,%s):%.20Lf", left, right,  (anc->timePUnits - u->timePUnits) )<0)
                    return NULL;
                //snprintf(newickString, size, "(%s,%s):%10.9lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate );
                free(left);
                left = NULL;
                free(right);
                right=NULL;
                //free(outgroup);
                //outgroup=NULL;
            }
            else if (left!=NULL  &&  right!=NULL  && p->parent == NULL)
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
                //free(outgroup);
                //outgroup=NULL;
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
        //if (isnan(temp) || isinf(temp))
        //   fprintf (stderr, "\n isNan temp LogDensityTime, Delta: %Lf, T:%.20Lf \n",popI->delta ,popI->timeOriginSTD);
        result = result + temp;
        // fprintf (stderr, "\n Product log Density Time   = %lf after pop order %d,  popI->LogDensityTime( popI->timeOriginSTD) %lf,  popI->timeOriginSTD: %lf, popI->delta: %lf\n", product, popI->order, popI->LogDensityTime( popI->timeOriginSTD), popI->timeOriginSTD,popI->delta );
    }
    return result;
}

long double Chain::SumLogProbFatherPopulations(long double K) {
    
    Population *fatherPop;
    Population *popJ,*popK;
    long double  numerator=0;
    long double  denominator=1;
    long double  result=0;
    long double temp;
    for ( unsigned int j = 0; j < numClones-1 ; j++)
    {
        popJ = populations[j];
        for ( unsigned int k = popJ->order+1; k < numClones-1 ; k++){
            popK = populations[k];
            
            denominator = denominator + log( popK->x);
            temp=popJ->timeOriginSTD * popJ->x / popK->x;
            temp=Population::LogCalculateH(temp, popK->timeOriginSTD, popK->delta, K);
            denominator = denominator  + temp;
        }
        
        if(popJ->order < numClones-1 && popJ->FatherPop == NULL)
                   fprintf (stderr, "\n The pop with order %d has father pop = null\n",popJ->order );
        if (popJ->FatherPop !=NULL && popJ->FatherPop->order > popJ->order )
        {
            fatherPop = popJ->FatherPop;
            // result = result + log( fatherPop->popSize);//no popSize size available
            numerator = numerator + log( fatherPop->x);
            temp=popJ->timeOriginSTD * popJ->x / fatherPop->x;
            numerator=Population::LogCalculateH(temp, fatherPop->timeOriginSTD, fatherPop->delta, K);
            numerator = numerator  + temp;
            // fprintf (stderr, "\n Product calculate H    = %lf after pop order %d \n", product, popJ->order);
        }
        else{
            fprintf (stderr, "\n The pop with order %d has father pop with incorrect order!\n",popJ->order );
            
           }
       
    }
    result = numerator /denominator;
    return(result);
}

/************************ LogConditionalLikelihoodTree ***********************/
/*  LogConditionalLikelihoodTree */
long double  Chain::LogConditionalLikelihoodTree( ProgramOptions &programOptions, MCMCoptions &mcmcOptions  )
{
    
    long double result;
    result =SumLogDensitiesTimeOriginSTDPopulations();
    result= result + SumLogProbFatherPopulations(programOptions.K);
    // fprintf (stderr, "\n Product after  = %lf \n", result);
    result = result + LogDensityCoalescentTimesForPopulation2(programOptions.K);
    
    return result;
}
/************************ LogDensityCoalescentTimesForPopulation ***********************/
/*  LogDensityCoalescentTimesForPopulation */
long double  Chain::LogDensityCoalescentTimesForPopulation(long double K)
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
                
                temp =  Population::LogLambda(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta, K);
                if (isnan(temp) || isinf(temp) )
                    fprintf (stderr, "\n isNan product\n");
                result= result + temp;
                termOnlyAfterFirstCoalEvent =(currentCoalescentEvent == 0)?0:Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent-1], popI->timeOriginSTD, popI->delta, K);
                temp =  (numberAliveCells/ 2.0)* (numberAliveCells-1)*(Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta, K)-termOnlyAfterFirstCoalEvent);
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
                temp= popI->LogProbNoCoalescentEventBetweenTimes(lastEventTimeBeforeMigration,popI->immigrantsPopOrderedByModelTime[currentMigrationEvent].first, numberAliveCells, K );
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
            temp =  Population::LogLambda(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta, K);
            
            if (isnan(temp) || isinf(temp) )
                fprintf (stderr, "\n isNan temp\n");
            
            temp =  Population::LogLambda(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta, K);
            result= result + temp;
            
            termOnlyAfterFirstCoalEvent =(currentCoalescentEvent == 0)?0:Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent-1], popI->timeOriginSTD, popI->delta, K);
            temp=( numberAliveCells/ 2.0)* (numberAliveCells-1.0)*(Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta, K)-termOnlyAfterFirstCoalEvent);
            if (isnan(temp) || isinf(temp) )
                fprintf (stderr, "\n isNan temp\n");
            temp=( numberAliveCells/ 2.0)* (numberAliveCells-1.0)*(Population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta, K)-termOnlyAfterFirstCoalEvent);
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
long double  Chain::LogDensityCoalescentTimesForPopulation2(long double K){
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
            else// if(std::binary_search(immigrantsTimes.begin(),immigrantsTimes.end(), timeCurrentEvent) )//is a migration event
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
long double Chain::LogDensityCoalTimes(std::vector<long double> allEventsSorted, std::vector<long double> coalEventTimes, std::vector<long double> immigrantsTimes, long double timeOriginSTD, long double delta, int sampleSize, long double K){
    long double result =0.0;
    long double timeCurrentEvent;
    long int numberAliveCells;
    long double temp;
    long int currentCoalescentEventInThisEpoch=0;
    long int currentMigrationEvent=0;
    long double termOnlyAfterFirstCoalEvent;
    long double lastEventTimeBeforeMigration=0;
    
    numberAliveCells= sampleSize;
    currentCoalescentEventInThisEpoch=0;
    currentMigrationEvent=0;
    lastEventTimeBeforeMigration=0;
    
    for ( unsigned j = 0; j < allEventsSorted.size(); j++)
    {
        timeCurrentEvent= allEventsSorted.at(j);
        if (timeCurrentEvent== timeOriginSTD)
            break;
        if(std::binary_search(coalEventTimes.begin(), coalEventTimes.end(), timeCurrentEvent) )//is a coalescent event
        {
            if (numberAliveCells > 1)
            {
                temp=log(numberAliveCells * (numberAliveCells-1.0)/2.0);
                result= result + temp;
                
                temp =  Population::LogLambda(timeCurrentEvent,timeOriginSTD, delta, K);
                result= result + temp;
                termOnlyAfterFirstCoalEvent =(currentCoalescentEventInThisEpoch == 0)?Population::FmodelTstandard(lastEventTimeBeforeMigration, timeOriginSTD, delta, K):Population::FmodelTstandard(coalEventTimes[currentCoalescentEventInThisEpoch-1], timeOriginSTD, delta, K);//if no coalescent
                temp =  (numberAliveCells / 2.0)* (numberAliveCells - 1.0)*(Population::FmodelTstandard(timeCurrentEvent,timeOriginSTD, delta, K)-termOnlyAfterFirstCoalEvent);
                result= result - temp;
                lastEventTimeBeforeMigration = timeCurrentEvent;
                currentCoalescentEventInThisEpoch++;
                numberAliveCells--;
            }
            
        }
        else //if(std::binary_search(immigrantsTimes.begin(),immigrantsTimes.end(), timeCurrentEvent) )//is a migration event
        {
            if (numberAliveCells > 1)
            {
                temp= Population::LogProbNoCoalescentEventBetweenTimes(lastEventTimeBeforeMigration, timeCurrentEvent,numberAliveCells, timeOriginSTD, delta, K);
                result= result+ temp;
            }
            lastEventTimeBeforeMigration=timeCurrentEvent;
            currentMigrationEvent++;
            numberAliveCells++;
            currentCoalescentEventInThisEpoch=0;//new epoch is starting
        }
    }
    return result;
    
}
long double  Chain::LogConditionalLikelihoodSequences(pll_msa_t * msa, char* NewickString, ProgramOptions &programOptions, long double seqError,long double dropoutError)
{
    
    //pll_state_t pll_map_gt10_2[256];
    if (NewickString == NULL)
    {
        fprintf (stderr, "\nERROR: The newick representation of the tree cannot be empty\n\n");
        Output::PrintUsage();
        return 0;
    }
    
    unsigned int i;
    unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
    pll_partition_t * partition;
    
    pll_utree_t *unrootedTree = pll_utree_parse_newick_string_unroot(NewickString);
    
    //pll_rtree_t *rootedTree = pll_rtree_parse_newick_string(NewickString);
    //pll_utree_t * unrootedTree = pll_rtree_unroot(rootedTree);
    
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
        //        Output::PrintUsage();
        //        return 0;
    }
    
//
  
    
    pllmod_subst_model_t * model = pllmod_util_model_info_genotype(JC_MODEL);
    //pllmod_subst_model_t * model = pllmod_util_model_info_genotype(GT_MODEL);
    
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
   // unsigned int attributes = PLL_ATTRIB_ARCH_CPU | PLL_ATTRIB_ARCH_SSE | PLL_ATTRIB_ARCH_AVX;
    unsigned int attributes = PLL_ATTRIB_ARCH_AVX;
    partition = pll_partition_create(tip_nodes_count,
                                     inner_nodes_count,//unsigned  int  clv_buffers
                                     model->states,
                                     (unsigned int)(msa->length),//unsigned  int  sites
                                     1,//unsigned  int  rate_matrices
                                     branch_count, // unsigned int prob_matrices
                                     RATE_CATS, // unsigned  int  rate_cats
                                     inner_nodes_count, // unsigned  int   scale_buffers
                                     attributes );
   
    //pll_set_asc_bias_type(partition, 0);
    //pll_set_asc_bias_type(partition, PLL_ATTRIB_AB_LEWIS);// also could be PLL_ATTRIB_AB_FELSENSTEIN or PLL_ATTRIB_AB_STAMATAKIS
   // pll_update_invariant_sites_proportion(partition, 0, 0.5);//"Invariant sites are not compatible with asc bias correction"
    
     //set_partition_tips_costum( partition, msa, programOptions,  seqError, dropoutError);
    
    
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
    
    set_partition_tips_costum( partition, msa, programOptions,  seqError, dropoutError);
    
    double user_freqs[10] = { 0.254504, 0.157238, 0.073667, 0.287452, 0.027336,
              0.058670, 0.041628, 0.024572, 0.064522, 0.010412 };
   
    //double  *  empirical_frequencies = pllmod_msa_empirical_frequencies(partition);
    
    if (model->freqs)
        pll_set_frequencies(partition, 0,  model->freqs);
    else{
        //pll_set_frequencies(partition, 0,  empirical_frequencies);
        pll_set_frequencies(partition, 0,  user_freqs);
    }
    unsigned int * weight = pll_compress_site_patterns(msa->sequence,
                                                         pll_map_gt10,
                                                         tip_nodes_count,
                                                         &(msa->length));
    
    double  * empirical_subst_rates = pllmod_msa_empirical_subst_rates( partition);
 
    
    //  printf("Number of unique site patterns: %d\n\n", msa->length);
    
    /* initialize the array of base frequencies  AA CC GG TT AC/CA AG/GA AT/TA CG/GC CT/TC GT/TG  */
//        double user_freqs[10] = { 0.254504, 0.157238, 0.073667, 0.287452, 0.027336,
//            0.058670, 0.041628, 0.024572, 0.064522, 0.010412 };
    
    
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
//    if (model->freqs)
//        pll_set_frequencies(partition, 0,  model->freqs);
//    else{
//        //pll_set_frequencies(partition, 0,  empirical_frequencies);
//        pll_set_frequencies(partition, 0,  user_freqs);
//    }
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
            
            //pll_set_tip_clv(partition, tip_clv_index, tipCLV, PLL_TRUE);
            
            
        }
    }
    
}


int Chain::set_tipclv_custom_error_model(pll_partition_t * partition,
                                         unsigned int tip_index,
                                         const pll_state_t * state,
                                         const char * sequence,
                                         long double  seqErrorRate,
                                         long double  dropoutRate
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
        if ((c = state[(int)sequence[index]]) == 0)
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
                        tipclv[j] = 1. - seqErrorRate + 0.5 * seqErrorRate * dropoutRate;
                    else
                        tipclv[j] =  (1. - dropoutRate ) * (1. - seqErrorRate) + one_12 * seqErrorRate * dropoutRate;
                }
                else if (mut_dist[state_id][j] == 1)
                {
                    /* 1 letter away */
                    if (HOMO(j))
                    {
                        tipclv[j] = one_12 * seqErrorRate * dropoutRate +
                        one_3  * (1. - dropoutRate) * seqErrorRate;
                    }
                    else
                    {
                        if (HOMO(state_id))
                        {
                            tipclv[j] = 0.5 * dropoutRate + one_6 * seqErrorRate -
                            three_8 * seqErrorRate * dropoutRate;
                        }
                        else
                        {
                            tipclv[j]= one_6 * seqErrorRate -
                            one_8 * seqErrorRate * dropoutRate;
                        }
                    }
                }
                else
                {
                    /* 2 letters away */
                    if (HOMO(state_id))
                        tipclv[j] = one_12 * seqErrorRate * dropoutRate;
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

void Chain::initMutationRate( MCMCoptions &mcmcOptions, ProgramOptions &programOptions, const gsl_rng * randomGenerator) {
    
    long double logPriorDensity;
    long double from =0.0;
    if (programOptions.doUsefixedMutationRate){
        theta = programOptions.mutationRate;
        logPriorDensity=0.0;
    }
    else
    {
        if(mcmcOptions.priorsType==0)
        {
            theta = Random::RandomLogUniform(mcmcOptions.MutRatefrom, mcmcOptions.MutRateto, mcmcOptions.useGSLRandomGenerator, randomGenerator, NULL);
            
            logPriorDensity=Distributions::LogUniformDensity(theta,mcmcOptions.MutRatefrom, mcmcOptions.MutRateto);
            
        }
        else if(mcmcOptions.priorsType==1)
        {
            theta = Random::RandomExponentialStartingFrom( mcmcOptions.lambdaExponentialPriorMutationRate, from, mcmcOptions.useGSLRandomGenerator, randomGenerator, NULL);
            logPriorDensity=Distributions::LogExponentialDensity(theta,from, mcmcOptions.lambdaExponentialPriorMutationRate);
        }
        
        else{
            
            theta = Random::RandomPowerLawDistribution(mcmcOptions.parameterPowerLawDistributionMutationRate, from, mcmcOptions.useGSLRandomGenerator, randomGenerator, NULL);
            logPriorDensity=Distributions::LogPowerLawDistibutionDensity(theta, mcmcOptions.parameterPowerLawDistributionMutationRate,from);
            
        }
    }
    thetaPar->setParameterValue(theta);
    
    //printf ( "initial  scaled mutation rate: %Lf with logPrior: %Lf  for chain %d \n",theta, thetaPar->getParameter()->getLogPriorCurrentValue(), chainNumber);
    
}

void Chain::initLogLikelihoods(MCMCoptions &mcmcOptions, pll_msa_t *msa, ProgramOptions &programOptions, Partition *partition) {
    
    //if (!mcmcOptions.noData)
    // {
    long double result=0;
    currentlogSumDensitiesTimeOriginSTDPopulations =SumLogDensitiesTimeOriginSTDPopulations();
    result=currentlogSumDensitiesTimeOriginSTDPopulations;
    currentlogProbFatherPopulations=SumLogProbFatherPopulations(programOptions.K);
    result= result + currentlogProbFatherPopulations;
    currentlogDensityCoalescentTimesForPopulation= LogDensityCoalescentTimesForPopulation2(programOptions.K);
    result = result + currentlogDensityCoalescentTimesForPopulation;
    currentlogConditionalLikelihoodTree= result;
    if (mcmcOptions.verbose>=0)
        std::cout << "\nInitial log likelihood of the tree of chain " << chainNumber << " is " << currentlogConditionalLikelihoodTree << std::endl;
    // if (mcmcOptions.useSequencesLikelihood ==1)
    //  {
    //char *  rootedNewick2 = pll_rtree_export_newick( initialRootedTree->root, NULL);
    //initialRootedTree = pll_rtree_parse_newick_string(rootedNewick3);
    
    //currentlogConditionalLikelihoodSequences = pll_utils::LogConditionalLikelihoodSequences( msa,  rootedNewick2, programOptions, seqErrorRate, dropoutRate);
    gtErrorModel= new GenotypeErrorModel("GT20", seqErrorRate, dropoutRate, 16);
    
     // treeLik=  new TreeLikelihood(*partition, RootedTree(initialRootedTree),  msa, gtErrorModel);
//     treeLik= new TreeLikelihood(
//                                                initialRootedTree->tip_count,
//                                                initialRootedTree->inner_count,
//                                                16,
//                                                (unsigned int)(msa->length),//int numberSites,
//                                                1,
//                                                initialRootedTree->edge_count,//int probMatrices,
//                                                1,//,int numberRateCats,
//                                                initialRootedTree->inner_count,
//                                                0,
//                                                false,
//                                                false,
//                                                false,
//                                                false,
//                                                false,
//                                                false,initialRootedTree ,  msa,  gtErrorModel);
//
    
  //  currentlogConditionalLikelihoodSequences = treeLik->computeRootLogLikelihood();
  // currentlogConditionalLikelihoodSequences =   pll_utils::LogConditionalLikelihoodSequencesRootedTree( msa,  initialRootedTree, programOptions, seqErrorRate, dropoutRate );
    
   // if (mcmcOptions.verbose>=0)
    std::cout << "\nInitial log likelihood of the sequences of chain " << chainNumber << " is " << currentlogConditionalLikelihoodSequences << std::endl;
    
   // free(rootedNewick2);
    //rootedNewick2=NULL;
    // }
  
    if (mcmcOptions.verbose>=2)
        std::cout<<  "\nInitial log likelihood of the sequences of chain "<< chainNumber << " is = " << currentlogConditionalLikelihoodSequences << std::endl;
    // }
   
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
void Chain::InitPopulationSampleSizes(vector<int> &sampleSizes)
{
    for (unsigned int i = 0; i < numClones; ++i)
    {
        auto popI =  populations[i];
        popI->sampleSize = sampleSizes.at(i);
        
    }
}
void Chain::InitPopulationGametes()
{
    for (unsigned int i = 0; i < numClones; ++i)
    {
        auto popI =  populations[i];
        
        popI->InitIdsGametes(numClones);
        popI->InitIdsActiveGametes();
        
    }
}
void Chain::InitPopulationRTips()
{
    for (unsigned int i = 0; i < numClones; ++i)
    {
        auto popI =  populations[i];
        
        popI->InitRTips();
        
    }
}

void Chain::initChainTree( std::string &healthyTipLabel, MCMCoptions &mcmcOptions, pll_msa_t *msa, ProgramOptions &programOptions ) {
    
    // assert(initialRootedTree !=NULL);
    //  if (initialRootedTree!=NULL)
    //      initialUnrootedTree = pll_rtree_unroot(initialRootedTree);
    
    
    rootRootedTree = initialRootedTree->root;
    //root =initialUnrootedTree->nodes[initialUnrootedTree->tip_count + initialUnrootedTree->inner_count - 1];
    
    //mcmcOptions.numberTumorCells =initialUnrootedTree->tip_count-1;
    mcmcOptions.numberTumorCells =initialRootedTree->tip_count-1;
    
    pll_rtree_reset_template_indices(initialRootedTree->root, initialRootedTree->tip_count);
    
    // pll_utree_traverse_apply(chain->initialUnrootedTree->vroot, NULL, Chain::computeNumberTipsSubTree, chain->root->data);
    
    //std::transform(alpha, alpha + chain->numClones , alpha,std::bind2nd(std::divides<double>(),totalSampleSize));
    // std::transform(alpha, alpha + chain->numClones , alpha,[totalSampleSize](long double  a) {return a /totalSampleSize; } );
}
void Chain::initTreeBranchLengths(std::string &healthyTipLabel)
{
    
    initBranches(healthyTipLabel, edgeLengths, edges);
}
void Chain::initTreeTips(std::string &healthyTipLabel, MCMCoptions &mcmcOptions,ProgramOptions &programOptions ){
    
    //initPopulationsTipsFromTree(initialUnrootedTree, NO, programOptions.healthyTipLabel);
    initPopulationsTipsFromRootedTree(initialRootedTree, NO, programOptions.healthyTipLabel);
    //initNodeDataFromTree();
    
    
    
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

Chain *Chain::initializeChain( int chainNumber,   ProgramOptions &programOptions,  MCMCoptions &mcmcOptions, vector<int> &sampleSizes, const gsl_rng * randomGenerator, std::vector<std::vector<int> > &ObservedData,char* ObservedCellNames[], pll_msa_t *msa, pll_rtree_t * inputRootedTree, StructuredCoalescentTree *structCoalTree,  string& healthyTipLabel, FilePaths &filePaths, Partition *partition)
{
    
    
    vector<long double> alpha;
    
    Chain *chain=new Chain(chainNumber, programOptions.numClones, 1, programOptions.mutationRate, programOptions.seqErrorRate, programOptions.dropoutRate, mcmcOptions);
    
    chain->PrepareFiles(filePaths, programOptions, chain->files, chainNumber);
    
    if (programOptions.numberClonesKnown)
    {
        chain->numNodes = programOptions.numNodes;
    }
    chain->InitChainPopulations( programOptions.noisy, programOptions.TotalNumSequences , mcmcOptions);
    
    chain->initMutationRate( mcmcOptions, programOptions, randomGenerator);
    if (mcmcOptions.verbose>=2)
        std::cout <<  "\n initial theta" << chain->theta << std::endl;
    chain->samplePopulationGrowthRateFromPriors(mcmcOptions, randomGenerator );
    
    chain->SetPopulationsBirthRate(mcmcOptions.fixedLambda);
    
    if (programOptions.doUseFixedTree == NO)
    {
        //        Population* oldestPop=chain->getPopulationbyIndex(chain->numClones -1);
        //        chain->InitPopulationSampleSizes(sampleSizes);
        //        chain->InitPopulationGametes();
        //        chain->InitPopulationRTips();
        //        chain->initProportionsVectorFromSampleSizes(sampleSizes);
        //        chain->initPopulationsThetaDelta();
        //        oldestPop->timeOriginSTD  =   Random::RandomDensityModelTimeOrigin (oldestPop->delta, mcmcOptions.useGSLRandomGenerator, 0, randomGenerator );
        //
        //        chain->InitListPossibleMigrations();
        //        chain->InitPopulationsCoalescentEvents();
        //        chain->resetMRCAMap();
        //        chain->initEdgesRootedTree(programOptions.TotalNumSequences);
        
        pll_rtree_t  * trueTree=structCoalTree->getTree();
        
        //char *  rootedNewick3 = pll_rtree_export_newick( tree->root,  Utils::cb_serialize);
        char *  rootedNewick3 = pll_rtree_export_newick( trueTree->root, NULL);
        
        chain->initialRootedTree = pll_rtree_parse_newick_string(rootedNewick3);
        if (!(chain->initialRootedTree))
            fprintf(stderr, "Error in pll_rtree_parse_newick %s ", rootedNewick3);
        
        //std::memcpy(chain->initialRootedTree,tree, sizePllRTree);
        chain->rootRootedTree =chain->initialRootedTree->root;
        //chain->initMRCAMap();
        
        chain->initChainTree( healthyTipLabel, mcmcOptions, msa, programOptions);
        
        chain->saveTrueTreeInfo(trueTree, rootedNewick3,  programOptions);
        chain->closeTrueTreeFiles( programOptions);
        //TreeNode *u= (TreeNode *)(oldestPop->rMRCA->data);
        //oldestPop->setLowerBoundTimeOriginInput(u->timeInputTreeUnits);
        
        //chain->ListClonesAccordingTimeToOrigin(chain->populations);
        // chain->initTimeOriginSTDYoungerPopulations(mcmcOptions);
        
        //chain->filterSortPopulationsCoalescentEvents();
        
        // chain->initLogLikelihoods(mcmcOptions, msa, programOptions);
    }
    else
    {
        
        char *  rootedNewick2 = pll_rtree_export_newick( inputRootedTree->root,  NULL);
        
        chain->initialRootedTree = pll_rtree_parse_newick_string(rootedNewick2);
        if (!(chain->initialRootedTree))
            fprintf(stderr, "Error in pll_rtree_parse_newick %s ", rootedNewick2);
        
        int success = chain->setInitialTreeRootedTree(inputRootedTree);
        
        if (success!=0)
        {
            fprintf(stderr, "\n There was an error reading the tree file! \n" );
            exit(-1);
        }
        chain->saveTrueTreeInfo(inputRootedTree,rootedNewick2, programOptions);
        chain->closeTrueTreeFiles( programOptions);
    }//else
    chain->rnodes = Utils::filterHealthyTip(chain->initialRootedTree,chain->initialRootedTree->nodes,chain->initialRootedTree->inner_count + chain->initialRootedTree->tip_count, healthyTipLabel );
    
    //write to a file initial tree(simulated or input)

    // chain->rnodes = Utils::vectorFromDoublePointer(chain->initialRootedTree->nodes,chain->initialRootedTree->inner_count + chain->initialRootedTree->tip_count );
    
    // chain->rnodes.erase(std::remove_if(chain->rnodes.begin(), chain->rnodes.end(),
    //                           [&](std::shared_ptr<pll_rnode_t>&  node) { return std::string(node->label).compare(healthyTipLabel)==0 ; }), chain->rnodes.end());
    Population* oldestPop=chain->getPopulationbyIndex(chain->numClones -1);
    
    mcmcOptions.numberTumorCells =chain->initialRootedTree->tip_count-1;
    
    chain->initTreeTips(healthyTipLabel, mcmcOptions,programOptions );
    chain->rescaleNodeDataFromRootedTree(chain->theta );
    chain->initNumberTipsSubTree( chain->rootRootedTree);
    chain->initTreeBranchLengths(healthyTipLabel);
    
    chain->rMRCAPopulation = chain->initTimeOfOriginsOnRootedTree( chain->edgeLengths, programOptions.numClones -1,  healthyTipLabel, mcmcOptions, randomGenerator);
    alpha=chain->initVectorSampleSizes(healthyTipLabel, mcmcOptions, programOptions);
    Random::randomDirichletFromVector (alpha, chain->proportionsVector,true, randomGenerator, NULL);
    chain->initPopulationsThetaDelta();
    
    TreeNode *u= (TreeNode *)(oldestPop->rMRCA->data);
    oldestPop->setLowerBoundTimeOriginInput(u->timeInputTreeUnits);
    
    oldestPop->timeOriginSTD =   Random::RandomDensityModelTimeOrigin (oldestPop->delta, mcmcOptions.useGSLRandomGenerator, oldestPop->lowerBoundTimeOriginInput/ (chain->theta * oldestPop->x), randomGenerator, NULL);
    oldestPop->TimeOriginSTD->setParameterValue(oldestPop->timeOriginSTD);
    
    //oldestPop->setTimeOriginSTD(Random::RandomDensityModelTimeOrigin (oldestPop->delta, mcmcOptions.useGSLRandomGenerator, oldestPop->lowerBoundTimeOriginInput/ (chain->theta * oldestPop->x), randomGenerator))
    
    chain->initOriginTimeOldestPopulation(healthyTipLabel, mcmcOptions);
    
    chain->ListClonesAccordingTimeToOrigin(chain->populations);
    chain->initTimeOriginSTDYoungerPopulations(mcmcOptions);
    chain->InitListPossibleMigrations();//after setting the timeSTD
    chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, healthyTipLabel);
    chain->filterSortPopulationsCoalescentEvents();
    
    // totalTreeLength = SumBranches2(rootRootedTree, mutationRate);
    
    partition = new Partition(chain->initialRootedTree->tip_count,// numberTips
     chain->initialRootedTree->inner_count,//unsigned  int  clvBuffers
     16,// model->states,//numberStates
      msa->length,//unsigned  int  sites
     1,//unsigned  int numberRateMatrices
     chain->initialRootedTree->edge_count, // unsigned int probMatrices
     RATE_CATS,//RATE_CATS, // unsigned  int  numberRateCats
      chain->initialRootedTree->inner_count, // unsigned  int numberScaleBuffers
     0, //int statesPadded
    false, false, false, false, false, false);
    
    chain->initLogLikelihoods(mcmcOptions, msa, programOptions, partition);
    
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
void  Chain::initListMoves(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    moves.clear();
    Population *popI;
    
    vector<MCMCParameterWithKernel *> vec;
    vector<IMCMCParameter<long double> *> affectedParameters;
    MCMCVectorParameterWithKernel *vectorParamKernel;
    
    vector<long double> dummy;
    if (numClones > 1){
        
        vector<long double> sampleSizes;
        for (unsigned int i = 0; i < populations.size(); ++i)
        {
            sampleSizes[i]=populations[i]->sampleSize;
        }
        
        long double logPriorDensity = Distributions::LogDirichletDensity(proportionsVector,  sampleSizes, dummy);
        vectorParamKernel = new MCMCVectorParameterWithKernel("ProportionsVector", proportionsVector, sampleSizes,logPriorDensity);
        vec ={thetaPar,mutationRatePar };
        
        for (unsigned int i = 0; i < populations.size(); ++i)
        {
            IMCMCParameter<long double> * deltaTPop = populations[i]->DeltaT->getParameter();
            IMCMCParameter<long double> * thetaPop = populations[i]->Theta->getParameter();
            affectedParameters.push_back(deltaTPop);
            affectedParameters.push_back(thetaPop);
        }
        
        //vector<long double> &nitialValue, vector<long double> &kernelParameter
        NewProportionsVectorMove *newProportionsVector= new NewProportionsVectorMove(chainNumber, "Proportions Vector Move", vec, affectedParameters, *vectorParamKernel, currentlogConditionalLikelihoodTree, currentlogConditionalLikelihoodSequences );
        // NewProportionsVectorMove(Chain *chain, std::string nameMove, vector<MCMCParameterWithKernel *> &mcmcParKernel, vector<MCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam);
        moves.push_back(newProportionsVector );
    }
    NewGlobalScaledMutationRateMove *newGlobalScaledMutationRateMove;
    NewGlobalScaledGrowthRateForPopulationMove *newGlobalScaledGrowthRateMove;
    NewTimeOriginOnTreeforPopulationMove *newTimeOriginOnTreeforPopulationMove ;
    
    //Pair of theta and r move
    // NewGlobalScaledMutationRateMove(Chain *chain, std::string nameMove, Population *pop, vector<MCMCParameterWithKernel *> &mcmcParKernel,  vector<MCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam)
    vec.clear();
    vec.push_back(thetaPar);
    affectedParameters.clear();
    for (unsigned int i = 0; i < populations.size(); ++i)
    {
        IMCMCParameter<long double> * thetaPop = populations[i]->Theta->getParameter();
        affectedParameters.push_back(thetaPop);
        
        if (!mcmcOptions.splitThetaDeltaTmoves){
            vec.push_back(populations[i]->DeltaT);
        }
    }
    
    string name;
    if (!mcmcOptions.splitThetaDeltaTmoves)
        name = "Pair Global scaled mutation rate(theta) and scaled growth rate(deltaT) Move";
    else
        name = "Global scaled mutation rate(theta)";
    
    newGlobalScaledMutationRateMove= new NewGlobalScaledMutationRateMove(chainNumber, name,  NULL, vec, affectedParameters, *vectorParamKernel, currentlogConditionalLikelihoodTree, currentlogConditionalLikelihoodSequences);
    
    //  NewGlobalScaledMutationRateMove(Chain *chain, std::string nameMove, Population *pop, vector<MCMCParameterWithKernel *> &mcmcParKernel,  vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam);
    moves.push_back(newGlobalScaledMutationRateMove);
    for( int i = 0 ; i < numClones; i++)
    {
        popI=populations[i];
        
        vec.clear();
        //IMCMCParameter<long double> * deltaTPop = populations[i]->DeltaT->getParameter();
        
        if (mcmcOptions.splitThetaDeltaTmoves){
            
            vec.push_back(populations[i]->DeltaT);
            affectedParameters.clear();
            affectedParameters.push_back(popI->TimeOriginSTD->getParameter());
            
            newGlobalScaledGrowthRateMove =  new NewGlobalScaledGrowthRateForPopulationMove(chainNumber, "Global scaled growth rate(deltaT)", popI, vec, affectedParameters, *vectorParamKernel, currentlogConditionalLikelihoodTree, currentlogConditionalLikelihoodSequences );
            
            moves.push_back(newGlobalScaledGrowthRateMove);
        }
        
        
        
        vec.clear();
        vec.push_back(popI->TimeOriginSTD);
        
        newTimeOriginOnTreeforPopulationMove= new NewTimeOriginOnTreeforPopulationMove(chainNumber, "new Time of Origin Move on Tree for population", popI, vec, affectedParameters, *vectorParamKernel, currentlogConditionalLikelihoodTree, currentlogConditionalLikelihoodSequences);
        moves.push_back(newTimeOriginOnTreeforPopulationMove);
    }
}

void Chain::stepAllMoves(   MCMCoptions &mcmcOptions, const gsl_rng *rngGsl,    ProgramOptions &programOptions){
    
    std::random_device rng;
    std::mt19937 urng(rng());
    std::shuffle ( moves.begin(), moves.end(), urng );
    vector<MCMCmove*>::iterator it;
    MCMCmove* currentMove;
    for ( it=moves.begin(); it!=moves.end(); ++it)
    {
        currentMove=  *it;
        //fprintf (stderr, "\n>> started  move %s ", currentMove->name().c_str());
        currentMove->move(programOptions, mcmcOptions, rngGsl,NULL, this);
        // fprintf (stderr, "\n>> finished  move %s \n", currentMove->name().c_str());
    }
    //moves.clear();
    //   }
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
    int sampleSize= totalSampleSize();
    
    if ( utree->tip_count - 1 != sampleSize )
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
                if (node->label != 0 && node->label)
                    std::copy(node->label, node->label + strlen(node->label), u1->observedCellName);
                u1->scaledByThetaTimeInputTreeUnits= time ;
                u1->timePUnits= time ;
                u1->timeInputTreeUnits = time_input_tree_units;
                u1->initNumberTipsVector(numClones);
                // printf("updated data for node %d and time input %Lf and scaled time %Lf by theta %Lf \n", node->node_index, u1->timeInputTreeUnits, u1->scaledByThetaTimeInputTreeUnits, scale );
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
    Population *pop;
    TreeNode *u;
    for (unsigned int i = 0; i < initialRootedTree->tip_count+ initialRootedTree->inner_count; ++i)
    {
        auto node =initialRootedTree->nodes[i];
        u=((TreeNode *) node->data);
        
        assert(newTheta >0);
        u->scaledByThetaTimeInputTreeUnits =  u->timeInputTreeUnits / newTheta;
        
        pop= getPopulationbyIndex(u->indexCurrentClone);
        if (pop!=NULL)
            u->timePUnits =  u->timeInputTreeUnits / (newTheta * (pop->x));
        else
            printf("population is null for node  %d  and label %s, observed: %s \n", node->node_index, u->cellName, u->observedCellName );
    }
}
void Chain::saveTreeNodeCurrentTimePUnits(){
    
    TreeNode *u;
    for (unsigned int i = 0; i < initialRootedTree->tip_count+ initialRootedTree->inner_count; ++i)
    {
        auto node =initialRootedTree->nodes[i];
        u = (TreeNode *) node->data;
        u->oldScaledByThetaTimeInputTreeUnits = u->scaledByThetaTimeInputTreeUnits;
        u->oldtimePUnits = u->timePUnits;
    }
}

void Chain::rollbackTreeNodeCurrentTimePUnits(){
    
    TreeNode *u;
    for (unsigned int i = 0; i < initialRootedTree->tip_count+ initialRootedTree->inner_count; ++i)
    {
        auto node =initialRootedTree->nodes[i];
        u = (TreeNode *) node->data;
        assert(!(u->oldtimePUnits != u->oldtimePUnits));
        u->timePUnits = u->oldtimePUnits;
        u->scaledByThetaTimeInputTreeUnits=u->oldScaledByThetaTimeInputTreeUnits ;
    }
}
std::map<pll_unode_t*, Population*>  Chain::chooseTimeOfOriginsOnTree( const gsl_rng * randomGenerator)
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
std::map<pll_rnode_t*, vector<Population*> >  Chain::initTimeOfOriginsOnRootedTree( vector<pair<double, pll_tree_edge_t *> > &edgeLengths, int numberPoints, string &healthyCellLabel, MCMCoptions &mcmcOptions, const gsl_rng * randomGenerator)
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
            
            //            if (assignedPop ==0)
            //                nextEvent =15;
            //            else if (assignedPop ==1)
            //                nextEvent =15;
            
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
            popI->TimeOriginInput->setParameterValue(proposedTime);
            popI->scaledtimeOriginInput = proposedTime / theta;
            
            if (mcmcOptions.verbose>=1)
                printf( "\n MRCA node id %d with input time %Lf and scaled %Lf and parent with input time %Lf and scaled %Lf was assigned to the oldest pop %d and time of origin %Lf, scaled %Lf \n", MRCA->node_index,u->timeInputTreeUnits, u->scaledByThetaTimeInputTreeUnits, v->timeInputTreeUnits,  v->scaledByThetaTimeInputTreeUnits, popI->index, popI->timeOriginInput, popI->scaledtimeOriginInput  );
            
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
        printf( "\n MRCA node id %d with input time %.10Lf and scaled %.10Lf and parent with input time %.10Lf and scaled %.10Lf was assigned to the oldest pop %d and time of origin %.10Lf, scaled %.10Lf std %.10Lf \n", oldestPop->rMRCA->node_index,u->timeInputTreeUnits, u->scaledByThetaTimeInputTreeUnits, v->timeInputTreeUnits,  v->scaledByThetaTimeInputTreeUnits, oldestPop->index, oldestPop->timeOriginInput, oldestPop->scaledtimeOriginInput,  oldestPop->timeOriginSTD  );
    //update the father node of the MRCA if needed
    if (v->timeInputTreeUnits < oldestPop->timeOriginInput){
        //printf( "\n updated node from %.10Lf to %.10Lf", v->timeInputTreeUnits, oldestPop->timeOriginSTD);
        v->timeInputTreeUnits = oldestPop->timeOriginInput;
        v->scaledByThetaTimeInputTreeUnits =  v->timeInputTreeUnits / theta;
        assert(oldestPop->x >0);
        v->timePUnits =  v->scaledByThetaTimeInputTreeUnits / oldestPop->x;
        oldestPop->rMRCA->parent->length =  v->timeInputTreeUnits - u->timeInputTreeUnits ;
        oldestPop->rMRCA->length= v->timeInputTreeUnits - u->timeInputTreeUnits;
    }
    
}
void Chain::initPopulationCoalescentAndMigrationEventsFromRootNodeOnTree(pll_unode_t *p, Population *currentPopulation, std::map<pll_unode_t*, Population*> &mrcaOfPopulation )
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
bool Chain::initPopulationsSampleSizes( std::map<pll_rnode_t*, vector<Population*>>  &rmrcaOfPopulation, string &healthyTipLabel){
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
void Chain::initPopulationsCoalescentAndMigrationEventsFromRootedTree(std::map<pll_rnode_t*, vector<Population*>> &rmrcaOfPopulation, string& healthyTipLabel){
    
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
        if (p->left == NULL &&  p->right == NULL  )//tip
        {
            assert(p->label !=NULL);
            if (std::string(p->label).compare(healthyTipLabel)!=0)//different from healthy tip
                currentPopulation->rtips.push_back(p);
            else
                return;//nothing to do with the healthy cell
        }
        else
        {
            //initialize currentPopulation
            if (p->parent == NULL)//the root ot the tree
                // && (std::string(p->right->label).compare(healthyTipLabel)==0) && rmrcaOfPopulation.count(p->left)>0)
            {
                if (p->right!=NULL && p->right->label!=NULL && (std::string(p->right->label).compare(healthyTipLabel)==0))
                {
                    
                    if (rmrcaOfPopulation.find(p->left)!= rmrcaOfPopulation.end())
                    {   updateFatherPopOnSameEdge(rmrcaOfPopulation, p->left,NULL);
                        currentPopulation = rmrcaOfPopulation[p->left].front();
                        initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(p->left, currentPopulation, rmrcaOfPopulation , healthyTipLabel);
                        return;
                    }
                }
                else if(p->left!=NULL && p->left->label!=NULL && (std::string(p->left->label).compare(healthyTipLabel)==0)){
                    
                    if (rmrcaOfPopulation.find(p->right)!= rmrcaOfPopulation.end())
                    {
                        updateFatherPopOnSameEdge(rmrcaOfPopulation, p->right,NULL);
                        currentPopulation = rmrcaOfPopulation[p->right].front();
                        initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(p->right, currentPopulation, rmrcaOfPopulation , healthyTipLabel);
                        return;
                    }
                }
            }
            
            //if is not the root ot the tree
            TreeNode * u= (TreeNode *)p->data;
            if (p->parent != NULL &&  currentPopulation != NULL )
            {
                //if is not the root of the tree(the parent of the MRCA of the oldest population
                //u->timePUnits=u->scaledByThetaTimeInputTreeUnits / (currentPopulation->x);
                assert(currentPopulation->x >0);
                assert(u->timePUnits == u->scaledByThetaTimeInputTreeUnits / (currentPopulation->x) );
                
                currentPopulation->CoalescentEventTimes.push_back(u->scaledByThetaTimeInputTreeUnits / (currentPopulation->x) );
                //printf("\n coalescent event node index %d for pop %d, time input tree %lf, scaled time %lf, model time  %lf\n ", p->clv_index,currentPopulation->index, u->timeInputTreeUnits, u->scaledByThetaTimeInputTreeUnits, u->scaledByThetaTimeInputTreeUnits / (currentPopulation->effectPopSize));
                currentPopulation->numCompletedCoalescences= currentPopulation->numCompletedCoalescences +1;
            }
            if(p->left !=NULL && rmrcaOfPopulation.find(p->left) != rmrcaOfPopulation.end() )
            { // in this branch we have  migration
                
                updateFatherPopOnSameEdge(rmrcaOfPopulation, p->left,currentPopulation);
                
                initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(p->left, rmrcaOfPopulation[p->left].front(), rmrcaOfPopulation, healthyTipLabel);
            }
            else{
                initPopulationCoalescentAndMigrationEventsFromNodeOnRootedTree(p->left, currentPopulation, rmrcaOfPopulation, healthyTipLabel);
            }
            if(p->right !=NULL && rmrcaOfPopulation.find(p->right)  != rmrcaOfPopulation.end() )
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


void Chain::initPopulationSampleSizesFromNodeOnRootedTree(pll_rnode_t *p, Population *currentPopulation, std::map<pll_rnode_t*, vector<Population*>>  &rmrcaOfPopulation, string &healthyTipLabel)
{
    if (p!=NULL)
    {
        TreeNode * treeNode=(TreeNode *)(p->data);
        //initialize currentPopulation
        if (p->parent == NULL && (std::string(p->right->label).compare(healthyTipLabel)==0) && rmrcaOfPopulation.find(p->left) != rmrcaOfPopulation.end())
        {
            
            currentPopulation = getYoungestPopulationOnEdge(p->left, rmrcaOfPopulation);
            initPopulationSampleSizesFromNodeOnRootedTree(p->left, currentPopulation, rmrcaOfPopulation, healthyTipLabel);
            return;
            
        }
        else if (p->parent == NULL && (std::string(p->left->label).compare(healthyTipLabel)==0) && rmrcaOfPopulation.find(p->right) != rmrcaOfPopulation.end())
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
            if (rmrcaOfPopulation.find(p->left)!= rmrcaOfPopulation.end() )
            {
                
                initPopulationSampleSizesFromNodeOnRootedTree(p->left, ((vector<Population*>)rmrcaOfPopulation[p->left]).front(), rmrcaOfPopulation, healthyTipLabel );
            }
            else{
                initPopulationSampleSizesFromNodeOnRootedTree(p->left, currentPopulation, rmrcaOfPopulation, healthyTipLabel );
            }
            if (rmrcaOfPopulation.find(p->right)!= rmrcaOfPopulation.end() )
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
Population* Chain::getYoungestPopulationOnEdge(pll_rnode_t* p, std::map<pll_rnode_t*, vector<Population*> >  &rmrcaOfPopulation)
{
    Population* youngestPop=NULL;
    vector<Population*> vec;
    int numberOccurNode = rmrcaOfPopulation.count(p);
    int numberPoints;
    map<pll_rnode_t*, vector<Population*>>::iterator it ;
    if (p!=NULL &&  numberOccurNode > 0)
    {
        it=rmrcaOfPopulation.find(p);
        if(it != rmrcaOfPopulation.end())
        {
            vec = it->second;//rmrcaOfPopulation[p];
            numberPoints =  vec.size();
            youngestPop = vec.front();
        }
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
void Chain::samplePopulationDeltaFromPriors(MCMCoptions &mcmcOptions, long int *seed ,const gsl_rng * rngGsl)
{
    int i;
    Population *popI;
    long double  randomDelta;
    for( i = 0 ; i < numClones; i++)
    {
        popI=populations[i];
        randomDelta = Random::RandomLogUniform(mcmcOptions.Deltafrom, mcmcOptions.Deltato, mcmcOptions.useGSLRandomGenerator, rngGsl, NULL);
        //popI->delta = chain->proportionsVector[i] * randomDelta;
        popI->delta =  randomDelta;
        popI->growthRate =popI->delta  / popI->x;
        popI->deathRate= popI->birthRate - popI->growthRate;
    }
}
void Chain::samplePopulationGrowthRateFromPriors(MCMCoptions &mcmcOptions, const gsl_rng * rngGsl )
{
    int i;
    Population *popI;
    long double  randomGrowthRate;
    long double  logPriorDensity;
    long double from =0.0;
    for( i = 0 ; i < numClones; i++)
    {
        popI=populations[i];
        if (mcmcOptions.priorsType ==0)
        {
            randomGrowthRate = Random::RandomLogUniform(mcmcOptions.GrowthRatefrom, mcmcOptions.GrowthRateto, mcmcOptions.useGSLRandomGenerator, rngGsl, NULL);
            logPriorDensity=Distributions::LogUniformDensity(randomGrowthRate,mcmcOptions.GrowthRatefrom, mcmcOptions.GrowthRateto);
        }
        if (mcmcOptions.priorsType ==2)
        {
            randomGrowthRate =Random::RandomPowerLawDistribution(mcmcOptions.parameterPowerLawDistributionGrowthRate, from, mcmcOptions.useGSLRandomGenerator, rngGsl, NULL);
            logPriorDensity=Distributions::LogPowerLawDistibutionDensity(randomGrowthRate, mcmcOptions.parameterPowerLawDistributionGrowthRate,from);
            
        }
        if (mcmcOptions.priorsType ==1)
        {
            randomGrowthRate =Random::RandomExponentialStartingFrom(mcmcOptions.lambdaExponentialPriorGrowthRate,from, mcmcOptions.useGSLRandomGenerator, rngGsl, NULL);
            logPriorDensity=Distributions::LogExponentialDensity(randomGrowthRate,from,mcmcOptions.lambdaExponentialPriorGrowthRate);
            
        }
        
        popI->deltaT =  randomGrowthRate;
        popI->DeltaT->setParameterValue(randomGrowthRate);
        
        if (mcmcOptions.verbose>=2)
            std::cout << "\n initial population DeltaT "  << popI->deltaT << "with log prior value: " <<  popI->DeltaT->getParameter()->getLogPriorCurrentValue() << std::endl;
        
        
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
        else if(mrca->parent != NULL &&  mrca->parent->parent == NULL &&  currentMrcaOfPopulation[mrca].size() >=1 && pop->timeOriginInput == ((TreeNode *)mrca->parent->data)->timeInputTreeUnits )// if it is the edge from the tumor MRCA to the root and  we have  at least 1 event and is the oldest population(the time of origin is the root of the tree)
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

std::map<pll_rnode_t*, vector<Population*>> Chain::chooseAvailableEdgeOnRootedTreeForPopulation(Population *pop, std::map<pll_rnode_t*, vector<Population*>> &currentMrcaOfPopulation, string &healthyCellLabel, const gsl_rng* rngGsl)
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
        
        random =Random::randomUniformFromGsl2(rngGsl);
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
        printf( "\n MRCA node id %d with time %Lf and parent with time %Lf was assigned to pop %d with order %d and time of origin %Lf, scaled %Lf \n", MRCA->node_index, u->scaledByThetaTimeInputTreeUnits, v->scaledByThetaTimeInputTreeUnits, pop->index, pop->order,  pop->timeOriginInput, pop->scaledtimeOriginInput );
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
long double  Chain::sumAvailableBranchLengths(std::map<pll_rnode_t*, vector<Population*>> &currentMRCAPopulation)
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
long double  Chain::sumAdjacentEdges(std::map<pll_rnode_t*, vector<Population*>> &currentMRCAPopulation, pll_rnode_t* MRCA, Population *pop)
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
std::map<pll_rnode_t*, vector<Population*>> Chain::chooseNewTimeofOriginOnEdge(Population *pop, MCMCoptions &mcmcOptions, const gsl_rng* rngGsl, boost::mt19937* rngBoost)
{
    TreeNode * u, *v;
    double m;
    long double randomNumber;
    long double currentTransformedTime, newTransformedTime;
    std::map<pll_rnode_t*, vector<Population*>> copyMRCAOfPopulation;
    copyMRCAOfPopulation.clear();
    copyMRCAOfPopulation.insert(rMRCAPopulation.begin(), rMRCAPopulation.end());
    
    
    u= (TreeNode *)( pop->rMRCA->data);
    v= (TreeNode *)( pop->rMRCA->parent->data);
    
    // pop->setLowerBoundTimeOriginInput(u->timeInputTreeUnits);
    
    long double proposedTimeOriginInput=0.0;
    long double proposedScaledTimeOriginInput=0.0;
    long double proposedTimeOriginSTD=0.0;
    
    if(mcmcOptions.kernelType ==1)
    {
        proposedTimeOriginSTD= Random::randomNormalGreaterThan(pop->timeOriginInput, mcmcOptions.sigmaNormalKernelTimeofOrigin, u->timeInputTreeUnits, true, rngGsl,rngBoost);
        proposedScaledTimeOriginInput = proposedTimeOriginSTD *pop->x;
        proposedTimeOriginInput= proposedScaledTimeOriginInput *theta;
    }
    else if(mcmcOptions.kernelType ==0)
    {
        if(mcmcOptions.useGSLRandomGenerator)
            randomNumber= Random::randomUniformFromGsl2(rngGsl);
        else
            randomNumber= Random::randomUniformBoost(rngBoost);
        
        long double kernelParameter = pop->TimeOriginSTD->getKernelParameter();
        // m = exp(2 * log(mcmcOptions.paramMultiplierTimeOriginOldestPop) *(randomNumber -0.5));
        
        m = exp(2 * log(kernelParameter) *(randomNumber -0.5));
        
        if (mcmcOptions.doMCMCMoveTimeOriginInputOldestPop)
        {
            if (mcmcOptions.verbose >=2)
                std::cout << "\n move time origin input \n" << std::endl;
            currentTransformedTime = (mcmcOptions.upperBoundTimeOriginInputOldestPop - pop->timeOriginInput) /(pop->timeOriginInput - pop->lowerBoundTimeOriginInput);
            newTransformedTime = m * currentTransformedTime;
            proposedTimeOriginInput= ( mcmcOptions.upperBoundTimeOriginInputOldestPop + (newTransformedTime * pop->lowerBoundTimeOriginInput))  /( 1+newTransformedTime);
            proposedScaledTimeOriginInput= proposedTimeOriginInput / theta;
            proposedTimeOriginSTD = proposedScaledTimeOriginInput / pop->x;
        }
        else{
            if (mcmcOptions.verbose >=2)
                std::cout << "\n move time origin STD \n" << std::endl;
            pop->TimeOriginSTD->stepKernel(rngGsl);
            proposedTimeOriginSTD=pop->TimeOriginSTD->getParameter()->getCurrentValue();
            
            //proposedTimeOriginSTD = m*  pop->timeOriginSTD;
            proposedScaledTimeOriginInput = proposedTimeOriginSTD *pop->x;
            proposedTimeOriginInput= proposedScaledTimeOriginInput *theta;
            
        }
        // if (proposedTime * pop->theta < pop->lowerBoundTimeOriginInput)
        // {
        //printf( "\n reflected back from %.10Lf to  %.10Lf, mult by theta: %.10Lf, bound: %.10Lf \n", proposedTime,(pop->lowerBoundTimeOriginInput * pop->lowerBoundTimeOriginInput) / proposedTime, proposedTime * pop->theta, pop->lowerBoundTimeOriginInput );
        //     proposedTime = (pop->lowerBoundTimeOriginInput * pop->lowerBoundTimeOriginInput) / (pop->theta* pop->theta *proposedTime);
        
        // }
    }
    //printf( "\n proposed time of origin input tree units %Lf,  old %Lf \n", proposedTime,pop->timeOriginInput );
    
    
    pop->timeOriginSTD =proposedTimeOriginSTD;
    pop->scaledtimeOriginInput = proposedScaledTimeOriginInput;
    pop->timeOriginInput=proposedTimeOriginInput;
    
    
    //pop->TimeOriginInput->setParameterValue(proposedTimeOriginInput);
    
    if (mcmcOptions.verbose>=2)
        std::cout << "\n proposed new timeOriginSTD "<< pop->timeOriginSTD << " , " << pop->oldTimeOriginSTD << std::endl;
    
    //update the node if needed
    //if (!mcmcOptions.noData && v->timeInputTreeUnits < pop->timeOriginInput)
    if ( v->timeInputTreeUnits < pop->timeOriginInput)
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
    
    if (programOptions.doUseFixedTree == YES)
    {
        sprintf(File,"%s/Chain%d", filePaths.resultsDir,  chainNumber );
        mkdir(File,S_IRWXU);
        
    }
    //strcpy (resultsDir, dir);
    if (programOptions.doPrintTrees == YES)
    {
        if (programOptions.doUseFixedTree == YES)
        {
            sprintf(File,"%s/Chain%d/%s", filePaths.resultsDir,  chainNumber ,filePaths.treeDir);
            mkdir(File,S_IRWXU);
            sprintf(File,"%s/Chain%d/%s/%s.trees", filePaths.resultsDir, chainNumber, filePaths.treeDir,filePaths.treeFile);
            sprintf(files.fpTrees->path,"%s/Chain%d/%s/%s.trees", filePaths.resultsDir, chainNumber, filePaths.treeDir,filePaths.treeFile);
        }
        else
        {
            sprintf(File,"%s/Chain%d_%s.trees", filePaths.resultsDir, chainNumber, filePaths.treeFile);
            sprintf(files.fpTrees->path,"%s/Chain%d_%s.trees", filePaths.resultsDir, chainNumber, filePaths.treeFile);
        }
        
        
        if (openFile(&files.fpTrees->f, File) == -1)
        {
            fprintf (stderr, "Can't open \"%s\"\n", File);
            exit(-1);
        }
        if (programOptions.doUseFixedTree == YES)
        {
            sprintf(File,"%s/Chain%d/%s/%s_2.trees", filePaths.resultsDir,chainNumber,  filePaths.treeDir, filePaths.treeFile);
            sprintf(files.fpTrees2->path,"%s/Chain%d/%s/%s_2.trees", filePaths.resultsDir,chainNumber,  filePaths.treeDir, filePaths.treeFile);
        }
        else{
            
            sprintf(File,"%s/Chain%d_%s_2.trees", filePaths.resultsDir,chainNumber,   filePaths.treeFile);
            sprintf(files.fpTrees2->path,"%s/Chain%d_%s_2.trees", filePaths.resultsDir,chainNumber,   filePaths.treeFile);
            
        }
        
        if (openFile(&files.fpTrees2->f, File) == -1)
        {
            fprintf(stderr, "Can't open %s.\n", File);
            exit(-1);
        }
    }
    if (programOptions.doPrintTimes == YES)
    {
        if (programOptions.doUseFixedTree == YES)
        {
            sprintf(File,"%s/Chain%d/%s", filePaths.resultsDir,chainNumber, filePaths.timesDir );
            mkdir(File,S_IRWXU);
            sprintf(File,"%s/Chain%d/%s/%s.txt", filePaths.resultsDir, chainNumber, filePaths.timesDir ,filePaths.timesFile);
            sprintf(files.fpTimes->path,"%s/Chain%d/%s/%s.txt", filePaths.resultsDir, chainNumber, filePaths.timesDir ,filePaths.timesFile);
        }
        else{
            sprintf(File,"%s/Chain%d_%s.txt", filePaths.resultsDir, chainNumber, filePaths.timesFile);
            sprintf(files.fpTimes->path,"%s/Chain%d_%s.txt", filePaths.resultsDir, chainNumber, filePaths.timesFile);
        }
        
        if (openFile(&files.fpTimes->f, File) == -1)
        {
            fprintf (stderr, "Can't open \"%s\"\n", File);
            exit(-1);
        }
        if (programOptions.doUseFixedTree == YES){
            //sprintf(File,"%s/%s", filePaths.resultsDir, filePaths.timesDir);
            sprintf(File,"%s/Chain%d/%s/%s_2.txt", filePaths.resultsDir,chainNumber, filePaths.timesDir, filePaths.timesFile);
            sprintf(files.fpTimes2->path,"%s/Chain%d/%s/%s_2.txt", filePaths.resultsDir,chainNumber, filePaths.timesDir, filePaths.timesFile);
        }
        else {
            sprintf(File,"%s/Chain%d_%s_2.txt", filePaths.resultsDir,chainNumber,  filePaths.timesFile);
            sprintf(files.fpTimes2->path,"%s/Chain%d_%s_2.txt", filePaths.resultsDir,chainNumber,  filePaths.timesFile);
        }
        
        if (openFile(&files.fpTimes2->f, File) == -1)
        {
            fprintf(stderr, "Can't open %s.\n", File);
            exit(-1);
        }
    }
    //log file
    
    if (programOptions.doUseFixedTree == YES)
    {
        sprintf(File,"%s/Chain%d/%s.log", filePaths.resultsDir, chainNumber,  filePaths.logFile);
        sprintf(files.fplog->path,"%s/Chain%d/%s.log", filePaths.resultsDir, chainNumber,  filePaths.logFile);
    }
    else{
        sprintf(File,"%s/Chain%d_%s.log", filePaths.resultsDir, chainNumber,  filePaths.logFile);
        sprintf(files.fplog->path,"%s/Chain%d_%s.log", filePaths.resultsDir, chainNumber,  filePaths.logFile);
    }
    
    if (openFile(&files.fplog->f, File) == -1)
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
    
    fprintf (files.fplog->f, "%d\t", currentIteration);
    fprintf (files.fplog->f, "%.20Lf\t", currentlogConditionalLikelihoodTree);
    fprintf (files.fplog->f, "%.20Lf\t", currentlogConditionalLikelihoodSequences);
    //
    fprintf (files.fplog->f, "%.20Lf\t", theta);
    fprintf (files.fplog->f, "%.20Lf\t", exp(thetaPar->getParameter()->getLogPriorCurrentValue()));
    for (unsigned int i = 0; i < populations.size(); ++i){
        popI =populations[i];
        //fprintf (files.fplog, "%.8Lf\t", pop->effectPopSize);
        fprintf (files.fplog->f, "%.20Lf\t", popI->delta);
        //fprintf (files.fplog, "%.40Lf\t", pop->growthRate);
        fprintf (files.fplog->f, "%.10Lf\t", popI->x);
        fprintf (files.fplog->f, "%.20Lf\t", popI->deltaT);
        fprintf (files.fplog->f, "%.20Lf\t", exp(popI->DeltaT->getParameter()->getLogPriorCurrentValue()));
        //if (currentIteration<= mcmcOptions.iterationToComputeThinnig)
        //    posteriorDeltaT.push_back(popI->deltaT);
        fprintf (files.fplog->f, "%.20Lf\t", popI->theta);
        // fprintf (files.fplog, "%.40Lf\t", mutationRate);
        fprintf (files.fplog->f, "%d\t", popI->sampleSize);
        fprintf (files.fplog->f, "%.20Lf\t", popI->timeOriginInput);//multiplied by theta_i
        fprintf (files.fplog->f, "%.20Lf\t", popI->timeOriginSTD);//T
        //if (currentIteration<= mcmcOptions.iterationToComputeThinnig)
        //  posteriorT.push_back(popI->timeOriginSTD);
        fprintf (files.fplog->f, "%.20Lf\t", popI->scaledtimeOriginInput);//divided by theta_i
        fprintf (files.fplog->f, "%.20Lf\t", popI->delta*popI->timeOriginSTD);//divided by theta_i
        // fprintf (files.fplog, "%.20Lf\t", popI->timeOriginInput / mutationRate);//physical time(if lambda=1 is the time in generations)
    }
    
    if (currentIteration <= mcmcOptions.Niterations -1)
        fprintf (files.fplog->f, "\n");
}
void Chain::saveMCMCState( int  currentIteration, ProgramOptions &programOptions,  MCMCoptions &mcmcOptions)
{
    //    fprintf (files.fplog, "%d  %lf %lf %lf effect_pop_size(pop1) effect_pop_size(pop2) effect_pop_size(pop3) sample_size(pop1)  sample_size(pop2) sample_size(pop3) time_origin(pop1)  time_origin(pop2) time_origin(pop3) \n", pop);
    string paramName;
    Population *popI;
    
    posteriorTheta.push_back(theta);
    thetaPar->saveCurrentValue();
    thetaPar->saveCurrentValueToList();
    
    //vector<long double > row(stored.cols());
    //int idx = getIndexMCMCParameter(thetaPar);
    
    //assert(idx >=0 );
    //row.at(idx)=theta;
    //stored(currentIteration, idx)  = theta;
    
    for (size_t i = 0; i < populations.size(); ++i)
    {
        popI =populations[i];
        popI->savePosteriorValues();
        
        //idx = getIndexMCMCParameter(popI->DeltaT);
        //assert(idx >=0 & idx < stored.cols());
        //row.at(idx)=popI->DeltaT->getParameter()->getCurrentValue();
        //  stored(currentIteration, idx)  = popI->DeltaT->getParameter()->getCurrentValue();
        
        //  idx = getIndexMCMCParameter(popI->TimeOriginSTD);
        // assert(idx >=0);
        //row.at(idx)=popI->TimeOriginSTD->getParameter()->getCurrentValue();
        // stored(currentIteration, idx)  =popI->TimeOriginSTD->getParameter()->getCurrentValue();
    }
    
}


void Chain::resizeStoredMCMCparameters(int rows, int cols){
    //    assert(rows>0);
    //    assert(cols>0);
    //    if (stored.rows()==0)
    //      stored.resize(rows, cols);
    //    else
    //      stored.conservativeResize(rows, cols);
    
}
void Chain::writeHeaderOutputChain(const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files, const pll_rtree_t* trueTree, const long double trueTheta,
                                   const vector<long double> trueDeltaTs,
                                   const vector<long double> trueTs,
                                   const vector<int> trueSampleSizes, StructuredCoalescentTree* tree)
{
    
    //fout.open(files.fplog->path, ios::out | ios::app);
    //fout << std::setprecision(15);
    string paramName;
    
    
    std::time_t timeNow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    
    assert(trueDeltaTs.size()==numClones);
    assert(trueTs.size()==numClones);
    assert(trueSampleSizes.size()==numClones);
    
    //    fout <<  "#Tumor_coal" << '\n';
    //    fout <<  "# Generated %s " << std::ctime(&timeNow) << '\n';
    //    fout << "# true values: theta=" << trueTheta ;
    //    for (size_t  i = 0; i < numClones; ++i){
    //        fout <<" deltaT(pop" << std::to_string(i) << ")=" << trueDeltaTs.at(i);
    //        fout <<" time_origin_std(pop" << std::to_string(i) << ")" << trueTs.at(i);
    //        fout <<" prop_effect_pop_size(pop" << std::to_string(i) << ")" << trueSampleSizes.at(i);
    //    }
    //    fout <<  "\n";
    //    fout << "state\t";
    //    fout << "loglikTree\t";
    //    fout << "loglikSeq\t";
    //    paramName =  "theta";
    //    fout << "%s\t";
    //    for (size_t  i = 0; i < populations.size(); ++i){
    //        // paramName =  "effect_pop_size(pop" + std::to_string(i) + ")";
    //        // fprintf (files.fplog, "%s\t", paramName.c_str());
    //        paramName =  "delta(pop" + std::to_string(i) + ")";
    //        fout <<  paramName.c_str() << "%s\t";
    //        // paramName =  "growth_rate(pop" + std::to_string(i) + ")";
    //        //fprintf (files.fplog, "%s\t", paramName.c_str());
    //        paramName =  "prop_effect_pop_size(pop" + std::to_string(i) + ")";
    //        fout <<  paramName.c_str() << "%s\t";
    //        paramName =  "deltaT(pop" + std::to_string(i) + ")";
    //        fout <<  paramName.c_str() << "%s\t";
    //        paramName =  "theta(pop" + std::to_string(i) + ")";
    //        fout <<  paramName.c_str() << "%s\t";
    //        //paramName =  "mut_rate(pop" + std::to_string(i) + ")";
    //        //fprintf (files.fplog, "%s\t", paramName.c_str());
    //        paramName =  "sample_size(pop" + std::to_string(i) + ")";
    //        fout <<  paramName.c_str() << "%s\t";
    //        paramName =  "time_origin_input(pop" + std::to_string(i) + ")";
    //        fout <<  paramName.c_str() << "%s\t";
    //        paramName =  "time_origin_std(pop" + std::to_string(i) + ")";
    //        fout <<  paramName.c_str() << "%s\t";
    //        paramName =  "scaled_time_origin(pop" + std::to_string(i) + ")";
    //        fout <<  paramName.c_str() << "%s\t";
    //        // paramName =  "physical_time(pop" + std::to_string(i) + ")";
    //        // fprintf (files.fplog, "%s\t", paramName.c_str());
    //    }
    //    fout <<  "\n";
    //    fout.close();
    
    fprintf (files.fplog->f, "#Tumor_coal\t");
    fprintf (files.fplog->f, "#Generated %s\t", std::ctime(&timeNow));
    
    fprintf (files.fplog->f, "#true values  theta=%.15Lf, ", trueTheta);
    fprintf (files.fplog->f, "loglikTree=%.15Lf, ", tree->getLogLikelihoodTree());
    for (int  i = 0; i < numClones; ++i){
        fprintf (files.fplog->f, "deltaT(pop%d)=%.15Lf, ",i, trueDeltaTs.at(i));
        fprintf (files.fplog->f, "time_origin_std(pop%d)=%.15Lf, ",i, trueTs.at(i));
        fprintf (files.fplog->f, "prop_effect_pop_size(pop%d)=%d, ",i, trueSampleSizes.at(i) );
        fprintf (files.fplog->f, "Delta*T(pop%d)=%.15Lf\n",i, trueDeltaTs.at(i)* trueTs.at(i));
    }
    
    fprintf (files.fplog->f, "state\t");
    fprintf (files.fplog->f, "loglikTree\t");
    fprintf (files.fplog->f, "loglikSeq\t");
    
    paramName =  "theta";
    fprintf (files.fplog->f, "%s\t", paramName.c_str());
    paramName =  "prior_theta";
    fprintf (files.fplog->f, "%s\t", paramName.c_str());
    for (size_t  i = 0; i < populations.size(); ++i)
    {
        // paramName =  "effect_pop_size(pop" + std::to_string(i) + ")";
        // fprintf (files.fplog->f, "%s\t", paramName.c_str());
        paramName =  "delta(pop" + std::to_string(i) + ")";
        fprintf (files.fplog->f, "%s\t", paramName.c_str());
        // paramName =  "growth_rate(pop" + std::to_string(i) + ")";
        //fprintf (files.fplog->f, "%s\t", paramName.c_str());
        paramName =  "prop_effect_pop_size(pop" + std::to_string(i) + ")";
        fprintf (files.fplog->f, "%s\t", paramName.c_str());
        paramName =  "deltaT(pop" + std::to_string(i) + ")";
        fprintf (files.fplog->f, "%s\t", paramName.c_str());
        paramName =  "prior_deltaT(pop" + std::to_string(i) + ")";
        fprintf (files.fplog->f, "%s\t", paramName.c_str());
        paramName =  "theta(pop" + std::to_string(i) + ")";
        fprintf (files.fplog->f, "%s\t", paramName.c_str());
        //paramName =  "mut_rate(pop" + std::to_string(i) + ")";
        //fprintf (files.fplog->f, "%s\t", paramName.c_str());
        paramName =  "sample_size(pop" + std::to_string(i) + ")";
        fprintf (files.fplog->f, "%s\t", paramName.c_str());
        paramName =  "time_origin_input(pop" + std::to_string(i) + ")";
        fprintf (files.fplog->f, "%s\t", paramName.c_str());
        paramName =  "time_origin_std(pop" + std::to_string(i) + ")";
        fprintf (files.fplog->f, "%s\t", paramName.c_str());
        paramName =  "scaled_time_origin(pop" + std::to_string(i) + ")";
        fprintf (files.fplog->f, "%s\t", paramName.c_str());
        paramName =  "Delta*T(pop" + std::to_string(i) + ")";
        fprintf (files.fplog->f, "%s\t", paramName.c_str());
        // paramName =  "physical_time(pop" + std::to_string(i) + ")";
        // fprintf (files.fplog->f, "%s\t", paramName.c_str());
    }
    fprintf (files.fplog->f, "\n");
}
void Chain::saveTrueTreeInfo(pll_rtree_t * trueTree, char *  rootedNewick3,  ProgramOptions& programOptions)
{
    if (programOptions.doPrintTrees == YES)
    {
        Output::PrintTrees(0, trueTree->root, files.fpTrees->f, programOptions.mutationRate, programOptions.doUseObservedCellNames);
        Output::PrintTrees2(0, trueTree->root, files.fpTrees2->f, programOptions.mutationRate, NULL, NO);
    }
    if (programOptions.doPrintTimes == YES)
    {
        vector<pll_rnode_t *> trueTreeNodes=Utils::filterHealthyTip(trueTree,trueTree->nodes,trueTree->inner_count + trueTree->tip_count, programOptions.healthyTipLabel );
        //  vector<pll_rnode_t *> trueTreeNodes=Utils::vectorFromDoublePointer(trueTree->nodes,trueTree->inner_count + trueTree->tip_count);
        Output::PrintTimes(0, files.fpTimes->f, programOptions.mutationRate, trueTreeNodes, NO);
        Output::PrintTimes2(0, files.fpTimes2->f, programOptions.mutationRate, trueTreeNodes, NO);
    }
}
//following https://www.csee.usf.edu/~kchriste/tools/autoc.c
long double  Chain::autoCorrelation(int lag, vector<long double> &values, long double meanP, long double varianceP)
{
    long double  result=0;
    long double  autoCovariance=0;
    //long double  variance=0;
    if (values.size() >0)
    {
        //long double  mean = std::accumulate( values.begin(), values.end(), 0.0)/ values.size();
        //long double mean1 = mcmcPar->getCurrentMean();
        // assert(abs(meanP-mean)<0.0000000001);
        // long double variance1 = mcmcPar->getCurrentVariance();
        assert(values.size() -lag >0);
        for (size_t  i = 0; i < (values.size() -lag) ; ++i)
        {
            autoCovariance += (values.at(i)- meanP)*(values.at(i+lag)- meanP);
            //variance += pow(values.at(i)- meanP,2);
        }
        autoCovariance=(1.0 / (values.size() - lag)) * autoCovariance;
        //        for (size_t j = values.size() -lag ; j < values.size()  ; ++j)
        //        {
        //            variance += pow(values.at(j)- mean,2);
        //        }
        //variance = variance / values.size();
        // assert(abs(varianceP- variance)<0.0001);
        result = autoCovariance/ varianceP;
    }
    return result;
}
long double  Chain::autoCorrelation(vector<int> &lagsVector, vector<long double> &values, long double correlationThreshold, int &indexLag, long double meanP, long double varianceP)
{
    long double  result=2;
    
    int k=0;
    
    if (values.size() >0)
    {
        if(lagsVector.at(k) < values.size()){
            //result= autoCorrelation(lagsVector.at(k), values);
            while (abs(result) > correlationThreshold && k < lagsVector.size())
            {
                
                result=autoCorrelation(lagsVector.at(k), values, meanP,  varianceP);
                k++;
            }
            
            if (k < lagsVector.size() )
                indexLag= k -1 ;
            else
                indexLag =lagsVector.size()-1;
            
        }
    }
    return result;
}
int  Chain::getLagWithCorrelationCloseToZero(vector<long double>& correlationVector, long double correlationThreshold)
{
    long double  result=-1;
    long double sumConsecutive = 0.0;
    
    for (size_t i = 0; i < (correlationVector.size()-1) ; i+=2)
    {
        sumConsecutive = correlationVector[i] + correlationVector[i+1];
        if ( sumConsecutive <=0)
        {
            result = i;
            return result;
        }
    }
    
    assert(result >=0);
    return result;
}
vector<int> Chain::vectorLags(int valuesSize){
    int lag=min(100, valuesSize );
    vector<int> result;
    while(lag < valuesSize){
        result.push_back(lag);
        lag+=100;
        
    }
    return result;
}
long double  Chain::ESS(int lag, vector<long double> &values)
{
    long double  result=0;
    long double  sumAutoCorrelations =0;
    for (size_t i = 0; i < values.size()  ; ++i)
    {
        sumAutoCorrelations +=autoCorrelation(i, values, 0.0, 0.0);
    }
    result = values.size() / (1 + 2 * sumAutoCorrelations);
    return result;
}
bool Chain::checkMigrationsOrder()
{
    bool migrationsOrderCorrect = true;
    Population *popI;
    long double lastMigrationTime;
    for(size_t  i = 0 ; i < numClones; i++)
    {
        popI = populations[i];
        lastMigrationTime= popI->immigrantsPopOrderedByModelTime.at(popI->immigrantsPopOrderedByModelTime.size()-1).first ;
        if (lastMigrationTime != popI->timeOriginSTD)//check the  last position
        {
            migrationsOrderCorrect=false;
            fprintf(stderr, "\n the last migration time  is not the time of origin \n");
            break;
        }
        for( size_t i = 0 ; i < popI->immigrantsPopOrderedByModelTime.size()-1; i++)
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
    fclose(files.fplog->f);
}
void Chain::updateNodeInfoOldestPopulation(Population * oldestPop, double newMRCAtimeInputTreeUnits){
    TreeNode *u, *v;
    u= (TreeNode *)( oldestPop->rMRCA->data);
    v= (TreeNode *)( oldestPop->rMRCA->parent->data);
    //printf( "\n updated node from %.10Lf to %.10lf", v->timeInputTreeUnits, newMRCAtimeInputTreeUnits);
    v->timeInputTreeUnits = newMRCAtimeInputTreeUnits;
    
    v->scaledByThetaTimeInputTreeUnits =  v->timeInputTreeUnits / theta;
    assert(oldestPop->x >0);
    v->timePUnits  = v->scaledByThetaTimeInputTreeUnits / oldestPop->x;
    oldestPop->rMRCA->parent->length =  v->timeInputTreeUnits - u->timeInputTreeUnits ;
    oldestPop->rMRCA->length =  v->timeInputTreeUnits - u->timeInputTreeUnits ;
    
}
void Chain::drawModelTimeOriginFromConditionalDensity(Population * oldestPop, MCMCoptions &mcmcOptions, const gsl_rng * randomGenerator){
    
    long double modelTimeOriginOldestPop=Random::RandomDensityModelTimeOrigin(oldestPop->delta, mcmcOptions.useGSLRandomGenerator,  oldestPop->lowerBoundTimeOriginInput, randomGenerator, NULL );
    
    oldestPop->timeOriginSTD = modelTimeOriginOldestPop;
    
    oldestPop->scaledtimeOriginInput=oldestPop->timeOriginSTD*(oldestPop->x);
    
    oldestPop->timeOriginInput =oldestPop->scaledtimeOriginInput * theta;
    
    
    updateNodeInfoOldestPopulation(oldestPop, oldestPop->timeOriginInput);
    
}
void Chain::printMovesSummary(ProgramOptions &programOptions, MCMCoptions &mcmcOptions){
    
    vector<MCMCmove*>::iterator it;
    MCMCmove* currentMove;
    string name;
    double percentAcceptance=0.0;
    int totalNumberMoves=0;
    long double currentKernelParameterValue;
    long double updatedKernelParameterValue;
    for ( it=moves.begin(); it!=moves.end(); ++it)
    {
        currentMove=  *it;
        name =currentMove->name();
        totalNumberMoves = ((int)(currentMove->numberAccepted()) + (int)(currentMove->numberRejected()) );
        percentAcceptance = (double)(currentMove->numberAccepted())/ totalNumberMoves ;
        
        currentKernelParameterValue= currentMove->getKernelValue();
        
        updatedKernelParameterValue= Utils::tune(currentKernelParameterValue,percentAcceptance);
        
        //        if (percentAcceptance <0.40){
        //            currentMove->decreaseParameter(programOptions, mcmcOptions);
        //            std::cout << "\n increased parameter for move " << name.c_str();
        //        }
        //        else if (percentAcceptance >0.50){
        //            currentMove->increaseParameter(programOptions, mcmcOptions);
        //            std::cout << "\n decreased parameter for move " << name.c_str();
        //        }
        
        currentMove->setKernelValue(updatedKernelParameterValue);
        std::cout << "\n Number accepted moves of type " << name.c_str() << ": "<< currentMove->numberAccepted() <<"(" << percentAcceptance*100 << " %) out of " << totalNumberMoves << " moves " << std::endl;
        if (currentKernelParameterValue !=updatedKernelParameterValue)
            std::cout << "\n The lenth of the interval of the multiplier for move " << name.c_str() << " changed from  "<< currentKernelParameterValue <<"  to " <<  updatedKernelParameterValue  << std::endl;
        
        
        currentMove->resetNumberAccepted();
        currentMove->resetNumberRejected();
    }
    
}
void Chain::printLastMovesSummary(){
    vector<MCMCmove*>::iterator it;
    MCMCmove* currentMove;
    string name;
    double percentAcceptance=0.0;
    int totalNumberMoves=0;
    for ( it=moves.begin(); it!=moves.end(); ++it)
    {
        currentMove=  *it;
        name =currentMove->name();
        totalNumberMoves = ((int)(currentMove->numberAccepted()) + (int)(currentMove->numberRejected()) );
        percentAcceptance = (double)(currentMove->numberAccepted()) *100.00/ totalNumberMoves ;
        
        
        std::cout << "\n Number accepted moves of type " << name.c_str() << ": "<< currentMove->numberAccepted() <<"(" << percentAcceptance << " %%) out of " << totalNumberMoves << " moves " << endl;
    }
}
void  Chain::InitListPossibleMigrations()
{
    sort(populations.begin(), populations.end(), comparePopulationsByTimeOrigin);
    Population *p;
    
    for ( size_t i = 0; i < numClones; i++) {
        p = populations[i];
        p->InitListPossibleMigrations(i);
    }
}
void Chain::InitPopulationsCoalescentEvents( ) {
    
    Population *popI;
    for( size_t i = 0 ; i < numClones; i++)
    {
        popI = populations[i];
        popI->InitCoalescentEvents(numClones);
    }
}
void Chain::deleteRNodes()
{
    
    for (std::vector<pll_rnode_t *>::iterator i = rnodes.begin(); i != rnodes.end(); ++i)
    {
        delete *i;
    }
    delete healthyTip;
    rnodes.clear();
    rtreeTips.clear();
    
}
void Chain::addEdgeFromNode(pll_rnode_t *node ){
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
        edgeLengths[currentNumberEdgesRootedTree]=make_pair(edge->length, edge);
        edges[currentNumberEdgesRootedTree]=edge;
        currentNumberEdgesRootedTree++;
    }
    
}
void Chain::resetMRCAMap()
{
    //    Population *popI;
    //    for(std::map<pll_rnode_t*, std::vector<Population*>>::iterator itr = rMRCAPopulation.begin(); itr != rMRCAPopulation.end(); itr++)
    //    {
    //        delete itr->first;
    //        for(vector<Population *>::iterator it = itr->second.begin(); it != itr->second.end(); ++it) {
    //            popI = *it;
    //            delete popI;
    //        }
    //        (itr->second).clear();
    //    }
    rMRCAPopulation.clear();
}
void Chain::initEdgesRootedTree(int numberTips)
{
    
    assert(numberTips >0);
    edges.clear();
    edgeLengths.clear();
    
    //edgeLengths (2*numberTips-2, std::make_pair(0.0, NULL));
    for(size_t i = 0 ; i < 2*numberTips-2; i++)
    {
        edges.push_back(NULL);
        edgeLengths.push_back(std::make_pair(0.0, nullptr));
    }
}
void Chain::insertMRCAMap(pll_rnode_t   *r, Population *population)
{
    assert(r!=NULL);
    if (rMRCAPopulation.find(r)!= rMRCAPopulation.end())
    {
        //int numberPop =((vector<Population*>)rMRCAPopulation[r]).size();
        // assert(numberPop>0);
        vector<Population*> populationsWithSameMRCA =(vector<Population*>)rMRCAPopulation[r];
        if (std::find(populationsWithSameMRCA.begin(), populationsWithSameMRCA.end(), population)== populationsWithSameMRCA.end())
        {
            rMRCAPopulation[r].push_back(population);
        }
    }
    else{
        
        rMRCAPopulation.insert(std::pair<pll_rnode_t*, std::vector<Population*>>(r, {population}) );
    }
}
void Chain::initMRCAMap(){
    
    for (size_t i = 0; i < numClones; ++i)
    {
        auto popI =  populations[i];
        
        insertMRCAMap(popI->rMRCA, popI);
    }
}
void Chain::closeFiles( ProgramOptions &programOptions){
    fclose(files.fplog->f);
}
void Chain::closeTrueTreeFiles( ProgramOptions &programOptions)
{
    if (programOptions.doPrintTrees == YES && programOptions.doUseFixedTree==NO)
    {
        fclose(files.fpTrees->f);
        fclose(files.fpTrees2->f);
    }
    if (programOptions.doPrintTimes == YES && programOptions.doUseFixedTree==NO)
    {
        fclose(files.fpTimes->f);
        fclose(files.fpTimes2->f);
    }
}
void Chain::computeThinnig(MCMCoptions & mcmcOptions){
    
    long double autocorrPosteriorT=0.0;
    long double autocorrPosteriorDeltaT=0.0;
    long double autocorrPosteriorProportion=0.0;
    long double autocorrPosteriorTheta=0.0;
    int k=0;
    
    // vector<int> lags= vectorLags(posteriorTheta.size());
    //std::sort (lags.begin(), lags.end());
    
    vector<long double > vec =thetaPar->getParameter()->getChainValues();
    
    assert(std::equal(posteriorTheta.begin(), posteriorTheta.end(), vec.begin()));
    //    int indexLagTheta=0;
    //    int indexLagDelta=0;
    //    int indexLagProportion=0;
    //    int indexLagTimeOrigin=0;
    
    std::cout << "\n computing auto correlation of  "<< vec.size() << " values of theta";
    //autocorrPosteriorTheta=autoCorrelation(lags, posteriorTheta, mcmcOptions.thresholdAutoCorrelation,  indexLagTheta, thetaPar->getParameter()->getCurrentMean(), thetaPar->getParameter()->getCurrentVariance());
    
    vector<long double> result;
    
    if ( std::adjacent_find( vec.begin(), vec.end(), std::not_equal_to<>() ) == vec.end() )//all values equal
    {
        autocorrPosteriorTheta = 1.0;
    }
    else{
        mcmc_utils::autocorrelationVector(vec,result);
        autocorrPosteriorTheta = getLagWithCorrelationCloseToZero(result,  mcmcOptions.thresholdAutoCorrelation);
    }
    
    if (mcmcOptions.verbose>=2)
        std::cout << "\n auto correlation for theta "<< autocorrPosteriorTheta << std::endl ;
    
    //    int maxIndexLagDelta=lags.size();
    //    int maxIndexLagTimeOrigin=lags.size();
    //    int maxIndexLagProportion=lags.size();
    for (size_t i = 0; i < numClones; ++i)
    {
        auto popI =  populations[i];
        std::cout << "\n computing auto correlation of  "<< popI->posteriorTimeOriginSTD.size() << " values of time origin STD";
        vec = popI->TimeOriginSTD->getParameter()->getChainValues();
        assert(std::equal(popI->posteriorTimeOriginSTD.begin(), popI->posteriorTimeOriginSTD.end(), vec.begin()));
        
        // autocorrPosteriorT = autoCorrelation(lags, popI->posteriorTimeOriginSTD, mcmcOptions.thresholdAutoCorrelation,indexLagTimeOrigin, popI->TimeOriginSTD->getParameter()->getCurrentMean(), popI->TimeOriginSTD->getParameter()->getCurrentVariance());
        
        result.clear();
        if ( std::adjacent_find( vec.begin(), vec.end(), std::not_equal_to<>() ) == vec.end() )//all values equal
        {
            autocorrPosteriorT = 1.0;
        }
        else{
            mcmc_utils::autocorrelationVector(vec,result);
            autocorrPosteriorT = getLagWithCorrelationCloseToZero(result,  mcmcOptions.thresholdAutoCorrelation);
        }
        
        vec = popI->DeltaT->getParameter()->getChainValues();
        assert(std::equal(popI->posteriorDeltaT.begin(), popI->posteriorDeltaT.end(), vec.begin()));
        std::cout << "\n computing auto correlation of  "<< popI->posteriorDeltaT.size() << " values of deltaT";
        
        //autocorrPosteriorDeltaT= autoCorrelation(lags, popI->posteriorDeltaT, mcmcOptions.thresholdAutoCorrelation, indexLagDelta, popI->DeltaT->getParameter()->getCurrentMean(), popI->DeltaT->getParameter()->getCurrentVariance()  );
        
        result.clear();
        if ( std::adjacent_find( vec.begin(), vec.end(), std::not_equal_to<>() ) == vec.end() )//all values equal
        {
            autocorrPosteriorDeltaT = 1.0;
        }
        else{
            mcmc_utils::autocorrelationVector(vec,result);
            autocorrPosteriorDeltaT = getLagWithCorrelationCloseToZero(result,  mcmcOptions.thresholdAutoCorrelation);
        }
        result.clear();
        
        if (numClones>1){
            result.clear();
            mcmc_utils::autocorrelationVector(popI->posteriorProportion,result);
            autocorrPosteriorProportion= getLagWithCorrelationCloseToZero(result,  mcmcOptions.thresholdAutoCorrelation);
            //           autocorrPosteriorProportion=autoCorrelation(lags, popI->posteriorProportion, mcmcOptions.thresholdAutoCorrelation, indexLagProportion, 0.0, 0.0);
            //            if (indexLagProportion < maxIndexLagProportion)
            //                maxIndexLagProportion = indexLagProportion;
        }
        //popI->resetPosteriorValues();
        
        //        if (indexLagDelta < maxIndexLagDelta)
        //            maxIndexLagDelta = indexLagDelta;
        //        if (indexLagTimeOrigin < maxIndexLagTimeOrigin)
        //            maxIndexLagTimeOrigin = indexLagTimeOrigin;
    }
    if (mcmcOptions.verbose>=2){
        //std::cout << "\n auto correlation for DeltaT  and lag  " << lags.at(maxIndexLagDelta);
        std::cout << "\n Chain " << chainNumber << "  auto correlation for DeltaT  and lag  " << autocorrPosteriorDeltaT << std::endl;
    }
    if (mcmcOptions.verbose>=2){
        //std::cout << "\n auto correlation for TOriginSTD  and lag " << lags.at(maxIndexLagTimeOrigin);
        std::cout << "\n Chain " << chainNumber << " auto correlation for TOriginSTD  and lag  " << autocorrPosteriorT;
    }
    if (mcmcOptions.verbose>=2){
        std::cout << "\n Chain  " << chainNumber << " auto correlation for theta  and lag  " << autocorrPosteriorTheta << std::endl ;
    }
    
    if (numClones >1){
        // k = max({indexLagTheta,maxIndexLagTimeOrigin, maxIndexLagDelta, maxIndexLagProportion} );
        k = max({autocorrPosteriorTheta,autocorrPosteriorT, autocorrPosteriorDeltaT, autocorrPosteriorProportion} );
    }
    else{
        // k = max({indexLagTheta,maxIndexLagTimeOrigin, maxIndexLagDelta} );
        k = max({autocorrPosteriorTheta,autocorrPosteriorT, autocorrPosteriorDeltaT} );
    }
    
    // thinning =lags.at(k);
    thinning= k;
    if (thinning == -1)
    {
        std::cout << "\n Chain " << chainNumber << ", not MCMC moves accepted! " << std::endl;
        // here we have to do something to re initialize the parameter values for the current group
        // of chains
    }
    std::cout << "\n Chain " << chainNumber << ", thinning " << thinning << std::endl;
}
void Chain::resetPosteriorValues()
{
    posteriorTheta.clear();
    thetaPar->resetParameterValues();
    for (size_t i = 0; i < numClones; ++i)
    {
        auto popI =  populations[i];
        popI->resetPosteriorValues();
    }
}
int Chain::numParams(){
    return indexMCMCParameters.size();
}
int Chain::index(const std::string& name) const {
    int index = -1;
    for (int i = 0; i < paramNames.size(); i++)
        if (paramNames[i] == name)
            return i;
    return index;
}
std::vector<std::vector<long double>> Chain::getSamples(){
    
    std::vector<std::vector<long double>> result;
    result.push_back(thetaPar->getParameter()->getChainValues());
    
    //add sampled entries of vector x for multiple clones
    for (size_t i = 0; i < numClones; ++i)
    {
        auto popI =  populations[i];
        result.push_back(popI->DeltaT->getParameter()->getChainValues());
        result.push_back(popI->TimeOriginSTD->getParameter()->getChainValues());
        result.push_back(popI->X->getChainValues());
    }
    return(result);
}
//Eigen::MatrixXd Chain::getStoredSamples() const{
//
//    return(stored);
//
//}
std::vector<long double> Chain::getSampleByName(const std::string& name){
    
    std::vector<long double> result;
    MCMCParameterWithKernel *param;
    //get the mcmcParameter by name
    int idx=-1;
    for (size_t i = 0; i < indexMCMCParameters.size(); ++i){
        param=indexMCMCParameters[i].first;
        if (param->getParameter()->getName() == name){
            idx=i;
            break;
        }
    }
    if (idx>0)
        result = param->getParameter()->getChainValues();
    return(result);
}
int Chain::getIndexMCMCParameter(MCMCParameterWithKernel *parameter){
    
    std::vector<long double> result;
    MCMCParameterWithKernel *param;
    //get the mcmcParameter by name
    int idx=-1;
    for (size_t i = 0; i < indexMCMCParameters.size(); ++i){
        param=indexMCMCParameters[i].first;
        if (param->getParameter()->getName() == parameter->getParameter()->getName()){
            idx=i;
            break;
        }
    }
    return idx;
}
std::vector<long double> Chain::getSampleByIndex(const int idx) const{
    std::vector<long double> result;
    MCMCParameterWithKernel *param;
    assert(idx >=0 && idx < indexMCMCParameters.size());
    param=indexMCMCParameters[idx].first;
    result = param->getParameter()->getChainValues();
    return(result);
}
int Chain::getSampleSizeByIndex(const int idx) const{
    int result;
    MCMCParameterWithKernel *param;
    assert(idx >=0 && idx < indexMCMCParameters.size());
    param=indexMCMCParameters[idx].first;
    result = param->getParameter()->getChainValues().size();
    return(result);
}
std::string Chain::getParameterNameByIndex(const int idx) const{
    MCMCParameterWithKernel *param;
    assert(idx >=0 && idx < indexMCMCParameters.size());
    param=indexMCMCParameters[idx].first;
    return(param->getParameter()->getName());
}
