//
//  utils.cpp
//  tumor_coal
//
//  Created by Fausto Fabian Crespo Fernandez on 5/10/19.
//

#include "utils.hpp"
#include "data_utils.hpp"
#include "definitions.h"
#include "kseq.h"
#include <zlib.h>
#include <unistd.h>
#include "Population.hpp"
#include "Chain.hpp"
#include "libpll/pll_optimize.h"
#include "libpll/pll_tree.h"
#include "libpll/pllmod_algorithm.h"
#include "libpll/pll_msa.h"
#include "libpll/pllmod_common.h"


#include <stdarg.h>
#include <search.h>
#include <time.h>

#define GT_MODEL "GTGTR4"
//#define GT_MODEL "GTJC"
#define RATE_CATS 1
#define BRLEN_MIN 1e-6
#define BRLEN_MAX 1e+2
KSEQ_INIT(int, read);

/************************ ReadParametersFromFastaFile ***********************/
/*  ReadParametersFromFastaFile */
void ReadParametersFromFastaFile(char *fileName, ProgramOptions *programOptions){
    //read fasta
    FILE *fastaFile;
    kseq_t *seq;
    int l1;
    int current=0;
    char *currentSeq;
    char *currentQual;
    int  seqlength=0;
    char *t;
    int index,i;
    int max_length=0.0;
    int numberSeq;
  
    
    if ((fastaFile = freopen(fileName, "r", stdin)) != NULL)
    {
        seq = kseq_init(fileno(fastaFile));
        numberSeq=0;
        while ((l1 = kseq_read(seq)) >= 0 )
        {
            //printf("name: %s\n", seq->name.s);
            //if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
            //printf("seq: %s\n", seq->seq.s);
            numberSeq=numberSeq+1;
            if ( l1 > max_length)
            {
                max_length=l1;
            }
        }
        //*numSites=max_length;
        programOptions->numSites=max_length;
        if (numberSeq >=1){
            // *numCells =numberSeq;//not counting the healthy cell
            programOptions->numCells=numberSeq;
        }
        else{
            //   *numCells =0;
            programOptions->numCells=numberSeq;
            programOptions->TotalNumSequences=numberSeq;
        }
        kseq_destroy(seq);
        fclose(fastaFile);
    }
    else{
        
        fprintf (stderr, "\nERROR: Can't read parameters file.");
        PrintUsage();
        
    }
}
/************************ ReadFastaFile ***********************/
/*  ReadFastaFile */
void ReadFastaFile(char *fileName, int** ObservedData,  char **ObservedCellNames, ProgramOptions *programOptions){
    FILE *fastaFile;
    kseq_t *seq;
    int l1;
    int current=0;
    char *currentSeq;
    char *currentQual;
    int  seqlength=0;
    char *t;
    int index,i;
    int max_length=0.0;
    int numberSeq;
    
    if ((fastaFile = freopen(fileName, "r", stdin)) != NULL){
        seq = kseq_init(fileno(fastaFile));
        while ((l1 = kseq_read(seq)) >= 0 ) {
            //printf("name: %s\n", seq->name.s);
            //strcpy( cellNames[current] , seq->name.s);
            ObservedCellNames[current] = (char*) malloc(MAX_NAME);
            if (ObservedCellNames[current] != NULL)
            {
                strcpy(ObservedCellNames[current], seq->name.s);
                //memcpy(ObservedCellNames[current], seq->name.s, sizeof(seq->name.s));
                
            }
            
            //if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
            //printf("seq: %s\n", seq->seq.s);
            seqlength=0;
            //if(current < *numSites){
            currentSeq=seq->seq.s;
            // for ( t= currentSeq; *t != '\0'; t++) {
            for ( index= 0; index < l1; index++) {
                t= currentSeq+index;
                if (programOptions->doUseGenotypes == NO) // use sequences
                    ObservedData[current][index]= WhichNucChar(*t);
                else // use genotypypes
                    ObservedData[current][index]= WhichGenotypeChar(*t);
                
                seqlength++;
                //   }
                //seqlength= sizeof( seq->seq.s) / sizeof(char);
                //dataFromFile[current++]= WhichNucChar(*temp);
            }
            current++;
            if (seq->qual.l){
                // printf("qual: %s\n", seq->qual.s);
                currentQual=seq->qual.s;
            }
        }
        //printf("return value: %d\n", l1);
        kseq_destroy(seq);
        //  gzclose(fastaFile);
        fclose(fastaFile);
    }
}
/************************ LogConditionalLikelihoodTree ***********************/
/*  LogConditionalLikelihoodTree */
double LogConditionalLikelihoodTree(pll_unode_t  *tree, pll_unode_t *nodes, population **populations, int numClones)
{
    population* popI;
    population* popJ;
    population* fatherPop;
    double product=0;
    int i, j;
    double temp;
    for ( i = 0; i < numClones; i++)
    {
        popI=*(populations + i );
        product = product + log( DensityTime(popI->delta, popI->timeOriginSTD));
    }
    for ( j = 0; j < numClones - 1; j++)
    {
        popJ = *(populations + j );
        product = product + log( popJ->popSize);
        fatherPop = popJ -> FatherPop;
        temp=popJ->timeOriginSTD * popJ->effectPopSize / fatherPop->effectPopSize;
        temp=population::CalculateH(popJ->timeOriginSTD * popJ->effectPopSize / fatherPop->effectPopSize, fatherPop->timeOriginSTD, fatherPop->delta);
        product = product  + log( temp);
    }
    //for ( i = 0; i < numClones; i++)
    //  {
    //     popI=*(populations + i );
    product = product + LogDensityCoalescentTimesForPopulation(tree, nodes, populations, numClones);
    //}
    return product;
}
/************************ LogDensityCoalescentTimesForPopulation ***********************/
/*  LogDensityCoalescentTimesForPopulation */
double LogDensityCoalescentTimesForPopulation(pll_unode_t  *tree, pll_unode_t *nodes,  population **populations, int numClones)
{
    double result =0;
    int i, k;
    population *popI;
    int numberLeftCoalescences;
    int numberLeftMigrations;
    int numberAliveCells;
    int currentCoalescentEvent=0;
    int currentMigrationEvent=0;
    double temp;
    for ( i = 0; i < numClones; i++){
        popI = *(populations + i );
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
            while(popI->CoalescentEventTimes[currentCoalescentEvent] < popI->migrationTimes[currentMigrationEvent] && numberAliveCells > 1)
            {
                temp=log(numberAliveCells * (numberAliveCells-1)/2);
                result= result + temp;
                temp = log(1/population::CalculateH(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta));
                result= result + temp;
                temp =  (numberAliveCells/2)* (numberAliveCells-1)*(population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta)-population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent+1], popI->timeOriginSTD, popI->delta));
                result= result -temp;
                currentCoalescentEvent++;
                numberLeftCoalescences--;
                numberAliveCells--;
            }
            //if (numberLeftMigrations > 0 && numberAliveCells > 1)// if there are migrations
            if (numberLeftMigrations > 0 )// if there are migrations
            {  temp= LogProbNoCoalescentEventBetweenTimes(popI,popI->CoalescentEventTimes[currentCoalescentEvent],popI-> migrationTimes[currentMigrationEvent], numberAliveCells );
                result= result+ temp;
                numberLeftMigrations--;
                currentMigrationEvent++;
            }
        }
        //here there are only coalescents events left(al least one event)
        while(numberLeftCoalescences > 0 && numberAliveCells > 1)
        {   temp = log(numberAliveCells * (numberAliveCells-1)/2);
            result= result + temp;
            temp = log(1/population::CalculateH(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta));
            result= result + temp;
            temp=( numberAliveCells/2)* (numberAliveCells-1)*(population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent],popI->timeOriginSTD, popI->delta)-population::FmodelTstandard(popI->CoalescentEventTimes[currentCoalescentEvent+1], popI->timeOriginSTD, popI->delta));
            result= result -  temp;
            currentCoalescentEvent++;
            numberLeftCoalescences--;
            numberAliveCells--;
        }
    }
    
    return result;
}
double DensityTime(double delta, double u){
    double term1=delta * exp(-1*delta*u);
    double term2=1-exp(-1*delta*u);
    return delta * term1 * exp(-1*term1/term2) /(term2 * term2);
}
double LogProbNoCoalescentEventBetweenTimes(population *popI,double from, double to, int numberActiveInd)
{   int j=numberActiveInd;
    double result=0.0;

    result=  -1 * j* (j-1)*(population::FmodelTstandard(to,popI->timeOriginSTD, popI->delta)-population::FmodelTstandard(from, popI->timeOriginSTD, popI->delta))/2;
    return result;
}

void InitializeChains( chain **chains,  ProgramOptions *programOptions,  MCMCoptions *mcmcOptions,int * sampleSizes, long int *seed, char* ObservedCellNames[], pll_msa_t * msa, pll_utree_t * initialTree)
{
    int chainNumber;
    chain *currentChain=NULL;
    double totalTreeLength;
    int        numCA, numMIG;
    double      cumNumCA, meanNumCA, cumNumMIG, meanNumMIG, numEventsTot;
    char *newickString2;
    for(chainNumber=0; chainNumber< mcmcOptions->numChains;chainNumber++)
    {
        
        chains[chainNumber] = (chain*)malloc(sizeof(chain));
        
        chains[chainNumber]->chainNumber = chainNumber;
        chains[chainNumber]->currentNumberIerations=0;
        chains[chainNumber]->numClones = programOptions->numClones;
        chains[chainNumber]->mutationRate =programOptions->mutationRate;
        chains[chainNumber]->gammaParam=1;
        chains[chainNumber]->seqErrorRate= programOptions->seqErrorRate;
        chains[chainNumber]->dropoutRate =programOptions->dropoutRate;
        if (programOptions->numberClonesKnown)
        {
            chains[chainNumber]->numNodes = programOptions->numNodes;
            
            
        }
        
        chains[chainNumber]->proportionsVector = (double *) calloc((programOptions->numClones), (long) sizeof(double));
        if (!(chains[chainNumber]->proportionsVector))
        {
            fprintf (stderr, "Could not allocate proportions vector (%lu bytes)\n", (programOptions->numClones ) * (long) sizeof(double));
            exit (-1);
        }
        chains[chainNumber]->oldproportionsVector = (double *) calloc((programOptions->numClones), (long) sizeof(double));
        if (!(chains[chainNumber]->oldproportionsVector))
        {
            fprintf (stderr, "Could not allocate proportions vector (%lu bytes)\n", (programOptions->numClones ) * (long) sizeof(double));
            exit (-1);
        }
        chains[chainNumber]->populations = (population**)malloc (sizeof(struct Population*)  * programOptions->numClones);
        if (!(chains[chainNumber]->populations))
        {
            fprintf (stderr, "Could not allocate populations (%lu bytes)\n", (programOptions->numClones)  * (long) sizeof(Population*));
            exit (1);
        }
        
        
        
        chains[chainNumber]->treeTips =(pll_unode_t**)  malloc (programOptions->TotalNumSequences* sizeof(pll_unode_t*));
        if (!(chains[chainNumber]->treeTips))
        {
            fprintf (stderr, "Could not allocate the treeTips array\n");
            exit (-1);
        }
        chains[chainNumber]->oldtreeTips =(pll_unode_t**) malloc (programOptions->TotalNumSequences* sizeof(pll_unode_t*));
        if (!(chains[chainNumber]->oldtreeTips))
        {
            fprintf (stderr, "Could not allocate the treeTips array\n");
            exit (-1);
        }
        chains[chainNumber]->InitChainPopulations( programOptions->noisy, programOptions->TotalNumSequences );
       chains[chainNumber]->FillChainPopulationsFromPriors( programOptions,mcmcOptions, sampleSizes, seed );
        
        if (programOptions-> doUseFixedTree)
            chains[chainNumber]->initialTree=initialTree;
        else
        {
                                 
            chains[chainNumber]->MakeCoalescenceTree (seed,
                                 &(chains[chainNumber]->numNodes),
                                 chains[chainNumber]->numClones,
                                 programOptions,
                                 cumNumCA,
                                 meanNumCA,
                                 cumNumMIG,
                                 meanNumMIG,
                                 &numMIG,
                                 &numCA,
                                 &numEventsTot,
                                 ObservedCellNames, sampleSizes
                                 ) ;
        }
        //        cumNumCA += numCA;
        //        cumNumMIG += numMIG;
        //        countTMRCA = treeRootInit[0]->timePUnits;
        //
        //        varTimeGMRCA[currentIteration] = countTMRCA;
        
        // totalTreeLength = SumBranches(chains[chainNumber]->root, chains[chainNumber]->mutationRate);
        //        cumNumMUperTree=0;
        
        newickString2=NULL;
//        newickString2 = toNewickString2 ( chains[chainNumber]->root, chains[chainNumber]->mutationRate,     programOptions->doUseObservedCellNames);
//        printf("\n newick = %s  \n", newickString2);
        
        
        chains[chainNumber]->currentlogConditionalLikelihoodTree= LogConditionalLikelihoodTree(chains[chainNumber]->root, chains[chainNumber]->nodes, chains[chainNumber]->populations,  chains[chainNumber]->numClones);
        printf ( "Initial likelihood of the tree of chain %d is:  %lf \n",chainNumber, chains[chainNumber]->currentlogConditionalLikelihoodTree );
        
        chains[chainNumber]->currentlogConditionalLikelihoodSequences= LogConditionalLikelihoodSequences( msa,  newickString2, programOptions,  chains[chainNumber]->seqErrorRate,
                                                                    chains[chainNumber]->dropoutRate);
        
        fprintf (stderr, "Initial likelihood of the sequences of chain %d  is = %lf  \n", chainNumber,chains[chainNumber]->currentlogConditionalLikelihoodSequences );
        
        free(newickString2);
        newickString2=NULL;
        
    }
}
double  LogConditionalLikelihoodSequences(pll_msa_t * msa, char* NewickString, ProgramOptions *programOptions, double seqError,double dropoutError)
{
    
    FILE    *outputShell;
    char script[80];
    char buf[1000];
    
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
    
    pllmod_treeinfo_t * treeinfo = pllmod_treeinfo_create(unrootedTree->vroot,
                                                          tip_nodes_count,
                                                          1,
                                                          PLLMOD_COMMON_BRLEN_LINKED);
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
        fprintf (stderr, "\nERROR: The multiple sequence alignment is empty\n\n");
        PrintUsage();
        return 0;
    }
    
    
    
    pllmod_subst_model_t * model = pllmod_util_model_info_genotype(GT_MODEL);
    
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
    
    set_partition_tips( partition, msa, programOptions,  seqError, dropoutError);
    
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
    
    double * empirical_frequencies = pllmod_msa_empirical_frequencies(partition);
    
    unsigned int * weight = pll_compress_site_patterns(msa->sequence,
                                                       pll_map_gt10,
                                                       tip_nodes_count,
                                                       &(msa->length));
    printf("Number of unique site patterns: %d\n\n", msa->length);
    
    /* initialize the array of base frequencies  AA CC GG TT AC/CA AG/GA AT/TA CG/GC CT/TC GT/TG  */
    double user_freqs[10] = { 0.254504, 0.157238, 0.073667, 0.287452, 0.027336,
        0.058670, 0.041628, 0.024572, 0.064522, 0.010412 };
    
    //computeGenotypesFreq(user_freqs,  msa);
    
    
    /* substitution rates: for GTR4 model those are 6 "regular" DNA susbt. rates + 1 rate
     * for "unlikely" double substitutions (eg A/A -> C/T) */
    double unique_subst_rates[7] = { 0.001000, 0.101223, 0.001000, 0.001000, 1.000000,
        0.001000, 0.447050 };
    
    /* get full above-diagonal half-matrix */
    double * user_subst_rates = expand_uniq_rates(model->states, unique_subst_rates,
                                                  model->rate_sym);
    
    double rate_cats[RATE_CATS] = {0};
    
    /* compute the discretized category rates from a gamma distribution
     with alpha shape 1 and store them in rate_cats  */
    pll_compute_gamma_cats(1, RATE_CATS, rate_cats, PLL_GAMMA_RATES_MEAN);
    
    /* set frequencies at model with index 0 (we currently have only one model) */
    //pll_set_frequencies(partition, 0, model->freqs ? model->freqs : user_freqs);
    pll_set_frequencies(partition, 0, model->freqs ? model->freqs : empirical_frequencies);
    
    /* set substitution parameters at model with index 0 */
    pll_set_subst_params(partition, 0, model->rates ? model->rates : user_subst_rates);
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
void set_partition_tips( pll_partition_t * partition, pll_msa_t * msa, ProgramOptions *programOptions, double seqError, double dropoutError)
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
        
        if (programOptions->doUseGenotypes == NO)
        {
            pll_set_tip_states(partition, tip_clv_index, pll_map_gt10, msa->sequence[i]);
        }
        else
        {
            
            //            for ( currentState = 0; currentState < states; currentState++)
            //            {
            //                from=1;
            //                to=1;
            set_tipclv1( partition,
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
/********************* expand_uniq_rates **********************/
/* expand_uniq_rates(Function that uses code from genotype.c from pll-modules/examples) */
double * expand_uniq_rates(int states, const double * uniq_rates, const int * rate_sym)
{
    unsigned int i;
    unsigned int num_rates = states * (states-1) / 2;
    double * subst_rates =(double *) calloc(num_rates, sizeof(double));
    for (i = 0; i < num_rates; ++i)
        subst_rates[i] = rate_sym ? uniq_rates[rate_sym[i]] : uniq_rates[i];
    return subst_rates;
}
int set_tipclv1(pll_partition_t * partition,
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
void  destroyTree(pll_utree_t * tree, void (*cb_destroy)(void *))
{
    
    unsigned int i;
    
    /* deallocate tip nodes */
    for (i = 0; i < tree->tip_count; ++i)
    {
        dealloc_data(tree->nodes[i], cb_destroy);
        
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
            dealloc_data(node, cb_destroy);
            free(node);
            node = next;
        }
        while(node && node != first);
    }
    
    /* deallocate tree structure */
    free(tree->nodes);
    free(tree);
    
}
void dealloc_data(pll_unode_t * node, void (*cb_destroy)(void *))
{
    if (node->data)
    {
        if (cb_destroy)
            cb_destroy(node->data);
    }
}



double RandomLogUniform( double from, double to, long int *seed){
    
    return(exp(from + RandomUniform(seed)*(to -from)));
}
void  RandomDirichlet (double s, int vectorSize, double **outputVector, long int *seed)
{   int i;
    double sum=0.0;
    double current;
    // *outputVector = malloc(vectorSize * sizeof(double));
    // if (*outputVector == NULL)
    //     return;
    for (i=0; i < vectorSize; i++){
        current = RandomGamma(s, seed);
        (*outputVector)[i] = current;
        //*(*outputVector + i)=current;
        sum=sum+current;
    }
    for (i=0; i < vectorSize; i++){
        (*outputVector)[i] =  (*outputVector)[i] / sum;
        // *(*outputVector + i)= *(*outputVector + i)/sum;
    }
}
void InitPopulationSampleSizes(population **populations, int TotalSampleSize, int numClones, double *proportionsVector, long int *seed)
{
    int i,j;
    population *popI, *popJ;
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
        popJ = *(populations + j - 1);
        popJ->sampleSize=0;
    }
    for (i = 0; i < TotalSampleSize; i++)
    {
        rand=RandomUniform(seed);
        for (j = 1; j <= numClones;j ++)
        {
            popJ = *(populations + j - 1);
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
