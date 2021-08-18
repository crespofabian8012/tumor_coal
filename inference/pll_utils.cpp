//
//  pll_utils.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 15/03/2021.
//

#include "pll_utils.hpp"


#include "constants.hpp"
#include "output_functions.hpp"

//#include "genotype_error_model.cpp"
#include <stdlib.h>


//long double  pll_utils::LogConditionalLikelihoodSequences(pll_msa_t * msa, char* NewickString, ProgramOptions &programOptions, long double seqError,long double dropoutError)
//{
//
//    if (NewickString == NULL)
//    {
//        fprintf (stderr, "\nERROR: The newick representation of the tree cannot be empty\n\n");
//        Output::PrintUsage();
//        return 0;
//    }
//    
//    pll_utree_t *rootedTree = pll_utree_parse_newick_string(NewickString);
//    return LogConditionalLikelihoodSequences( msa, rootedTree,  programOptions, seqError,dropoutError);
//}
long double  pll_utils::LogConditionalLikelihoodSequences(pll_msa_t * msa,
                                                          pll_rtree_t *rootedTree,
                                                          ProgramOptions &programOptions,
                                                          long double seqError,
                                                          long double dropoutError){
    
    char*  rootedNewick = pll_rtree_export_newick( rootedTree->root, NULL);
    return LogConditionalLikelihoodSequences(msa, rootedNewick, programOptions, seqError, dropoutError);
    
}
long double  pll_utils::LogConditionalLikelihoodSequencesRootedTree(pll_msa_t * msa,
                                                                    pll_rtree_t *rootedTree,
                                                                    ProgramOptions &programOptions,
                                                                    long double seqError,
                                                                    long double dropoutError){
    
    //char*  rootedNewick = pll_rtree_export_newick( rootedTree->root, NULL);
    unsigned int i;
    unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
    pll_partition_t * partition;
    
    pll_utree_t * unrootedTree = pll_rtree_unroot(rootedTree);
    
    
    /* compute node count information */
    tip_nodes_count = unrootedTree->tip_count;
    inner_nodes_count = unrootedTree->inner_count;
    nodes_count = inner_nodes_count + tip_nodes_count;
    branch_count = unrootedTree->edge_count;
    
    
    
    pll_utree_reset_template_indices((pll_unode_s*)unrootedTree->vroot, tip_nodes_count);
    
    pllmod_treeinfo_t * treeinfo = pllmod_treeinfo_create((pll_unode_s*)unrootedTree->vroot,
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
#ifdef __APPLE__
        entry.key = xstrdup(tipnodes[i]->label);
#else
        entry.key = tipnodes[i]->label;
#endif
        entry.key = tipnodes[i]->label;
        entry.data = (void *)(data+i);
        hsearch(entry, ENTER);
    }
    
    if (msa ==NULL)
    {
        std::cout << "\nERROR: The multiple sequence alignment is empty\n\n"<< std::endl;
        
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
    
    unsigned int attributes = PLL_ATTRIB_ARCH_AVX | PLL_ATTRIB_AB_LEWIS;
    //PLL_ATTRIB_ASC_BIAS_LEWIS or PLL_ATTRIB_ASC_BIAS_FELSENSTEIN or PLL_ATRIB_ASC_BIAS_STAMATAKIS
    partition = pll_partition_create(tip_nodes_count,
                                     inner_nodes_count,
                                     model->states,
                                     (unsigned int)(msa->length),
                                     1,
                                     branch_count,
                                     RATE_CATS,
                                     inner_nodes_count,
                                     attributes);
    
    //pll_set_asc_bias_type(partition, PLL_ATTRIB_AB_LEWIS);
    
    //static unsigned int invar_weights[STATES] = { 50, 40, 60, 20 };
    //pll_set_asc_state_weights(partition, invar_weights);
    
    //setPartitionTipsCostum( partition, msa, programOptions,  seqError, dropoutError);
    
    
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
    
    
    setPartitionTipsCostum( partition, msa, programOptions,  seqError, dropoutError);
    double  * empirical_frequencies;
    empirical_frequencies = pllmod_msa_empirical_frequencies(partition);
    
    double user_freqs[10] = { 0.254504, 0.157238, 0.073667, 0.287452, 0.027336,
        0.058670, 0.041628, 0.024572, 0.064522, 0.010412 };
    
    if (model->freqs)
        pll_set_frequencies(partition, 0,  model->freqs);
    else{
        pll_set_frequencies(partition, 0,  empirical_frequencies);
        // pll_set_frequencies(partition, 0,  user_freqs);
    }
    
    unsigned int * weight = pll_compress_site_patterns(msa->sequence,
                                                       pll_map_gt10,
                                                       tip_nodes_count,
                                                       &(msa->length));
    
    
    // double  * empirical_subst_rates = pllmod_msa_empirical_subst_rates( partition);
    
    /* initialize the array of base frequencies  AA CC GG TT AC/CA AG/GA AT/TA CG/GC CT/TC GT/TG  */
    //        double user_freqs[10] = { 0.254504, 0.157238, 0.073667, 0.287452, 0.027336,
    //            0.058670, 0.041628, 0.024572, 0.064522, 0.010412 };
    
    
    /* substitution rates: for GTR4 model those are 6 "regular" DNA susbt. rates + 1 rate
     * for "unlikely" double substitutions (eg A/A -> C/T) */
    double  unique_subst_rates[7] = { 0.001000, 0.101223, 0.001000, 0.001000, 1.000000,
        0.001000, 0.447050 };
    
    /* get full above-diagonal half-matrix */
    double  * user_subst_rates = expandUniqRates(model->states, unique_subst_rates,
                                                 model->rate_sym);
    
    double  rate_cats[RATE_CATS] = {0};
    
    /* compute the discretized category rates from a gamma distribution
     with alpha shape 1 and store them in rate_cats  */
    pll_compute_gamma_cats(1, RATE_CATS, rate_cats, PLL_GAMMA_RATES_MEAN);
    
    if (model->freqs)
        pll_set_frequencies(partition, 0,  model->freqs);
    else{
        //pll_set_frequencies(partition, 0,  empirical_frequencies);
        pll_set_frequencies(partition, 0,  user_freqs);
    }
    /* set substitution parameters at model with index 0 */
    //pll_set_subst_params(partition, 0, model->rates ? model->rates : user_subst_rates);
    //pll_set_subst_params(partition, 0, model->rates ? model->rates : empirical_subst_rates);
    
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
    
    pll_partition_destroy(partition);
    //pll_utree_destroy(unrootedTree, NULL);
    destroyUTree(unrootedTree, NULL);
    pllmod_util_model_destroy(model);
    return loglh;
    
}
long double  pll_utils::LogConditionalLikelihoodSequences(pll_msa_t * msa, char* NewickString, ProgramOptions &programOptions, long double seqError,long double dropoutError )
{
    unsigned int i;
    unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
    pll_partition_t * partition;
    
    pll_utree_t *unrootedTree = pll_utree_parse_newick_string_unroot(NewickString);
    //pll_utree_t * unrootedTree = pll_rtree_unroot(rootedTree);
    
    /* compute node count information */
    tip_nodes_count = unrootedTree->tip_count;
    inner_nodes_count = unrootedTree->inner_count;
    nodes_count = inner_nodes_count + tip_nodes_count;
    branch_count = unrootedTree->edge_count;
    
    
    pllmod_treeinfo_t * treeinfo = pllmod_treeinfo_create((pll_unode_s*)unrootedTree->vroot,
                                                          tip_nodes_count,
                                                          1, PLLMOD_COMMON_BRLEN_LINKED);
    
    
    pll_unode_t ** tipnodes = unrootedTree->nodes;
    
    /* create a libc hash table of size tip_nodes_count */
    hcreate(tip_nodes_count);
    
    /* populate a libc hash table with tree tip labels */
    unsigned int * data = (unsigned int *)malloc(tip_nodes_count *
                                                 sizeof(unsigned int));
    //char *label;
    for (i = 0; i < tip_nodes_count; ++i)
    {
        data[i] = tipnodes[i]->clv_index;
        ENTRY entry;
        //label= tipnodes[i]->label;
#ifdef __APPLE__
        entry.key = xstrdup(tipnodes[i]->label);
#else
        entry.key = tipnodes[i]->label;
#endif
        entry.key = tipnodes[i]->label;
        entry.data = (void *)(data+i);
        hsearch(entry, ENTER);
    }
    
    if (msa ==NULL)
    {
        std::cout << "\nERROR: The multiple sequence alignment is empty\n\n"<< std::endl;
        
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
    
    unsigned int attributes = PLL_ATTRIB_ARCH_AVX | PLL_ATTRIB_AB_LEWIS;
    //PLL_ATTRIB_ASC_BIAS_LEWIS or PLL_ATTRIB_ASC_BIAS_FELSENSTEIN or PLL_ATRIB_ASC_BIAS_STAMATAKIS
    partition = pll_partition_create(tip_nodes_count,
                                     inner_nodes_count,
                                     model->states,
                                     (unsigned int)(msa->length),
                                     1,
                                     branch_count,
                                     RATE_CATS,
                                     inner_nodes_count,
                                     attributes);
    
    //pll_set_asc_bias_type(partition, PLL_ATTRIB_AB_LEWIS);
    
    //static unsigned int invar_weights[STATES] = { 50, 40, 60, 20 };
    //pll_set_asc_state_weights(partition, invar_weights);
    
    //setPartitionTipsCostum( partition, msa, programOptions,  seqError, dropoutError);
    
    
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
    
    
    setPartitionTipsCostum( partition, msa, programOptions,  seqError, dropoutError);
    double  * empirical_frequencies;
    empirical_frequencies = pllmod_msa_empirical_frequencies(partition);
    
    double user_freqs[10] = { 0.254504, 0.157238, 0.073667, 0.287452, 0.027336,
        0.058670, 0.041628, 0.024572, 0.064522, 0.010412 };
    
    if (model->freqs)
        pll_set_frequencies(partition, 0,  model->freqs);
    else{
        pll_set_frequencies(partition, 0,  empirical_frequencies);
        // pll_set_frequencies(partition, 0,  user_freqs);
    }
    
    unsigned int * weight = pll_compress_site_patterns(msa->sequence,
                                                       pll_map_gt10,
                                                       tip_nodes_count,
                                                       &(msa->length));
    
    
    // double  * empirical_subst_rates = pllmod_msa_empirical_subst_rates( partition);
    
    /* initialize the array of base frequencies  AA CC GG TT AC/CA AG/GA AT/TA CG/GC CT/TC GT/TG  */
    //        double user_freqs[10] = { 0.254504, 0.157238, 0.073667, 0.287452, 0.027336,
    //            0.058670, 0.041628, 0.024572, 0.064522, 0.010412 };
    
    
    /* substitution rates: for GTR4 model those are 6 "regular" DNA susbt. rates + 1 rate
     * for "unlikely" double substitutions (eg A/A -> C/T) */
    double  unique_subst_rates[7] = { 0.001000, 0.101223, 0.001000, 0.001000, 1.000000,
        0.001000, 0.447050 };
    
    /* get full above-diagonal half-matrix */
    double  * user_subst_rates = expandUniqRates(model->states, unique_subst_rates,
                                                 model->rate_sym);
    
    double  rate_cats[RATE_CATS] = {0};
    
    /* compute the discretized category rates from a gamma distribution
     with alpha shape 1 and store them in rate_cats  */
    pll_compute_gamma_cats(1, RATE_CATS, rate_cats, PLL_GAMMA_RATES_MEAN);
    
    if (model->freqs)
        pll_set_frequencies(partition, 0,  model->freqs);
    else{
        //pll_set_frequencies(partition, 0,  empirical_frequencies);
        pll_set_frequencies(partition, 0,  user_freqs);
    }
    /* set substitution parameters at model with index 0 */
    //pll_set_subst_params(partition, 0, model->rates ? model->rates : user_subst_rates);
    //pll_set_subst_params(partition, 0, model->rates ? model->rates : empirical_subst_rates);
    
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
    
    pll_partition_destroy(partition);
    //pll_utree_destroy(unrootedTree, NULL);
    destroyUTree(unrootedTree, NULL);
    pllmod_util_model_destroy(model);
    return loglh;
}
void pll_utils::setPartitionTipsCostum( pll_partition_t * partition, pll_msa_t * msa, ProgramOptions &programOptions, long double  seqError, long double  dropoutError)
{
    
    int i;
    //GenotypeErrorModel gtErrorModel= GenotypeErrorModel("PT19", seqError, dropoutError);
    unsigned int tip_clv_index;
    
    /* find sequences in hash table and link them with the corresponding taxa */
    for (i = 0; i < msa->count; ++i)
    {
        ENTRY query;
        query.key = msa->label[i];
        ENTRY * found = NULL;
        
        found = hsearch(query,FIND);
        
        if (!found)
            fprintf(stderr,"Sequence with header %s does not appear in the tree", msa->label[i]);
        
        tip_clv_index = *((unsigned int *)(found->data));
        //std::cout << "clv index for tip with label " <<msa->label[i] << " is "<< tip_clv_index << std::endl;
        
        if (programOptions.doUseGenotypes == NO)
        {
            if (pll_set_tip_states(partition, tip_clv_index, pll_map_gt10, msa->sequence[i]) != PLL_SUCCESS)
            {
                fprintf(stderr, "PLL error %d: %s\n", pll_errno, pll_errmsg);
                exit(pll_errno);
            }
        }
        else
        {
            if (setTipClvCustomErrorModel( partition,
                                          tip_clv_index,
                                          pll_map_gt10,
                                          msa->sequence[i],seqError,dropoutError)!= PLL_SUCCESS)
            {
                fprintf(stderr, "PLL error %d: %s\n", pll_errno, pll_errmsg);
                exit(pll_errno);
            }
        }
    }
}
int pll_utils::setTipClvCustomErrorModel(pll_partition_t * partition,
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
double  * pll_utils::expandUniqRates(int states, const  double  * uniq_rates, const int * rate_sym)
{
    unsigned int i;
    
    unsigned int num_rates = states * (states-1) / 2;
    //long double  * subst_rates = calloc(num_rates, sizeof(double));
    double  * subst_rates =(double  *) calloc(num_rates, sizeof(double));
    for (i = 0; i < num_rates; ++i)
        subst_rates[i] = rate_sym ? uniq_rates[rate_sym[i]] : uniq_rates[i];
    
    return subst_rates;
}
void  pll_utils::destroyUTree(pll_utree_t * tree, void (*cb_destroy)(void *))
{
    
    unsigned int i;
    
    /* deallocate tip nodes */
    for (i = 0; i < tree->tip_count; ++i)
    {
        deallocUDataCostum(tree->nodes[i], cb_destroy);
        
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
            deallocUDataCostum(node, cb_destroy);
            free(node);
            node = next;
        }
        while(node && node != first);
    }
    
    /* deallocate tree structure */
    free(tree->nodes);
    free(tree);
    
}
void pll_utils::destroyRTree(pll_rtree_t * tree,
                             void (*cb_destroy)(void *))
{
    unsigned int i;
    pll_rnode_t * node;
    
    /* deallocate all nodes */
    for (i = 0; i < tree->tip_count + tree->inner_count; ++i)
    {
        node = tree->nodes[i];
        deallocRDataCostum(node, cb_destroy);
        
        if (node->label)
            free(node->label);
        
        free(node);
    }
    /* deallocate tree structure */
    free(tree->nodes);
    free(tree);
}
void pll_utils::deallocUDataCostum(pll_unode_t * node, void (*cb_destroy)(void *))
{
    if (node->data)
    {
        if (cb_destroy)
            cb_destroy(node->data);
    }
}
void pll_utils::deallocRDataCostum(pll_rnode_t * node, void (*cb_destroy)(void *))
{
    if (node->data)
    {
        if (cb_destroy)
            cb_destroy(node->data);
    }
}
pll_rnode_t * pll_utils::cloneRNode(const pll_rnode_t * node)
{
    if (node!=NULL){
        pll_rnode_t * new_node = (pll_rnode_t *)malloc(sizeof(pll_rnode_t));
        memcpy(new_node, node, sizeof(pll_rnode_t));
        
        new_node->clv_index=node->clv_index;
        new_node->length=node->length;
        new_node->pmatrix_index=node->pmatrix_index;
        new_node->scaler_index = node->scaler_index;
        new_node->clv_index=node->clv_index;
        
        if (node->label)
        {
            new_node->label = (char *)malloc(strlen(node->label)+1);
            strcpy(new_node->label,node->label);
        }
        if (node->data){
            TreeNode *currentTreeNode = (TreeNode*)(node->data);
            //TreeNode * newTreeNode= new TreeNode(currentTreeNode->genotypeSequence.size());
            TreeNode * newTreeNode= new TreeNode(0);//no need to store IUPAC code of the sequences
            
            std::copy(currentTreeNode->observedCellName, currentTreeNode->observedCellName + strlen(currentTreeNode->observedCellName), newTreeNode->observedCellName);
            //or for char* do newTreeNode->observedCellName= strdup(currentTreeNode->observedCellName)
            newTreeNode->indexSequenceMSA= currentTreeNode->indexSequenceMSA;
            newTreeNode->time = currentTreeNode->time;//model time
            newTreeNode->timePUnits = currentTreeNode->timePUnits;// model time * theta * x
            //newTreeNode->genotypeSequence = currentTreeNode->genotypeSequence;
            new_node->data = newTreeNode;
        }
        
        new_node->left = cloneRNode(node->left);
        new_node->right = cloneRNode(node->right);
        return new_node;
    }
    else
        return NULL;
}
pll_rtree_t * pll_utils::cloneRTree(const pll_rtree_t * tree){
    
    if (tree ==NULL){
        return NULL;
    }
    else{
        pll_rtree_t * new_tree = (pll_rtree_t *)malloc(sizeof(pll_rtree_t));
        
        if (tree->root){
            new_tree->root = cloneRNode(tree->root);
        }
        return new_tree;
    }
}
char * pll_utils::xstrdup(const char * s)
{
    size_t len = strlen(s);
    //char * p = (char *)xmalloc(len+1);
    char * p = (char *)malloc(len+1);
    return strcpy(p,s);
}
const pll_partition_t *pll_utils::createReferencePartition(pll_msa_t * msa) {
    
    
    const unsigned int rate_category_count = 4;
    double rate_categories[4] = {0, 0, 0, 0};
    pll_compute_gamma_cats(1, 4, rate_categories, PLL_GAMMA_RATES_MEAN);
    
    const unsigned int subst_model_count = 1;
    const double subst_params[6] = {1, 1, 1, 1, 1, 1};
    
    const unsigned int nucleotide_states = 4;
    const double nucleotide_frequencies[4] = {0.25, 0.25, 0.25, 0.25};
    
    pll_partition *partition = pll_partition_create(
                                                    msa->count,
                                                    0, // Don't allocate any inner CLV's.
                                                    nucleotide_states, msa->length, subst_model_count,
                                                    0, // Don't allocate any pmatrices.
                                                    rate_category_count,
                                                    0, // Don't allocate any scale buffers.
                                                    PLL_ATTRIB_ARCH_SSE);
    
    assert(partition);
    pll_set_frequencies(partition, 0, nucleotide_frequencies);
    pll_set_category_rates(partition, rate_categories);
    pll_set_subst_params(partition, 0, subst_params);
    
    for (unsigned int i = 0; i < msa->count; i++) {
        
        if (pll_set_tip_states(partition, i, pll_map_nt, msa->sequence[i]) != PLL_SUCCESS)
        {
            fprintf(stderr, "PLL error %d: %s\n", pll_errno, pll_errmsg);
            exit(pll_errno);
        }
    }
    
    pll_update_eigen(partition, 0);
    
    return partition;
}
const pll_partition_t *pll_utils::createGTReferencePartition(pll_msa_t * msa) {
    
    
    const unsigned int rate_category_count = 1;
    unsigned int number_states=16;
    // double rate_categories[4] = {0, 0, 0, 0};
    //pll_compute_gamma_cats(1, 4, rate_categories, PLL_GAMMA_RATES_MEAN);
    
    const unsigned int subst_model_count = 1;
    //const double subst_params[6] = {1, 1, 1, 1, 1, 1};
    
    const unsigned int nucleotide_states = 16;
    
    
    pll_partition *partition = pll_partition_create(
                                                    msa->count,
                                                    0, // Don't allocate any inner CLV's.
                                                    nucleotide_states, msa->length, subst_model_count,
                                                    0, // Don't allocate any pmatrices.
                                                    rate_category_count,
                                                    0, // Don't allocate any scale buffers.
                                                    PLL_ATTRIB_ARCH_SSE);
    
    assert(partition);
    
    double subst_params[ number_states*(number_states-1) /2] ;
    std::fill_n (subst_params,number_states*(number_states-1) /2, 1);
    
    double  rate_categories[number_states];
    std::fill_n (rate_categories, number_states*(number_states-1) /2, 1);
    
    double  category_weights[number_states];
    std::fill_n (category_weights, number_states*(number_states-1) /2, 1);
    
    double nucleotide_frequencies[number_states];
    std::fill_n (nucleotide_frequencies, number_states, 1.0/(number_states));
    
    pll_compute_gamma_cats(1, RATE_CATS, rate_categories, PLL_GAMMA_RATES_MEAN);
    
    pll_set_frequencies(partition, 0, nucleotide_frequencies);
    pll_set_category_rates(partition, rate_categories);
    pll_set_subst_params(partition, 0, subst_params);
    
    for (unsigned int i = 0; i < msa->count; i++) {
        
        if (pll_set_tip_states(partition, i, pll_map_gt10, msa->sequence[i]) != PLL_SUCCESS)
        {
            fprintf(stderr, "PLL error %d: %s\n", pll_errno, pll_errmsg);
            exit(pll_errno);
        }
    }
    
    pll_update_eigen(partition, 0);
    
    return partition;
}

double pll_utils::computeLogLikelihood(double *clv, unsigned int *scale_buffer,
                                       const pll_partition_t *p) {
    const unsigned int parameter_indices[4] = {0};
    
    
    return pll_core_root_loglikelihood(
                                       p->states, p->sites, p->rate_cats,
                                       clv, scale_buffer,
                                       p->frequencies, p->rate_weights, p->pattern_weights, p->prop_invar,
                                       p->invariant, parameter_indices, NULL, p->attributes);
}
int pll_utils::cb_rfull_traversal(pll_rnode_t * node)
{
    return 1;
}
double pll_utils::meanHeterozygosity(long seq_count, long site_count, pll_msa_t* msa)
{
    /* This calculates the mean heterozygosity for sequences in the alignment,
     assuming that all sequences are from the same species and each sequence is
     diploid unphased.
     The code needs to be changed if there are 2 or more species.
     */
    long j,h;
    double H = 0;
    
    assert(seq_count > 0 && site_count > 0);
    for (j = 0; j < seq_count; j++)
    {
        for (h = 0; h < site_count; h++)
            if (msa->sequence[j][h] != 1 &&
                msa->sequence[j][h] != 2 &&
                msa->sequence[j][h] != 4 &&
                msa->sequence[j][h] != 8)
                H++;
    }
    H /= seq_count * site_count;
    return(H);
}
/*** Ziheng 2020-10-17
 * Calculation of average sequence distances (p-distance) over loci and sequnce pairs
 * for simulated data.  Sequences have to be fulled resolved/phased.
 ***/
double pll_utils::meanDistance(long seq_count, long site_count, pll_msa_t* msa)
{
    long j1,j2,h;
    double md = 0;
    
    assert(seq_count > 0 && site_count > 0);
    for (j1 = 0; j1 < seq_count; j1++)
    {
        for (j2 = 0; j2 < j1; j2++)
        {
            for (h = 0; h < site_count; h++)
                if (msa->sequence[j1][h] != msa->sequence[j2][h]) md++;
        }
    }
    md /= seq_count * (seq_count - 1.0) / 2.0 * site_count;
    return(md);
}
