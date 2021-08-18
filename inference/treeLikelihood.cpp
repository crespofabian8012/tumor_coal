//
//  treeLikelihood.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 26/06/2021.
//
#include <numeric>
#include "treeLikelihood.hpp"

#include "pll_utils.hpp"

#include "search.h"

TreeLikelihood::TreeLikelihood(Partition *reference_partition, pll_rtree_t * rooted_tree, pll_msa_t * msa, GenotypeErrorModel *gtErrorModel):
rooted_tree(rooted_tree),reference_partition(reference_partition),msa(msa),gtErrorModel(gtErrorModel)
{
    unsigned int leafCount, internalNodeCount, nodeCount, branchCount, patternCount;
    
    
    leafCount = rooted_tree->tip_count;
    internalNodeCount = rooted_tree->inner_count;
    nodeCount = internalNodeCount + leafCount;
    branchCount = rooted_tree->edge_count;
    patternCount=msa->length;
    
    number_genotype_states= 10;
    number_phased_genotype_states= 16;
    number_states = number_phased_genotype_states;
    
    branch_lengths = (double *)malloc(branchCount * sizeof(double));
    matrix_indices = (unsigned int *)malloc(branchCount * sizeof(unsigned int));
    operations = (pll_operation_t *)malloc(internalNodeCount *
                                           sizeof(pll_operation_t));
    
    recomputeTipClvs();
    
    double subst_params[ number_states*(number_states-1) /2] ;
    std::fill_n (subst_params,number_states*(number_states-1) /2, 1);
    
    double  rate_categories[number_states];
    std::fill_n (rate_categories, number_states*(number_states-1) /2, 1);
    
    double  category_weights[number_states];
    std::fill_n (category_weights, number_states*(number_states-1) /2, 1);
    
    double nucleotide_frequencies[number_states];
    std::fill_n (nucleotide_frequencies, number_states, 1.0/(number_states));
    
    pll_compute_gamma_cats(1, RATE_CATS, rate_categories, PLL_GAMMA_RATES_MEAN);
    reference_partition->setCategoryRates(rate_categories);
    
    reference_partition->setCategoryWeights(category_weights);
    reference_partition->setFrequencies(0, nucleotide_frequencies);
    reference_partition->setSubstParams(0,subst_params);
    reference_partition->updateEigen(0);
    
    initOperations();
};

TreeLikelihood::TreeLikelihood(int numberTips,
                               int clvBuffers,
                               int numberStates,
                               int numberSites,
                               int numberRateMatrices,
                               int probMatrices,
                               int numberRateCats,
                               int numberScaleBuffers,
                               int statesPadded,
                               bool sse,
                               bool avx,
                               bool  avx2,
                               bool avx512,
                               bool  asc,
                               bool tipPatternCompression, pll_rtree_t * rooted_tree, pll_msa_t * msa, GenotypeErrorModel *gtErrorModel ):rooted_tree(rooted_tree),msa(msa), gtErrorModel(gtErrorModel){
    
    
    unsigned int leafCount, internalNodeCount, nodeCount, branchCount, patternCount;
    
    leafCount = rooted_tree->tip_count;
    internalNodeCount = rooted_tree->inner_count;
    nodeCount = internalNodeCount + leafCount;
    branchCount = rooted_tree->edge_count;
    patternCount=msa->length;
    
    branch_lengths = (double *)malloc(branchCount * sizeof(double));
    matrix_indices = (unsigned int *)malloc(branchCount * sizeof(unsigned int));
    operations = (pll_operation_t *)malloc(internalNodeCount *
                                           sizeof(pll_operation_t));
    
    //pllmod_subst_model_t * model = pllmod_util_model_info_genotype(JC_MODEL);
    
    number_genotype_states= 10;
    number_phased_genotype_states= 16;
    number_states = number_phased_genotype_states;
    
    
    reference_partition = new Partition(leafCount,// numberTips
                                        internalNodeCount,//unsigned  int  clvBuffers
                                        //3* internalNodeCount,//unsigned  int  clvBuffers
                                        number_states,// model->states,//numberStates
                                        patternCount,//unsigned  int  sites
                                        1,//unsigned  int numberRateMatrices
                                        branchCount, // unsigned int probMatrices
                                        RATE_CATS,//RATE_CATS, // unsigned  int  numberRateCats
                                        internalNodeCount, // unsigned  int numberScaleBuffers
                                        0, //int statesPadded
                                        sse, avx, avx2, avx512, asc, tipPatternCompression);
    
    // unsigned int * patterns_weights = pll_compress_site_patterns(msa->sequence,
    //                  pll_map_gt10,
    //                   rooted_tree->tip_count,
    //                   &(msa->length));
    
    //   reference_partition->setPatternWeights(patterns_weights);
    
    if (!(reference_partition->getPartition())){
        
        fprintf(stderr, "PLL error %d: %s\n", pll_errno, pll_errmsg);
        exit(pll_errno);
        
    }
    //pll_repeats_t * repeats = reference_partition->getPartition()->repeats;
    // int use_repeats = pll_repeats_enabled(reference_partition->getPartition());
    
    // std::cout<< "use repeats " <<use_repeats << std::endl;
    //    if (use_repeats)
    //        pll_repeats_t * repeats = reference_partition->getPartition()->repeats;
    //
    // createHashTable();
    
    recomputeTipClvs();
    
    double subst_params[ number_states*(number_states-1) /2] ;
    std::fill_n (subst_params,number_states*(number_states-1) /2, 1);
    
    double  rate_categories[number_states];
    std::fill_n (rate_categories, number_states*(number_states-1) /2, 1);
    
    double  category_weights[number_states];
    std::fill_n (category_weights, number_states*(number_states-1) /2, 1);
    
    double nucleotide_frequencies[number_states];
    std::fill_n (nucleotide_frequencies, number_states, 1.0/(number_states));
    
    pll_compute_gamma_cats(1, RATE_CATS, rate_categories, PLL_GAMMA_RATES_MEAN);
    reference_partition->setCategoryRates(rate_categories);
    
    reference_partition->setCategoryWeights(category_weights);
    reference_partition->setFrequencies(0, nucleotide_frequencies);
    reference_partition->setSubstParams(0,subst_params);
    reference_partition->updateEigen(0);
    
    initOperations();
    
};

double TreeLikelihood::computeRootLogLikelihood(){
    
    double logLik=0.0;
    unsigned int paramsIndexes[RATE_CATS] = {0};
    int ntips= rooted_tree->tip_count;
    std::vector<unsigned int> a(ntips);
    std::iota(a.begin(), a.end(), 0);
    
    double * persite_lnl = (double *) malloc(reference_partition->numberSites * sizeof(double));
    logLik = reference_partition->computeRootLogLikelihood(rooted_tree->root->clv_index,
                                                           rooted_tree->root->scaler_index,
                                                           paramsIndexes,
                                                           persite_lnl);
    free(persite_lnl);
    
    return(logLik);
}
void TreeLikelihood::createHashTable(){
    int leafsCount=rooted_tree->tip_count;
    pll_rnode_t ** tipnodes = rooted_tree->nodes;
    
    hcreate(leafsCount);
    
    /* populate a libc hash table with tree tip labels */
    unsigned int * data = (unsigned int *)malloc(leafsCount *
                                                 sizeof(unsigned int));
    //char *label;
    for (size_t i = 0; i < leafsCount; ++i)
    {
        
        data[i] = tipnodes[i]->clv_index;
        ENTRY entry;
        //label= tipnodes[i]->label;
#ifdef __APPLE__
        entry.key = pll_utils::xstrdup(tipnodes[i]->label);
#else
        entry.key = tipnodes[i]->label;
#endif
        entry.key = tipnodes[i]->label;
        entry.data = (void *)(data+i);
        hsearch(entry, ENTER);
    }
    
    
    hdestroy();
    free(data);
}
int TreeLikelihood::setPartitionTips(bool doUseGenotypes){
    
    int leafsCount=rooted_tree->tip_count;
    pll_rnode_t ** tipnodes = rooted_tree->nodes;
    
    hcreate(leafsCount);
    
    /* populate a libc hash table with tree tip labels */
    unsigned int * data = (unsigned int *)malloc(leafsCount *
                                                 sizeof(unsigned int));
    char *label;
    for (size_t i = 0; i < leafsCount; ++i)
    {
        
        data[i] = tipnodes[i]->clv_index;
        ENTRY entry;
        label= tipnodes[i]->label;
#ifdef __APPLE__
        entry.key = pll_utils::xstrdup(tipnodes[i]->label);
#else
        entry.key = tipnodes[i]->label;
#endif
        entry.key = tipnodes[i]->label;
        entry.data = (void *)(data+i);
        hsearch(entry, ENTER);
    }
    
    unsigned int tip_clv_index;
    
    /* find sequences in hash table and link them with the corresponding taxa */
    for (unsigned i = 0; i < msa->count; ++i)
    {
        ENTRY query;
        query.key = msa->label[i];
        ENTRY * found = NULL;
        
        found = hsearch(query,FIND);
        
        if (!found)
            fprintf(stderr,"Sequence with header %s does not appear in the tree", msa->label[i]);
        
        tip_clv_index = *((unsigned int *)(found->data));
        // std::cout << "clv index for tip with label " <<msa->label[i] << " is "<< tip_clv_index << std::endl;
        
        if (reference_partition->setTipStates(tip_clv_index, msa->sequence[i]) != PLL_SUCCESS)
        {
            fprintf(stderr, "PLL error %d: %s\n", pll_errno, pll_errmsg);
            exit(pll_errno);
        }
    }
    hdestroy();
    
    return PLL_SUCCESS;
}
int TreeLikelihood::initOperations(){
    
    int tip_nodes_count = rooted_tree->tip_count;
    int   inner_nodes_count = tip_nodes_count-1;
    int   nodes_count = inner_nodes_count + tip_nodes_count;
    int   branch_count = nodes_count - 1;
    unsigned int matrix_count;
    
    pll_rnode_t ** travbuffer = (pll_rnode_t **)malloc(nodes_count * sizeof(pll_rnode_t *));
    
    
    /* get inner nodes */
    pll_rnode_t ** inner_nodes_list = (pll_rnode_t **)malloc(rooted_tree->tip_count *
                                                             sizeof(pll_rnode_t *));
    memcpy(inner_nodes_list,
           rooted_tree->nodes+tip_nodes_count,
           inner_nodes_count*sizeof(pll_rnode_t *));
    
    unsigned int traversal_size;
    
    /* compute a partial traversal starting from the randomly selected
     inner node */
    
    if (!pll_rtree_traverse(rooted_tree->root,
                            PLL_TREE_TRAVERSE_POSTORDER,
                            pll_utils::cb_rfull_traversal,
                            travbuffer,
                            &traversal_size)){
        fprintf(stderr, "Function pll_rtree_traverse() root node as parameter %d: %s\n", pll_errno, pll_errmsg);
        exit(pll_errno);
        
    }
    
    /* given the computed traversal descriptor, generate the operations
     structure, and the corresponding probability matrix indices that
     may need recomputing */
    pll_rtree_create_operations(travbuffer,
                                traversal_size,
                                branch_lengths,
                                matrix_indices,
                                operations,
                                &matrix_count,
                                &ops_count);
    
    unsigned int paramsIndexes[RATE_CATS] = {0};
    // int ntips= rooted_tree->tip_count;
    // std::vector<unsigned int> a(ntips);
    std::vector<unsigned int> a(branch_count);
    std::iota(a.begin(), a.end(), 0);
    
    reference_partition->updateProbMatrices(paramsIndexes, &a[0], branch_lengths, rooted_tree->edge_count);
    
    reference_partition->updatePartials(operations, ops_count);
    //    unsigned  int precision =6;
    //    for (unsigned int i = 0; i < branch_count; ++i)
    //       {
    //         printf ("P-matrix (%d) for branch length %f\n", a[i], branch_lengths[i]);
    //         printf ("P-matrix (%d) for branch length %f\n", i, branch_lengths[i]);
    //
    //         pll_show_pmatrix(reference_partition->getPartition(), i,precision);
    //         printf ("\n");
    //       }
    free(inner_nodes_list);
    free(travbuffer);
    
    return PLL_SUCCESS;
}
void TreeLikelihood::fillTipClv(unsigned int tip_id, std::vector<double> &clv) const
{
    
    auto clv_size = msa->length * number_states;
    if( clv.size()!=clv_size)
        clv.resize(clv_size);
    
    auto clvp = clv.begin();
    auto seq = msa->sequence[tip_id];
    //auto charmap = _model.charmap();
    auto charmap = pll_map_gt10;
    
    for (size_t j = 0; j < msa->length; ++j)
    {
        auto charstate = (pll_state_t) seq[j];
        pll_state_t state = charmap ? charmap[(int) charstate] : charstate;
        
        
        gtErrorModel->computeStateErrorProbPT17(state, clvp);
        
        //if (j == 0 && 0)
        //{
        //      printf("state: %llu ", state);
        //      for (size_t k = 0; k < number_states; ++k)
        //        printf("%lf ", clvp[k]);
        //      printf("\n");
        // }
        
        clvp += number_states;
    }
    
    assert(clvp == clv.end());
}
void TreeLikelihood::fillTipClv(unsigned int tip_id, std::vector<double> &clv, pll_msa_t *msaP, GenotypeErrorModel &gtError,unsigned int number_statesP)
{
    
    
    auto clv_size = msaP->length * number_statesP;
    if( clv.size()!=clv_size)
        clv.resize(clv_size);
    
    auto clvp = clv.begin();
    //auto seq = msa->at(tip_id);
    auto seq = msaP->sequence[tip_id];
    //auto charmap = _model.charmap();
    auto charmap = pll_map_gt10;
    
    for (size_t j = 0; j < msaP->length; ++j)
    {
        auto charstate = (pll_state_t) seq[j];
        pll_state_t state = charmap ? charmap[(int) charstate] : charstate;
        
        
        gtError.computeStateErrorProbPT17(state, clvp);
        
        //if (j == 0 && 0)
        //{
        //      printf("state: %llu ", state);
        //      for (size_t k = 0; k < number_states; ++k)
        //        printf("%lf ", clvp[k]);
        //      printf("\n");
        // }
        
        clvp += number_statesP;
    }
    
    assert(clvp == clv.end());
}
void TreeLikelihood::recomputeTipClvs()
{
    //assert(!(reference_partition->getPartition()->attributes & PLL_ATTRIB_PATTERN_TIP));
    
    std::vector<double> tmp_clv;
    for (size_t i = 0; i < msa->count; ++i)
    {
        fillTipClv(i, tmp_clv);
        reference_partition->setTipCLV(i, tmp_clv.data());
    }
}

void TreeLikelihood::destroyHashTable(){
    
    hdestroy();
    
}
void TreeLikelihood::changeGenotypeErrorModel(GenotypeErrorModel *newGtErrorModel){
    
    *gtErrorModel=*newGtErrorModel;
    //update clv tips in the  partitions
    
    recomputeTipClvs();
    
    double subst_params[ number_states*(number_states-1) /2] ;
    std::fill_n (subst_params,number_states*(number_states-1) /2, 1);
    
    double  rate_categories[number_states];
    std::fill_n (rate_categories, number_states*(number_states-1) /2, 1);
    
    double  category_weights[number_states];
    std::fill_n (category_weights, number_states*(number_states-1) /2, 1);
    
    double nucleotide_frequencies[number_states];
    std::fill_n (nucleotide_frequencies, number_states, 1.0/(number_states));
    
    pll_compute_gamma_cats(1, RATE_CATS, rate_categories, PLL_GAMMA_RATES_MEAN);
    reference_partition->setCategoryRates(rate_categories);
    
    reference_partition->setCategoryWeights(category_weights);
    reference_partition->setFrequencies(0, nucleotide_frequencies);
    reference_partition->setSubstParams(0,subst_params);
    reference_partition->updateEigen(0);
    
    initOperations();
}
