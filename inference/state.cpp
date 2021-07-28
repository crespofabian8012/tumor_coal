//
//  state.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 16/03/2021.
//


#include "state.hpp"
#include "pll_utils.hpp"
#include <memory>
#include "poset_smc_params.hpp"


State::State(PosetSMCParams &smcParams):
pll_buffer_manager(smcParams.pll_buffer_manager),partition(smcParams.partition), num_sites(smcParams.num_sites)
{
    height = 0.0;
    
    initForest(  smcParams.sampleSize, smcParams.msa, smcParams.positions,  smcParams.getProgramOptions());
    
    populationSet = new PopulationSet(smcParams.getProgramOptions().numClones);
    
    //auto deltaTs  = populationSet->samplePopulationGrowthRateFromPriors(mcmcOptions, rngGsl );
       

       populationSet->initPopulationSampleSizes(smcParams.sampleSizes);
       
       populationSet->initPopulationGametes();
       populationSet->initPopulationRTips();
       populationSet->initProportionsVectorFromSampleSizes(smcParams.sampleSizes);
    ;
}

void State::initForest( int sampleSize, pll_msa_t *msa, std::vector<int> &positions, ProgramOptions &programOptions){
    
    assert(sampleSize==positions.size());
    
    const  pll_partition_t *reference_partition = partition;

    //pll_rtree_t * tree;
    std::shared_ptr<PartialTreeNode>  p;
    //TreeNode *t;
    for(size_t i=0; i< positions.size();++i){
        
        unsigned int sites_alloc =
            reference_partition->asc_bias_alloc ? reference_partition->sites + reference_partition->states : reference_partition->sites;

        unsigned int clv_size =
            sites_alloc * reference_partition->states_padded * reference_partition->rate_cats * sizeof(double);
        unsigned int scaler_size = (reference_partition->attributes & PLL_ATTRIB_RATE_SCALERS)
                                       ? sites_alloc * reference_partition->rate_cats
                                       : sites_alloc;

        
        p=std::make_shared<PartialTreeNode>(pll_buffer_manager, nullptr, nullptr, msa->label[i], height, clv_size,
        scaler_size);

        delete (p->clv);
        delete (p->scale_buffer);
        p->clv = reference_partition->clv[i];
        p->scale_buffer = nullptr;
        p->ln_likelihood = compute_ln_likelihood(p->clv, nullptr, reference_partition);

        roots.push_back(p);
    }
    height = 0.0;
}

//std::shared_ptr<State> State::makeDeepCopy(){
//    
//    std::shared_ptr<State> newState = std::make_shared<State>(height);
//    std::vector<pll_rnode_t *> forestCopy;
//    pll_rnode_t *tree;
//    pll_rnode_t *treeCopy;
//    
////    for(size_t i=0; i< roots.size();++i){
////
////        char * newick=pll_rtree_export_newick(tree, NULL);
////        //pll_utree_graph_clone
////
////       // treeCopy =pll_rtree_parse_newick_string(newick);
////        //treeCopy = pll_utils::cloneRTree(tree);
////
////        forestCopy.push_back(treeCopy);
////    }
////    newState->setForest(forestCopy);
//    return newState;
//}

State &State::operator=(const State &original){
    
    if (this == &original)
       return *this;

     partition = original.partition;
      height = original.height;
     roots = original.roots;
     return *this;
}

std::vector<std::shared_ptr<PartialTreeNode>> State::propose(int i, int j, double height_delta){
    
    std::vector<std::shared_ptr<PartialTreeNode>> result;
    
    assert(roots.size() > 1 && "Expected more than one root");
    assert(i != j && "Cannot connect, this would make a loop");
    assert(height_delta >= 0 && "Height change can't be negative");
    assert(i >= 0 && i < roots.size() && j >= 0 && j < roots.size() &&
           "Index out of bounds");

    const pll_partition_t *p =partition;

    height += height_delta;

    std::shared_ptr<PartialTreeNode> child_left = roots[i];
    std::shared_ptr<PartialTreeNode> child_right = roots[j];

    double left_length =height - child_left->height;
    double right_length = height - child_right->height;

    unsigned int pmatrix_size =
        p->states * p->states_padded * p->rate_cats * sizeof(double);

    std::shared_ptr<PartialTreeEdge> edge_left = std::make_shared<PartialTreeEdge>(
        pll_buffer_manager, child_left, left_length, pmatrix_size);

    std::shared_ptr<PartialTreeEdge> edge_right = std::make_shared<PartialTreeEdge>(
        pll_buffer_manager, child_right, right_length, pmatrix_size);

    unsigned int sites_alloc =
        p->asc_bias_alloc ? p->sites + p->states : p->sites;

    unsigned int clv_size =
        sites_alloc * p->states_padded * p->rate_cats * sizeof(double);
    unsigned int scaler_size = (p->attributes & PLL_ATTRIB_RATE_SCALERS)
                                   ? sites_alloc * p->rate_cats
                                   : sites_alloc;
    unsigned int scale_buffer_size = scaler_size * sizeof(double);

    std::shared_ptr<PartialTreeNode> parent = std::make_shared<PartialTreeNode>(
        pll_buffer_manager, edge_left, edge_right, "", height, clv_size,
        scale_buffer_size);

    const unsigned int matrix_indices[1] = {0};
    const unsigned int param_indices[4] = {0, 0, 0, 0};

    int left_edge_pmatrix_result = pll_core_update_pmatrix(
        &parent->edge_l->pmatrix, p->states, p->rate_cats, p->rates,
        &parent->edge_l->length, matrix_indices, param_indices, p->prop_invar,
        p->eigenvals, p->eigenvecs, p->inv_eigenvecs, 1, p->attributes);
    assert(left_edge_pmatrix_result == PLL_SUCCESS);

    int right_edge_pmatrix_result = pll_core_update_pmatrix(
        &parent->edge_r->pmatrix, p->states, p->rate_cats, p->rates,
        &parent->edge_r->length, matrix_indices, param_indices, p->prop_invar,
        p->eigenvals, p->eigenvecs, p->inv_eigenvecs, 1, p->attributes);
    assert(right_edge_pmatrix_result == PLL_SUCCESS);

    pll_core_update_partial_ii(p->states, (p->asc_bias_alloc ? p->sites + p->states : p->sites),
        p->rate_cats, parent->clv, parent->scale_buffer, child_left->clv,
        child_right->clv, edge_left->pmatrix, edge_right->pmatrix,
        child_left->scale_buffer, child_right->scale_buffer, p->attributes);

    parent->ln_likelihood =
        compute_ln_likelihood(parent->clv, parent->scale_buffer, p);

    assert(parent->ln_likelihood <= 0 && "Likelihood can't be more than 100%");

    // Remove children
    remove_roots(i, j);
    // Add new internal node
    roots.push_back(parent);
    
    return result;
}

double State::likelihood_factor(std::shared_ptr<PartialTreeNode> root){
    
    double result=0.0;
    assert(root->edge_l && root->edge_r && "Root cannot be a leaf");

     std::shared_ptr<PartialTreeNode> left = root->edge_l->child;
     std::shared_ptr<PartialTreeNode> right = root->edge_r->child;

     double ln_m = root->ln_likelihood;
     double ln_l = left->ln_likelihood;
     double ln_r = right->ln_likelihood;

     assert(ln_m <= 0 && ln_l <= 0 && ln_r <= 0 &&
            "Likelihood can't be more than 100%");

     return ln_m - (ln_l + ln_r);
    
    return result;
    
}

void State::remove_roots(int i, int j){
    assert(i != j && "Expected different indices");
    assert(i >= 0 && i < roots.size() && j >= 0 && j < roots.size() &&
           "Index out of bounds");

    if (i < j) {
      roots.erase(roots.begin() + j, roots.begin() + j + 1);
      roots.erase(roots.begin() + i, roots.begin() + i + 1);
    } else {
      roots.erase(roots.begin() + i, roots.begin() + i + 1);
      roots.erase(roots.begin() + j, roots.begin() + j + 1);
    }
}

double State::compute_ln_likelihood(double *clv, unsigned int *scale_buffer,
                             const pll_partition_t *p) {
  const unsigned int parameter_indices[4] = {0, 0, 0, 0};

  return pll_core_root_loglikelihood(
      p->states, p->sites, p->rate_cats,

      clv, scale_buffer,

      p->frequencies, p->rate_weights, p->pattern_weights, p->prop_invar,
      p->invariant, parameter_indices, NULL, p->attributes);
}

State::~State(){
    

}
