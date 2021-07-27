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

std::vector<std::shared_ptr<PartialTreeNode>> State::propose(int i, int j, double height){
    
    std::vector<std::shared_ptr<PartialTreeNode>> result;
    
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
