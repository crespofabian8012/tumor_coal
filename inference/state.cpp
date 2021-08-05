//
//  state.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 16/03/2021.
//


#include "state.hpp"
#include "pll_utils.hpp"
#include "treeLikelihood.hpp"
#include <memory>
#include "poset_smc_params.hpp"


State::State(PosetSMCParams &smcParams, gsl_rng *random):
pll_buffer_manager(smcParams.pll_buffer_manager),partition(smcParams.partition), num_sites(smcParams.num_sites)
{
    height = 0.0;

    //TODO put random values on Genotype errors model
    gtError =  new GenotypeErrorModel("GT20", smcParams.getProgramOptions().meanGenotypingError,  1.0 - sqrt (1.0 - smcParams.getProgramOptions().fixedADOrate), 16);;
    populationSet = new PopulationSet(smcParams.getProgramOptions().numClones);
    
    populationSet->initPopulationSampleSizes(smcParams.sampleSizes);
    
    populationSet->initPopulationGametes();
    populationSet->initPopulationRTips();
    populationSet->initProportionsVectorFromSampleSizes(smcParams.sampleSizes);
    
    populationSet->initDeltaThetaFromPriors(random);
    populationSet->setPopulationsToriginConditionalDelta(random);
    populationSet->sortPopulationsByTorigin();
    
    populationSet->initListPossibleMigrations();
    populationSet->initPopulationsCoalescentEvents();
    
    initForest(smcParams.sampleSize, smcParams.msa, smcParams.positions,  smcParams.getProgramOptions());
    nextAvailable = smcParams.sampleSize+1;

}
State::State(const State &original):pll_buffer_manager(original.pll_buffer_manager){

    height=original.height;
    num_sites=original.num_sites;
    populationSet=original.populationSet;
    gtError=original.gtError;
    roots=original.roots;
    partition= original.partition;
    nextAvailable = original.nextAvailable;
   
}
void State::initForest( int sampleSize, pll_msa_t *msa, std::vector<int> &positions, ProgramOptions &programOptions){
    
    
    const Partition *reference_partition = partition;
    
    std::shared_ptr<PartialTreeNode>  p;
    std::vector<int> CumSumNodes(programOptions.numClones+1);
    CumSumNodes[0] = 0;
    Population *pop;
    int currentPopIndex;
    int  cumIndivid = 0;
    
    for (size_t i = 1; i <= programOptions.numClones; i++)
    {
        pop = populationSet->getPopulationbyIndex(i-1);
        pop->resetGametesCounters();
        CumSumNodes[i] = CumSumNodes[i - 1] + pop->sampleSize;
    }
    
    
    for(size_t i=0; i< positions.size();++i){
        
        cumIndivid++;
        unsigned int sites_alloc =
        reference_partition->ascBiasCorrection() ? reference_partition->numberSites + reference_partition->numberStates : reference_partition->numberSites;
        
        unsigned int clv_size =
        sites_alloc * reference_partition->getStatesPadded() * reference_partition->numberRateCats * sizeof(double);
        unsigned int scaler_size = (reference_partition->attributes & PLL_ATTRIB_RATE_SCALERS)
        ? sites_alloc * reference_partition->numberRateCats
        : sites_alloc;
        
        
        p=std::make_shared<PartialTreeNode>(pll_buffer_manager, nullptr, nullptr, msa->label[i], height, clv_size,
                                            scaler_size);
        
        delete (p->clv);
        delete (p->scale_buffer);
        
        std::vector<double> tmp_clv;
     
        reference_partition->computeCLV(i, msa, gtError, tmp_clv);
   
       // reference_partition->initTipCLV(i, tmp_clv.data());
        
        //p->clv = reference_partition->getCLV(i);
        p->clv =tmp_clv.data();
        p->scale_buffer = nullptr;
        p->ln_likelihood = compute_ln_likelihood(p->clv, nullptr, reference_partition);
        for (size_t j = 1; j <= programOptions.numClones; j++)
        {
            pop =populationSet->getPopulationbyIndex(j-1);
            currentPopIndex = pop->index;
            // Identify to which clone belongs this node
            if (cumIndivid <= CumSumNodes[j] && cumIndivid > CumSumNodes[j - 1]){
                pop->idsActiveGametes[pop->numActiveGametes]=i;
                pop->numActiveGametes=pop->numActiveGametes+1;
                pop->idsGametes[pop->numGametes]=i;
                pop->numGametes=pop->numGametes+1;
                p->index_population=currentPopIndex;
                break;
            }
        }
        
        
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
    num_sites=original.num_sites;
    roots = original.roots;
    gtError=original.gtError;
    populationSet = original.populationSet;
    nextAvailable = original.nextAvailable;
    return *this;
 
}


std::shared_ptr<PartialTreeNode> State::connect(int i, int j, double height_delta){
    
    std::vector<std::shared_ptr<PartialTreeNode>> result;
    
    assert(roots.size() > 1 && "Expected more than one root");
    assert(i != j && "Cannot connect, this would make a loop");
    assert(height_delta >= 0 && "Height change can't be negative");
    assert(i >= 0 && i < roots.size() && j >= 0 && j < roots.size() &&
           "Index out of bounds");
    
    const Partition *p =partition;
    
    height += height_delta;
    
    std::shared_ptr<PartialTreeNode> child_left = roots[i];
    std::shared_ptr<PartialTreeNode> child_right = roots[j];
    
    double left_length =height - child_left->height;
    double right_length = height - child_right->height;
    
    unsigned int pmatrix_size =
    p->numberStates * p->getStatesPadded() * p->numberRateCats * sizeof(double);
    
    std::shared_ptr<PartialTreeEdge> edge_left = std::make_shared<PartialTreeEdge>(
                                                                                   pll_buffer_manager, child_left, left_length, pmatrix_size);
    
    std::shared_ptr<PartialTreeEdge> edge_right = std::make_shared<PartialTreeEdge>(
                                                                                    pll_buffer_manager, child_right, right_length, pmatrix_size);
    
    unsigned int sites_alloc =
    p->ascBiasCorrection() ? p->numberSites + p->numberStates : p->numberSites;
    
    unsigned int clv_size =
    sites_alloc * p->getStatesPadded() * p->numberRateCats * sizeof(double);
    unsigned int scaler_size = (p->attributes & PLL_ATTRIB_RATE_SCALERS)
    ? sites_alloc * p->numberRateCats
    : sites_alloc;
    unsigned int scale_buffer_size = scaler_size * sizeof(double);
    
    std::shared_ptr<PartialTreeNode> parent = std::make_shared<PartialTreeNode>(
                                                                                pll_buffer_manager, edge_left, edge_right, "", height, clv_size,
                                                                                scale_buffer_size);
    
    const unsigned int matrix_indices[1] = {0};
    const unsigned int param_indices[4] = {0, 0, 0, 0};
    
    int left_edge_pmatrix_result = pll_core_update_pmatrix(
                                                           &parent->edge_l->pmatrix, p->numberStates, p->numberRateCats, p->rates(),
                                                           &parent->edge_l->length, matrix_indices, param_indices, p->propInvar(),
                                                           p->eigenVals(), p->eigenVecs(), p->invEigenVecs(), 1, p->attributes);
    assert(left_edge_pmatrix_result == PLL_SUCCESS);
    
    int right_edge_pmatrix_result = pll_core_update_pmatrix(
                                                            &parent->edge_r->pmatrix, p->numberStates, p->numberRateCats, p->rates(),
                                                            &parent->edge_r->length, matrix_indices, param_indices, p->propInvar(),
                                                            p->eigenVals(), p->invEigenVecs(), p->invEigenVecs(), 1, p->attributes);
    assert(right_edge_pmatrix_result == PLL_SUCCESS);
    
    pll_core_update_partial_ii(p->numberStates, (p->ascBiasCorrection() ? p->numberSites + p->numberStates : p->numberSites),
                               p->numberRateCats, parent->clv, parent->scale_buffer, child_left->clv,
                               child_right->clv, edge_left->pmatrix, edge_right->pmatrix,
                               child_left->scale_buffer, child_right->scale_buffer, p->attributes);
    
    parent->ln_likelihood =
    compute_ln_likelihood(parent->clv, parent->scale_buffer, p);
    
    assert(parent->ln_likelihood <= 0 && "Likelihood can't be more than 100%");
    
    // Remove children
    remove_roots(i, j);
    // Add new internal node
    roots.push_back(parent);
    
    return parent;
}

double State::likelihood_factor(std::shared_ptr<PartialTreeNode> root)const{
    
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
                                    const Partition *p) {
    const unsigned int parameter_indices[4] = {0, 0, 0, 0};
   // const unsigned int parameter_indices[16] = {0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    return pll_core_root_loglikelihood(
                                       p->numberStates, p->numberSites, p->numberRateCats,
                                       
                                       clv, scale_buffer,
                                       
                                       p->frequencies(), p->rateWeights(), p->patternWeights(), p->propInvar(),
                                       p->invariant(), parameter_indices, NULL, p->attributes);
}

State::~State(){
    
    
}
