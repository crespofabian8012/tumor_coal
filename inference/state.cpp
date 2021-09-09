//
//  state.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 16/03/2021.
//


#include "state.hpp"
#include "utils.hpp"
#include "pll_utils.hpp"
#include "core_likelihood.hpp"
#include "treeLikelihood.hpp"
#include <memory>
#include "poset_smc_params.hpp"


State::State(PosetSMCParams &smcParams, gsl_rng *random):
pll_buffer_manager(smcParams.pll_buffer_manager),partition(smcParams.partition), num_sites(smcParams.num_sites)
{
    topHeightModelTime =0.0;
    topHeightScaledByTheta = 0.0;
    logWeight = 0.0;
    theta = *smcParams.theta;
    initIdsNextCoalEvents(smcParams.getProgramOptions().numClones);
    initIdsInmigrationEvents(smcParams.getProgramOptions().numClones);
    gtError =  smcParams.gtErrorModel;
    
    
    populationSet = new PopulationSet(smcParams.getProgramOptions().numClones);
    
    populationSet->initPopulationSampleSizes(smcParams.sampleSizes);
    populationSet->initPopulationGametes();
    populationSet->initPopulationRTips();
    populationSet->initProportionsVectorFromSampleSizes(smcParams.sampleSizes);
    populationSet->initPopulationDeltas( smcParams.populationDeltaTs);
    populationSet->initPopulationTOriginSTD( smcParams.populationToriginSTDs);
    populationSet->initTheta(theta);
    
    // populationSet->setPopulationsToriginConditionalDelta(random);
    populationSet->sortPopulationsByTorigin();
    
    populationSet->initListPossibleMigrations();
    populationSet->initPopulationsCoalescentEvents();
  
    
    
    initForest(smcParams.sampleSize, smcParams.msa, smcParams.positions,  smcParams.getProgramOptions());
    nextAvailable = smcParams.sampleSize-1;
    
//    double K= smcParams.getProgramOptions().K;
//       populationSet->sampleEventTimesScaledByProportion( random, K,smcParams.getProgramOptions().noisy);
//      populationSet->resetNumActiveGametesCounter();
    assert(smcParams.positions.size()>0);
    postorder.reserve(2*smcParams.positions.size()-1);
    coalEventTimesScaledByTheta.reserve(smcParams.positions.size()-1);
}
State::State(const State &original):pll_buffer_manager(original.pll_buffer_manager){
    
    topHeightModelTime = original.topHeightModelTime;
    topHeightScaledByTheta = original.topHeightScaledByTheta;
    num_sites=original.num_sites;
    
    populationSet=new PopulationSet(*(original.populationSet));
    gtError = new GenotypeErrorModel(*(original.gtError));
    
    idsNextCoalEvents = original.idsNextCoalEvents;
    idsNextInmigrationEvents = original.idsNextInmigrationEvents;
    roots = original.roots;
    postorder  = original.postorder;
    coalEventTimesScaledByTheta = original.coalEventTimesScaledByTheta;
    //    roots.reserve(original.root_count());
    //    for (auto const& fptr : original.getRoots())
    //         roots.emplace_back(fptr->Clone());
    
    
    partition = original.partition;
    nextAvailable = original.nextAvailable;
    theta = original.theta;
    logWeight = original.logWeight;
    
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
        
        unsigned int clv_num_elements= sites_alloc * reference_partition->numberStates * reference_partition->numberRateCats;
        // unsigned int clv_size =
        // clv_num_elements * sizeof(double);
        unsigned int scaler_size = (reference_partition->attributes & PLL_ATTRIB_RATE_SCALERS)
        ? sites_alloc * reference_partition->numberRateCats
        : sites_alloc;
        
        p=std::make_shared<PartialTreeNode>(pll_buffer_manager, nullptr, nullptr, msa->label[i], 0.0, clv_num_elements,
                                            scaler_size, reference_partition->alignment(), i);
        
        //delete (p->pclv);
        //delete (p->pscale_buffer);
        // auto clv_elem = reference_partition->numberSites * reference_partition->numberStates;
        
        p->buildCLV(i,reference_partition->numberStates,  msa, gtError,  false);
        
        std::vector<unsigned  int> invSites(reference_partition->numberStates, 0);
        
        pll_count_invariant_sites(reference_partition->getPartition(),  invSites.data());
        
        p->ln_likelihood = compute_ln_likelihood(p->pclv, nullptr, reference_partition->getPartition());
        
        logWeight +=p->ln_likelihood;
        
        // std::cout << " tip " << i << " loglik " <<p->ln_likelihood << std::endl;
        if (p->ln_likelihood <= -10000000 ){
            std::cout << " tip " << i << " loglik " <<p->ln_likelihood << std::endl;
            p->showpClV(reference_partition->numberStates, reference_partition->numberRateCats, reference_partition->getStatesPadded(), reference_partition->numberSites, 3);
        }
        
        for (size_t j = 1; j <= programOptions.numClones; j++)
        {
            pop =populationSet->getPopulationbyIndex(j-1);
            currentPopIndex = pop->index;
            // Identify to which clone belongs this node
            if (cumIndivid <= CumSumNodes[j] && cumIndivid > CumSumNodes[j - 1]){
                pop->idsActiveGametes[pop->numActiveGametes]=p->index;
                pop->numActiveGametes=pop->numActiveGametes+1;
                pop->idsGametes[pop->numGametes]=i;
                pop->numGametes=pop->numGametes+1;
                p->index_population=currentPopIndex;
                break;
            }
        }
        
        roots.push_back(p);
    }
    if (programOptions.numClones==1)
        assert(populationSet->getPopulationbyIndex(0)->numActiveGametes == positions.size());
    
}

void State::resetNumActiveGametesCounter(){
    
    populationSet->resetNumActiveGametesCounter();
}


State &State::operator=(const State &original){
    
    if (this == &original)
        return *this;
    
    partition = original.partition;
    topHeightScaledByTheta = original.topHeightScaledByTheta;
    topHeightModelTime = original.topHeightModelTime;
    theta = original.theta;
    num_sites=original.num_sites;
    roots = original.roots;
    postorder  = original.postorder;
    coalEventTimesScaledByTheta  = original.coalEventTimesScaledByTheta;
    gtError=original.gtError;
    populationSet = original.populationSet;
    nextAvailable = original.nextAvailable;
    idsNextCoalEvents = original.idsNextCoalEvents;
    idsNextInmigrationEvents = original.idsNextInmigrationEvents;
    return *this;
    
}


std::shared_ptr<PartialTreeNode> State::connect(int firstId, int secondId, size_t index_pop_new_node, double newNodeHeight, double logLikNewNode ){
    
    // std::vector<std::shared_ptr<PartialTreeNode>> result;
    //unsigned int float_precision = 3;
    
    assert(roots.size() > 1 && "Expected more than one root");
    
    
    const Partition *p =partition;
    
    //heightScaledByTheta = heightScaledByTheta+ height_delta;
    int idxFirstID = getNodeIdxById(firstId);
    int idxSecondId = getNodeIdxById(secondId);
    
    assert(idxFirstID>=0 &&  idxSecondId >=0 && "Indexes must be positive");
    // assert(height_delta >= 0.0 && "Height change can't be negative");
    assert(idxFirstID != idxSecondId && "Cannot connect, this would make a loop");
    // assert(height_delta >= 0.0 && "Height change can't be negative");
    assert(idxFirstID >= 0 && idxFirstID < roots.size() && idxSecondId >= 0 && idxSecondId < roots.size() &&
           "Index out of bounds");
    
    std::shared_ptr<PartialTreeNode> child_left = roots[idxFirstID];
    std::shared_ptr<PartialTreeNode> child_right = roots[idxSecondId];
    
    //std::cout<< "roots first idx "<< idxFirstID <<  " idxSecondId " << idxSecondId<< std::endl;
    
    double left_length = newNodeHeight - child_left->height;
    double right_length =  newNodeHeight - child_right->height;
    
    unsigned int pmatrix_elements = p->numberStates * p->numberStates  * p->numberRateCats;
    //unsigned int pmatrix_size = pmatrix_elements * sizeof(double);
    
    std::shared_ptr<PartialTreeEdge> edge_left = std::make_shared<PartialTreeEdge>(
                                                                                   pll_buffer_manager, child_left, left_length, pmatrix_elements, partition->alignment());
    
    std::shared_ptr<PartialTreeEdge> edge_right = std::make_shared<PartialTreeEdge>(
                                                                                    pll_buffer_manager, child_right, right_length, pmatrix_elements, partition->alignment());
    
    unsigned int sites_alloc =
    p->ascBiasCorrection() ? p->numberSites + p->numberStates : p->numberSites;
    
    unsigned int clv_elements = sites_alloc  * p->numberStates *  p->numberRateCats;
    
    // unsigned int clv_size = clv_elements * sizeof(double);
    //sites_alloc * p->getStatesPadded() * p->numberRateCats * sizeof(double);
    
    unsigned int scaler_size = (p->attributes & PLL_ATTRIB_RATE_SCALERS)
    ? sites_alloc * p->numberRateCats
    : sites_alloc;
    // unsigned int scale_buffer_size = scaler_size * sizeof(double);
    
    
    
    std::shared_ptr<PartialTreeNode> parent = std::make_shared<PartialTreeNode>(
                                                                                pll_buffer_manager, edge_left, edge_right, "", newNodeHeight, clv_elements,
                                                                                scaler_size, partition->alignment(), nextAvailable);
    
    parent->index_population=index_pop_new_node;
    nextAvailable++;
    
    const unsigned int matrix_indices[1] = {0};
    //Array params_indices holds the indices of rate matrices that will be used for each rate category. This array must be of size rate_cats
    
    //const unsigned int param_indices[4] = {0, 0, 0, 0};
    const unsigned int param_indices[1] = {0};
    
    
    //    int left_edge_pmatrix_result= pll_core_update_pmatrix2( &parent->edge_l->pmatrix, p->numberStates,p->numberRateCats, p->rates(),  &parent->edge_l->length,matrix_indices, param_indices, p->eigenVals(), p->eigenVecs(),  p->invEigenVecs(), 1,  p->attributes);
    
    //    int left_edge_pmatrix_result= pll_core_update_pmatrix_16x16_jc69( &parent->edge_l->pmatrix, p->numberStates,p->numberRateCats, p->rates(),  &parent->edge_l->length,matrix_indices, param_indices,
    //                                                                     1,
    //                                                                     p->attributes);
    int left_edge_pmatrix_result = pll_core_update_pmatrix(
                                                           &parent->edge_l->pmatrix,//double **pmatrix,
                                                           p->numberStates,// unsigned int states
                                                           p->numberRateCats,//unsigned int rate_cats
                                                           p->rates(),//const double *rates
                                                           &parent->edge_l->length,//const double *branch_lengths
                                                           matrix_indices,//const unsigned int *matrix_indices
                                                           param_indices,//const unsigned int *param_indices
                                                           p->propInvar(),//const double *prop_invar
                                                           p->eigenVals(),//double *const *eigenvals
                                                           p->eigenVecs(),//double *const *eigenvecs
                                                           p->invEigenVecs(),//double *const *inv_eigenvecs
                                                           1,//unsigned int count
                                                           p->attributes);
    assert(left_edge_pmatrix_result == PLL_SUCCESS);
    
    //       int right_edge_pmatrix_result= pll_core_update_pmatrix2( &parent->edge_r->pmatrix, p->numberStates,p->numberRateCats, p->rates(),  &parent->edge_r->length,matrix_indices, param_indices, p->eigenVals(), p->eigenVecs(),  p->invEigenVecs(), 1,  p->attributes);
    //
    //    int right_edge_pmatrix_result= pll_core_update_pmatrix_16x16_jc69( &parent->edge_r->pmatrix, p->numberStates,p->numberRateCats, p->rates(),  &parent->edge_r->length,matrix_indices, param_indices,
    //                                                                      1,
    //                                                                      p->attributes);
    
    int right_edge_pmatrix_result = pll_core_update_pmatrix(
                                                            &parent->edge_r->pmatrix, p->numberStates, p->numberRateCats, p->rates(),
                                                            &parent->edge_r->length, matrix_indices, param_indices, p->propInvar(),
                                                            p->eigenVals(), p->eigenVecs(), p->invEigenVecs(), 1, p->attributes);
    assert(right_edge_pmatrix_result == PLL_SUCCESS);
    
    
    unsigned int sites= (p->ascBiasCorrection() ? p->numberSites + p->numberStates : p->numberSites);
    
    
    pll_core_update_partial_ii(p->numberStates, sites,
                               p->numberRateCats, parent->pclv,
                               parent->pscale_buffer,
                               child_left->pclv,
                               child_right->pclv,
                               edge_left->pmatrix,
                               edge_right->pmatrix,
                               child_left->pscale_buffer,
                               child_right->pscale_buffer,
                               p->attributes);
    
    parent->ln_likelihood =
    compute_ln_likelihood(parent->pclv,
                          parent->pscale_buffer,
                          p->getPartition());
    
    
    // std::cout << "log lik of parent node "<<  " is "<< parent->ln_likelihood<< std::endl;
    
    if(isnan(parent->ln_likelihood) || isinf(parent->ln_likelihood)){
        
        std::cout << "Error: log lik  "<<  " is "<< parent->ln_likelihood<< std::endl;
        
    }
    
    // std::cout << "log lik of parent node "<<  " is "<< parent->ln_likelihood<< std::endl;
    
    assert(!isnan(parent->ln_likelihood) && !isinf(parent->ln_likelihood));
    
    parent->ln_coal_likelihood = (child_left->height >  child_right->height)? child_left->ln_coal_likelihood: child_right->ln_coal_likelihood;
      
    parent->ln_coal_likelihood+= logLikNewNode;
    
    // assert(parent->ln_likelihood <= 0 && "Likelihood can't be more than 100%");
    
    // Remove children
    remove_roots(idxFirstID, idxSecondId);
    // Add new internal node
    addRoot(parent);
    //std::cout<< "root count "<< roots.size()<< std::endl;
    
    return parent;
}
std::shared_ptr<PartialTreeNode> State::proposeNewNode(int firstId, int secondId, size_t index_pop_new_node, double newNodeHeight, double logLikNewNode){
    
    assert(roots.size() > 1 && "Expected more than one root");
    
    const Partition *p =partition;
    
    //heightScaledByTheta = heightScaledByTheta+ height_delta;
    int idxFirstID = getNodeIdxById(firstId);
    int idxSecondId = getNodeIdxById(secondId);
    
    assert(idxFirstID>=0 &&  idxSecondId >=0 && "Indexes must be positive");
    // assert(height_delta >= 0.0 && "Height change can't be negative");
    assert(idxFirstID != idxSecondId && "Cannot connect, this would make a loop");
    // assert(height_delta >= 0.0 && "Height change can't be negative");
    assert(idxFirstID >= 0 && idxFirstID < roots.size() && idxSecondId >= 0 && idxSecondId < roots.size() &&
           "Index out of bounds");
    
    std::shared_ptr<PartialTreeNode> child_left = roots[idxFirstID];
    std::shared_ptr<PartialTreeNode> child_right = roots[idxSecondId];
    
    //  std::cout<< "roots first idx "<< idxFirstID <<  " idxSecondId " << idxSecondId<< std::endl;
    
    //  if (!(child_left->label.empty()) && !(child_right->label.empty()))
    //      std::cout<< "child_left->label "<< child_left->label <<  " child_right->label " << child_right->label<< std::endl;
    
    double left_length = newNodeHeight - child_left->height;
    double right_length = newNodeHeight - child_right->height;
    
    unsigned int pmatrix_elements = p->numberStates * p->numberStates  * p->numberRateCats;
    //unsigned int pmatrix_size = pmatrix_elements * sizeof(double);
    
    //  std::cout<< "left_length "<< left_length <<  " right_length " << right_length<< std::endl;
    
    std::shared_ptr<PartialTreeEdge> edge_left = std::make_shared<PartialTreeEdge>(
                                                                                   pll_buffer_manager, child_left, left_length, pmatrix_elements, partition->alignment());
    
    std::shared_ptr<PartialTreeEdge> edge_right = std::make_shared<PartialTreeEdge>(
                                                                                    pll_buffer_manager, child_right, right_length, pmatrix_elements, partition->alignment());
    
    unsigned int sites_alloc =
    p->ascBiasCorrection() ? p->numberSites + p->numberStates : p->numberSites;
    
    unsigned int clv_elements = sites_alloc  * p->numberStates *  p->numberRateCats;
    
    // unsigned int clv_size = clv_elements * sizeof(double);
    //sites_alloc * p->getStatesPadded() * p->numberRateCats * sizeof(double);
    
    unsigned int scaler_size = (p->attributes & PLL_ATTRIB_RATE_SCALERS)
    ? sites_alloc * p->numberRateCats
    : sites_alloc;
    // unsigned int scale_buffer_size = scaler_size * sizeof(double);
    
    std::shared_ptr<PartialTreeNode> parent = std::make_shared<PartialTreeNode>(
                                                                                pll_buffer_manager, edge_left, edge_right, "", newNodeHeight, clv_elements,
                                                                                scaler_size, partition->alignment(), nextAvailable);
    
    parent->index_population=index_pop_new_node;
    
    const unsigned int matrix_indices[1] = {0};
    //Array params_indices holds the indices of rate matrices that will be used for each rate category. This array must be of size rate_cats
    
    const unsigned int param_indices[1] = {0};
    
    
    //    int left_edge_pmatrix_result= pll_core_update_pmatrix_16x16_jc69( &parent->edge_l->pmatrix, p->numberStates,p->numberRateCats, p->rates(),  &parent->edge_l->length,matrix_indices, param_indices,
    //                                                                     1,
    //                                                                    p->attributes);
    
    
    int left_edge_pmatrix_result = pll_core_update_pmatrix(
                                                           &parent->edge_l->pmatrix,//double **pmatrix,
                                                           p->numberStates,// unsigned int states
                                                           p->numberRateCats,//unsigned int rate_cats
                                                           p->rates(),//const double *rates
                                                           &parent->edge_l->length,//const double *branch_lengths
                                                           matrix_indices,//const unsigned int *matrix_indices
                                                           param_indices,//const unsigned int *param_indices
                                                           p->propInvar(),//const double *prop_invar
                                                           p->eigenVals(),//double *const *eigenvals
                                                           p->eigenVecs(),//double *const *eigenvecs
                                                           p->invEigenVecs(),//double *const *inv_eigenvecs
                                                           1,//unsigned int count
                                                           p->attributes);
    
    assert(left_edge_pmatrix_result == PLL_SUCCESS);
    
    //    int right_edge_pmatrix_result= pll_core_update_pmatrix_16x16_jc69( &parent->edge_r->pmatrix, p->numberStates,p->numberRateCats, p->rates(),  &parent->edge_r->length,matrix_indices, param_indices,
    //                                                                      1,
    //                                                                      p->attributes);
    
    int right_edge_pmatrix_result = pll_core_update_pmatrix(&parent->edge_r->pmatrix, p->numberStates,     p->numberRateCats, p->rates(),
                                                            &parent->edge_r->length, matrix_indices, param_indices, p->propInvar(),
                                                            p->eigenVals(), p->eigenVecs(), p->invEigenVecs(), 1, p->attributes);
    assert(right_edge_pmatrix_result == PLL_SUCCESS);
    
    
    unsigned int sites= (p->ascBiasCorrection() ? p->numberSites + p->numberStates : p->numberSites);
    
    
    pll_core_update_partial_ii(p->numberStates, sites,
                               p->numberRateCats, parent->pclv,
                               parent->pscale_buffer,
                               child_left->pclv,
                               child_right->pclv,
                               edge_left->pmatrix,
                               edge_right->pmatrix,
                               child_left->pscale_buffer,
                               child_right->pscale_buffer,
                               p->attributes);
    //    if (1)
    //        parent-> showpClV(p->numberStates, p->numberRateCats, p->statesPadded, sites, 3);
    
    parent->ln_likelihood =
    compute_ln_likelihood(parent->pclv,
                          parent->pscale_buffer,
                          p->getPartition());
    if(isnan(parent->ln_likelihood) || isinf(parent->ln_likelihood)){
          
          std::cout << "Error: log lik  "<<  " is "<< parent->ln_likelihood<< std::endl;
      }
    parent->ln_coal_likelihood = (child_left->height >  child_right->height)? child_left->ln_coal_likelihood: child_right->ln_coal_likelihood;
   
    parent->ln_coal_likelihood+= logLikNewNode;
    

    assert(!isnan(parent->ln_likelihood) && !isinf(parent->ln_likelihood));
    
    //assert(parent->ln_likelihood <= 0 && "Likelihood can't be more than 100%");
    
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
    
    result= ln_m - (ln_l + ln_r);
    
    return result;
    
}

void State::remove_roots(int i, int j){
    assert(i != j && "Expected different indices");
    assert(i >= 0 && i < roots.size() && j >= 0 && j < roots.size() &&
           "Index out of bounds");
    int previous_size= roots.size();
    
    addNodeToPostorderByIndex(i);
    addNodeToPostorderByIndex(j);
    if (i < j) {
        roots.erase(roots.begin() + j, roots.begin() + j + 1);
        roots.erase(roots.begin() + i, roots.begin() + i + 1);
    } else {
        roots.erase(roots.begin() + i, roots.begin() + i + 1);
        roots.erase(roots.begin() + j, roots.begin() + j + 1);
    }
    assert(roots.size()==(previous_size-2));
}

void State::updateIndexesActiveGametes(int idxFirstId, int idxSecondId, size_t idxPopFirst,size_t idxPopSecond, size_t newNodeId, size_t idxPopNewNode){
    
    assert(idxPopFirst <=populationSet->numClones);
    assert(idxPopSecond <=populationSet->numClones);
    Population *popFirst, *popSecond;
    
    if (idxPopFirst == idxPopSecond)//coalescent event
    {
        popFirst = populationSet->getPopulationbyIndex(idxPopFirst);
        popFirst->idsActiveGametes[idxFirstId] = newNodeId;
        popFirst->idsActiveGametes[idxSecondId] = popFirst->idsActiveGametes[ popFirst->numActiveGametes-1];
        
    }
    else//migration
    {
        std::cout<< "migration" << std::endl;
        popFirst = populationSet->getPopulationbyIndex(idxPopFirst);
        popSecond = populationSet->getPopulationbyIndex(idxSecondId);
        
        if (popFirst->numActiveGametes == 1){
            
            popSecond->idsActiveGametes[idxSecondId] = newNodeId;
            
            popFirst->numActiveGametes= popFirst->numActiveGametes-1;
            popSecond->numActiveGametes= popSecond->numActiveGametes+1;
        }
        else{ //(popSecond->numActiveGametes == 1)
            
            popFirst->idsActiveGametes[idxPopFirst] = newNodeId;
            
            popSecond->numActiveGametes= popSecond->numActiveGametes-1;
            popFirst->numActiveGametes= popFirst->numActiveGametes+1;
        }
        
    }
}

double State::compute_ln_likelihood(std::vector<double> &clv, std::vector<unsigned int>  &scale_buffer,
                                    const Partition *p) {
    const unsigned int parameter_indices[1] = {0};
    
    //    double result = pll_core_root_loglikelihood2(  p->numberStates, p->numberSites, p->numberRateCats,
    //                                                 clv.data(), scale_buffer.data(),
    //                                                 p->frequencies(), p->rateWeights(), p->patternWeights(), parameter_indices,  nullptr, p->attributes);
    double result=  pll_core_root_loglikelihood(
                                                p->numberStates, p->numberSites, p->numberRateCats,
                                                
                                                clv.data(), scale_buffer.data(),
                                                
                                                p->frequencies(), p->rateWeights(), p->patternWeights(), p->propInvar(),
                                                p->invariant(), parameter_indices, nullptr, p->attributes);
    return result;
}
double State::compute_ln_likelihood(double *pclv, unsigned int *pscale_buffer,
                                    const Partition *p) {
    const unsigned int parameter_indices[1] = {0};
    //    double result = pll_core_root_loglikelihood2(  p->numberStates, p->numberSites, p->numberRateCats,
    //                                                 pclv, pscale_buffer,
    //                                                 p->frequencies(), p->rateWeights(), p->patternWeights(), parameter_indices,  nullptr, p->attributes);
    double result=  pll_core_root_loglikelihood(
                                                p->numberStates, p->numberSites, p->numberRateCats,
                                                
                                                pclv, pscale_buffer,
                                                
                                                p->frequencies(), p->rateWeights(), p->patternWeights(), p->propInvar(),
                                                p->invariant(), parameter_indices, nullptr, p->attributes);
    return result;
}

double State::compute_ln_likelihood(std::vector<double> &clv, std::vector<unsigned int>  &scale_buffer,
                                    const pll_partition_t *p) {
    
    const unsigned int parameter_indices[1] = {0};
    
    double result = pll_core_root_loglikelihood2(  p->states, p->sites, p->rate_cats,
                                                 clv.data(), scale_buffer.data(),
                                                 p->frequencies, p->rate_weights, p->pattern_weights, parameter_indices,  nullptr, p->attributes);
    
    //        double result= pll_core_root_loglikelihood(
    //                                                   p->states, p->sites, p->rate_cats,
    //                                                   clv, scale_buffer,
    //                                                   p->frequencies, p->rate_weights, p->pattern_weights, p->prop_invar,
    //                                                   p->invariant, parameter_indices, nullptr, p->attributes);
    return result;
}

double State::compute_ln_likelihood(double *pclv, unsigned int* pscale_buffer,
                                    const pll_partition_t *p) {
    // const unsigned int parameter_indices[4] = {0, 0, 0, 0};
    const unsigned int parameter_indices[1] = {0};
    
    
    //    double result = pll_core_root_loglikelihood2(  p->states, p->sites, p->rate_cats,
    //                                                 pclv, pscale_buffer,
    //                                                 p->frequencies, p->rate_weights, p->pattern_weights, parameter_indices,  nullptr, p->attributes);
    
    double result= pll_core_root_loglikelihood(
                                               p->states, p->sites, p->rate_cats,
                                               pclv, pscale_buffer,
                                               p->frequencies, p->rate_weights, p->pattern_weights, p->prop_invar,
                                               p->invariant, parameter_indices, nullptr, p->attributes);
    return result;
}
int State::getNodeIdxById(size_t id){
    
    if (roots.size()==0)
        return -1;
    
    unsigned int i = 0;
    while(i<roots.size() && roots[i]->index!=id)
        i++;
    if (i==roots.size())
        return -1;
    return i;
}
void State::printTree(std::shared_ptr<PartialTreeNode> root, std::ostream &stream) {
    if (root->edge_l && root->edge_r) {
        stream << "(";
        printTree(root->edge_l->child, stream);
        stream << ":" << root->edge_l->length;
        stream << ", ";
        printTree(root->edge_r->child, stream);
        stream << ":" << root->edge_r->length;
        stream << ")";
    } else {
        std::string label(root->label);
        stream << label;
    }
}
std::string State::getNewick(std::shared_ptr<PartialTreeNode> root) {
    std::string newick = "";
    
    newick = getNewickRecursive(root, newick);
    newick.append(";");
    return newick;
}
std::string State::getNewickRecursive(std::shared_ptr<PartialTreeNode> root, std::string &newick) {
    
    if (root->edge_l && root->edge_r) {
        newick.append("(");
        getNewickRecursive(root->edge_l->child, newick);
        newick.append(":");
        newick.append(std::to_string(root->edge_l->length));
        newick.append(",");
        getNewickRecursive(root->edge_r->child, newick);
        newick.append(":");
        newick.append(std::to_string(root->edge_r->length));
        newick.append(")");
    } else {
        std::string label(root->label);
        newick.append(label);
    }
    return newick;
}
void State::printTreeChronologicalOrder(std::shared_ptr<PartialTreeNode> root, std::ostream &stream) {
    
    for(size_t i =0; i < postorder.size(); i++){
        
        stream << postorder[i]->index <<"("<< postorder[i]->number_nodes_cluster<< ")";
        if (i<(postorder.size()-1))
            stream <<", ";
    }
}
void State::initIdsNextCoalEvents(int numClones){
    for(size_t i =0; i< numClones; i++ ){
        
        idsNextCoalEvents.push_back(0);
    }
    
}
void State::initIdsInmigrationEvents(int numClones){
    for(size_t i =0; i< numClones; i++ ){
        
        idsNextInmigrationEvents.push_back(0);
    }
    
}
int State::getIdNextCoalEventForPopulation(int i){
    return idsNextCoalEvents[i];
}
int State::getIdNextInmigrationEventForPopulation(int i){
    return idsNextInmigrationEvents[i];
}
void State::moveNextIdEventForPopulation(int i, bool thereIsInmigration){
    
    idsNextCoalEvents[i]= idsNextCoalEvents[i]+1;
    if (thereIsInmigration)
        idsNextInmigrationEvents[i]= idsNextInmigrationEvents[i]+1;
}

void State::addRoot(std::shared_ptr<PartialTreeNode> node){
    
    roots.push_back(node);
    addNodeToPostorderByIndex(roots.size()-1);
    
}
void State::addNodeToPostorderByIndex(int idx){
    
    postorder.push_back(roots[idx]);
    
}
double State::getNextCoalTime(gsl_rng *random, int& idxLeftNodePop, int& idxRightNodePop, double &logLik, double K){
    
    double result = populationSet->proposeNextCoalEventTime(  random,  idxLeftNodePop,  idxRightNodePop, logLik,  K);
    assert(result >0);
    return result;
    
}
double   State::proposalPriorPost(gsl_rng * random, Population *leftNodePop,Population *rightNodePop,  double newNodeHeight, double logLikNewHeight){
    
    if (leftNodePop->index == rightNodePop->index)
        return proposalCoalNodePriorPost( random, leftNodePop,  newNodeHeight,  logLikNewHeight);
    else{
        
        Population *inmigrantPop = (leftNodePop->numActiveGametes == 1)?leftNodePop:rightNodePop;
        Population *receiverPop = (leftNodePop->numActiveGametes == 1)?rightNodePop:leftNodePop;
        assert(inmigrantPop->numActiveGametes==1);
        return proposalCoalMRCANodePriorPost(random, inmigrantPop, receiverPop,  newNodeHeight,  logLikNewHeight);
    
    }
}
double   State::proposalCoalMRCANodePriorPost(gsl_rng * random, Population *inmigrantPop,Population *receiverPop, double newNodeHeight, double logLikNewHeight){
    
    int idxFirst, idxSecond, idxFirstRoot, idxSecondRoot;
    double logWeightDiff;
    int numberProposals= receiverPop->idsActiveGametes.size();
    std::vector<double> logWeights(numberProposals, 0.0);
    std::vector<double> normWeights(numberProposals, 0.0);
    
    std::vector<std::shared_ptr<PartialTreeNode>> nodeProposals(numberProposals);
    
    double maxlogWeight = DOUBLE_NEG_INF;
    size_t posMax;
    
    for(size_t i=0; i <  receiverPop->idsActiveGametes.size(); i++){
        
        idxFirstRoot = inmigrantPop->idsActiveGametes[0];
        idxSecondRoot = receiverPop->idsActiveGametes[i];
        
        // std::cout << "idx first" << idxFirst << " idx second "<<idxSecond << std::endl;
        assert(idxFirstRoot != idxSecondRoot);
        nodeProposals[i] = proposeNewNode( idxFirstRoot,  idxSecondRoot, receiverPop->index, newNodeHeight, logLikNewHeight );
        logWeights[i] = nodeProposals[i]->ln_likelihood;
       
        if (logWeights[i]>maxlogWeight)
        {
            maxlogWeight = logWeights[i];
            posMax = i;
        }
    }
    unsigned int pos;
    Utils::normalize(logWeights, normWeights);
    pos = Random::randomDiscreteFromProbabilityVector(random, &normWeights[0], normWeights.size());
    
    
   
    logWeightDiff = likelihood_factor(nodeProposals[pos]);
       assert(!isnan(logWeightDiff) && !isinf(logWeightDiff));
       logWeight = logWeight+ logWeightDiff;
    
    idxFirst =inmigrantPop->idsActiveGametes[0];//first idx is the inmigrant root
    idxSecond = receiverPop->idsActiveGametes[pos];//second idx inside the list of active gametes of the chosen pop
    
    int posFirstRoot= getNodeIdxById(idxFirstRoot);
    int posSecondRoot = getNodeIdxById(idxSecondRoot);
    
    
    remove_roots(posFirstRoot, posSecondRoot);
    addRoot(nodeProposals[pos]);
    
    
    increaseNextAvailable();
    updateIndexesActiveGametes( idxFirst, idxSecond,  inmigrantPop->index, receiverPop->index,  nodeProposals[pos]->index, nodeProposals[pos]->index_population);
    
     receiverPop->numActiveGametes= receiverPop->numActiveGametes-1;
    setTopHeightScaledByTheta(newNodeHeight);
    return logWeight;
}
double   State::proposalCoalNodePriorPost(gsl_rng * random, Population *chosenPop, double newNodeHeight, double logLikNewHeight){
    
    int idxFirst, idxSecond, idxFirstRoot, idxSecondRoot;
    double logWeightDiff = 0.0;
    std::vector<std::pair<int, int>> allCoalPairs = Utils::allCombinations(chosenPop->numActiveGametes, 2);
    int numberProposals = allCoalPairs.size();
    
    std::vector<double> logLik(numberProposals, 0.0);
    std::vector<double> normlogLik(numberProposals, 0.0);
    std::vector<double> logRatio(numberProposals, 0.0);
      std::vector<double> normlogRatio(numberProposals, 0.0);
    
    std::vector<std::shared_ptr<PartialTreeNode>> nodeProposals(numberProposals);
    
    double maxlogLik = DOUBLE_NEG_INF;
    double maxlogRatio = DOUBLE_NEG_INF;
    size_t posMaxLik = 0, posMaxRatio = 0;
    for(size_t i=0; i <  allCoalPairs.size(); i++){
        idxFirst = allCoalPairs[i]. first;
        idxSecond = allCoalPairs[i]. second;
        
        idxFirstRoot = chosenPop->idsActiveGametes[idxFirst];
        idxSecondRoot = chosenPop->idsActiveGametes[idxSecond];
        
        // std::cout << "idx first" << idxFirst << " idx second "<<idxSecond << std::endl;
        assert(idxFirstRoot != idxSecondRoot);
        nodeProposals[i] = proposeNewNode( idxFirstRoot,  idxSecondRoot, chosenPop->index, newNodeHeight, logLikNewHeight );
        
        logLik[i] = nodeProposals[i]->ln_likelihood;
        logRatio[i] = likelihood_factor(nodeProposals[i]);
    
        if (logLik[i]>maxlogLik)
        {
            maxlogLik = logLik[i];
            posMaxLik = i;
        }
        if (logRatio[i]>maxlogRatio)
               {
                   maxlogRatio = logRatio[i];
                   posMaxRatio = i;
               }
    }
    unsigned int posRatios, posLik, pos;
    double logSumNormLik = Utils::normalize(logLik, normlogLik);
    double logSumNormRatios = Utils::normalize(logRatio, normlogRatio);
    posLik = Random::randomDiscreteFromProbabilityVector(random, &normlogLik[0], normlogLik.size());
    posRatios = Random::randomDiscreteFromProbabilityVector(random, &normlogRatio[0], normlogRatio.size());
    
    if (posLik != posRatios){
  
        std::cout<< "posLik " << posLik<< " norm logLik " << normlogLik[posLik] << std::endl;
         std::cout<<"parent "<< nodeProposals[posLik]->index<<  " left " << nodeProposals[posLik]->getIndexLeftChild() << " right " << nodeProposals[posLik]->getIndexRightChild() <<std::endl;
        std::cout<< "max loglik pair " <<  nodeProposals[posMaxLik]->getIndexLeftChild() << " second " << nodeProposals[posMaxLik]->getIndexRightChild() <<std::endl;
        
    
        std::cout<< "posRatios " << posRatios<< " norm logRatio " <<normlogRatio[posRatios] << std::endl;
        std::cout<<"parent "<< nodeProposals[posRatios]->index<< " left " << nodeProposals[posRatios]->getIndexLeftChild() << " right " << nodeProposals[posRatios]->getIndexRightChild() <<std::endl;
               std::cout<< "max logRatio pair " << nodeProposals[posMaxRatio]->getIndexLeftChild() << " second " << nodeProposals[posMaxRatio]->getIndexRightChild() <<std::endl;
        std::cout<< "\n";
    }
        
    pos = posRatios;
    logWeightDiff = likelihood_factor(nodeProposals[pos]);
    logWeightDiff = logSumNormRatios;
    assert(!isnan(logWeightDiff) && !isinf(logWeightDiff));
    logWeight = logWeight+ logWeightDiff;
    

    idxFirst = allCoalPairs[pos]. first;//first idx inside the list of active gametes of the chosen pop
    idxSecond = allCoalPairs[pos]. second;//second idx inside the list of active gametes of the chosen pop
    if (true){
         std::cout<<"parent "<< nodeProposals[pos]->index<< " left " << nodeProposals[pos]->getIndexLeftChild() << " right " << nodeProposals[pos]->getIndexRightChild() <<std::endl;

    }
    
    // std::cout<< "Max weight: first " << allCoalPairs[posMax]. first << "second " << allCoalPairs[posMax]. second << " weight " << weight << " max weight " << logWeights[posMax] << std::endl;
    idxFirstRoot = chosenPop->idsActiveGametes[idxFirst];
    idxSecondRoot = chosenPop->idsActiveGametes[idxSecond];
    
    assert(idxFirstRoot != idxSecondRoot);
    
    int posFirstRoot= getNodeIdxById(idxFirstRoot);
    int posSecondRoot = getNodeIdxById(idxSecondRoot);
    
    remove_roots(posFirstRoot, posSecondRoot);
    assert(nodeProposals[pos]->number_leaves_cluster>1);
    addRoot(nodeProposals[pos]);
    
    increaseNextAvailable();
    updateIndexesActiveGametes( idxFirst, idxSecond,  chosenPop->index, chosenPop->index,  nodeProposals[pos]->index, nodeProposals[pos]->index_population);
    
    chosenPop->numActiveGametes= chosenPop->numActiveGametes-1;
    setTopHeightScaledByTheta(newNodeHeight);
    return logWeight;
}
double  State::proposalPriorPrior(gsl_rng * random, Population *leftNodePop,Population *rightNodePop, double newNodeHeight, double logLikNewHeight){
    
    if (leftNodePop->index == rightNodePop->index)
        return proposalCoalNodePriorPrior( random, leftNodePop,  newNodeHeight,  logLikNewHeight);
    else{
        
        Population *inmigrantPop = (leftNodePop->numActiveGametes == 1)?leftNodePop:rightNodePop;
        Population *receiverPop = (leftNodePop->numActiveGametes == 1)?rightNodePop:leftNodePop;
        assert(inmigrantPop->numActiveGametes==1);
        return proposalCoalMRCANodePriorPrior(random, inmigrantPop, receiverPop,  newNodeHeight,  logLikNewHeight);
    
    }
    
    
}
double  State::proposalPostPost(gsl_rng * random,  double &newHeight , double& logLikNewHeight, int K){
    
    double result = 0.0;
    
    
    
    
    return result;
}
double State::proposalCoalNodePostPost(gsl_rng * random, Population *chosenPop, double newNodeHeight, double logLikNewHeight){
    
     
       int idxFirst, idxSecond, idxFirstRoot, idxSecondRoot;
       double logWeightDiff = 0.0;
       std::vector<std::pair<int, int>> allCoalPairs = Utils::allCombinations(chosenPop->numActiveGametes, 2);
       int numberProposals = allCoalPairs.size();
       
       std::vector<double> logWeights(numberProposals, 0.0);
       std::vector<double> normWeights(numberProposals, 0.0);
       
       std::vector<std::shared_ptr<PartialTreeNode>> nodeProposals(numberProposals);
       
       double maxlogWeight = DOUBLE_NEG_INF;
       size_t posMax;
       for(size_t i=0; i <  allCoalPairs.size(); i++){
           idxFirst = allCoalPairs[i]. first;
           idxSecond = allCoalPairs[i]. second;
           
           idxFirstRoot = chosenPop->idsActiveGametes[idxFirst];
           idxSecondRoot = chosenPop->idsActiveGametes[idxSecond];
           
           // std::cout << "idx first" << idxFirst << " idx second "<<idxSecond << std::endl;
           assert(idxFirstRoot != idxSecondRoot);
           nodeProposals[i] = proposeNewNode( idxFirstRoot,  idxSecondRoot, chosenPop->index, newNodeHeight, logLikNewHeight );
           
           logWeights[i] = nodeProposals[i]->ln_likelihood;
       
           if (logWeights[i]>maxlogWeight)
           {
               maxlogWeight = logWeights[i];
               posMax = i;
           }
       }
       unsigned int pos;
       Utils::normalize(logWeights, normWeights);
       pos = Random::randomDiscreteFromProbabilityVector(random, &normWeights[0], normWeights.size());
       
      logWeightDiff = likelihood_factor(nodeProposals[pos]);
       assert(!isnan(logWeightDiff) && !isinf(logWeightDiff));
       logWeight = logWeight+ logWeightDiff;
       
    
       
       idxFirst = allCoalPairs[pos]. first;//first idx inside the list of active gametes of the chosen pop
       idxSecond = allCoalPairs[pos]. second;//second idx inside the list of active gametes of the chosen pop
       if (0){
       std::cout<< "first " << idxFirst << " second " << idxSecond << " time "<< newNodeHeight<< " weight " << logWeight << " norm weight " << normWeights[pos] << std::endl;
       }
       
       // std::cout<< "Max weight: first " << allCoalPairs[posMax]. first << "second " << allCoalPairs[posMax]. second << " weight " << weight << " max weight " << logWeights[posMax] << std::endl;
       idxFirstRoot = chosenPop->idsActiveGametes[idxFirst];
       idxSecondRoot = chosenPop->idsActiveGametes[idxSecond];
       
       assert(idxFirstRoot != idxSecondRoot);
       
       int posFirstRoot= getNodeIdxById(idxFirstRoot);
       int posSecondRoot = getNodeIdxById(idxSecondRoot);
       
       remove_roots(posFirstRoot, posSecondRoot);
       assert(nodeProposals[pos]->number_leaves_cluster>1);
       addRoot(nodeProposals[pos]);
       
       increaseNextAvailable();
       updateIndexesActiveGametes( idxFirst, idxSecond,  chosenPop->index, chosenPop->index,  nodeProposals[pos]->index, nodeProposals[pos]->index_population);
       
       chosenPop->numActiveGametes= chosenPop->numActiveGametes-1;
       setTopHeightScaledByTheta(newNodeHeight);
       return logWeight;
}
  double   proposalCoalMRCANodePostPost(gsl_rng * random, Population *inmigrantPop,Population *receiverPop, double newNodeHeight, double logLikNewHeight){
      
      double result = 0.0;
      
      
      return result;
  }
double  State::proposalCoalNodePriorPrior(gsl_rng * random, Population *chosenPop, double newNodeHeight,  double logLikNewHeight){
    int idxFirst, idxSecond, idxFirstRoot, idxSecondRoot;
    double logWeightDiff;
    int choosePairIndividuals = YES;
    if (getNumberPopulations()==1)
        assert(chosenPop->numActiveGametes == root_count());
    
    chosenPop->ChooseRandomIndividual(&idxFirst, getNumberPopulations(),   &idxSecond, random, choosePairIndividuals);
    
    idxFirstRoot = chosenPop->idsActiveGametes[idxFirst];
    idxSecondRoot = chosenPop->idsActiveGametes[idxSecond];
    

    assert(idxFirstRoot != idxSecondRoot);
    
    setTopHeightScaledByTheta(newNodeHeight);
     
    std::shared_ptr<PartialTreeNode> node = connect(idxFirstRoot, idxSecondRoot, chosenPop->index, newNodeHeight, logLikNewHeight );
    
    updateIndexesActiveGametes( idxFirst, idxSecond,  chosenPop->index, chosenPop->index,  node->index, node->index_population);
    
    chosenPop->numActiveGametes= chosenPop->numActiveGametes-1;
    
    logWeightDiff = likelihood_factor(node);
    assert(!isnan(logWeightDiff) && !isinf(logWeightDiff));
    logWeight = logWeight+ logWeightDiff;
    

     std::cout<< "first " << idxFirstRoot << " second " << idxSecondRoot << " time "<< newNodeHeight<< " weight " << logWeight  << std::endl;
    return logWeight;
}
double  State::proposalCoalMRCANodePriorPrior(gsl_rng * random, Population *inmigrantPop,Population *receiverPop, double newNodeHeight, double logLikNewHeight){
    int idxFirst, idxSecond, idxFirstRoot, idxSecondRoot;
    double logWeightDiff = 0.0;
    int choosePairIndividuals = NO;
 
    receiverPop->ChooseRandomIndividual(&idxFirst, getNumberPopulations(),   &idxSecond, random, choosePairIndividuals);
    
    idxFirstRoot = receiverPop->idsActiveGametes[idxFirst];
    idxSecondRoot = inmigrantPop->idsActiveGametes[0];
    
    setTopHeightScaledByTheta(newNodeHeight);
    assert(idxFirstRoot != idxSecondRoot);
    
    std::shared_ptr<PartialTreeNode> node = connect(idxFirstRoot, idxSecondRoot, receiverPop->index, newNodeHeight, logLikNewHeight );
    
    updateIndexesActiveGametes( idxFirst, idxSecond,  receiverPop->index, receiverPop->index,  node->index, node->index_population);
    
    receiverPop->numActiveGametes= receiverPop->numActiveGametes-1;
    
    logWeightDiff = likelihood_factor(node);
    assert(!isnan(logWeightDiff) && !isinf(logWeightDiff));
    logWeight = logWeight+ logWeightDiff;

  
    return logWeight;
}
unsigned int State::numberNonTrivialTrees(){
    unsigned int result = 0;
    for( int i=0; i< roots.size();++i){
        if (roots[i]->number_leaves_cluster >1)
            result++;
    }
    return(result);
}
State::~State(){
    
    
}
