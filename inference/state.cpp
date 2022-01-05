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
#include <Eigen/Core>

#include "poset_smc_params.hpp"
#include "poset_smc.hpp"
#include "ccars.hpp"
#include "gnu_plotter.hpp"

State::State(PosetSMCParams &smcParams, gsl_rng *random):
pll_buffer_manager(smcParams.pll_buffer_manager),partition(smcParams.partition), num_sites(smcParams.num_sites)
{
    
    
    topHeightModelTime =0.0;
    topHeightScaledByTheta = 0.0;
    logWeight = 0.0;
    theta = smcParams.theta;
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
    
    
    initForest(smcParams.sampleSize, smcParams.msa->getNonConstRawPtr(), smcParams.positions,  smcParams.getProgramOptions());
    nextAvailable = smcParams.sampleSize-1;
    
    if (smcParams.usePriorInSMC1){
        pairMinModelTime = populationSet->initPairProposalsNextEventTimePerPopulation(random, smcParams.getProgramOptions().K);
    }
    else{
        pairMinModelTime = initPairProposalsNextEventTime(random, smcParams.getProgramOptions().normalizeClv, smcParams.getProgramOptions().K);
    }
    
    
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
    num_sites= original.num_sites;
    
    populationSet=new PopulationSet(*(original.populationSet));
    gtError = new GenotypeErrorModel(*(original.gtError));
    
    //populationSet= std::move(original.populationSet);
    // gtError = std::move(original.gtError);
    
    //    idsNextCoalEvents = original.idsNextCoalEvents;
    //    idsNextInmigrationEvents = original.idsNextInmigrationEvents;
    
    //idsNextCoalEvents = std::move(original.idsNextCoalEvents);
    //idsNextInmigrationEvents = std::move( original.idsNextInmigrationEvents);
    
    roots = original.roots;
    postorder  = original.postorder;
    coalEventTimesScaledByTheta = original.coalEventTimesScaledByTheta;
    
    //roots = std::move(original.roots);
    //postorder  = std::move(original.postorder);
    coalEventTimesScaledByTheta = std::move(original.coalEventTimesScaledByTheta);
    //    roots.reserve(original.root_count());
    //    for (auto const& fptr : original.getRoots())
    //         roots.emplace_back(fptr->Clone());
    
    
    partition = original.partition;
    nextAvailable = original.nextAvailable;
    theta = original.theta;
    logWeight = original.logWeight;
    pairMinModelTime =  original.pairMinModelTime;
    
}
State::State(State &&rhs):
pll_buffer_manager(rhs.pll_buffer_manager){
    
    topHeightModelTime = std::move(rhs.topHeightModelTime);
    topHeightScaledByTheta = std::move(rhs.topHeightScaledByTheta);
    num_sites= std::move(rhs.num_sites);
    
    populationSet= std::move(rhs.populationSet);
    gtError = std::move(rhs.gtError);
    
    //idsNextCoalEvents = std::move(rhs.idsNextCoalEvents);
    //idsNextInmigrationEvents = std::move(rhs.idsNextInmigrationEvents);
    roots = std::move(rhs.roots);
    postorder  = std::move(rhs.postorder);
    coalEventTimesScaledByTheta =std::move( rhs.coalEventTimesScaledByTheta);
    
    
    partition = std::move(rhs.partition);
    nextAvailable = std::move(rhs.nextAvailable);
    theta = std::move(rhs.theta);
    logWeight = std::move(rhs.logWeight);
    pairMinModelTime =  std::move(rhs.pairMinModelTime);
    
}
void State::initForest( int sampleSize, pll_msa_t *msa, std::vector<int> &positions, ProgramOptions &programOptions){
    
    //TODO receive a map from tip to population and use it in the constructor of
    //PartialTreeNode instead of 0
    const Partition *reference_partition = partition;
    
    // std::shared_ptr<PartialTreeNode>  p;
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
    
    auto delLeavePartialTreeNode = [](PartialTreeNode* pTreeNode)
    {
        PLLBufferManager *manager = pTreeNode->getManager();
        manager->clv_buffer.push(pTreeNode->pclv);
        manager->scale_buffer_buffer.push(pTreeNode->pscale_buffer);
        
        pTreeNode->pclv = nullptr;
        pTreeNode->pscale_buffer = nullptr;
        
        delete pTreeNode;
    };
    
    
    for(size_t i=0; i< positions.size();++i)
    {
        
        cumIndivid++;
        
        unsigned int sites_alloc =
        reference_partition->ascBiasCorrection() ? reference_partition->numberSites + reference_partition->numberStates : reference_partition->numberSites;
        
        unsigned int clv_num_elements= sites_alloc * reference_partition->numberStates * reference_partition->numberRateCats;
        // unsigned int clv_size =
        // clv_num_elements * sizeof(double);
        unsigned int scaler_size = (reference_partition->attributes & PLL_ATTRIB_RATE_SCALERS)
        ? sites_alloc * reference_partition->numberRateCats
        : sites_alloc;
        std::shared_ptr<PartialTreeNode> p( new PartialTreeNode(pll_buffer_manager,
                                                                nullptr, nullptr,0.0, 0.0, 0, msa->label[i], 0.0, clv_num_elements,
                                                                scaler_size, reference_partition->alignment(), i, 0, nullptr, nullptr), delLeavePartialTreeNode);
        
        //delete (p->pclv);
        //delete (p->pscale_buffer);
        // auto clv_elem = reference_partition->numberSites * reference_partition->numberStates;
        p->buildCLV(i,reference_partition->numberStates,  msa, gtError,  programOptions.normalizeLeavesClv);
        
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
        roots.push_back(std::move(p));
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
    //idsNextCoalEvents = original.idsNextCoalEvents;
    //idsNextInmigrationEvents = original.idsNextInmigrationEvents;
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
                                                                                pll_buffer_manager, child_left, child_right,
                                                                                left_length, right_length,
                                                                                pmatrix_elements,
                                                                                "", newNodeHeight, clv_elements,
                                                                                scaler_size, partition->alignment(), nextAvailable,
                                                                                index_pop_new_node, nullptr, nullptr);
    
    
    
    
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
                               parent->edge_l->pmatrix,
                               parent->edge_r->pmatrix,
                               child_left->pscale_buffer,
                               child_right->pscale_buffer,
                               p->attributes);
    
    parent->ln_likelihood =
    compute_ln_likelihood(parent->pclv,
                          parent->pscale_buffer,
                          p->getPartition());
    
    
     if (1)
           parent-> showpClV(p->numberStates, p->numberRateCats, p->statesPadded, sites, 3);
    
    if(isnan(parent->ln_likelihood) || isinf(-1*parent->ln_likelihood)){
        
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
std::shared_ptr<PartialTreeNode> State::proposeNewNode(int firstId, int secondId, size_t index_pop_new_node, double newNodeHeight, double logLikNewNode, double* left_pmatrix, double *right_pmatrix, bool normalize){
    
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
    
    
    double left_length = newNodeHeight - child_left->height;
    double right_length = newNodeHeight - child_right->height;
    
    unsigned int pmatrix_elements = p->numberStates * p->numberStates  * p->numberRateCats;
    
    unsigned int sites_alloc =
    p->ascBiasCorrection() ? p->numberSites + p->numberStates : p->numberSites;
    
    unsigned int clv_elements = sites_alloc  * p->numberStates *  p->numberRateCats;
    
    // unsigned int clv_size = clv_elements * sizeof(double);
    //sites_alloc * p->getStatesPadded() * p->numberRateCats * sizeof(double);
    
    unsigned int scaler_size = (p->attributes & PLL_ATTRIB_RATE_SCALERS)
    ? sites_alloc * p->numberRateCats
    : sites_alloc;
    // unsigned int scale_buffer_size = scaler_size * sizeof(double);
    
    auto delPartialTreeNode = [](PartialTreeNode* pTreeNode)
    {
        PLLBufferManager *manager = pTreeNode->getManager();
        manager->clv_buffer.push(pTreeNode->pclv);
        manager->scale_buffer_buffer.push(pTreeNode->pscale_buffer);
        
        pTreeNode->pclv = nullptr;
        pTreeNode->pscale_buffer = nullptr;
        
        manager->pmatrix_buffer.push(pTreeNode->edge_l->pmatrix);
        manager->pmatrix_buffer.push(pTreeNode->edge_r->pmatrix);
        
        // pll_aligned_free(pmatrix);
        pTreeNode->edge_l->pmatrix = nullptr;
        pTreeNode->edge_r->pmatrix = nullptr;
        delete pTreeNode;
    };
    
    std::shared_ptr<PartialTreeNode> parent( new PartialTreeNode(
                                                                 pll_buffer_manager, child_left, child_right,
                                                                 left_length, right_length,
                                                                 pmatrix_elements,
                                                                 "", newNodeHeight, clv_elements,
                                                                 scaler_size, partition->alignment(), nextAvailable,
                                                                 index_pop_new_node,
                                                                 left_pmatrix,
                                                                 right_pmatrix), delPartialTreeNode);
    
    const unsigned int matrix_indices[1] = {0};
    //Array params_indices holds the indices of rate matrices that will be used for each rate category. This array must be of size rate_cats
    
    const unsigned int param_indices[1] = {0};
    
    int left_edge_pmatrix_result;
    
 
    if (left_pmatrix == nullptr){
  
        left_edge_pmatrix_result = pll_core_update_pmatrix(
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
        

//           parent->edge_l->showPMatrix( p->numberStates, p->numberRateCats, p->statesPadded,  3);
                          
    }
    if (right_pmatrix == nullptr){
        int right_edge_pmatrix_result;
      
        right_edge_pmatrix_result = pll_core_update_pmatrix(&parent->edge_r->pmatrix, p->numberStates,     p->numberRateCats, p->rates(),
                                                            &parent->edge_r->length, matrix_indices, param_indices, p->propInvar(),
                                                            p->eigenVals(), p->eigenVecs(), p->invEigenVecs(), 1, p->attributes);
        //   }
        assert(right_edge_pmatrix_result == PLL_SUCCESS);
        // parent->edge_r->showPMatrix( p->numberStates, p->numberRateCats, p->statesPadded,  3);
    }
    
    unsigned int sites= (p->ascBiasCorrection() ? p->numberSites + p->numberStates : p->numberSites);
    
    if (0){
        std::cout << "left clv" << std::endl;
        child_left-> showpClV(p->numberStates, p->numberRateCats, p->statesPadded, sites, 3);
        std::cout << "right clv" << std::endl;
        child_right-> showpClV(p->numberStates, p->numberRateCats, p->statesPadded, sites, 3);
    }
    
    pll_core_update_partial_ii_norm(p->numberStates, sites,
                               p->numberRateCats, parent->pclv,
                               parent->pscale_buffer,
                               child_left->pclv,
                               child_right->pclv,
                               parent->edge_l->pmatrix,
                               parent->edge_r->pmatrix,
                               child_left->pscale_buffer,
                               child_right->pscale_buffer,
                               p->attributes,
                                    normalize);
       if (0)
           parent-> showpClV(p->numberStates, p->numberRateCats, p->statesPadded, sites, 3);
    
    parent->ln_likelihood =
    compute_ln_likelihood(parent->pclv,
                          parent->pscale_buffer,
                          p->getPartition());
    
   
    
    if (0)
        parent-> showpClV(p->numberStates, p->numberRateCats, p->statesPadded, sites, 3);
    
    std::vector<double> result(num_sites);
    
    if(isnan(parent->ln_likelihood) || isinf(-1*parent->ln_likelihood)){
        
        std::cout << "Error: log lik  "<<  " is "<< parent->ln_likelihood<< std::endl;
    }
    parent->ln_coal_likelihood = (child_left->height >  child_right->height)? child_left->ln_coal_likelihood: child_right->ln_coal_likelihood;
    
    parent->ln_coal_likelihood+= logLikNewNode;
    
    
    assert(!isnan(parent->ln_likelihood) && !isinf(parent->ln_likelihood));
    
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
double State::likelihood_factor(std::unique_ptr<PartialTreeNode> root)const{
    
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
    
        double result = pll_core_root_loglikelihood2(  p->numberStates, p->numberSites, p->numberRateCats,
                                                     clv.data(), scale_buffer.data(),
                                                     p->frequencies(), p->rateWeights(), p->patternWeights(), parameter_indices,  nullptr, p->attributes);
//    double result=  pll_core_root_loglikelihood(
//                                                p->numberStates, p->numberSites, p->numberRateCats,
//
//                                                clv.data(), scale_buffer.data(),
//
//                                                p->frequencies(), p->rateWeights(), p->patternWeights(), p->propInvar(),
//                                                p->invariant(), parameter_indices, nullptr, p->attributes);
    return result;
}
double State::compute_ln_likelihood(double *pclv, unsigned int *pscale_buffer,
                                    const Partition *p) {
    const unsigned int parameter_indices[1] = {0};
        double result = pll_core_root_loglikelihood2(  p->numberStates, p->numberSites, p->numberRateCats,
                                                     pclv, pscale_buffer,
                                                     p->frequencies(), p->rateWeights(), p->patternWeights(), parameter_indices,  nullptr, p->attributes);
//    double result=  pll_core_root_loglikelihood(
//                                                p->numberStates, p->numberSites, p->numberRateCats,
//
//                                                pclv, pscale_buffer,
//
//                                                p->frequencies(), p->rateWeights(), p->patternWeights(), p->propInvar(),
//                                                p->invariant(), parameter_indices, nullptr, p->attributes);
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
    
    
    double result = pll_core_root_loglikelihood2(  p->states, p->sites, p->rate_cats,
                                                 pclv, pscale_buffer,
                                                 p->frequencies, p->rate_weights, p->pattern_weights, parameter_indices,  nullptr, p->attributes);
    
    //    double result= pll_core_root_loglikelihood(
    //                                               p->states, p->sites, p->rate_cats,
    //                                               pclv, pscale_buffer,
    //                                               p->frequencies, p->rate_weights, p->pattern_weights, p->prop_invar,
    //                                               p->invariant, parameter_indices, nullptr, p->attributes);
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
double State::getNextCoalTime(gsl_rng *random, int& idxReceiverPop, int& idxIncomingPop, double &logLik, double K){
    
    int idxIncomingNode;
    double proposalCoalTimeScaledByProportion = populationSet->proposeNextCoalEventTime(  random,  idxReceiverPop,  idxIncomingPop, idxIncomingNode, logLik,  K);
    assert(proposalCoalTimeScaledByProportion >0);
    if (idxReceiverPop!=idxIncomingPop){
        Population *receiverPop= populationSet->getPopulationbyIndex(idxReceiverPop);
        receiverPop->idsActiveGametes[receiverPop->numActiveGametes] = idxIncomingNode;
        receiverPop->numActiveGametes = receiverPop->numActiveGametes + 1;
    }
    return proposalCoalTimeScaledByProportion;
    
}
void State::acceptNextCoalTimeProposal(double nextCoaTimeScaledByProportion, gsl_rng *random, int idxReceiverPop, int idxIncomingPop,  double K)
{
    double proportion= populationSet->getPopulationbyIndex( idxReceiverPop)->x;
    populationSet->acceptNextCoalEventTime(random, nextCoaTimeScaledByProportion/proportion,  idxReceiverPop,  idxIncomingPop, K);
    
    
    
}
Population* State::getCurrentPopulation(){
    
    Population* result = populationSet->getCurrentPopulation();
    return result;
    
}
double   State::proposalPriorPost(gsl_rng * random, Population *leftNodePop,Population *rightNodePop,  double newNodeHeight, double logLikNewHeight, size_t iteration){
    
    assert(logLikNewHeight!= 0);
    
    
    
    if (leftNodePop->index == rightNodePop->index){
        //const Partition *p =partition;
        //        unsigned int pmatrix_elements = p->numberStates * p->numberStates  * p->numberRateCats;
        //        const unsigned int matrix_indices[2] = {0,1};
        //        const unsigned int param_indices[2] = {0,1};
        //        double * left_pmatrix;
        //        double * right_pmatrix;
        if (iteration ==1){//the same pmatrix for every pair
            
            // left_pmatrix = (double *)pll_aligned_alloc(pmatrix_elements * sizeof(double), p->alignment()); //(double *)std::malloc(pmatrix_elements*sizeof(double));
            // right_pmatrix = (double *)pll_aligned_alloc(pmatrix_elements * sizeof(double), p->alignment()); //(double *)std::malloc(pmatrix_elements*sizeof(double));
            //                   double * matrices[2] = { left_pmatrix, right_pmatrix };
            //                   double branch_lengths[2] = {newNodeHeight, newNodeHeight};
            //            int pmatrix_result;
            //            // if (p->getModelName() == "GT16JC")
            //               pmatrix_result= pll_core_update_pmatrix_16x16_jc69( matrices, p->numberStates,p->numberRateCats, p->rates(), branch_lengths,matrix_indices, param_indices,
            //                1,
            //                p->attributes);
            //            // else{
            ////                 pmatrix_result = pll_core_update_pmatrix(
            ////                                                                        matrices, p->numberStates, p->numberRateCats, p->rates(),
            ////                                                                        branch_lengths, matrix_indices, param_indices, p->propInvar(),
            ////                                                                        p->eigenVals(), p->eigenVecs(), p->invEigenVecs(), 1, p->attributes);
            //             //    }
            double result = proposalCoalNodePriorPost( random, leftNodePop,  newNodeHeight,  logLikNewHeight, iteration,  nullptr);
            //            pll_aligned_free(left_pmatrix);
            //            pll_aligned_free(right_pmatrix);
            return result;
        }
        else
            return proposalCoalNodePriorPost( random, leftNodePop,  newNodeHeight,  logLikNewHeight, iteration,  nullptr);
    }
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
        nodeProposals[i] = proposeNewNode( idxFirstRoot,  idxSecondRoot, receiverPop->index, newNodeHeight, logLikNewHeight, nullptr, nullptr, false);
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
    // logWeight = logWeight+ logWeightDiff;
    
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
double   State::proposalCoalNodePriorPost(gsl_rng * random, Population *chosenPop, double newNodeHeight, double logLikNewHeight,size_t iteration, double** pmatrices){
    
    int idxFirst, idxSecond, idxFirstRoot, idxSecondRoot;
    double logWeightDiff = 0.0;
    
    std::vector<std::pair<int, int>> allCoalPairs = PosetSMC::getCombinationsRandomOrder(chosenPop->numActiveGametes);
    int numberProposals = allCoalPairs.size();
    
    std::vector<double> logLik(numberProposals, 0.0);
    std::vector<double> normlogLik(numberProposals, 0.0);
    std::vector<double> logRatio(numberProposals, 0.0);
    std::vector<double> normlogRatio(numberProposals, 0.0);
    
    std::vector<std::shared_ptr<PartialTreeNode>> nodeProposals(numberProposals);
    
    double maxlogLik = DOUBLE_NEG_INF;
    double maxlogRatio = DOUBLE_NEG_INF;
    size_t posMaxLik = 0, posMaxRatio = 0;
    //for(size_t i=0; i <  allCoalPairs.size(); i++){
    size_t i = 0;
    std::vector<pairs> allPairs(allCoalPairs.size());
    for (auto &coalPair : allCoalPairs)
    {
        idxFirst = coalPair. first;
        idxSecond = coalPair. second;
        
        idxFirstRoot = chosenPop->idsActiveGametes[idxFirst];
        idxSecondRoot = chosenPop->idsActiveGametes[idxSecond];
        
        allPairs[i] = coalPair;
        // std::cout << "idx first" << idxFirst << " idx second "<<idxSecond << std::endl;
        assert(idxFirstRoot != idxSecondRoot);
        if (pmatrices!=nullptr)
            nodeProposals[i] = proposeNewNode( idxFirstRoot,  idxSecondRoot, chosenPop->index, newNodeHeight, logLikNewHeight, pmatrices[0], pmatrices[1] , false );
        else
            nodeProposals[i] = proposeNewNode( idxFirstRoot,  idxSecondRoot, chosenPop->index, newNodeHeight, logLikNewHeight, nullptr, nullptr, false );
        logLik[i] = nodeProposals[i]->ln_likelihood;
        logRatio[i] = nodeProposals[i]->likelihood_factor();
        
        
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
        i++;
    }
    unsigned int posRatios;
    //unsigned int posLik;
    unsigned int pos;
    //double logSumNormLik = Utils::normalize(logLik, normlogLik);
    double logSumNormRatios = Utils::normalize(logRatio, normlogRatio);
    // posLik = Random::randomDiscreteFromProbabilityVector(random, &normlogLik[0], normlogLik.size());
    posRatios = Random::randomDiscreteFromProbabilityVector(random, &normlogRatio[0], normlogRatio.size());
    
    // if (posLik != posRatios){
    
    //        std::cout<< "posLik " << posLik<< " norm logLik " << normlogLik[posLik] << std::endl;
    //        std::cout<<"parent "<< nodeProposals[posLik]->index<<  " left " << nodeProposals[posLik]->getIndexLeftChild() << " right " << nodeProposals[posLik]->getIndexRightChild() <<std::endl;
    //        std::cout<< "max loglik pair " <<  nodeProposals[posMaxLik]->getIndexLeftChild() << " second " << nodeProposals[posMaxLik]->getIndexRightChild() <<std::endl;
    //
    //
    //        std::cout<< "posRatios " << posRatios<< " norm logRatio " <<normlogRatio[posRatios] << std::endl;
    //        std::cout<<"parent "<< nodeProposals[posRatios]->index<< " left " << nodeProposals[posRatios]->getIndexLeftChild() << " right " << nodeProposals[posRatios]->getIndexRightChild() <<std::endl;
    //        std::cout<< "max logRatio pair " << nodeProposals[posMaxRatio]->getIndexLeftChild() << " second " << nodeProposals[posMaxRatio]->getIndexRightChild() <<std::endl;
    //        std::cout<< "\n";
    // }
    
    pos = posRatios;
    
    logWeightDiff = logSumNormRatios;
    assert(!isnan(logWeightDiff) && !isinf(logWeightDiff));
    //logWeight = logWeight+ logWeightDiff;
    
    
    idxFirst = allPairs[pos]. first;//first idx inside the list of active gametes of the chosen pop
    idxSecond = allPairs[pos]. second;//second idx inside the list of active gametes of the chosen pop
    if (false){
        std::cout<<"logRatio final proposal: parent "<< nodeProposals[pos]->index<< " left " << nodeProposals[pos]->getIndexLeftChild() << " right " << nodeProposals[pos]->getIndexRightChild() <<std::endl;
        std::cout<< "\n";
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
    return logWeightDiff;
}
double   State::proposalPostPost1(gsl_rng * random,  double &newHeight , double& logLikNewHeight,size_t numIncrementsPOSTPOST, double K){
    
    int idxFirst, idxSecond, idxFirstRoot, idxSecondRoot;
    double logWeightDiff = 0.0;
    int idxReceiverPop,  idxIncomingPop, idxIncomingNode;
    
    Population * currPop = populationSet->getCurrentPopulation();
    std::set<std::pair<int, int>> allCoalPairs = PosetSMC::getCombinations(currPop->numActiveGametes);
    int numberPairs = allCoalPairs.size();
    int numberProposals = numberPairs *numIncrementsPOSTPOST;
    std::vector<double> logWeights(numberProposals);
    std::vector<double> normWeights(numberProposals);
    std::vector<double> logLik(numberProposals);
    std::vector<double> normlogLik(numberProposals);
    std::vector<double> logRatio(numberProposals);
    std::vector<double> normlogRatio(numberPairs);
    
    std::vector<double> heightProposals(numIncrementsPOSTPOST);
    std::vector<double> logPriorHeightProposals(numIncrementsPOSTPOST);
    std::vector<double> logRatioForCurrentIncrement(numberPairs);
    std::vector<double> normLogRatioForCurrentIncrement(numberPairs);
    std::vector<std::shared_ptr<PartialTreeNode>> nodeProposals;
    std::vector<pairs> allPairs(numberProposals);
    
    std::vector<double> logWeightIncrements(numIncrementsPOSTPOST);
    std::vector<double> normLogWeightIncrements(numIncrementsPOSTPOST);
    double logPrior;
    
    newHeight = 0.0;
    logLikNewHeight = 0.0;
    double proposalCoalTimeScaledByProportion = 0.0;
    double  newNodeHeight = 0.0;
    double logLikFactor = 0.0;
    size_t i = 0;
    double normIncrement = 0.0;
    for( size_t j = 0; j< numIncrementsPOSTPOST; ++j){
        
        proposalCoalTimeScaledByProportion = populationSet->proposeNextCoalEventTime(  random,  idxReceiverPop,  idxIncomingPop, idxIncomingNode, logPrior,  K);
        assert(proposalCoalTimeScaledByProportion >0);
        
        if (idxReceiverPop != idxIncomingPop){
            std::cout<< "migration"<< std::endl;
            // do something with theidxIncomingNode
            //idxFirstRoot = currPop->idsActiveGametes[idxIncomingPop];
        }
        newNodeHeight = proposalCoalTimeScaledByProportion* theta;
        heightProposals[j]=proposalCoalTimeScaledByProportion;
        logPriorHeightProposals[j]=logPrior;
        
        //logRatioForCurrentIncrement.clear();
        //normLogRatioForCurrentIncrement.clear();
        i = 0;
        
        
        for (auto &coalPair : allCoalPairs)
        {
            idxFirst = coalPair. first;
            idxSecond = coalPair. second;
            
            idxFirstRoot = currPop->idsActiveGametes[idxFirst];
            idxSecondRoot = currPop->idsActiveGametes[idxSecond];
            
            assert(idxFirstRoot != idxSecondRoot);
            size_t posProposal = j*numIncrementsPOSTPOST+i;
            nodeProposals[posProposal] = proposeNewNode( idxFirstRoot,  idxSecondRoot, currPop->index, newNodeHeight, logPrior, nullptr, nullptr, false );
            
            
            logLik[posProposal] =nodeProposals[posProposal] ->ln_likelihood;
            logLikFactor = likelihood_factor(nodeProposals[posProposal] );
            logRatio[posProposal] =logLikFactor;
            logRatioForCurrentIncrement[i]= logLikFactor;
            
            allPairs[posProposal] = coalPair;
            
            i++;
        }
        logWeightIncrements[j] = Utils::norm(logRatioForCurrentIncrement);
        normIncrement = Utils::logSumExp(logRatioForCurrentIncrement);
    }
    
    unsigned int posPair, pos, posIncrement;
    //select increment
    double logSumWeightIncrements = Utils::normalize(logWeightIncrements, normLogWeightIncrements);
    posIncrement = Random::randomDiscreteFromProbabilityVector(random, &normLogWeightIncrements[0], normLogWeightIncrements.size());
    newHeight = heightProposals[posIncrement] * theta;
    logLikNewHeight = logPriorHeightProposals[posIncrement];
    //now select a pair
    std::vector<double>::const_iterator from = logRatio.begin() + posIncrement*numberPairs;
    std::vector<double>::const_iterator to = logRatio.begin() + (posIncrement+1)*numberPairs-1;
    std::vector<double> logRatioChosenIncrement(from, to);
    
    
    //double logSumNormRatios = Utils::normalize(logRatioChosenIncrement, normlogRatio);
    //assert(logSumNormRatios== logWeightIncrements[posIncrement] );
    //posLik = Random::randomDiscreteFromProbabilityVector(random, &normlogLik[0], normlogLik.size());
    posPair = Random::randomDiscreteFromProbabilityVector(random, &normlogRatio[0], normlogRatio.size());
    
    if (true){
        
        std::cout<< "posPair " << posPair<< " norm logRatio " <<normlogRatio[posPair] << std::endl;
        std::cout<<"parent "<< nodeProposals[posIncrement*numberPairs+posPair]->index<< " left " << nodeProposals[posIncrement*numberPairs+posPair]->getIndexLeftChild() << " right " << nodeProposals[posIncrement*numberPairs+posPair]->getIndexRightChild() << " height " << newHeight << std::endl;
        //        std::cout<< "max logRatio pair " << nodeProposals[posMaxRatio]->getIndexLeftChild() << " second " << nodeProposals[posMaxRatio]->getIndexRightChild() <<std::endl;
        std::cout<< "\n";
    }
    
    pos = posIncrement*numberPairs+posPair;
    assert(pos < nodeProposals.size());
    logWeightDiff = logSumWeightIncrements;
    assert(!isnan(logWeightDiff) && !isinf(logWeightDiff));
    //logWeight = logWeight+ logWeightDiff;
    
    
    idxFirst = allPairs[pos]. first;//first idx inside the list of active gametes of the chosen pop
    idxSecond = allPairs[pos]. second;//second idx inside the list of active gametes of the chosen pop
    if (false){
        std::cout<<"logRatio final proposal: parent "<< nodeProposals[pos]->index<< " left " << nodeProposals[pos]->getIndexLeftChild() << " right " << nodeProposals[pos]->getIndexRightChild() <<std::endl;
        std::cout<< "\n";
    }
    
    idxFirstRoot = currPop->idsActiveGametes[idxFirst];
    idxSecondRoot = currPop->idsActiveGametes[idxSecond];
    
    assert(idxFirstRoot != idxSecondRoot);
    
    int posFirstRoot= getNodeIdxById(idxFirstRoot);
    int posSecondRoot = getNodeIdxById(idxSecondRoot);
    
    remove_roots(posFirstRoot, posSecondRoot);
    assert(nodeProposals[pos]->number_leaves_cluster>1);
    addRoot(nodeProposals[pos]);
    
    increaseNextAvailable();
    updateIndexesActiveGametes( idxFirst, idxSecond,  currPop->index, currPop->index,  nodeProposals[pos]->index, nodeProposals[pos]->index_population);
    
    currPop->numActiveGametes= currPop->numActiveGametes-1;
    setTopHeightScaledByTheta(newHeight);
    return logWeightDiff;
}
double   State::proposalTSMC1(int iter, gsl_rng * random,  double &newHeight , double& logLikNewHeight, double K, bool usePriorInSMC1, bool normalizedCLV){
    
    int idxFirst, idxSecond, idxFirstRoot, idxSecondRoot;
    double logWeightDiff = 0.0;
    double pairModelTimeEntry, pairLifeModelTime;
    double  newNodeHeight = 0.0, newNodeHeightModelTime= 0.0;
    double logPrior = 0.0;

    Population * currPop = populationSet->getCurrentPopulation();
    
    idxFirst = pairMinModelTime.first.first;
    idxSecond = pairMinModelTime.first.second;
    
    idxFirstRoot = currPop->idsActiveGametes[idxFirst];
    idxSecondRoot = currPop->idsActiveGametes[idxSecond];
    
    assert(idxFirstRoot != idxSecondRoot);
    
    ProposalDistribInfo pairModelTimes = pairMinModelTime.second;
    newNodeHeightModelTime  = pairModelTimes.timeProposal  * currPop->x ;
    newNodeHeight = newNodeHeightModelTime * currPop->theta;
    
    pairModelTimeEntry = pairModelTimes.creationTime * currPop->x;
    pairLifeModelTime = newNodeHeightModelTime- pairModelTimeEntry;
    std::shared_ptr<PartialTreeNode> newNode = proposeNewNode( idxFirstRoot,  idxSecondRoot, currPop->index, newNodeHeight, logPrior, nullptr, nullptr, normalizedCLV );
    
    double logLikFactor = newNode->likelihood_factor();
    assert(logLikFactor!= DOUBLE_NEG_INF);
    
    logWeightDiff+= logLikFactor;
    
    if (1){
        
        assert(idxFirstRoot != idxSecondRoot);
                   
        int idxFirstID = getNodeIdxById(idxFirstRoot);
        int idxSecondId = getNodeIdxById(idxSecondRoot);
        
        CCLogDensity *log_density =
             new GenotypeJCPairLogProposal(partition->numberSites, partition->numberStates, currPop->timeOriginSTD, currPop->delta, currPop->theta, pairModelTimeEntry,
                                           roots[idxFirstID]->height,
                                           roots[idxSecondId]->height, roots[idxFirstID]->pclv, roots[idxSecondId]->pclv, normalizedCLV);
          
        double proposal_density =log_density->h_concave_Felsenstein(newNodeHeightModelTime)+ log_density->h_convex_Felsenstein(newNodeHeightModelTime);
        
          std::cout << "winner pair, left " << idxFirstRoot<< " right "<< idxSecondRoot << std::endl;
         std::cout << "proposed time " << newNodeHeight<< std::endl;
          std::cout << "log weight diff  after log Fels lik " << logWeightDiff<< std::endl;
          
      }
    
    double diffPreviousCoalModelTime= newNodeHeightModelTime - currPop->currentModelTime*currPop->x  ;
    diffPreviousCoalModelTime= newNodeHeightModelTime - pairModelTimeEntry  ;
    
    //denominator: log proposal distribution
    if (usePriorInSMC1){
        logWeightDiff-= -Population::FmodelTstandard(newNodeHeightModelTime, currPop->timeOriginSTD, currPop->delta,K)+ Population::FmodelTstandard(pairModelTimeEntry, currPop->timeOriginSTD, currPop->delta,K);
        
    }
    else{
            
        logWeightDiff -= log(exp(pairModelTimes.logLikAtProposal-pairModelTimes.logIntegral));
        //std::cout << "denominator  " << pairModelTimes.logLikAtProposal-pairModelTimes.logIntegral<< std::endl;
        
        //std::cout << "other denominator  " << -Population::FmodelTstandard(newNodeHeightModelTime, currPop->timeOriginSTD, currPop->delta,K)+ Population::FmodelTstandard(pairModelTimeEntry, currPop->timeOriginSTD, currPop->delta,K)<< std::endl;
    }
    
    if (0){
           std::cout << "log weight diff  after  " << logWeightDiff<< std::endl;
       }
    //numerator
//    size_t numActiveGametes = currPop->numActiveGametes;
//    logWeightDiff+=  -log(currPop->theta)+log(numActiveGametes*(numActiveGametes-1)/ 2.0) +
//    Population::LogLambda(newNodeHeightModelTime, currPop->timeOriginSTD, currPop->delta,K)-
//    (numActiveGametes*(numActiveGametes-1)/ 2.0)*Population::FmodelTstandard(newNodeHeightModelTime, currPop->timeOriginSTD, currPop->delta,K)  ;
    
    logWeightDiff+=  -log(currPop->theta) +
       Population::LogLambda(newNodeHeightModelTime, currPop->timeOriginSTD, currPop->delta,K)-
      Population::FmodelTstandard(newNodeHeightModelTime, currPop->timeOriginSTD, currPop->delta,K);
//    
    std::pair<Pair,ProposalDistribInfo>  copyPairMinModelTime = std::make_pair(pairMinModelTime.first, pairMinModelTime.second);
    
    currPop->pairCurrentProposals.erase(pairMinModelTime.first);
    
    
    assert(!isnan(logWeightDiff));
    assert(!isinf(-1*logWeightDiff));
    //update pairs with their proposals
    if (0){
        
        
        std::cout << "log weight diff  before " << logWeightDiff<< std::endl;
        
    }
    pairMinModelTime = updatePairCurrentProposalsMap(iter, idxFirst, newNodeHeightModelTime, K, logWeightDiff,  copyPairMinModelTime,  random,  usePriorInSMC1, normalizedCLV, currPop);
    
    increaseNextAvailable();
    updateIndexesActiveGametes( idxFirst, idxSecond,  currPop->index, currPop->index,  newNode->index, newNode->index_population);
    
    int posFirstRoot= getNodeIdxById(idxFirstRoot);
    int posSecondRoot = getNodeIdxById(idxSecondRoot);
    
    remove_roots(posFirstRoot, posSecondRoot);
    
    assert(newNode->number_leaves_cluster>1);
    addRoot(newNode);
    
    currPop->numActiveGametes= currPop->numActiveGametes-1;
    newHeight = newNodeHeight;
    setTopHeightScaledByTheta(newHeight);
    
    if (0){
        std::cout << "pair left " << newNode->getIndexLeftChild()<< " right " << newNode->getIndexRightChild() << std::endl;
        std::cout << "new node scaled height " << newHeight<< std::endl;

        std::cout << "log weight diff after " << logWeightDiff<< std::endl;
        
    }
    if ( logWeightDiff< -1e10)
       std::cout << "log weight diff after " << logWeightDiff<< std::endl;
    
    return logWeightDiff;
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
double  State::proposalPostPost(gsl_rng * random,  double &newHeight , double& logLikNewHeight,size_t numIncrementsPOSTPOST, double K){
    
    
    int idxFirst, idxSecond, idxFirstRoot, idxSecondRoot;
    int idxReceiverPop,  idxIncomingPop;
    int idxIncomingNode;
    double logLikNewCoalTime = 0.0;
    double logWeightDiff = 0.0;
    
    Population * currPop = populationSet->getCurrentPopulation();
    std::set<std::pair<int, int>> allCoalPairs = PosetSMC::getCombinations(currPop->numActiveGametes);
    int numberProposals = allCoalPairs.size() *numIncrementsPOSTPOST;
    std::vector<double> logWeights(numberProposals, 0.0);
    std::vector<double> normWeights(numberProposals, 0.0);
    std::vector<double> logLik(numberProposals, 0.0);
    std::vector<double> normlogLik(numberProposals, 0.0);
    std::vector<double> logRatio(numberProposals, 0.0);
    std::vector<double> normlogRatio(numberProposals, 0.0);
    std::vector<double> heightProposals(numberProposals, 0.0);
    
    std::vector<std::shared_ptr<PartialTreeNode>> nodeProposals(numberProposals);
    std::vector<pairs> allPairs(numberProposals);
    
    int i=0;
    for( auto coalPair : allCoalPairs){
        
        
        idxFirst = coalPair.first;
        idxSecond = coalPair.second;
        
        assert(idxFirst!=idxSecond);
        idxFirstRoot = currPop->idsActiveGametes[idxFirst];
        idxSecondRoot = currPop->idsActiveGametes[idxSecond];
        
        assert(idxFirstRoot != idxSecondRoot);
        
        for( size_t j = 0; j!= numIncrementsPOSTPOST; ++j){
            
            double proposalCoalTimeScaledByProportion = populationSet->proposeNextCoalEventTime(  random,  idxReceiverPop,  idxIncomingPop, idxIncomingNode, logLikNewCoalTime,  K);
            assert(proposalCoalTimeScaledByProportion >0);
            if (idxReceiverPop != idxIncomingPop){
                std::cout<< "migration"<< std::endl;
                // do something with theidxIncomingNode
                //idxFirstRoot = currPop->idsActiveGametes[idxIncomingPop];
            }
            
            double  newNodeHeight = proposalCoalTimeScaledByProportion* theta;
            heightProposals[i*numIncrementsPOSTPOST+j] =proposalCoalTimeScaledByProportion;
            nodeProposals[i*numIncrementsPOSTPOST+j] = proposeNewNode( idxFirstRoot,  idxSecondRoot, currPop->index, newNodeHeight, logLikNewCoalTime, nullptr, nullptr, false );
            
            
            logLik[i*numIncrementsPOSTPOST+j] = nodeProposals[i*numIncrementsPOSTPOST+j]->ln_likelihood;
            logRatio[i*numIncrementsPOSTPOST+j] = likelihood_factor(nodeProposals[i*numIncrementsPOSTPOST+j]);
            
            allPairs[i*numIncrementsPOSTPOST+j] = std::make_pair(coalPair.first, coalPair.second);
        }
        i++;
        
    }
    unsigned int posRatios, posLik, pos;
    posLik = 0;
    //double logSumNormLik = Utils::normalize(logLik, normlogLik);
    double logSumNormRatios = Utils::normalize(logRatio, normlogRatio);
    //posLik = Random::randomDiscreteFromProbabilityVector(random, &normlogLik[0], normlogLik.size());
    posRatios = Random::randomDiscreteFromProbabilityVector(random, &normlogRatio[0], normlogRatio.size());
    pos = posRatios;
    double proportion= currPop->x;
    populationSet->acceptNextCoalEventTime(random, heightProposals[pos]/proportion,  currPop->index,  currPop->index,  K);
    
    logWeightDiff = likelihood_factor(nodeProposals[pos]);
    logWeightDiff = logSumNormRatios;
    assert(!isnan(logWeightDiff) && !isinf(logWeightDiff));
    //logWeight = logWeight+ logWeightDiff;
    
    
    
    if (true){
        std::cout<<"parent "<< nodeProposals[pos]->index<< " left " << nodeProposals[pos]->getIndexLeftChild() << " right " << nodeProposals[pos]->getIndexRightChild() << " height "<< heightProposals[pos]*theta<< std::endl;
        
    }
    
    idxFirst = allPairs[pos].first;//first idx inside the list of active gametes of the chosen pop
    idxSecond = allPairs[pos].second;//second idx inside the list of active gametes of the chosen pop
    assert(idxFirst != idxSecond);
    
    // std::cout<< "Max weight: first " << allCoalPairs[posMax]. first << "second " << allCoalPairs[posMax]. second << " weight " << weight << " max weight " << logWeights[posMax] << std::endl;
    
    idxFirstRoot = currPop->idsActiveGametes[idxFirst];
    idxSecondRoot = currPop->idsActiveGametes[idxSecond];
    
    assert(idxFirstRoot != idxSecondRoot);
    
    int posFirstRoot= getNodeIdxById(idxFirstRoot);
    int posSecondRoot = getNodeIdxById(idxSecondRoot);
    
    remove_roots(posFirstRoot, posSecondRoot);
    assert(nodeProposals[pos]->number_leaves_cluster>1);
    addRoot(nodeProposals[pos]);
    
    increaseNextAvailable();
    updateIndexesActiveGametes( idxFirst, idxSecond,  currPop->index, currPop->index,  nodeProposals[pos]->index, nodeProposals[pos]->index_population);
    
    currPop->numActiveGametes= currPop->numActiveGametes-1;
    newHeight = heightProposals[pos];
    setTopHeightScaledByTheta(heightProposals[pos]*theta);
    return logWeightDiff;
    
}
//double State::proposalCoalNodePostPost(gsl_rng * random, Population *chosenPop, double newNodeHeight, double logLikNewHeight){
//
//
//    int idxFirst, idxSecond, idxFirstRoot, idxSecondRoot;
//    double logWeightDiff = 0.0;
//    std::vector<std::pair<int, int>> allCoalPairs = Utils::allCombinations(chosenPop->numActiveGametes, 2);
//    int numberProposals = allCoalPairs.size();
//
//    std::vector<double> logWeights(numberProposals, 0.0);
//    std::vector<double> normWeights(numberProposals, 0.0);
//
//    std::vector<std::shared_ptr<PartialTreeNode>> nodeProposals(numberProposals);
//
//    double maxlogWeight = DOUBLE_NEG_INF;
//    size_t posMax;
//    for(size_t i=0; i <  allCoalPairs.size(); i++){
//        idxFirst = allCoalPairs[i]. first;
//        idxSecond = allCoalPairs[i]. second;
//
//        idxFirstRoot = chosenPop->idsActiveGametes[idxFirst];
//        idxSecondRoot = chosenPop->idsActiveGametes[idxSecond];
//
//        // std::cout << "idx first" << idxFirst << " idx second "<<idxSecond << std::endl;
//        assert(idxFirstRoot != idxSecondRoot);
//        nodeProposals[i] = proposeNewNode( idxFirstRoot,  idxSecondRoot, chosenPop->index, newNodeHeight, logLikNewHeight );
//
//        logWeights[i] = nodeProposals[i]->ln_likelihood;
//
//        if (logWeights[i]>maxlogWeight)
//        {
//            maxlogWeight = logWeights[i];
//            posMax = i;
//        }
//    }
//    unsigned int pos;
//    Utils::normalize(logWeights, normWeights);
//    pos = Random::randomDiscreteFromProbabilityVector(random, &normWeights[0], normWeights.size());
//
//    logWeightDiff = likelihood_factor(nodeProposals[pos]);
//    assert(!isnan(logWeightDiff) && !isinf(logWeightDiff));
//    //logWeight = logWeight+ logWeightDiff;
//
//
//
//    idxFirst = allCoalPairs[pos]. first;//first idx inside the list of active gametes of the chosen pop
//    idxSecond = allCoalPairs[pos]. second;//second idx inside the list of active gametes of the chosen pop
//    if (0){
//        std::cout<< "first " << idxFirst << " second " << idxSecond << " time "<< newNodeHeight<< " weight " << logWeight << " norm weight " << normWeights[pos] << std::endl;
//    }
//
//    // std::cout<< "Max weight: first " << allCoalPairs[posMax]. first << "second " << allCoalPairs[posMax]. second << " weight " << weight << " max weight " << logWeights[posMax] << std::endl;
//    idxFirstRoot = chosenPop->idsActiveGametes[idxFirst];
//    idxSecondRoot = chosenPop->idsActiveGametes[idxSecond];
//
//    assert(idxFirstRoot != idxSecondRoot);
//
//    int posFirstRoot= getNodeIdxById(idxFirstRoot);
//    int posSecondRoot = getNodeIdxById(idxSecondRoot);
//
//    remove_roots(posFirstRoot, posSecondRoot);
//    assert(nodeProposals[pos]->number_leaves_cluster>1);
//    addRoot(nodeProposals[pos]);
//
//    increaseNextAvailable();
//    updateIndexesActiveGametes( idxFirst, idxSecond,  chosenPop->index, chosenPop->index,  nodeProposals[pos]->index, nodeProposals[pos]->index_population);
//
//    chosenPop->numActiveGametes= chosenPop->numActiveGametes-1;
//    setTopHeightScaledByTheta(newNodeHeight);
//    return logWeightDiff;
//}
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
    //logWeight = logWeight+ logWeightDiff;
    
    
    std::cout<< "first " << idxFirstRoot << " second " << idxSecondRoot << " time "<< newNodeHeight<< " weight " << logWeight  << std::endl;
    return logWeightDiff;
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
    //logWeight = logWeight+ logWeightDiff;
    
    
    return logWeightDiff;
}
unsigned int State::numberNonTrivialTrees(){
    unsigned int result = 0;
    for( int i=0; i< roots.size();++i){
        if (roots[i]->number_leaves_cluster >1)
            result++;
    }
    return(result);
}
Eigen::ArrayXXd State::getLogLikelihoodGrid( double from, double to, size_t n){
    
    int idxFirst, idxSecond, idxFirstRoot, idxSecondRoot;
    
    
    Population *currPop = getCurrentPopulation();
    double newNodeHeight;
    double logLikNewHeight = 0.0;
    std::set<std::pair<int, int>> allCoalPairs = PosetSMC::getCombinations(currPop->numActiveGametes);
    int numberProposals = allCoalPairs.size()*(n+1);
    Eigen::ArrayXXd logLikWeightsPerPair(n+1,allCoalPairs.size());
    
    std::vector<double> logLikRatioGrid(numberProposals);
    std::vector<double> logLik(numberProposals, 0.0);
    std::vector<double> normlogLik(numberProposals, 0.0);
    std::vector<double> logRatio(numberProposals, 0.0);
    std::vector<double> normlogRatio(numberProposals, 0.0);
    
    std::vector<std::shared_ptr<PartialTreeNode>> nodeProposals(numberProposals);
    
    size_t i = 0;
    std::vector<pairs> allPairs(allCoalPairs.size());
    
    i=0;
    for (auto coalPair : allCoalPairs)
    {
        idxFirst = coalPair. first;
        idxSecond = coalPair. second;
        
        idxFirstRoot = currPop->idsActiveGametes[idxFirst];
        idxSecondRoot = currPop->idsActiveGametes[idxSecond];
        
        allPairs[i] = coalPair;
        // std::cout << "idx first" << idxFirst << " idx second "<<idxSecond << std::endl;
        assert(idxFirstRoot != idxSecondRoot);
        for( size_t j = 0; j<= n; ++j){
            
            newNodeHeight = from + j*(to-from)/(1.0*n);
            
            nodeProposals[i*(n+1)+j] = proposeNewNode( idxFirstRoot,  idxSecondRoot, currPop->index, newNodeHeight, logLikNewHeight, nullptr, nullptr, false );
            
            logLikWeightsPerPair(j, i) = likelihood_factor(nodeProposals[i*(n+1)+j]);
        }
        i++;
    }
    return logLikWeightsPerPair;
}
using ListDouble = std::vector<double>;
using ListListDouble = std::vector<std::vector<double>>;
using PairListListDouble = std::pair<ListListDouble,ListListDouble>;
PairListListDouble State::evalLogLikRatioGrid( double from, double to, size_t n, std::vector<std::string> &labels){
    
    int idxFirst, idxSecond, idxFirstRoot, idxSecondRoot;
    
    Population *currPop = getCurrentPopulation();
    double newNodeHeight;
    double logLikNewHeight = 0.0;
    std::set<std::pair<int, int>> allCoalPairs = PosetSMC::getCombinations(currPop->numActiveGametes);
    int numberProposals = allCoalPairs.size()*(n+1);
    std::vector<std::vector<double>> logLikWeightsPerPair(allCoalPairs.size());
    std::vector<std::vector<double>> incrementsByPair(allCoalPairs.size());
    
    std::vector<double> logLikRatioGrid(numberProposals);
    std::vector<double> logLik(numberProposals, 0.0);
    std::vector<double> normlogLik(numberProposals, 0.0);
    std::vector<double> logRatio(numberProposals, 0.0);
    std::vector<double> normlogRatio(numberProposals, 0.0);
    
    std::vector<std::shared_ptr<PartialTreeNode>> nodeProposals(numberProposals);
    
    
    
    double maxHeightPair = 0;
    size_t i = 0;
    std::vector<pairs> allPairs(allCoalPairs.size());
    std::string lab;
    i=0;
    for (auto coalPair : allCoalPairs)
    {
        idxFirst = coalPair. first;
        idxSecond = coalPair. second;
        
        idxFirstRoot = currPop->idsActiveGametes[idxFirst];
        idxSecondRoot = currPop->idsActiveGametes[idxSecond];
        
        allPairs[i] = coalPair;
        // std::cout << "idx first" << idxFirst << " idx second "<<idxSecond << std::endl;
        assert(idxFirstRoot != idxSecondRoot);
        int idxFirstID = getNodeIdxById(idxFirstRoot);
        int idxSecondId = getNodeIdxById(idxSecondRoot);
        for( size_t j = 0; j<= n; ++j){
            
            newNodeHeight = from + j*(to-from)/(1.0*n);
            
            maxHeightPair = max(roots[idxFirstID]->height,roots[idxSecondId]->height);
            if (newNodeHeight >maxHeightPair){
                
                nodeProposals[i*(n+1)+j] = proposeNewNode( idxFirstRoot,  idxSecondRoot, currPop->index, newNodeHeight, logLikNewHeight, nullptr, nullptr, false );
                
                lab = "(" + std::to_string(nodeProposals[i*(n+1)+j]->getIndexLeftChild())+", "+ std::to_string(nodeProposals[i*(n+1)+j]->getIndexRightChild())+")";
                logLikWeightsPerPair[i].push_back(likelihood_factor(nodeProposals[i*(n+1)+j]));
                incrementsByPair[i].push_back(newNodeHeight);
            }
        }
        labels.push_back(lab);
        i++;
    }
    return std::make_pair(incrementsByPair,logLikWeightsPerPair) ;
}
using PairListDouble = std::pair<ListDouble,ListDouble>;
PairListDouble State::evalLogLikRatioGridPerIncrement( double from, double to, size_t n, std::vector<std::string> &labels){
    
    int idxFirst, idxSecond, idxFirstRoot, idxSecondRoot;
    double logLikFactor;
    double normIncrement;
    Population *currPop = getCurrentPopulation();
    double newIncrement;
    double logLikNewHeight = 0.0;
    double newNodeHeight, maxHeightPair;
    std::set<std::pair<int, int>> allCoalPairs = PosetSMC::getCombinations(currPop->numActiveGametes);
    size_t numPairs = allCoalPairs.size();
    int numberProposals  = numPairs*(n+1);
    std::vector<std::vector<double>> logLikWeightsPerPair(numPairs);
    std::vector<double> logRatioForCurrentIncrement(numPairs);
    
    std::vector<double> logLikWeightIncrements;
    std::vector<double> increments;
    
    std::vector<double> logLikRatioGrid(numberProposals);
    std::vector<double> logLik(numberProposals, 0.0);
    std::vector<double> normlogLik(numberProposals, 0.0);
    std::vector<double> logRatio(numberProposals, 0.0);
    std::vector<double> normlogRatio(numberProposals, 0.0);
    
    
    size_t i = 0;
    //std::vector<pairs> allPairs(numberProposals);
    std::string lab;
    i=0;
    for( size_t j = 0; j<= n; ++j){
        
        newIncrement = from + j*(to-from)/(1.0*n);
        increments.push_back(newIncrement);
        i=0;
        for (auto coalPair : allCoalPairs)
        {
            idxFirst = coalPair. first;
            idxSecond = coalPair. second;
            
            idxFirstRoot = currPop->idsActiveGametes[idxFirst];
            idxSecondRoot = currPop->idsActiveGametes[idxSecond];
            
            //allPairs[j*(n+1)+i] = coalPair;
            int idxFirstID = getNodeIdxById(idxFirstRoot);
            int idxSecondId = getNodeIdxById(idxSecondRoot);
            // std::cout << "idx first" << idxFirst << " idx second "<<idxSecond << std::endl;
            assert(idxFirstRoot != idxSecondRoot);
            
            maxHeightPair = max(roots[idxFirstID]->height,roots[idxSecondId]->height);
            newNodeHeight = currPop->getCurrentModelTime()*currPop->x* theta +newIncrement;
            auto newNode =  uniquePtrNewNode( idxFirstRoot,  idxSecondRoot, currPop->index, newNodeHeight, logLikNewHeight, false );
            //nodeProposals[j*numPairs+i] = proposeNewNode(idxFirstRoot,  idxSecondRoot, currPop->index, newNodeHeight, logLikNewHeight, nullptr, nullptr );
            logLikFactor = newNode->likelihood_factor();
            logRatioForCurrentIncrement[i]= logLikFactor;
            
            i++;
        }
        logLikWeightIncrements.push_back(Utils::norm(logRatioForCurrentIncrement));
        
        normIncrement = Utils::logSumExp(logRatioForCurrentIncrement);
        lab =  std::to_string(newNodeHeight);
        labels.push_back(lab);
    }
    
    return std::make_pair(increments,logLikWeightIncrements);
}
std::pair<Pair, ProposalDistribInfo> State::initPairProposalsNextEventTime(gsl_rng *rngGsl, bool normalizedCLV
                                                                           , double K){
    
    int idxFirst, idxSecond, idxFirstRoot, idxSecondRoot;
    Population *popI = populationSet->getCurrentPopulation();
   
    assert(popI->numActiveGametes > 1);
    
    std::vector<Pair > pairs= PosetSMC::getCombinationsRandomOrder(popI->numActiveGametes);
    
    auto rd = std::random_device {};
    auto rng = std::default_random_engine { rd() };
    
    std::shuffle(std::begin(pairs), std::end(pairs), rng);
    
    double minModeltime = DOUBLE_INF;
    double proposalModelTime;
    Pair currPair;
    std::pair<Pair, ProposalDistribInfo> currEntry;
    std::pair<Pair, ProposalDistribInfo> winnerModelTime;
    
    int max_internal_points = 50;
    
    double currModeltime = popI->currentModelTime;
    if (currModeltime== 0.0)
        currModeltime = 1e-10;
    double torigin = 0.999*popI->timeOriginSTD;
    double interPoint1;
    double interPoint2 ;

    
    GNUPlotter plotter;
    size_t pair_idx=0;

    std::vector<double> sampledIncrementsPerPair(pairs.size());
    double lastSampledX = 0.0;
    double lastSampledY = 0.0;
    double logComplArea = 0.0;
    bool doPlots = true;
    if (doPlots){
        
         GNUPlotter plotter;
        
    }
        
    for(std::vector<Pair>::iterator it = pairs.begin(), end = pairs.end(); it != end; ++it)
    {
        currPair = *it;
        idxFirst = currPair. first;
        idxSecond = currPair. second;
        
        idxFirstRoot = popI->idsActiveGametes[idxFirst];
        idxSecondRoot = popI->idsActiveGametes[idxSecond];
        
        int idxFirstID = getNodeIdxById(idxFirstRoot);
        int idxSecondId = getNodeIdxById(idxSecondRoot);
        
        //interPoint1 = currModeltime+ Random::randomUniformFromGsl2(rngGsl)*(torigin-currModeltime);
        
        interPoint1 = currModeltime + 0.05*(torigin-currModeltime);
        //interPoint2 =currModeltime+ Random::randomUniformFromGsl2(rngGsl)*(torigin-currModeltime);
        interPoint2 = currModeltime + 0.75*(torigin-currModeltime);
    
        interPoint1 = currModeltime + 0.01*(torigin-currModeltime);
        std::vector<double> x = {currModeltime, interPoint1, interPoint2, torigin};
        
        Eigen::VectorXd xEigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(x.data(), x.size());
        CCLogDensity *log_density =
        new GenotypeJCPairLogProposal(partition->numberSites, partition->numberStates, popI->timeOriginSTD, popI->delta, popI->theta, currModeltime,
                                      roots[idxFirstID]->height,
                                      roots[idxSecondId]->height, roots[idxFirstID]->pclv, roots[idxSecondId]->pclv, normalizedCLV);
        
        if (doPlots){
                         //GNUPlotter plotter;
                         size_t numPoints= 50;
                         log_density->plot(plotter, numPoints, 0, idxFirstID,  idxSecondId );
        }
       
        CCARS *ccars = new CCARS(log_density, xEigen, currModeltime, torigin,  max_internal_points);
        lastSampledX = 0.0;
        lastSampledY = 0.0;
        logComplArea = 0.0;
        
       // std::cout << "approximating integral for pair " << std::to_string(idxFirstID)+ ", "+ std::to_string(idxSecondId)<< std::endl;
        double logAreaIntegral = ccars->approximateLogIntegral(rngGsl, max_internal_points);
        
        if (doPlots){
                   //GNUPlotter plotter;
                   size_t numPoints= 50;
                  ccars->plot( plotter, numPoints, 0, idxFirstID,  idxSecondId);
        
               }
        lastSampledX = ccars->sample(rngGsl, lastSampledY, logComplArea);
        
        if(logComplArea >logAreaIntegral){
            std::cout << " Pair: ( " <<idxFirst << " ,"<< idxSecond << ")" << std::endl;
            std::cout << "log comp Area= " <<logComplArea << std::endl;
            std::cout << " logAreaIntegral= " << logAreaIntegral << std::endl;
            logAreaIntegral  = logComplArea;
        }
        sampledIncrementsPerPair[pair_idx] =lastSampledX;
        pair_idx++;
        
        proposalModelTime = lastSampledX;
        ProposalDistribInfo proposalInfo;
        proposalInfo.creationTime  = currModeltime;
        proposalInfo.logIntegral = logAreaIntegral;
        proposalInfo.timeProposal = lastSampledX;
        proposalInfo.logComplementArea = logComplArea;
        proposalInfo.logLikAtProposal = lastSampledY;
        
        currEntry = std::make_pair(currPair,proposalInfo);
        if (proposalModelTime<minModeltime){
            
            winnerModelTime = currEntry;
            minModeltime =proposalModelTime;
        }
        
        popI->pairCurrentProposals.insert(currEntry);
    }
    
    return winnerModelTime;
}
std::pair<Pair, ProposalDistribInfo> State::updatePairCurrentProposalsMap(int iter, int posNewNodeIdsGametes, double modelTimeNewNode,double K, double &logWeightDiff, std::pair<Pair, ProposalDistribInfo> copyPairMinModelTime,   gsl_rng *rngGsl, bool usePriorInSMC1, bool normalizedCLV, Population* currPop){
    
    double waitingTimeKingman, proposalKingmanTime, proposalModelTime;
    double minModeltime = DOUBLE_INF;
    double currKingmanCoalTime = Population::FmodelTstandard(modelTimeNewNode, currPop->timeOriginSTD, currPop->delta, K);
    std::pair<Pair, ProposalDistribInfo>  pairMinModelTime;
    int idxFirst, idxSecond, idxFirstRoot, idxSecondRoot;
    double logLikNewHeight = 0.0;
    int max_internal_points = 50;
    bool doPlots = false;
//    if (doPlots)
//        GNUPlotter plotter;
    
    double torigin  = 0.999*currPop->timeOriginSTD;
    assert(currPop->pairCurrentProposals.size()== (currPop->numActiveGametes)*(currPop->numActiveGametes-1)/ 2.0 -1);
    PairInfoMap::iterator it = currPop->pairCurrentProposals.begin();
    double lik_factor_newNode;
    
    
    while ( it != currPop->pairCurrentProposals.end()) {
        
        auto currPair = it->first;

        
        if (Utils::pairsIntersected(currPair, copyPairMinModelTime.first)){
            
            double dropoutNodeProposalModelTime =  it->second.timeProposal  ;
            double dropoutNodeModelTimeEntry =   it->second.creationTime ;
            //numerator
//            logWeightDiff +=   -log(theta)+log(currPop->numActiveGametes*(currPop->numActiveGametes-1)/ 2.0) +
//            Population::LogLambda(dropoutNodeProposalModelTime, currPop->timeOriginSTD, currPop->delta,K)-
//            (currPop->numActiveGametes*(currPop->numActiveGametes-1)/ 2.0)*Population::FmodelTstandard(dropoutNodeProposalModelTime, currPop->timeOriginSTD, currPop->delta,K) ;
//
            logWeightDiff +=   -log(theta) +
            Population::LogLambda(dropoutNodeProposalModelTime, currPop->timeOriginSTD, currPop->delta,K)-
            Population::FmodelTstandard(dropoutNodeProposalModelTime, currPop->timeOriginSTD, currPop->delta,K) ;
            
            //denominator
            if (usePriorInSMC1){
                logWeightDiff -= -Population::FmodelTstandard(dropoutNodeProposalModelTime, currPop->timeOriginSTD, currPop->delta,K)+ Population::FmodelTstandard(dropoutNodeModelTimeEntry, currPop->timeOriginSTD, currPop->delta,K);
            }
            else{
                logWeightDiff -=  log(exp(it->second.logComplementArea-it->second.logIntegral));
            }
            idxFirst = currPair.first;
            idxSecond = currPair.second;
            
            idxFirstRoot = currPop->idsActiveGametes[idxFirst];
            idxSecondRoot = currPop->idsActiveGametes[idxSecond];
            
            assert(idxFirstRoot != idxSecondRoot);
            
            int idxFirstID = getNodeIdxById(idxFirstRoot);
            int idxSecondId = getNodeIdxById(idxSecondRoot);
            
//            logLikNewHeight = 0.0;
//
//            auto newNode =  uniquePtrNewNode( idxFirstRoot,  idxSecondRoot, currPop->index, dropoutNodeProposalModelTime, logLikNewHeight, !usePriorInSMC1 );
//
//            lik_factor_newNode = newNode->likelihood_factor();
//
//            if (lik_factor_newNode > -1e10){
//                //logWeightDiff -= lik_factor_newNode;
//            }
//            else{
//                std::cout << "log lik factor -inf for pair "<< currPair.first << " , " << currPair.second << std::endl;
//                std::cout <<lik_factor_newNode << std::endl;
//            }
            //only erase the pairs related with the last position of idsActiveGametes(the other ones will be updated)
            if (currPair.first == copyPairMinModelTime.first.second  || currPair.second ==copyPairMinModelTime.first.second)
                it = currPop->pairCurrentProposals.erase(it);
            else{
                if (usePriorInSMC1){
                    waitingTimeKingman =  Random::RandomExponential(1.0,0,  true, rngGsl,NULL );
                    proposalKingmanTime = waitingTimeKingman + currKingmanCoalTime;
                    
                    proposalModelTime = Population::GstandardTmodel(proposalKingmanTime, currPop->timeOriginSTD, currPop->delta, K);
                    
                }
                else{
                    
                    double interPoint1 =modelTimeNewNode + 0.05*(torigin-modelTimeNewNode);
                    
                    double interPoint2 =modelTimeNewNode + 0.75*(torigin-modelTimeNewNode);
                    
                    interPoint2 =modelTimeNewNode + 0.01*(torigin-modelTimeNewNode);
                    
                    std::vector<double>  x = {modelTimeNewNode, interPoint2, torigin};
                    
                    Eigen::VectorXd xEigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(x.data(), x.size());
                    
                    CCLogDensity *log_density =
                    new GenotypeJCPairLogProposal(partition->numberSites, partition->numberStates, currPop->timeOriginSTD, currPop->delta,theta, modelTimeNewNode,
                                                  roots[idxFirstID]->height,
                                                  roots[idxSecondId]->height, roots[idxFirstID]->pclv, roots[idxSecondId]->pclv, normalizedCLV);
                    
                    if (doPlots){
                                      GNUPlotter plotter;
                                      size_t numPoints= 50;
                                      log_density->plot(plotter, numPoints, iter,  idxFirstID,  idxSecondId );
                     }
                    
    
                    CCARS *ccars = new CCARS(log_density, xEigen, modelTimeNewNode, torigin,  max_internal_points);
                    double lastSampledX = 0.0;
                    double lastSampledY = 0.0;
                    double complementArea = 0.0;
                    double logAreaIntegral = ccars->approximateLogIntegral(rngGsl, max_internal_points);
                    
                    if (doPlots){
                               GNUPlotter plotter;
                               size_t numPoints= 50;
                              ccars->plot( plotter, numPoints, iter,  idxFirstID,  idxSecondId);
                    
                           }
                    lastSampledX = ccars->sample(rngGsl, lastSampledY, complementArea);
                    it->second.creationTime    = modelTimeNewNode;//creation time
                    it->second.timeProposal  = lastSampledX;//new proposal time
            
                    it->second.logIntegral = logAreaIntegral;
                    it->second.timeProposal = lastSampledX;
                    it->second.logComplementArea = complementArea;
                    it->second.logLikAtProposal = lastSampledY;
                }
               
                if (it->second.timeProposal <minModeltime){
                    
                    pairMinModelTime = *it;
                    minModeltime = it->second.timeProposal ;
                }
                it++;
            }
        }
        else{
            //find the minimum proposal model time in the pairs that didnt dropout
            if (it->second.timeProposal < minModeltime){
                
                minModeltime = it->second.timeProposal;
                pairMinModelTime = *it;
            }
            it++;
        }
    }
    assert(currPop->pairCurrentProposals.size()== (currPop->numActiveGametes-1)*(currPop->numActiveGametes-2)/ 2.0);
    return pairMinModelTime;
}
State::~State(){
}
