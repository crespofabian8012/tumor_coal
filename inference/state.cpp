//
//  state.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 16/03/2021.
//


#include "state.hpp"
#include "pll_utils.hpp"
#include "core_likelihood.hpp"
#include "treeLikelihood.hpp"
#include <memory>
#include "poset_smc_params.hpp"


State::State(PosetSMCParams &smcParams, gsl_rng *random):
pll_buffer_manager(smcParams.pll_buffer_manager),partition(smcParams.partition), num_sites(smcParams.num_sites)
{
    heightModelTime = 0.0;
    heightScaledByTheta =0.0;
    initialLogWeight = 0.0;
    
    initIdsNextCoalEvents(smcParams.getProgramOptions().numClones);
    gtError =  smcParams.gtErrorModel;
    
    populationSet = new PopulationSet(smcParams.getProgramOptions().numClones);
    populationSet->initPopulationSampleSizes(smcParams.sampleSizes);
    populationSet->initPopulationGametes();
    populationSet->initPopulationRTips();
    populationSet->initProportionsVectorFromSampleSizes(smcParams.sampleSizes);
    
    //populationSet->initDeltaThetaFromPriors(random, theta);
    theta = *smcParams.theta;
    populationSet->initPopulationDeltas( smcParams.populationDeltaTs);
    populationSet->initPopulationTOriginSTD( smcParams.populationToriginSTDs);
    populationSet->initTheta(theta);
    
   // populationSet->setPopulationsToriginConditionalDelta(random);
    populationSet->sortPopulationsByTorigin();
    
    populationSet->initListPossibleMigrations();
    populationSet->initPopulationsCoalescentEvents();
    
    initForest(smcParams.sampleSize, smcParams.msa, smcParams.positions,  smcParams.getProgramOptions());
    nextAvailable = smcParams.sampleSize-1;
    
    assert(smcParams.positions.size()>0);
    postorder.reserve(2*smcParams.positions.size()-1);
}
State::State(const State &original):pll_buffer_manager(original.pll_buffer_manager){
    
    heightModelTime = original.heightModelTime;
    heightScaledByTheta = original.heightScaledByTheta;
    num_sites=original.num_sites;
    
    populationSet=new PopulationSet(*(original.populationSet));
    gtError = new GenotypeErrorModel(*(original.gtError));
    
    idsNextCoalEvents = original.idsNextCoalEvents;
     roots = original.roots;
     postorder  = original.postorder;
//    roots.reserve(original.root_count());
//    for (auto const& fptr : original.getRoots())
//         roots.emplace_back(fptr->Clone());
   
    
    partition = original.partition;
    nextAvailable = original.nextAvailable;
    theta = original.theta;
    initialLogWeight = original.initialLogWeight;

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

        p=std::make_shared<PartialTreeNode>(pll_buffer_manager, nullptr, nullptr, msa->label[i], heightScaledByTheta, clv_num_elements,
                                            scaler_size, reference_partition->alignment(), i);
        
        //delete (p->pclv);
       //delete (p->pscale_buffer);
        // auto clv_elem = reference_partition->numberSites * reference_partition->numberStates;
 
        p->buildCLV(i,reference_partition->numberStates,  msa, gtError,  false);
        
        std::vector<unsigned  int> invSites(reference_partition->numberStates, 0);
        
         pll_count_invariant_sites(reference_partition->getPartition(),  invSites.data());
        
        p->ln_likelihood = compute_ln_likelihood(p->pclv, nullptr, reference_partition->getPartition());
        
        initialLogWeight +=p->ln_likelihood;
       
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



State &State::operator=(const State &original){
    
    if (this == &original)
        return *this;
    
    partition = original.partition;
    heightScaledByTheta = original.heightScaledByTheta;
    heightModelTime = original.heightModelTime;
    theta = original.theta;
    num_sites=original.num_sites;
    roots = original.roots;
    postorder  = original.postorder;
    gtError=original.gtError;
    populationSet = original.populationSet;
    nextAvailable = original.nextAvailable;
    idsNextCoalEvents = original.idsNextCoalEvents;
    return *this;
    
}


std::shared_ptr<PartialTreeNode> State::connect(int firstId, int secondId, size_t index_pop_new_node ){
    
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
    
    double left_length =heightScaledByTheta - child_left->height;
    double right_length = heightScaledByTheta - child_right->height;
    
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
                                                                                pll_buffer_manager, edge_left, edge_right, "", heightScaledByTheta, clv_elements,
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
    
    // assert(parent->ln_likelihood <= 0 && "Likelihood can't be more than 100%");

    // Remove children
    remove_roots(idxFirstID, idxSecondId);
    // Add new internal node
    addRoot(parent);
    //std::cout<< "root count "<< roots.size()<< std::endl;

    return parent;
}
std::shared_ptr<PartialTreeNode> State::proposeNewNode(int firstId, int secondId, size_t index_pop_new_node ){
        
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
    
    double left_length =heightScaledByTheta - child_left->height;
    double right_length = heightScaledByTheta - child_right->height;
    
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
                                                                                pll_buffer_manager, edge_left, edge_right, "", heightScaledByTheta, clv_elements,
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
    

    parent->ln_likelihood =
    compute_ln_likelihood(parent->pclv,
                          parent->pscale_buffer,
                          p->getPartition());
   

    if(isnan(parent->ln_likelihood) || isinf(parent->ln_likelihood)){
        
        std::cout << "Error: log lik  "<<  " is "<< parent->ln_likelihood<< std::endl;
        
    }
    
    assert(!isnan(parent->ln_likelihood) && !isinf(parent->ln_likelihood));
    
    // assert(parent->ln_likelihood <= 0 && "Likelihood can't be more than 100%");
  
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
    
//    if (ln_l <= -10000000 || ln_r <= -10000000 || ln_m <= -10000000)
//       std::cout << " left loglik " <<ln_l << " right loglik " <<ln_r << " parent loglik " <<ln_m << std::endl;
//    
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
    // const unsigned int parameter_indices[4] = {0, 0, 0, 0};
    const unsigned int parameter_indices[1] = {0};
    // const unsigned int parameter_indices[16] = {0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    
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
    // const unsigned int parameter_indices[4] = {0, 0, 0, 0};
    const unsigned int parameter_indices[1] = {0};
    // const unsigned int parameter_indices[16] = {0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    
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
    // const unsigned int parameter_indices[4] = {0, 0, 0, 0};
    const unsigned int parameter_indices[1] = {0};
    
    
    double result = pll_core_root_loglikelihood2(  p->states, p->sites, p->rate_cats,
                                                 clv.data(), scale_buffer.data(),
                                                 p->frequencies, p->rate_weights, p->pattern_weights, parameter_indices,  nullptr, p->attributes);
    
    //    double result= pll_core_root_loglikelihood(
    //                                               p->states, p->sites, p->rate_cats,
    //                                               clv, scale_buffer,
    //                                               p->frequencies, p->rate_weights, p->pattern_weights, p->prop_invar,
    //                                               p->invariant, parameter_indices, nullptr, p->attributes);
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
int State::getIdNextCoalEventForPopulation(int i){
    return idsNextCoalEvents[i];
}
void State::moveNextIdEventForPopulation(int i){
    idsNextCoalEvents[i]= idsNextCoalEvents[i]+1;
}
void State::addRoot(std::shared_ptr<PartialTreeNode> node){
    
   roots.push_back(node);
   addNodeToPostorderByIndex(roots.size()-1);
   
}
void State::addNodeToPostorderByIndex(int idx){

    postorder.push_back(roots[idx]);
    
}

State::~State(){
    
    
}
