//
//  state.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 16/03/2021.
//

#ifndef state_hpp
#define state_hpp

//#include "poset_smc_params.hpp"
#include "partial_tree_node.hpp"
#include "core_likelihood.hpp"
#include "pll_buffer_manager.hpp"
#include "population.hpp"
#include "genotype_error_model.hpp"
#include "partition.hpp"



#include <Eigen/Core>
#include  <vector>
extern "C"
    {
#include "libpll/pll.h"
    }
class PosetSMCParams;

using Pair = std::pair<int, int>;

using ListDouble = std::vector<double>;
using PairListDouble = std::pair<ListDouble,ListDouble>;
using ListListDouble = std::vector<std::vector<double>>;
using PairListListDouble = std::pair<ListListDouble,ListListDouble>;


class State{
    
    long double topHeightScaledByTheta;
    long double topHeightModelTime;
    //std::shared_ptr<PopulationSet> populationSet;
    PopulationSet * populationSet;
    unsigned int nextAvailable;
    PLLBufferManager *const pll_buffer_manager;
    const Partition *partition;
    unsigned int num_sites;
    std::vector<std::shared_ptr<PartialTreeNode>> roots;
    std::vector<std::shared_ptr<PartialTreeNode>> postorder;
    std::vector<int> idsNextCoalEvents;
    std::vector<int> idsNextInmigrationEvents;
    std::vector<double> coalEventTimesScaledByTheta;
    GenotypeErrorModel *gtError;
    long double theta;
    double logWeight;
    std::pair<Pair, ProposalDistribInfo> pairMinModelTime;
    
public:
    //constructors
    State( PosetSMCParams &params, gsl_rng *random);
    
    State(const State &original);
    
    State &operator=(const State &original);
    
    State(State&& rhs);
    
    //getters
    std::vector<pll_rnode_t *> getForest() const;
    
    PopulationSet& getPopulationSet()const  {return *populationSet;};
    
    Population* getPopulationByIndex(size_t i) const {return populationSet->getPopulationbyIndex(i);};
    
    int getNumberPopulations() const {return populationSet->numClones;};
    
    long double getTopHeightScaledByTheta()const{return topHeightScaledByTheta;}
    long double getTopHeightModelTime()const{return topHeightModelTime;}
    long double getTheta()const{return theta;}
    
    std::shared_ptr<PartialTreeNode> getRootAt(int i) const {return roots[i];};
    
    std::vector<std::shared_ptr<PartialTreeNode>> getRoots() const{ return roots;};
    
    double getLogWeight() const {return logWeight;}
    void  setLogWeight(double newLogWeight){ logWeight = newLogWeight;}
    std::vector<double> getCoalEventTimesScaledBytheta() const{return coalEventTimesScaledByTheta;}
    
    int getNodeIdxById(size_t id);
    int getIdNextCoalEventForPopulation(int i);
    int getIdNextInmigrationEventForPopulation(int i);
    GenotypeErrorModel &getErrorModel() const{return *gtError;}
    unsigned int getNextAvailable()const{return nextAvailable;}
    
    //methods
    void initForest( int sampleSize, pll_msa_t *msa, std::vector<int> &positions, ProgramOptions &programOptions);
    
    std::shared_ptr<PartialTreeNode> connect(int i, int j,  size_t index_pop_new_node, double newNodeHeight, double  logLikNewNode);
    
    double likelihood_factor(std::shared_ptr<PartialTreeNode> root)const;
    
    void remove_roots(int i, int j);
    
    std::shared_ptr<PartialTreeNode> proposeNewNode(int firstId, int secondId, size_t index_pop_new_node, double newNodeHeight , double logLikNewNode, double* left_pmatrix, double* right_pmatrix, bool normalize);
    

    template<typename... Ts>
    auto uniquePtrNewNode(int firstId, int secondId, size_t index_pop_new_node, double newNodeHeight, double logLikNewNode, bool normalizeClv){
        
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
        
        std::unique_ptr<PartialTreeNode, decltype(delPartialTreeNode)> parent(nullptr, delPartialTreeNode);
    //    parent.reset(new PartialTreeNode(std::forward<Ts>(
    //                                                      pll_buffer_manager, child_left, child_right,
    //                                                      left_length, right_length,
    //                                                      pmatrix_elements,
    //                                                      "", newNodeHeight, clv_elements,
    //                                                      scaler_size, partition->alignment(), nextAvailable,
    //                                                      index_pop_new_node, nullptr, nullptr)...)) ;
        
        
        parent.reset(new PartialTreeNode(
                                                          pll_buffer_manager, child_left, child_right,
                                                          left_length, right_length,
                                                          pmatrix_elements,
                                                          "", newNodeHeight, clv_elements,
                                                          scaler_size, partition->alignment(), nextAvailable,
                                                          index_pop_new_node, nullptr, nullptr)) ;
        const unsigned int matrix_indices[1] = {0};
        //Array params_indices holds the indices of rate matrices that will be used for each rate category. This array must be of size rate_cats
        
        const unsigned int param_indices[1] = {0};
        
        int left_edge_pmatrix_result;
        if (normalizeClv){
            
            left_edge_pmatrix_result= pll_core_update_pmatrix_16x16_jc69_matrix_second_form( &parent->edge_l->pmatrix, p->numberStates,p->numberRateCats, p->rates(),  &parent->edge_l->length,matrix_indices, param_indices,
                                                                                       1,
                                                                                       p->attributes);
            
        }
        else{
//            left_edge_pmatrix_result= pll_core_update_pmatrix_16x16_jc69( &parent->edge_l->pmatrix, p->numberStates,p->numberRateCats, p->rates(),  &parent->edge_l->length,matrix_indices, param_indices,
//                                                                            1,
//                                                                            p->attributes);

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
        }
        
        assert(left_edge_pmatrix_result == PLL_SUCCESS);
        int right_edge_pmatrix_result;
        if (normalizeClv){
                   
                   right_edge_pmatrix_result= pll_core_update_pmatrix_16x16_jc69_matrix_second_form( &parent->edge_r->pmatrix, p->numberStates,p->numberRateCats, p->rates(),  &parent->edge_r->length,matrix_indices, param_indices,
                   1,
                   p->attributes);
                   
               }
        else{
//                   right_edge_pmatrix_result= pll_core_update_pmatrix_16x16_jc69( &parent->edge_r->pmatrix, p->numberStates,p->numberRateCats, p->rates(),  &parent->edge_r->length,matrix_indices, param_indices,
//                                                                              1,
//                                                                              p->attributes);
        
                   right_edge_pmatrix_result = pll_core_update_pmatrix(&parent->edge_r->pmatrix, p->numberStates,     p->numberRateCats, p->rates(),
                                                                &parent->edge_r->length, matrix_indices, param_indices, p->propInvar(),
                                                                p->eigenVals(), p->eigenVecs(), p->invEigenVecs(), 1, p->attributes);
        }
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
                                   p->attributes//,normalizeClv
                                   );
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
    
    int root_count() const{return roots.size();};
    
    double compute_ln_likelihood(std::vector<double> &clv, std::vector<unsigned int>  &scale_buffer,
                                 const Partition *p);
    
    double compute_ln_likelihood(double *pclv, unsigned int* pscale_buffer,
                                 const Partition *p);
    double compute_ln_likelihood(std::vector<double> &clv, std::vector<unsigned int>  &scale_buffer,
                                 const pll_partition_t *p);
    
    double compute_ln_likelihood(double *pclv, unsigned int* pscale_buffer,
                                 const pll_partition_t *p);
    
    void updateIndexesActiveGametes(int idxFirstId, int idxSecondId, size_t idxPopFirst,size_t idxPopSecond, size_t newNodeId, size_t idxPopNewNode);
    
    
    void printTree(std::shared_ptr<PartialTreeNode> root, std::ostream &stream);
    
    void printTreeChronologicalOrder(std::shared_ptr<PartialTreeNode> root, std::ostream &stream);
    
    void initIdsNextCoalEvents(int numClones);
    void initIdsInmigrationEvents(int numClones);
    
    void addRoot(std::shared_ptr<PartialTreeNode> node);
    
    void setTopHeightScaledByTheta(long double  newHeight)
    {
        assert(newHeight>=topHeightScaledByTheta);
        topHeightScaledByTheta= newHeight;
        topHeightModelTime = topHeightScaledByTheta/theta;
        coalEventTimesScaledByTheta.push_back(newHeight);
    }
    void insertNewCoalTime(double newHeight){
        assert(newHeight>0);
        coalEventTimesScaledByTheta.push_back(newHeight);
    }
    void increaseNextAvailable(){nextAvailable++;}
    

    void moveNextIdEventForPopulation(int i, bool isCoal);
    
    void addNodeToPostorderByIndex(int idx);
    
    void resetNumActiveGametesCounter();
    
    std::string getNewickRecursive(std::shared_ptr<PartialTreeNode> root, std::string &newick);
    std::string getNewick(std::shared_ptr<PartialTreeNode> root);
    
    double getNextCoalTime(gsl_rng *random, int& idxLeftNodePop, int& idxRightNodePop,double &logLik, double K);
    
    //PriorPrior
    double  proposalCoalNodePriorPrior(gsl_rng * random, Population *chosenPop, double newNodeTime, double logLikNewHeight);
    double  proposalPriorPrior(gsl_rng * random, Population *leftNodePop,Population *rightNodePop, double newNodeHeight, double logLikNewHeight);
    double  proposalCoalMRCANodePriorPrior(gsl_rng * random, Population *inmigrantPop,Population *receiverPop, double newNodeHeight, double logLikNewHeight);
    //PriorPost
    double proposalCoalNodePriorPost(gsl_rng * random, Population *chosenPop, double newNodeTime, double logLikNewHeight, size_t iteration, double** pmatrices);
    double   proposalPriorPost(gsl_rng * random, Population *leftNodePop,Population *rightNodePop,  double newNodeHeight, double logLikNewHeight, size_t iteration);
    double   proposalCoalMRCANodePriorPost(gsl_rng * random, Population *inmigrantPop,Population *receiverPop, double newNodeHeight, double logLikNewHeight);
    
    //PostPost
    double  proposalPostPost(gsl_rng * random,  double &newHeight , double& logLikNewHeight, size_t numIncrementsPOSTPOST, double K);
     double proposalCoalNodePostPost(gsl_rng * random, Population *chosenPop, double newNodeTime, double logLikNewHeight);
    double   proposalCoalMRCANodePostPost(gsl_rng * random, Population *inmigrantPop,Population *receiverPop, double newNodeHeight, double logLikNewHeight);
    unsigned int numberNonTrivialTrees();
    Population* getCurrentPopulation();
    
    void acceptNextCoalTimeProposal(double nextCoaTimeScaledByProportion, gsl_rng *random, int idxReceiverPop, int idxIncomingPop,  double K);
    
    double  proposalPostPost1(gsl_rng * random,  double &newHeight , double& logLikNewHeight,size_t numIncrementsPOSTPOST, double K);
    
    double   proposalTSMC1(gsl_rng * random,  double &newHeight , double& logLikNewHeight, double K,  bool normalize);
    
    Eigen::ArrayXXd getLogLikelihoodGrid( double from, double to, size_t n);
    
    PairListListDouble evalLogLikRatioGrid( double from, double to, size_t n, std::vector<std::string> &labels);
    
    PairListDouble evalLogLikRatioGridPerIncrement( double from, double to, size_t n, std::vector<std::string> &labels);
    
    double likelihood_factor(std::unique_ptr<PartialTreeNode> root)const;
    
    std::pair<Pair, ProposalDistribInfo> initPairProposalsNextEventTime(gsl_rng *rngGsl, double K);
    
    ~State();
};
#endif /* state_hpp */
