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
#include "pll_buffer_manager.hpp"
#include "population.hpp"
#include "genotype_error_model.hpp"
#include "partition.hpp"

#include  <vector>
extern "C"
    {
#include "libpll/pll.h"
    }
class PosetSMCParams;
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
public:
    //constructors
    State( PosetSMCParams &params, gsl_rng *random);
    
    State(const State &original);
    
    State &operator=(const State &original);
    
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
    
    double getInitialLogWeight() const {return logWeight;}
    
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
    
    std::shared_ptr<PartialTreeNode> proposeNewNode(int firstId, int secondId, size_t index_pop_new_node, double newNodeHeight , double logLikNewNode);
    
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
    double proposalCoalNodePriorPost(gsl_rng * random, Population *chosenPop, double newNodeTime, double logLikNewHeight);
    double   proposalPriorPost(gsl_rng * random, Population *leftNodePop,Population *rightNodePop,  double newNodeHeight, double logLikNewHeight);
    double   proposalCoalMRCANodePriorPost(gsl_rng * random, Population *inmigrantPop,Population *receiverPop, double newNodeHeight, double logLikNewHeight);
    
    //PostPost
    double  proposalPostPost(gsl_rng * random,  double &newHeight , double& logLikNewHeight, int K);
     double proposalCoalNodePostPost(gsl_rng * random, Population *chosenPop, double newNodeTime, double logLikNewHeight);
    double   proposalCoalMRCANodePostPost(gsl_rng * random, Population *inmigrantPop,Population *receiverPop, double newNodeHeight, double logLikNewHeight);
    unsigned int numberNonTrivialTrees();
    ~State();
};
#endif /* state_hpp */
