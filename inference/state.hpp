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
   
    long double heightScaledByTheta;
    long double heightModelTime;
    //std::shared_ptr<PopulationSet> populationSet;
    PopulationSet * populationSet;
    unsigned int nextAvailable;
    PLLBufferManager *const pll_buffer_manager;

    const Partition *partition;
    unsigned int num_sites;
    std::vector<std::shared_ptr<PartialTreeNode>> roots;
   
    std::vector<int> idsNextCoalEvents;
    GenotypeErrorModel *gtError;
    long double theta;
    double initialLogWeight;
public:

    State( PosetSMCParams &params, gsl_rng *random);
    
    State(const State &original);
    std::vector<pll_rnode_t *> getForest() const;

    void initForest( int sampleSize, pll_msa_t *msa, std::vector<int> &positions, ProgramOptions &programOptions);
    std::shared_ptr<PartialTreeNode> getRootAt(int i) const {return roots[i];};
 
    
    State &operator=(const State &original);
    
    std::shared_ptr<PartialTreeNode> connect(int i, int j,  size_t index_pop_new_node);
    
    std::vector<std::shared_ptr<PartialTreeNode>> getRoots() const{ return roots;};
    
    double likelihood_factor(std::shared_ptr<PartialTreeNode> root)const;
    
    void remove_roots(int i, int j);
    
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
    
    PopulationSet& getPopulationSet()const  {return *populationSet;};
    
    Population* getPopulationByIndex(size_t i) const {return populationSet->getPopulationbyIndex(i);};
    
    int getNumberPopulations() const {return populationSet->numClones;};
    
    long double getHeightScaledByTheta()const{return heightScaledByTheta;}
    long double getHeightModelTime()const{return heightModelTime;}
    long double getTheta()const{return theta;}
    
    void printTree(std::shared_ptr<PartialTreeNode> root, std::ostream &stream);
    void initIdsNextCoalEvents(int numClones);
    
    void setHeightScaledByTheta(long double  newHeight)
    {
        assert(newHeight>=heightScaledByTheta);
        heightScaledByTheta= newHeight;
        heightModelTime = heightScaledByTheta/theta;
        
    }
    unsigned int getNextAvailable()const{return nextAvailable;}
    void increaseNextAvailable(){nextAvailable++;}
    
    GenotypeErrorModel &getErrorModel() const{return *gtError;}
    double getInitialLogWeight() const {return initialLogWeight;}
    
    int getNodeIdxById(size_t id);
    int getIdNextCoalEventForPopulation(int i);
    void moveNextIdEventForPopulation(int i);
    
    ~State();
};
#endif /* state_hpp */
