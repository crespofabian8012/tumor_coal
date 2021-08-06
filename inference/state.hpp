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
   
    long double height;
    //std::shared_ptr<PopulationSet> populationSet;
    PopulationSet * populationSet;
    unsigned int nextAvailable;
    PLLBufferManager *const pll_buffer_manager;

    const Partition *partition;
    unsigned int num_sites;
    std::vector<std::shared_ptr<PartialTreeNode>> roots;
    //std::shared_ptr<GenotypeErrorModel> gtError;
    GenotypeErrorModel *gtError;
public:

    State( PosetSMCParams &params, gsl_rng *random);
    
    State(const State &original);
    std::vector<pll_rnode_t *> getForest() const;

    void initForest( int sampleSize, pll_msa_t *msa, std::vector<int> &positions, ProgramOptions &programOptions);
    std::shared_ptr<PartialTreeNode> getRootAt(int i) const {return roots[i];};
    //std::shared_ptr<State> makeDeepCopy();
    
    State &operator=(const State &original);
    
    std::shared_ptr<PartialTreeNode> connect(int i, int j, double height_delta, unsigned int index_population);
    
    std::vector<std::shared_ptr<PartialTreeNode>> getRoots() const{ return roots;};
    
    double likelihood_factor(std::shared_ptr<PartialTreeNode> root)const;
    
    void remove_roots(int i, int j);
    
    int root_count() const{return roots.size();};
    
    double compute_ln_likelihood(double *clv, unsigned int *scale_buffer,
                                 const Partition *p);
    
    double compute_ln_likelihood(double *clv, unsigned int *scale_buffer,
                                        const pll_partition_t *p);
    
    void updateIndexesActiveGametes();
    
    PopulationSet* getPopulationSet()const  {return populationSet;};
    
    Population* getPopulationByIndex(size_t i) const {return populationSet->getPopulationbyIndex(i);};
    
    int getNumberPopulations() const {return populationSet->numClones;};
    
    long double getHeight()const{return height;}
    
    void setHeight(long double  newHeight)
    {
        assert(newHeight>=height);
        height= newHeight;
        
    }
    unsigned int getNextAvailable()const{return nextAvailable;}
    
    void increaseNextAvailable(){nextAvailable++;}
    
    
    ~State();
};
#endif /* state_hpp */
