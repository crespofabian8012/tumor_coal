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

#include  <vector>
extern "C"
{
#include "libpll/pll.h"
}
class PosetSMCParams;
class State{
   
    long double height;
    PopulationSet * populations;

    PLLBufferManager *const pll_buffer_manager;

    const pll_partition_t *partition;
    unsigned int num_sites;
    std::vector<std::shared_ptr<PartialTreeNode>> roots;
public:

    State( PosetSMCParams &params);
    std::vector<pll_rnode_t *> getForest() const;

    void initForest( int sampleSize, pll_msa_t *msa, std::vector<int> &positions, ProgramOptions &programOptions);
    std::shared_ptr<PartialTreeNode> getRootAt(int i) const {return roots[i];};
    //std::shared_ptr<State> makeDeepCopy();
    
    State &operator=(const State &original);
    
    std::vector<std::shared_ptr<PartialTreeNode>> propose(int i, int j, double height);
    
    std::vector<std::shared_ptr<PartialTreeNode>> getRoots() const{ return roots;};
    
    double likelihood_factor(std::shared_ptr<PartialTreeNode> root);
    
    void remove_roots(int i, int j);
    
    int root_count() const{return roots.size();};
    
    double compute_ln_likelihood(double *clv, unsigned int *scale_buffer,
                                 const pll_partition_t *p);
    ~State();
};
#endif /* state_hpp */
