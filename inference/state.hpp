//
//  state.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 16/03/2021.
//

#ifndef state_hpp
#define state_hpp

#include "mcmc_parameter.hpp"

#include  <vector>
extern "C"
{
#include "libpll/pll.h"
}


class State{
   
    long double currentModelTime;
    PopulationSet * populations;
    std::vector<pll_rtree_t *> forest;
    const pll_partition_t *partition;
    unsigned int num_sites;
    std::vector<std::shared_ptr<pll_rtree_t*>> roots;
public:
    State(long double currentModelTime);
    State(unsigned int num_sites, const pll_partition_t *partition);
    State(unsigned int num_sites, const pll_partition_t *partition,  int sampleSize, pll_msa_t *msa, std::vector<int> &positions, ProgramOptions &programOptions );
    std::vector<pll_rtree_t *> getForest() const;
    void setForest(std::vector<pll_rtree_t *> forestPar);
    
    void initForest( int sampleSize, pll_msa_t *msa, std::vector<int> &positions, ProgramOptions &programOptions);
    pll_rtree_t * getTreeAt(int i) const;
    std::shared_ptr<State> makeDeepCopy();
    
    State &operator=(const State &original);
    
    std::vector<std::shared_ptr<pll_rtree_t*>> propose(int i, int j, double height);
    
    std::vector<std::shared_ptr<pll_rtree_t*>> getRoots() const{ return roots;};
    
    double likelihood_factor(std::shared_ptr<pll_rtree_t*> root);
    
    void remove_roots(int i, int j);
    ~State();
};
#endif /* state_hpp */
