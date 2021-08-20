//
//  poset_smc.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 16/03/2021.
//

#ifndef poset_smc_hpp
#define poset_smc_hpp




#include "spf.hpp"

#include "poset_smc_params.hpp"

class State;
class PosetSMCParams;

class PosetSMC : public ProblemSpecification<State, PosetSMCParams>
{
    size_t numClones;
    size_t num_iter;
public:
    PosetSMC(size_t numClones, size_t num_iter);
    unsigned long num_iterations() override;
    std::shared_ptr<State> propose_initial(gsl_rng *random, double &log_w, PosetSMCParams &params) override;
    std::shared_ptr<State> propose_next(gsl_rng *random, unsigned int t, const State &curr, double &log_w, PosetSMCParams &params) override;
   
   // static void generate_data(gsl_rng *random, size_t T, SMCOptions &params, std::vector<double> &latent, std::vector<double> &obs)override ;
 
    double log_weight(unsigned int t, const shared_ptr<ParticleGenealogy<State> > &genealogy, const PosetSMCParams &p) override;
    
    void set_particle_population(const vector<shared_ptr<State> > &particles) override;
    
    ~PosetSMC(){};
};


#endif /* poset_smc_hpp */



