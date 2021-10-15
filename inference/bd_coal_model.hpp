//
//  struct_coal_model.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 05/10/2021.
//

#ifndef struct_coal_model_hpp
#define struct_coal_model_hpp

#include "spf.hpp"
#include "poset_smc_params.hpp"
#include "bd_coal_state.hpp"
#include "utils.hpp"
#include "msa.hpp"
#include "bd_coal_model_params.hpp"
#include <stdio.h>

using namespace std;

class BDCoalModel : public ProblemSpecification<BDCoalState, BDCoalPriorParams>
{
   
    MSA *inputMSA;
public:
    BDCoalModel(MSA *inputMSA);
    unsigned long num_iterations();
    shared_ptr<BDCoalState> propose_initial(gsl_rng *random, double &log_w, BDCoalPriorParams &params);
    shared_ptr<BDCoalState> propose_next(gsl_rng *random, unsigned int t, const BDCoalState &curr, double &log_w, BDCoalPriorParams &params);
    double log_weight(unsigned int t, const shared_ptr<ParticleGenealogy<BDCoalState> > &genealogy, const BDCoalPriorParams &params);
    static void generate_data(gsl_rng *random, size_t T, BDCoalPriorParams &params, vector<double> &latent, vector<double> &obs);
    ~BDCoalModel(){};
};


#endif /* struct_coal_model_hpp */
