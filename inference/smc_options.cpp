//
//  smc_options.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 15/03/2021.
//

#include "smc_options.hpp"
#include "random.h"


void SMCOptions::init()
{
    if (main_random == 0) {
        main_random = Random::generateRandomObject(main_seed);
    }
    if (resampling_random == 0) {
        resampling_random = Random::generateRandomObject(resampling_seed);
    }
    
    for (size_t i = 0; i < num_particles; i++) {
        proposal_randoms.push_back(Random::generateRandomObject(gsl_rng_get(main_random)));
    }
}

SMCOptions::~SMCOptions()
{
    //delete main_random;
    //delete resampling_random;
    gsl_rng_free(main_random);
    gsl_rng_free(resampling_random);
    for (size_t i = 0; i < num_particles; i++) {
        gsl_rng_free(proposal_randoms[i]);
        //delete proposal_randoms[i];
    }
}
