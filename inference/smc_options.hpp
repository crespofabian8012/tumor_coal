//
//  smc_options.hpp
//  run
//


#ifndef smc_options_hpp
#define smc_options_hpp


#include <gsl/gsl_rng.h>

#include <vector>
class SMCOptions
{
public:
    unsigned int num_particles = 1000;
    unsigned int max_virtual_particles = 1000000;

    enum ResamplingScheme {
        MULTINOMIAL = 0,
        STRATIFIED = 1,
        SYSTEMATIC = 2
    };
    ResamplingScheme resampling_scheme = MULTINOMIAL;
    double ess_threshold = 0.5;
    bool track_population = false;
    bool resample_last_round = false;
    bool use_SPF = false;
    bool debug = false;
    bool csmc_set_partile_population = true;

    long main_seed = 1;
    long resampling_seed = 2;
    
    gsl_rng *main_random = 0;
    gsl_rng *resampling_random = 0;
    
    std::vector<gsl_rng *> proposal_randoms;

    unsigned int num_threads = 1;
    
    void init();
    ~SMCOptions();
};

#endif /* smc_options_hpp */
