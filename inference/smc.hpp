//
//  smc.hpp
//  run
//
//  Created by  Seong-Hwan Jun
//

#ifndef smc_hpp
#define smc_hpp


#include "smc_options.hpp"
#include "utils.hpp"
#include "resampling.hpp"
#include "state.hpp"
#include "particle_population.cpp"


template <class S, class P> class ProblemSpecification
{
public:
    virtual unsigned long num_iterations() = 0;
    // propose new sample and return it with log weight
    //virtual std::pair<S, double> *propose_initial(gsl_rng *random, P &params) = 0;
    // design decision to ponder upon:
    // instead of passing in const S &curr, pass S, which means that the original data
    // is not to be modified
    // this approach is nice if default copy constructor can be used (which might be the case for most instances) but otherwise, the user needs to implement a copy constructor
    // for the object
    // alternative is to pass const reference: const S &curr.
    // in this case, the user has to first make a copy of the state and modify the copied state
    // in the function body
    //virtual std::pair<S, double> *propose_next(gsl_rng *random, int t, S curr, P &params) = 0;

    // return the state, set the logw
    virtual std::shared_ptr<S> propose_initial(gsl_rng *random, double &logw, P &params) = 0;
    // current state is not to be modified -- template class S should provide const functions so that necessary components can be accessed
    virtual std::shared_ptr<S> propose_next(gsl_rng *random, unsigned int t, const S &curr, double &logw, P &params) = 0;
    //virtual double log_weight(unsigned int t, const shared_ptr<ParticleGenealogy<S> > &genealogy, const P &p) = 0;

    // call back function, called by CSMC: user can sharpen the proposal using the sampled particles and log_weights if desired
    // can be useful for sharpening the proposal distribution
//    virtual void set_particle_population(const vector<vector<shared_ptr<S> > > &particles, const vector<vector<double> > &log_weights, const vector<double> &log_norms) { }
    virtual void set_particle_population(const std::vector<std::shared_ptr<S> > &particles) =0;
    virtual ~ProblemSpecification() { }
};
template <class S, class P> class SMC
{
    SMCOptions &options;
    ProblemSpecification<S, P> &proposal; // pointer to proposal object to be passed into SMC constructor

    double log_marginal_likelihood = 0;
    ParticlePopulation<S> *propose(gsl_rng *random, ParticlePopulation<S> *pop, int iter, int num_proposals, P &params);
    ParticlePopulation<S> *resample(gsl_rng *random, SMCOptions::ResamplingScheme resampling_scheme, ParticlePopulation<S> *pop, int N);
    std::vector<ParticlePopulation<S> *> *populations = 0;

public:
    SMC(ProblemSpecification<S, P> &proposal, SMCOptions &options);
    void run_smc(P &params);
    double get_log_marginal_likelihood();
    std::vector<S*> sample(gsl_rng *random);
    inline ParticlePopulation<S>* get_curr_population() { return populations->back(); }
    inline ParticlePopulation<S> *get_population(size_t i) { return populations->at(i); }
    ~SMC();
};

template <class S, class P>
SMC<S,P>::SMC(ProblemSpecification<S, P> &proposal, SMCOptions &options) :
options(options), proposal(proposal)
{
    this->options.init();
    this->populations = new std::vector<ParticlePopulation<S> *>();
}

template <class S, class P>
void SMC<S,P>::run_smc(P &params)
{
    unsigned long R = proposal.num_iterations();

    gsl_rng *random = options.main_random;

    log_marginal_likelihood = 0.0;
    ParticlePopulation<S> *curr_pop = 0;
    for (int r = 0; r < R; r++)
    {
        if (options.debug)
            std::cout << "iter: " << r << std::endl;
        curr_pop = propose(random, curr_pop, r, options.num_particles, params);
        populations->push_back(curr_pop);
        if (!options.track_population && r > 0) {
            if (!(*populations)[r-1]){
                //make a copy of population instead ?
                delete (*populations)[r-1]; // delete the ParticlePopulation (particles, log_weights, and normalized_weights)
                
            }
                
        }
        log_marginal_likelihood += curr_pop->get_log_Z_ratio();
        if (r == R - 1 && options.resample_last_round == false) {
            break;
        }

        if (options.ess_threshold == 0.0) // no resampling (essneitally SIS)
            continue;
        if (options.ess_threshold >= 1.0 || curr_pop->get_ess() <= options.ess_threshold) {
            // resample
            curr_pop = resample(random, options.resampling_scheme, curr_pop, options.num_particles);
            delete (*populations)[r]; // delete pre-resampling population
            (*populations)[r] = curr_pop;
        }
    }
}

template <class S, class P>
ParticlePopulation<S> *SMC<S,P>::propose(gsl_rng *random, ParticlePopulation<S> *pop, int iter, int num_proposals, P &params)
{
    bool is_resampled = true;
    if (pop != 0) {
        is_resampled = pop->is_resampled();
    }

    std::vector<std::shared_ptr<S> > *new_particles = new std::vector<std::shared_ptr<S> >();
    std::vector<double> *new_log_weights = new std::vector<double>(options.num_particles);

    double log_norm = DOUBLE_NEG_INF;
    double log_Z_ratio = DOUBLE_NEG_INF;
    const double log_n_proposals = -log(num_proposals);
    
    for (size_t n = 0; n < num_proposals; n++)
    {
        if (pop == 0) {
            new_particles->push_back(proposal.propose_initial(random, (*new_log_weights)[n], params));
        } else {
            S &parent = pop->get_particle(n);
            new_particles->push_back(proposal.propose_next(random, iter, parent, (*new_log_weights)[n], params));
        }

        double log_prev_w = log_n_proposals;
        if (is_resampled) {
            log_Z_ratio = Utils::log_add(log_Z_ratio, (*new_log_weights)[n] + log_prev_w);
        } else {
            log_prev_w = pop->get_log_weight(n);
            log_Z_ratio = log_add(log_Z_ratio, (*new_log_weights)[n] + log_prev_w - pop->get_log_norm());
        }
        (*new_log_weights)[n] += log_prev_w;
        log_norm = Utils::log_add(log_norm, new_log_weights->at(n));
    }

    // note: log_norm is an approximation of Z_t/Z_{t-1}
    ParticlePopulation<S> *new_pop = new ParticlePopulation<S>(new_particles, new_log_weights, log_norm, log_Z_ratio);
    return new_pop;
}

template <class S, class P>
ParticlePopulation<S>* SMC<S,P>::resample(gsl_rng *random, SMCOptions::ResamplingScheme resampling_scheme, ParticlePopulation<S> *pop, int N)
{
    unsigned int *indices = new unsigned int[N];
    switch (resampling_scheme)
    {
        case SMCOptions::ResamplingScheme::MULTINOMIAL:
            multinomial_resampling(random, pop->get_normalized_weights(), N, indices);
            break;
        case SMCOptions::ResamplingScheme::STRATIFIED:
            stratified_resampling(random, pop->get_normalized_weights(), N, indices);
            break;
        case SMCOptions::ResamplingScheme::SYSTEMATIC:
            systematic_resampling(random, pop->get_normalized_weights(), N, indices);
            break;
    }

    std::vector<std::shared_ptr<S> > *prev_particles = pop->get_particles();
    std::vector<std::shared_ptr<S> > *particles = new std::vector<std::shared_ptr<S> >();
    for (int n = 0; n < N; n++)
    {
        particles->push_back(prev_particles->at(indices[n])); // move ownership to new ParticlePopulation
    }

    ParticlePopulation<S> *new_pop = new ParticlePopulation<S>(particles, pop->get_log_norm());

    delete [] indices;
    return new_pop;
}

template <class S, class P>
double SMC<S,P>::get_log_marginal_likelihood()
{
    return log_marginal_likelihood;
}

template <class S, class P>
std::vector<S*> SMC<S,P>::sample(gsl_rng *random)
{
    std::vector<S*> ret;
    return ret;
}

template <class S, class P>
SMC<S,P>::~SMC()
{
}

#endif /* smc_hpp */
