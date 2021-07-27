//
//  particle_population.hpp
//  run


#ifndef particle_population_hpp
#define particle_population_hpp

//#include "utils.hpp"
//#include "constants.hpp"
//#include "utils.hpp"

//#include <vector>
//#include <math.h>
//
//#include "utils.hpp"
//#include "constants.hpp"
//
//template <class P>//P is the class ParticleSMC
//class ParticlePopulation
//{
//    std::vector<std::shared_ptr<P> > *particles;
//    std::vector<double> *log_weights = 0;
//    std::vector<double> *normalized_weights = 0;
//    unsigned long num_particles;
//    double ess = std::numeric_limits<double>::quiet_NaN();
//    double sum_weights;
//    double sum_weights_squared;
//    double log_norm = std::numeric_limits<double>::quiet_NaN();
//    double log_Z_ratio = std::numeric_limits<double>::quiet_NaN();
//    bool resampled = false;
//
//public:
//    ParticlePopulation(std::vector<std::shared_ptr<P> > *particles);
//    ParticlePopulation(std::vector<std::shared_ptr<P> > *particles, double log_norm);
//    ParticlePopulation(std::vector<std::shared_ptr<P> > *particles, std::vector<double> *log_weights);
//    ParticlePopulation(std::vector<std::shared_ptr<P> > *particles, std::vector<double> *log_weights, double log_norm, double log_Z_ratio);
//    double get_log_weight(size_t k);
//    P &get_particle(size_t k);
//    std::vector<std::shared_ptr<P> > *get_particles();
//    std::vector<double> *get_normalized_weights();
//    double get_ess();
//    double get_log_norm();
//    double get_log_Z_ratio();
//    bool is_resampled();
//    size_t get_num_particles();
//    static ParticlePopulation<P>* construct_equally_weighted_population(int num_particles);
//    static ParticlePopulation<P>* construct_equally_weighted_population(std::vector<P> *particles);
//    ~ParticlePopulation();
//};
//
//template <class P>
//ParticlePopulation<P>::ParticlePopulation(std::vector<std::shared_ptr<P> > *particles, std::vector<double> *log_weights)
//{
//    this->particles = particles;
//    this->log_weights = log_weights;
//    num_particles = particles->size();
//}
//
//template <class P>
//ParticlePopulation<P>::ParticlePopulation(std::vector<std::shared_ptr<P> > *particles, std::vector<double> *log_weights, double log_norm, double log_Z_ratio) :
//ParticlePopulation(particles, log_weights)
//{
//    this->log_norm = log_norm;
//    this->log_Z_ratio = log_Z_ratio;
//}
//
//template <class P>
//ParticlePopulation<P>::ParticlePopulation(std::vector<std::shared_ptr<P> > *particles)
//{
//    this->particles = particles;
//    this->num_particles = particles->size();
//    this->resampled = true;
//}
//
//template <class P>
//ParticlePopulation<P>::ParticlePopulation(std::vector<std::shared_ptr<P> > *particles, double log_norm) :
//ParticlePopulation(particles)
//{
//    this->log_norm = log_norm;
//}
//
//template <class P>
//double ParticlePopulation<P>::get_log_weight(size_t k)
//{
//    return log_weights->at(k);
//}
//
//template <class P>
//std::vector<std::shared_ptr<P> > *ParticlePopulation<P>::get_particles()
//{
//    return particles;
//}
//
//template <class P>
//P &ParticlePopulation<P>::get_particle(size_t k)
//{
//    P &p = *particles->at(k).get();
//    return p;
//}
//
//template <class P>
//double ParticlePopulation<P>::get_ess()
//{
//    // compute ESS
//    ess = 0.0;
//    sum_weights = 0.0;
//    sum_weights_squared = 0.0;
//    if (normalized_weights == 0) {
//        get_normalized_weights();
//    }
//    for (int i = 0; i < num_particles; i++) {
//        ess += (*normalized_weights)[i] * (*normalized_weights)[i];
//    }
//    ess = 1.0/ess;
//    ess /= num_particles;
//    return ess;
//}
//
//template <class P>
//double ParticlePopulation<P>::get_log_norm()
//{
//    if (log_norm == std::numeric_limits<double>::quiet_NaN()) {
//        // compute ESS
//        get_ess();
//    }
//    return log_norm;
//}
//
//template <class P>
//double ParticlePopulation<P>::get_log_Z_ratio()
//{
//    return log_Z_ratio;
//}
//
//template <class P>
//std::vector<double> *ParticlePopulation<P>::get_normalized_weights()
//{
//    if (normalized_weights == 0) {
//        normalized_weights = new std::vector<double>(num_particles);
//        if (is_resampled()) {
//            double w = 1./num_particles;
//            for (size_t k = 0; k < num_particles; k++) {
//                (*normalized_weights)[k] = w;
//            }
//        } else {
//            if (log_norm == std::numeric_limits<double>::quiet_NaN()) {
//                // normalize the log weights, compute log_norm
//                log_norm = Utils::normalize(*log_weights, *normalized_weights);
//            } else {
//                normalize(*log_weights, *normalized_weights, log_norm);
//            }
//        }
//    }
//    return normalized_weights;
//}
//
//template <class P>
//bool ParticlePopulation<P>::is_resampled()
//{
//    return resampled;
//}
//
//template <class P>
//size_t ParticlePopulation<P>::get_num_particles()
//{
//    return num_particles;
//}
//
//template<class P>
//ParticlePopulation<P>* ParticlePopulation<P>::construct_equally_weighted_population(int num_particles)
//{
//    std::vector<P> *particles = new std::vector<P>();
//    std::vector<double> *log_weights = new std::vector<double>();
//    double log_w = -log(num_particles);
//    for (int n = 0; n < num_particles; n++)
//    {
//        particles->push_back(0);
//        log_weights->push_back(log_w);
//    }
//    ParticlePopulation<P> *pop = new ParticlePopulation<P>(particles, log_weights);
//    return pop;
//}
//
//template<class P>
//ParticlePopulation<P>* ParticlePopulation<P>::construct_equally_weighted_population(std::vector<P> *particles)
//{
//    unsigned int num_particles = particles->size();
//    std::vector<double> *log_weights = new std::vector<double>();
//    double log_w = - log(num_particles);
//    for (int n = 0; n < num_particles; n++)
//    {
//        log_weights->push_back(log_w);
//    }
//    ParticlePopulation<P> *pop = new ParticlePopulation<P>(particles, log_weights);
//    pop->resampled = true;
//    return pop;
//}
//
//template <class P>
//ParticlePopulation<P>::~ParticlePopulation()
//{
//    if (!particles)
//        delete particles;
//    if (!log_weights)
//        delete log_weights;
//    if (!normalized_weights)
//        delete normalized_weights;
//}
#endif /* particle_population_hpp */
