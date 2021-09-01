//
//  hmc_nuts.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 23/08/2021.
//

#ifndef hmc_nuts_hpp
#define hmc_nuts_hpp
#include <Eigen/Eigen>
#include <stdio.h>

#include <numeric>
#include <iostream>
#include <math.h>
#include <random>
#include <boost/random/uniform_01.hpp>



struct pq_point {
    
    Eigen::VectorXd q;
    Eigen::VectorXd p;
    Eigen::VectorXd g;
    
    explicit pq_point(int n): q(n), p(n) {}
    pq_point(const pq_point& z): q(z.q.size()), p(z.p.size()){
        q = z.q;
        p = z.p;
    }
    
    pq_point& operator= (const pq_point& z) {
        if (this == &z)
            return *this;
        
        q = z.q;
        p = z.p;
        
        return *this;
    }
    virtual inline void getParamNames(std::vector<std::string>& modelNames,
                                      std::vector<std::string>& names) {
        names.reserve(q.size() + p.size());
        for (int i = 0; i < q.size(); ++i)
            names.emplace_back(modelNames[i]);
        for (int i = 0; i < p.size(); ++i)
            names.emplace_back(std::string("p_") + modelNames[i]);
    }
    
    virtual inline void get_params(std::vector<double>& values) {
        values.reserve(q.size() + p.size());
        for (int i = 0; i < q.size(); ++i)
            values.push_back(q[i]);
        for (int i = 0; i < p.size(); ++i)
            values.push_back(p[i]);
    }
    
    };
    
    
    struct nuts_util {
        // Constants through each recursion
        double log_u; // uniform sample
        double H0;     // Hamiltonian of starting point?
        int sign;     // direction of the tree in a given iteration/recursion
        
        // Aggregators through each recursion
        int n_tree;
        double sum_prob;
        bool criterion;
        
        // just to guarantee bool initializes to valid value
        nuts_util() : criterion(false) { }
    };
    
    
    struct posterior_params {
        Eigen::VectorXi v_n;
        Eigen::MatrixXd W;
        Eigen::MatrixXd norms_W;
        double alpha;
        double beta;
    };

template <class Model, template <class, class> class Hamiltonian,
    template <class> class Integrator, class BaseRNG>
class NUTS{
protected:
     typename Hamiltonian<Model, BaseRNG>::PointType z_;
     Integrator<Hamiltonian<Model, BaseRNG> > integrator_;
     Hamiltonian<Model, BaseRNG> hamiltonian_;

     BaseRNG& rand_int_;

     // Uniform(0, 1) RNG
     boost::uniform_01<BaseRNG&> rand_uniform_;

     double nom_epsilon_;
     double epsilon_;
     double epsilon_jitter_;
    
    unsigned int random_seed;
    unsigned int chain;
    double init_radius;
    int num_warmup;
    int num_samples;
    int num_thin;
    bool save_warmup;
    int refresh;
    double stepsize;
    double stepsize_jitter;
    int depth_;
    int max_depth_;
    double max_deltaH_;
    int n_leapfrog_;
    bool divergent_;
    double energy_;
    
public:
    NUTS(const Model& model, BaseRNG& rng):
          z_(model.num_params_r()),
           integrator_(),
           hamiltonian_(model),
           rand_int_(rng),
           rand_uniform_(rand_int_),
           nom_epsilon_(0.1),
           epsilon_(nom_epsilon_),
           epsilon_jitter_(0.0),
           depth_(0),
           max_depth_(5),
           max_deltaH_(1000),
           n_leapfrog_(0),
           divergent_(false),
           energy_(0) {};

    NUTS(const Model& model, BaseRNG& rng, Eigen::VectorXd& inv_e_metric):
             z_(model.num_params_r()),
              integrator_(),
              hamiltonian_(model),
              rand_int_(rng),
              rand_uniform_(rand_int_),
              nom_epsilon_(0.1),
              epsilon_(nom_epsilon_),
              epsilon_jitter_(0.0),
              depth_(0),
              max_depth_(5),
              max_deltaH_(1000),
              n_leapfrog_(0),
              divergent_(false),
              energy_(0) {};
    
    void setMetric(Eigen::VectorXd &inv_metric);
    void setNominalStepsize(double stepsize);
    void setStepsizeJitter(double stepsize_jitter);
    void setMaxDepth(int max_depth);
    
    bool buildTree(int depth, pq_point& z_propose, Eigen::VectorXd& p_sharp_beg,
           Eigen::VectorXd& p_sharp_end, Eigen::VectorXd& rho,
           Eigen::VectorXd& p_beg, Eigen::VectorXd& p_end, double H0,
           double sign, int& n_leapfrog, double& log_sum_weight,
                         double& sum_metro_prob);
 ~NUTS() {}
};
    
double loglike_cpp(const Eigen::VectorXd& eta_n, const posterior_params& postparams);
    Eigen::VectorXd grad_loglike_cpp(const Eigen::VectorXd& eta_n, const Eigen::VectorXi& v_n,
                                     const Eigen::MatrixXd& W, const Eigen::VectorXd& norms_W,
                                     float alpha, float beta);
bool computeCriterion(Eigen::VectorXd& p_sharp_minus,
                           Eigen::VectorXd& p_sharp_plus,
                           Eigen::VectorXd& rho);
int BuildNutsTree(pq_point& z, pq_point& z_propose,
                      Eigen::VectorXd& p_sharp_left,
                      Eigen::VectorXd& p_sharp_right,
                      Eigen::VectorXd& rho,
                      nuts_util& util,
                      int depth, float epsilon,
                      posterior_params& postparams,
                      std::default_random_engine& generator);
double posterior_h_cpp(const Eigen::VectorXd& h_n, const Eigen::VectorXi& v_n, const Eigen::MatrixXd& W,
                           float alpha, float beta);
    
double posterior_eta_cpp(const Eigen::VectorXd& eta_n, const Eigen::VectorXi& v_n, const Eigen::MatrixXd& W,
                             float alpha, float beta);
Eigen::MatrixXd runNutsSampler(const Eigen::VectorXi v_n,
                                    const Eigen::MatrixXd& W, Eigen::VectorXd current_q,
                                    double alpha, double beta,
                                    float epsilon,
                                    int iter);
#endif /* hmc_nuts_hpp */
    
