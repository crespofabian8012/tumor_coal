//
//  hmc_nuts.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 23/08/2021.
//

#ifndef hmc_nuts_hpp
#define hmc_nuts_hpp


#include <Eigen/Eigen>


#include <numeric>
#include <iostream>
#include <math.h>
#include <random>
#include <boost/random/uniform_01.hpp>

#include "utils.hpp"
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
    typename Hamiltonian<Model, BaseRNG>::PointType z;
    Integrator<Hamiltonian<Model, BaseRNG> > integrator;
    Hamiltonian<Model, BaseRNG> hamiltonian;
    
    BaseRNG& rand_int;
    
    // Uniform(0, 1) RNG
    boost::uniform_01<BaseRNG&> rand_uniform;
    
    double nom_epsilon;
    double epsilon;
    double epsilon_jitter;
    
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
    int depth;
    int max_depth;
    double max_deltaH;
    int n_leapfrog;
    bool divergent;
    double energy;
    
public:
    NUTS(const Model& model, BaseRNG& rng):
    z(model.num_params_r()),
    integrator(),
    hamiltonian(model),
    rand_int(rng),
    rand_uniform(rand_int),
    nom_epsilon(0.1),
    epsilon(nom_epsilon),
    epsilon_jitter(0.0),
    depth(0),
    max_depth(5),
    max_deltaH(1000),
    n_leapfrog(0),
    divergent(false),
    energy(0) {};
    
    NUTS(const Model& model, BaseRNG& rng, Eigen::VectorXd& inv_e_metric):
    z(model.num_params_r()),
    integrator(),
    hamiltonian(model),
    rand_int(rng),
    rand_uniform(rand_int),
    nom_epsilon(0.1),
    epsilon(nom_epsilon),
    epsilon_jitter(0.0),
    depth(0),
    max_depth(5),
    max_deltaH(1000),
    n_leapfrog(0),
    divergent(false),
    energy(0) {};
    
    void setMetric(Eigen::VectorXd &inv_metric);
    void setNominalStepsize(double stepsize);
    void setStepsizeJitter(double stepsize_jitter);
    void setMaxDepth(int max_depth);

    bool computeCriterion(Eigen::VectorXd& p_sharp_minus,
                                 Eigen::VectorXd& p_sharp_plus,
                                 Eigen::VectorXd& rho) {
              return (p_sharp_plus.array() * rho.array()).sum()  > 0 && (p_sharp_minus.array() * rho.array()).sum() > 0;
        }
    bool buildTree(int depth, pq_point& z_propose, Eigen::VectorXd& p_sharp_beg,
                      Eigen::VectorXd& p_sharp_end, Eigen::VectorXd& rho,
                      Eigen::VectorXd& p_beg, Eigen::VectorXd& p_end, double H0,
                      double sign, int& n_leapfrog, double& log_sum_weight,
                      double& sum_metro_prob)
    {
        // Base case
        if (depth == 0) {
          integrator.evolve(z, hamiltonian,
                                   sign * epsilon);
          ++n_leapfrog;

          double h = hamiltonian.H(z);
          if (std::isnan(h))
            h = std::numeric_limits<double>::infinity();

          if ((h - H0) > max_deltaH)
            divergent = true;

          log_sum_weight = Utils::logSumExp(log_sum_weight, H0 - h);

          if (H0 - h > 0)
            sum_metro_prob += 1;
          else
            sum_metro_prob += std::exp(H0 - h);

          z_propose = z;

          p_sharp_beg = hamiltonian.dtau_dp(z);
          p_sharp_end = p_sharp_beg;

          rho += z.p;
          p_beg = z.p;
          p_end = p_beg;

          return !divergent;
        }
        // General recursion

        // Build the initial subtree
        double log_sum_weight_init = -std::numeric_limits<double>::infinity();

        // Momentum and sharp momentum at end of the initial subtree
        Eigen::VectorXd p_init_end(z.p.size());
        Eigen::VectorXd p_sharp_init_end(z.p.size());

        Eigen::VectorXd rho_init = Eigen::VectorXd::Zero(rho.size());

        bool valid_init
            = buildTree(depth - 1, z_propose, p_sharp_beg, p_sharp_init_end,
                         rho_init, p_beg, p_init_end, H0, sign, n_leapfrog,
                         log_sum_weight_init, sum_metro_prob);

        if (!valid_init)
          return false;

        // Build the final subtree
        pq_point z_propose_final(z);

        double log_sum_weight_final = -std::numeric_limits<double>::infinity();

        // Momentum and sharp momentum at beginning of the final subtree
        Eigen::VectorXd p_final_beg(z.p.size());
        Eigen::VectorXd p_sharp_final_beg(z.p.size());

        Eigen::VectorXd rho_final = Eigen::VectorXd::Zero(rho.size());

        bool valid_final
            = buildTree(depth - 1, z_propose_final, p_sharp_final_beg, p_sharp_end,
                         rho_final, p_final_beg, p_end, H0, sign, n_leapfrog,
                         log_sum_weight_final, sum_metro_prob);

        if (!valid_final)
          return false;

        // Multinomial sample from right subtree
        double log_sum_weight_subtree
            = Utils::logSumExp(log_sum_weight_init, log_sum_weight_final);
       log_sum_weight = Utils::logSumExp(log_sum_weight, log_sum_weight_subtree);

        if (log_sum_weight_final > log_sum_weight_subtree) {
          z_propose = z_propose_final;
        } else {
          double accept_prob
              = std::exp(log_sum_weight_final - log_sum_weight_subtree);
          if (this->rand_uniform_() < accept_prob)
            z_propose = z_propose_final;
        }

        Eigen::VectorXd rho_subtree = rho_init + rho_final;
        rho += rho_subtree;

        // Demand satisfaction around merged subtrees
        bool persist_criterion
            = computeCriterion(p_sharp_beg, p_sharp_end, rho_subtree);

        // Demand satisfaction between subtrees
        rho_subtree = rho_init + p_final_beg;
        persist_criterion
            &= computeCriterion(p_sharp_beg, p_sharp_final_beg, rho_subtree);

        rho_subtree = rho_final + p_init_end;
        persist_criterion
            &= computeCriterion(p_sharp_init_end, p_sharp_end, rho_subtree);

        return persist_criterion;
      }
    

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

