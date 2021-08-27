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

struct pq_point {
    
    Eigen::VectorXd q;
    Eigen::VectorXd p;
    
    explicit pq_point(int n): q(n), p(n) {}
    pq_point(const pq_point& z): q(z.q.size()), p(z.p.size()) {
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
    double loglike_cpp(const Eigen::VectorXd& eta_n, const posterior_params& postparams);
    Eigen::VectorXd grad_loglike_cpp(const Eigen::VectorXd& eta_n, const Eigen::VectorXi& v_n,
                                     const Eigen::MatrixXd& W, const Eigen::VectorXd& norms_W,
                                     float alpha, float beta);
    bool compute_criterion(Eigen::VectorXd& p_sharp_minus,
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
#endif /* hmc_nuts_hpp */
    
