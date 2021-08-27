//
//  hmc_nuts.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 23/08/2021.
//

#include "hmc_nuts.hpp"
//#include "utils.hpp"


//#include <unsupported/Eigen/MatrixFunctions>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/poisson.hpp>
using Eigen::Dynamic;
using namespace std;

Eigen::MatrixXd sample_nuts_cpp(const Eigen::VectorXi v_n,
                                const Eigen::MatrixXd& W, Eigen::VectorXd current_q,
                                double alpha = 1, double beta = 1,
                                float epsilon = 0.01,
                                int iter=100){
  
  int K = W.cols();
  int F = W.rows();
  int MAXDEPTH = 10;
  
  std::default_random_engine generator;
  std::uniform_real_distribution<double> unif01(0.0,1.0);
  std::normal_distribution<double> normal(0,1);
  std::exponential_distribution<double> exp1(1);

  //const Eigen::VectorXd norms_W = u_ * W;
  // Store fixed data and parameters
  posterior_params postparams;
  const Eigen::VectorXd u_(F);
  postparams.W           = W;
  postparams.v_n      = v_n;
  postparams.norms_W = u_ * W;
  postparams.alpha      = alpha;
  postparams.beta      = beta;

  Eigen::MatrixXd h_n_samples(K, iter);   // traces of p
  Eigen::VectorXd p0(K);                  // initial momentum

  current_q = Eigen::log(current_q.array());         // Transform to unrestricted space
  h_n_samples.col(1) = current_q;
  
  pq_point z(K);

  // Used to compute the NUTS generalized stopping criterion
  Eigen::VectorXd rho(K);

  // Transition
  for(int i=2; i<iter; i++){
    std::cout << "******************* Sample: " << i << std::endl;

    nuts_util util;

    // Sample new momentum (K independent standard normal variates)
    for(int k=0; k<K; k++){
      p0[k] = normal(generator);
    }

    // Initialize the path. Proposed sample,
    // and leftmost/rightmost position and momentum
    ////////////////////////
    z.q = current_q;
    z.p = p0;
    pq_point z_plus(z);
    pq_point z_minus(z);
    pq_point z_propose(z);
    
    // Utils o compute NUTS stop criterion
    Eigen::VectorXd p_sharp_plus = z.p;
    Eigen::VectorXd p_sharp_dummy = p_sharp_plus;
    Eigen::VectorXd p_sharp_minus = p_sharp_plus;
    Eigen::VectorXd rho(z.p);

    // Hamiltonian
    // Joint log probability of position q and momentum p
    float joint = loglike_cpp(current_q, postparams) - 0.5* (p0.array() * p0.array()).sum();

    // Slice variable
    ///////////////////////
    // Sample the slice variable: u ~ uniform([0, exp(joint)]).
    // Equivalent to: (log(u) - joint) ~ exponential(1).
    // logu = joint - exprnd(1);
    std::exponential_distribution<double> exp1(1);
    float random = exp1(generator);
    util.log_u = joint - random;

    int n_valid = 1;
    util.criterion = true;

    // Build a trajectory until the NUTS criterion is no longer satisfied
    int depth_ = 0;
    int divergent_ = 0;
    util.n_tree = 0;
    util.sum_prob = 0;


    // Build a balanced binary tree until the NUTS criterion fails
    while(util.criterion && (depth_ < MAXDEPTH)){
      std::cout << "*****depth : " << depth_  << std::endl;

      // Build a new subtree in the chosen direction
      // (Modifies z_propose, z_minus, z_plus)
      Eigen::VectorXd rho_subtree(K);
      rho_subtree= Eigen::VectorXd::Zero(K);

      // Build a new subtree in a random direction
      util.sign = 2 * (unif01(generator) < 0.5) - 1;
      int n_valid_subtree=0;
      if(util.sign == 1){
             z.pq_point::operator=(z_minus);
           n_valid_subtree = BuildNutsTree(z, z_propose, p_sharp_dummy, p_sharp_plus, rho_subtree, util, depth_, epsilon, postparams, generator);
           z_plus.pq_point::operator=(z);
      } else {
           z.pq_point::operator=(z_plus);
           n_valid_subtree = BuildNutsTree(z, z_propose, p_sharp_dummy, p_sharp_minus, rho_subtree, util, depth_, epsilon, postparams, generator);
           z_minus.pq_point::operator=(z);
      }
      //if(!valid_subtree) break;
      
      ++depth_;  // Increment depth.
       

      if(util.criterion){
        // Use Metropolis-Hastings to decide whether or not to move to a
        // point from the half-tree we just generated.
        double subtree_prob = std::min(1.0, static_cast<double>(n_valid_subtree)/n_valid);

        if(unif01(generator) < subtree_prob){
          current_q = z_propose.q; // Accept proposal (it will be THE new sample when s=0)
        }
      }

      // Update number of valid points we've seen.
      n_valid += n_valid_subtree;

      // Break when NUTS criterion is no longer satisfied
      rho += rho_subtree;
      util.criterion = util.criterion && compute_criterion(p_sharp_minus, p_sharp_plus, rho);
    } // end while
    

    h_n_samples.col(i) = current_q;
    
  } // end for

  h_n_samples = h_n_samples.transpose();
  
  return(Eigen::exp(h_n_samples.array()));
}
void leapfrog(pq_point &z, float epsilon, posterior_params& postparams){
  
  z.p += epsilon * 0.5 * grad_loglike_cpp(z.q,
                                          postparams.v_n,
                                          postparams.W,
                                          postparams.norms_W,
                                          postparams.alpha,
                                          postparams.beta);
  z.q += epsilon * z.p;
  z.p += epsilon * 0.5 * grad_loglike_cpp(z.q,
                                          postparams.v_n,
                                          postparams.W,
                                          postparams.norms_W,
                                          postparams.alpha,
                                          postparams.beta);
}

// U-Turn criterion in the generalized form applicable to Riemanian spaces
// See Betancourt's Conceptual HMC (page 58)
bool compute_criterion(Eigen::VectorXd& p_sharp_minus,
                               Eigen::VectorXd& p_sharp_plus,
                               Eigen::VectorXd& rho) {
  return (p_sharp_plus.array() * rho.array()).sum()  > 0 && (p_sharp_minus.array() * rho.array()).sum() > 0;
      }
int BuildNutsTree(pq_point& z, pq_point& z_propose,
              Eigen::VectorXd& p_sharp_left,
              Eigen::VectorXd& p_sharp_right,
              Eigen::VectorXd& rho,
              nuts_util& util,
              int depth, float epsilon,
              posterior_params& postparams,
              std::default_random_engine& generator){

  //std::cout << "\n Tree direction:" << util.sign << " Depth:" << depth << std::endl;

  int K = z.q.rows();

  //std::default_random_engine generator;
  std::uniform_real_distribution<double> unif01(0.0,1.0);
  int F = postparams.W.rows();
  float delta_max = 1000; // Recommended in the NUTS paper: 1000
  
  // Base case - take a single leapfrog step in the direction v
  if(depth == 0){
    leapfrog(z, util.sign * epsilon, postparams);
    float joint = loglike_cpp(z.q, postparams) - 0.5 * (z.p.array() * z.p.array()).sum();
    int valid_subtree = (util.log_u <= joint);    // Is the new point in the slice?
    util.criterion = util.log_u - joint < delta_max; // Is the simulation wildly inaccurate? // TODO: review
    util.n_tree += 1;

    z_propose = z;
    rho += z.p;
    p_sharp_left = z.p;  // p_sharp = inv(M)*p (Betancourt 58)
    p_sharp_right = p_sharp_left;

    return valid_subtree;
  }

  // General recursion
  Eigen::VectorXd p_sharp_dummy(K);

  // Build the left subtree
    Eigen::VectorXd rho_left(K);
    rho_left.setZero();
 
    
  int n1 = BuildNutsTree(z, z_propose, p_sharp_left, p_sharp_dummy, rho_left, util, depth-1, epsilon, postparams, generator);

  if (!util.criterion) return 0; // early stopping

  // Build the right subtree
  pq_point z_propose_right(z);
    Eigen::VectorXd rho_right(K);// = Eigen::VectorXd::Zero(K);
    rho_right.setZero();
    
  int n2 = BuildNutsTree(z, z_propose_right, p_sharp_dummy, p_sharp_right, rho_right, util, depth-1, epsilon, postparams, generator);

  // Choose which subtree to propagate a sample up from.
  //double accept_prob = static_cast<double>(n2) / static_cast<double>(n1 + n2);
  double accept_prob = static_cast<double>(n2) / std::max((n1 + n2), 1); // avoids 0/0;
  float rand01 = unif01(generator);
  if(util.criterion && (rand01 < accept_prob)){
    z_propose = z_propose_right;
  }

  // Break when NUTS criterion is no longer satisfied
  Eigen::VectorXd rho_subtree = rho_left + rho_right;
  rho += rho_subtree;
  util.criterion = compute_criterion(p_sharp_left, p_sharp_right, rho);

  int n_valid_subtree = n1 + n2;
  return(n_valid_subtree);
}
// Posterior wrapper
double loglike_cpp(const Eigen::VectorXd& eta_n, const posterior_params& postparams){
  return posterior_eta_cpp(eta_n,
                          postparams.v_n,
                          postparams.W,
                          postparams.alpha,
                          postparams.beta);
}

// Gradient of the log posterior
Eigen::VectorXd grad_loglike_cpp(const Eigen::VectorXd& eta_n, const Eigen::VectorXi& v_n,
                                 const Eigen::MatrixXd& W, const Eigen::VectorXd& norms_W,
                                 float alpha, float beta){
  
  int K = eta_n.size();
  Eigen::VectorXd dh = Eigen::VectorXd::Zero(K);
 Eigen::VectorXd temp =  Eigen::exp(eta_n.array());
  Eigen::VectorXd lambdas = W * temp;
  for(int k=0; k < K; k++){
   
     dh[k] = alpha - exp(eta_n[k])*(beta + norms_W[k]) + (W.col(k) * exp(eta_n[k]) * lambdas.cwiseInverse()).sum();
  }
  return dh;
}
double posterior_h_cpp(const Eigen::VectorXd& h_n, const Eigen::VectorXi& v_n, const Eigen::MatrixXd& W,
                        float alpha, float beta){
  double logp = 0;
  
  if((h_n.array() < 0.0).any()){
    return -std::numeric_limits<double>::infinity();
  }
  
  int K = W.cols();
  int F = W.rows();
  
  Eigen::VectorXd lambdas = W * h_n;
  
  for(int f = 0; f < F; f++){
    logp += v_n[f] * log(lambdas[f]) - lambdas[f]; //- std::lgamma(v_n[f]+1);
  }

  for(int k = 0; k < K; k++){
    logp += alpha * log(h_n[k]) - log(h_n[k]) - beta*h_n[k];
  }
  
  return logp;
}

// Transformed Gamma-Poisson log posterior
double posterior_eta_cpp(const Eigen::VectorXd& eta_n, const Eigen::VectorXi& v_n, const Eigen::MatrixXd& W,
                        float alpha, float beta){
   
  Eigen::VectorXd exp_eta = Eigen::exp(eta_n.array());
  double logp = posterior_h_cpp(exp_eta, v_n, W, alpha, beta) + eta_n.sum();
  return logp;
}
