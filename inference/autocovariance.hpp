//
//  autocovariance.hpp
//  tumor_coal
//
//  taken from github Stan
// https://github.com/stan-dev/stan/blob/develop/src/stan/analyze/mcmc/compute_effective_sample_size.hpp
// https://github.com/stan-dev/stan/blob/develop/src/stan/analyze/mcmc/autocovariance.hpp
// https://github.com/stan-dev/stan/blob/develop/src/stan/analyze/mcmc/compute_potential_scale_reduction.hpp



#ifndef autocovariance_hpp
#define autocovariance_hpp

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include <unsupported/Eigen/FFT>

//#include <complex>
//#include <vector>
//#include <algorithm>
//#include <cmath>
//
//#include <limits>
namespace mcmc_utils {

template <typename T, typename DerivedA, typename DerivedB>
void autocorrelation(const Eigen::MatrixBase<DerivedA>& y,
                     Eigen::MatrixBase<DerivedB>& ac, Eigen::FFT<T>& fft) {
  size_t N = y.size();


 int i = 1;
  while(i < N)
         i = i << 1;

  size_t M =   i;//next 2  power
  size_t Mt2 = 2 * M;

  // centered_signal = y-mean(y) followed by N zeros
  Eigen::Matrix<T, Eigen::Dynamic, 1> centered_signal(Mt2);
  centered_signal.setZero();
  centered_signal.head(N) = y.array() - y.mean();

  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> freqvec(Mt2);
  //fft.SetFlag(fft.HalfSpectrum);
  fft.fwd(freqvec, centered_signal);

  freqvec = freqvec.cwiseAbs2();

  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> ac_tmp(Mt2);
  fft.inv(ac_tmp, freqvec);
  //fft.ClearFlag(fft.HalfSpectrum);

  // use "biased" estimate as recommended by Geyer (1992)
  ac = ac_tmp.head(N).real().array() / (N * N * 2);
  ac /= ac(0);
}
template <typename T, typename DerivedA, typename DerivedB>
void autocovariance(const Eigen::MatrixBase<DerivedA>& y,
                    Eigen::MatrixBase<DerivedB>& acov) {
  Eigen::FFT<T> fft;
  autocorrelation(y, acov, fft);

  using boost::accumulators::accumulator_set;
  using boost::accumulators::stats;
  using boost::accumulators::tag::variance;

  accumulator_set<double, stats<variance>> acc;
  for (int n = 0; n < y.size(); ++n) {
    acc(y(n));
  }

  acov = acov.array() * boost::accumulators::variance(acc);

}

template <typename T>
void autocovariance(const std::vector<T>& y, std::vector<T>& acov) {
  size_t N = y.size();
    
  acov.resize(N);

  const Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>> y_map(&y[0], N);
  Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>> acov_map(&acov[0], N);
  autocovariance<T>(y_map, acov_map);
    
  //  std::cout << acov_map << std::endl;
}

inline long  double compute_effective_sample_size(std::vector<const double*> draws,
                                            std::vector<size_t> sizes) {
  int num_chains = sizes.size();
  size_t num_draws = sizes[0];
  for (int chain = 1; chain < num_chains; ++chain) {
    num_draws = std::min(num_draws, sizes[chain]);
  }

  if (num_draws < 4) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // check if chains are constant; all equal to first draw's value
  bool are_all_const = false;
  Eigen::VectorXd init_draw = Eigen::VectorXd::Zero(num_chains);

  for (int chain_idx = 0; chain_idx < num_chains; chain_idx++) {
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1>> draw(
        draws[chain_idx], sizes[chain_idx]);

    for (int n = 0; n < num_draws; n++) {
      if (!std::isfinite(draw(n))) {
        return std::numeric_limits<double>::quiet_NaN();
      }
    }

    init_draw(chain_idx) = draw(0);

    if (draw.isApproxToConstant(draw(0))) {
      are_all_const |= true;
    }
  }

  if (are_all_const) {
    // If all chains are constant then return NaN
    // if they all equal the same constant value
    if (init_draw.isApproxToConstant(init_draw(0))) {
      return std::numeric_limits<double>::quiet_NaN();
    }
  }
   long  double mean_of_means=0.0;
   long  double variance_of_means =0.0;
  Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, 1> acov(num_chains);
  Eigen::VectorXd chain_mean(num_chains);
  Eigen::VectorXd chain_var(num_chains);
  for (int chain = 0; chain < num_chains; ++chain) {
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1>> draw(
        draws[chain], sizes[chain]);
    autocovariance<double>(draw, acov(chain));
    chain_mean(chain) = draw.mean();
    chain_var(chain) = acov(chain)(0) * num_draws / (num_draws - 1);
  }

 double mean_var = chain_var.mean();
  double var_plus = mean_var * (num_draws - 1) / num_draws;
    if (num_chains > 1){
        mean_of_means = (chain_mean.array() / chain_mean.size()).sum();
        variance_of_means = ((chain_mean.array() - mean_of_means) / std::sqrt((chain_mean.size() - 1.0))).square().sum();
       var_plus += variance_of_means;
    }
  Eigen::VectorXd rho_hat_s(num_draws);
  rho_hat_s.setZero();
  Eigen::VectorXd acov_s(num_chains);
  for (int chain = 0; chain < num_chains; ++chain)
    acov_s(chain) = acov(chain)(1);
  double rho_hat_even = 1.0;
  rho_hat_s(0) = rho_hat_even;
  double rho_hat_odd = 1 - (mean_var - acov_s.mean()) / var_plus;
  rho_hat_s(1) = rho_hat_odd;

  // Convert raw autocovariance estimators into Geyer's initial
  // positive sequence. Loop only until num_draws - 4 to
  // leave the last pair of autocorrelations as a bias term that
  // reduces variance in the case of antithetical chains.
  size_t s = 1;
  while (s < (num_draws - 4) && (rho_hat_even + rho_hat_odd) > 0) {
    for (int chain = 0; chain < num_chains; ++chain)
      acov_s(chain) = acov(chain)(s + 1);
    rho_hat_even = 1 - (mean_var - acov_s.mean()) / var_plus;
    for (int chain = 0; chain < num_chains; ++chain)
      acov_s(chain) = acov(chain)(s + 2);
    rho_hat_odd = 1 - (mean_var - acov_s.mean()) / var_plus;
    if ((rho_hat_even + rho_hat_odd) >= 0) {
      rho_hat_s(s + 1) = rho_hat_even;
      rho_hat_s(s + 2) = rho_hat_odd;
    }
    s += 2;
  }

  int max_s = s;
  // this is used in the improved estimate, which reduces variance
  // in antithetic case -- see tau_hat below
  if (rho_hat_even > 0)
    rho_hat_s(max_s + 1) = rho_hat_even;

  // Convert Geyer's initial positive sequence into an initial
  // monotone sequence
  for (int s = 1; s <= max_s - 3; s += 2) {
    if (rho_hat_s(s + 1) + rho_hat_s(s + 2) > rho_hat_s(s - 1) + rho_hat_s(s)) {
      rho_hat_s(s + 1) = (rho_hat_s(s - 1) + rho_hat_s(s)) / 2;
      rho_hat_s(s + 2) = rho_hat_s(s + 1);
    }
  }

  double num_total_draws = num_chains * num_draws;
  // Geyer's truncated estimator for the asymptotic variance
  // Improved estimate reduces variance in antithetic case
  double tau_hat = -1 + 2 * rho_hat_s.head(max_s).sum() + rho_hat_s(max_s + 1);
  return std::min(num_total_draws / tau_hat,
                  num_total_draws * std::log10(num_total_draws));
}

inline double compute_effective_sample_size(std::vector<const double*> draws,
                                            size_t size) {
  int num_chains = draws.size();
  std::vector<size_t> sizes(num_chains, size);
  return compute_effective_sample_size(draws, sizes);
}


inline double compute_split_effective_sample_size(std::vector<const double*> draws, std::vector<size_t> sizes) {
  int num_chains = sizes.size();
  size_t num_draws = sizes[0];
  for (int chain = 1; chain < num_chains; ++chain) {
    num_draws = std::min(num_draws, sizes[chain]);
  }
    long  double half = num_draws / 2.0;
     // double temp= std::ceil(half);
    double temp=  std::ceil(half);
     
     int half_draws = temp;
       
     std::vector<const double*> split_draws(2 * num_chains);
     for (int n = 0; n < num_chains; ++n) {
       split_draws[2 * n] = &draws[n][0];
       split_draws[2 * n + 1] = &draws[n][half_draws];
     }

  std::vector<size_t> half_sizes(2 * num_chains, std::floor(half));

  return compute_effective_sample_size(split_draws, half_sizes);
}

inline double compute_split_effective_sample_size(
    std::vector<const double*> draws, size_t size) {
  int num_chains = draws.size();
  std::vector<size_t> sizes(num_chains, size);
  return compute_split_effective_sample_size(draws, sizes);
}

//inline std::vector<const double*> split_chains(std::vector<const double*> draws, size_t size) ;
//inline std::vector<const double*> split_chains(const std::vector<const double*>& draws, const std::vector<size_t>& sizes) ;

inline std::vector<const   double*> split_chains(const std::vector<const double*>& draws, const std::vector<size_t>& sizes) {
  int num_chains = sizes.size();
  size_t num_draws = sizes[0];
  for (int chain = 1; chain < num_chains; ++chain) {
    num_draws = std::min(num_draws, sizes[chain]);
  }

 long  double half = num_draws / 2.0;
    
  double temp= std::ceil(half);
  int half_draws = temp;
    
  std::vector<const   double*> split_draws(2 * num_chains);
  for (int n = 0; n < num_chains; ++n) {
    split_draws[2 * n] = &draws[n][0];
    split_draws[2 * n + 1] = &draws[n][half_draws];
  }
  return split_draws;
}
inline std::vector<const   double*> split_chains(std::vector<const double*> draws,
                                               size_t size) {
  int num_chains = draws.size();
  std::vector<size_t> sizes(num_chains, size);
  return split_chains(draws, sizes);
}


inline double variance_eigen(const Eigen::VectorXd& x) {
  double m = (x.array() / x.size()).sum();
  return ((x.array() - m) / std::sqrt((x.size() - 1.0))).square().sum();
}

inline long  double mean_eigen(const Eigen::VectorXd& x) {
  return (x.array() / x.size()).sum();
}

inline long  double sd_eigen(const Eigen::VectorXd& x) { return std::sqrt(variance_eigen(x)); }
//
//double split_potential_scale_reduction(
//    const Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, 1>& samples)  {
//  int chains = samples.size();
//  int n_samples = samples(0).size();
//  for (int chain = 1; chain < chains; chain++) {
//    n_samples = std::min(n_samples, static_cast<int>(samples(chain).size()));
//  }
//  if (n_samples % 2 == 1)
//    n_samples--;
//  int n = n_samples / 2;
//
//  Eigen::VectorXd split_chain_mean(2 * chains);
//  Eigen::VectorXd split_chain_var(2 * chains);
//
//  for (int chain = 0; chain < chains; chain++) {
//    split_chain_mean(2 * chain) = mean_eigen(samples(chain).topRows(n));
//    split_chain_mean(2 * chain + 1) = mean_eigen(samples(chain).bottomRows(n));
//
//    split_chain_var(2 * chain) = variance_eigen(samples(chain).topRows(n));
//    split_chain_var(2 * chain + 1) = variance_eigen(samples(chain).bottomRows(n));
//  }
//
//  double var_between = n * variance_eigen(split_chain_mean);
//  double var_within = mean_eigen(split_chain_var);
//
//  // rewrote [(n-1)*W/n + B/n]/W as (n-1+ B/W)/n
//  return sqrt((var_between / var_within + n - 1) / n);
//}

inline double compute_potential_scale_reduction(
    std::vector<const double*> draws, std::vector<size_t> sizes) {
  int num_chains = sizes.size();
  size_t num_draws = sizes[0];
  for (int chain = 1; chain < num_chains; ++chain) {
    num_draws = std::min(num_draws, sizes[chain]);
  }

  // check if chains are constant; all equal to first draw's value
  bool are_all_const = false;
  Eigen::VectorXd init_draw = Eigen::VectorXd::Zero(num_chains);

  for (int chain = 0; chain < num_chains; chain++) {
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1>> draw(
        draws[chain], sizes[chain]);

    for (int n = 0; n < num_draws; n++) {
      if (!std::isfinite(draw(n))) {
        return std::numeric_limits<double>::quiet_NaN();
      }
    }

    init_draw(chain) = draw(0);

    if (draw.isApproxToConstant(draw(0))) {
      are_all_const |= true;
    }
  }

  if (are_all_const) {
    // If all chains are constant then return NaN
    // if they all equal the same constant value
    if (init_draw.isApproxToConstant(init_draw(0))) {
      return std::numeric_limits<double>::quiet_NaN();
    }
  }

  using boost::accumulators::accumulator_set;
  using boost::accumulators::stats;
  using boost::accumulators::tag::mean;
  using boost::accumulators::tag::variance;

  Eigen::VectorXd chain_mean(num_chains);
  accumulator_set<double, stats<variance>> acc_chain_mean;
  Eigen::VectorXd chain_var(num_chains);
  double unbiased_var_scale = num_draws / (num_draws - 1.0);

  for (int chain = 0; chain < num_chains; ++chain) {
    accumulator_set<double, stats<mean, variance>> acc_draw;
    for (int n = 0; n < num_draws; ++n) {
      acc_draw(draws[chain][n]);
    }

    chain_mean(chain) = boost::accumulators::mean(acc_draw);
    acc_chain_mean(chain_mean(chain));
    chain_var(chain)
        = boost::accumulators::variance(acc_draw) * unbiased_var_scale;
  }

  double var_between = num_draws * boost::accumulators::variance(acc_chain_mean)
                       * num_chains / (num_chains - 1);
  double var_within = chain_var.mean();

  // rewrote [(n-1)*W/n + B/n]/W as (n-1+ B/W)/n
  return sqrt((var_between / var_within + num_draws - 1) / num_draws);
}
inline double compute_potential_scale_reduction(
    std::vector<const double*> draws, size_t size) {
  int num_chains = draws.size();
  std::vector<size_t> sizes(num_chains, size);
  return compute_potential_scale_reduction(draws, sizes);
}
/**
 * Computes the split potential scale reduction (Rhat) for the
 * specified parameter across all kept samples.  When the number of
 * total draws N is odd, the (N+1)/2th draw is ignored.
 *
 * See more details in Stan reference manual section "Potential
 * Scale Reduction". http://mc-stan.org/users/documentation
 *
 * Current implementation assumes draws are stored in contiguous
 * blocks of memory.  Chains are trimmed from the back to match the
 * length of the shortest chain.
 *
 * @param draws stores pointers to arrays of chains
 * @param sizes stores sizes of chains
 * @return potential scale reduction for the specified parameter
 */
inline double compute_split_potential_scale_reduction(
    std::vector<const double*> draws, std::vector<size_t> sizes) {
  int num_chains = sizes.size();
  size_t num_draws = sizes[0];
  for (int chain = 1; chain < num_chains; ++chain) {
    num_draws = std::min(num_draws, sizes[chain]);
  }

  std::vector<const double*> split_draws = split_chains(draws, sizes);

  double half = num_draws / 2.0;
  std::vector<size_t> half_sizes(2 * num_chains, std::floor(half));

  return compute_potential_scale_reduction(split_draws, half_sizes);
}
/**
* Computes the split potential scale reduction (Rhat) for the
* specified parameter across all kept samples.  When the number of
* total draws N is odd, the (N+1)/2th draw is ignored.
*
* See more details in Stan reference manual section "Potential
* Scale Reduction". http://mc-stan.org/users/documentation
*
* Current implementation assumes draws are stored in contiguous
* blocks of memory.  Chains are trimmed from the back to match the
* length of the shortest chain.  Argument size will be broadcast to
* same length as draws.
*
* @param draws stores pointers to arrays of chains
* @param sizes stores sizes of chains
* @return potential scale reduction for the specified parameter
*/
inline  double compute_split_potential_scale_reduction(
    std::vector<const double*> draws, size_t size) {
  int num_chains = draws.size();
  std::vector<size_t> sizes(num_chains, size);
  return compute_split_potential_scale_reduction(draws, sizes);
}
}  // namespace mcmc_utils

#endif /* autocovariance_hpp */
