//
//  autocorrelation.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 22/02/2021.
//

#ifndef autocorrelation_hpp
#define autocorrelation_hpp

#include <stdio.h>

#endif /* autocorrelation_hpp */

#include <unsupported/Eigen/FFT>
#include <vector>

namespace mcmc_utils {
inline size_t fft_next_good_size(size_t N) {
  if (N <= 2) {
    return 2;
  }
  while (true) {
    size_t m = N;
    while ((m % 2) == 0) {
      m /= 2;
    }
    while ((m % 3) == 0) {
      m /= 3;
    }
    while ((m % 5) == 0) {
      m /= 5;
    }
    if (m <= 1) {
      return N;
    }
    N++;
  }
}

template <typename T>
T generic_mean(const std::vector<T> &vec){
    
    size_t sz = vec.size();
    
    if (sz == 0)
        return 0.0;
    if (sz == 1){
        long double elem = vec[0];
        return elem;
    }

 long double accum=0.0;
 for (auto it= vec.begin(); it != vec.end(); ++it)
     accum += *it;

   long double mean = accum / sz;
    
    return mean;
}
template <typename T>
void autocorrelationVector(const std::vector<T>& y, std::vector<T>& ac,
                     Eigen::FFT<T>& fft) {
  using std::complex;
  using std::vector;

  size_t N = y.size();
  size_t M = fft_next_good_size(N);
  size_t Mt2 = 2 * M;

  vector<complex<T> > freqvec;

  // centered_signal = y-mean(y) followed by N zeroes
  vector<T> centered_signal(y);
  centered_signal.insert(centered_signal.end(), Mt2 - N, 0.0);
  T mean = mean(y);
  for (size_t i = 0; i < N; i++) {
    centered_signal[i] -= mean;
  }

  fft.fwd(freqvec, centered_signal);
  for (size_t i = 0; i < Mt2; ++i) {
    freqvec[i] = complex<T>(norm(freqvec[i]), 0.0);
  }

  fft.inv(ac, freqvec);
  ac.resize(N);

  for (size_t i = 0; i < N; ++i) {
    ac[i] /= (N - i);
  }
  T var = ac[0];
  for (size_t i = 0; i < N; ++i) {
    ac[i] /= var;
  }
}

template <typename T>
void autocorrelationVector(const std::vector<T>& y, std::vector<T>& ac) {
    
  using std::complex;
  using std::vector;
  Eigen::FFT<T> fft;
    
  size_t N = y.size();
  size_t M = fft_next_good_size(N);
  size_t Mt2 = 2 * M;

  vector<complex<T> > freqvec;

  // centered_signal = y-mean(y) followed by N zeroes
  vector<T> centered_signal(y);
  centered_signal.insert(centered_signal.end(), Mt2 - N, 0.0);
  T mean = generic_mean(y);
  for (size_t i = 0; i < N; i++) {
    centered_signal[i] -= mean;
  }

  fft.fwd(freqvec, centered_signal);
  for (size_t i = 0; i < Mt2; ++i) {
    freqvec[i] = complex<T>(norm(freqvec[i]), 0.0);
  }

  fft.inv(ac, freqvec);
  ac.resize(N);

  for (size_t i = 0; i < N; ++i) {
    ac[i] /= (N - i);
  }
  T var = ac[0];
  for (size_t i = 0; i < N; ++i) {
    ac[i] /= var;
  }
}

template <typename T>
void autocovarianceVector(const std::vector<T>& y, std::vector<T>& acov,
                    Eigen::FFT<T>& fft) {
  autocorrelationVector(y, acov, fft);

  T var = variance(y) * (y.size() - 1) / y.size();
  for (size_t i = 0; i < y.size(); i++) {
    acov[i] *= var;
  }
}

template <typename T>
void autocovarianceVector(const std::vector<T>& y, std::vector<T>& acov) {
  
  Eigen::FFT<T> fft;
  autocorrelationVector(y, acov, fft);

  T var = variance(y) * (y.size() - 1) / y.size();
  for (size_t i = 0; i < y.size(); i++) {
    acov[i] *= var;
  }
}
}
