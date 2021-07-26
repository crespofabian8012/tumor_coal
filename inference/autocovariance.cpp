//
//  autocovariance.cpp
//  tumor_coal
//
//  Created by Fausto Fabian Crespo Fernandez on 28/01/2021.
//

#include "autocovariance.hpp"
namespace mcmc_utils {
//   inline std::vector<const double*> split_chains(std::vector<const double*> draws,
//                                                  size_t size) {
//     int num_chains = draws.size();
//     std::vector<size_t> sizes(num_chains, size);
//     return split_chains(draws, sizes);
//   }
//   inline std::vector<const double*> split_chains(const std::vector<const double*>& draws, const std::vector<size_t>& sizes) {
//     int num_chains = sizes.size();
//     size_t num_draws = sizes[0];
//     for (int chain = 1; chain < num_chains; ++chain) {
//       num_draws = std::min(num_draws, sizes[chain]);
//     }
//
//     double half = num_draws / 2.0;
//       
//    double temp= std::ceil(half);
//     int half_draws = temp;
//       
//     std::vector<const double*> split_draws(2 * num_chains);
//     for (int n = 0; n < num_chains; ++n) {
//       split_draws[2 * n] = &draws[n][0];
//       split_draws[2 * n + 1] = &draws[n][half_draws];
//     }
//     return split_draws;
//   }
}
