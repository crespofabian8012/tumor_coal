//
//  core_likelihood.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 11/08/2021.
//

#ifndef core_likelihood_hpp
#define core_likelihood_hpp

#define PLL_SCALE_FACTOR 115792089237316195423570985008687907853269984665640564039457584007913129639936.0  /*  2**256 (exactly)  */
#define PLL_SCALE_THRESHOLD (1.0/PLL_SCALE_FACTOR)
#define PLL_SCALE_FACTOR_SQRT 340282366920938463463374607431768211456.0 /* 2**128 */

double pll_core_root_loglikelihood2(unsigned int states,
                                    unsigned int sites,
                                    unsigned int rate_cats,
                                    const double * clv,
                                    const unsigned int * scaler,
                                    double * const * frequencies,
                                    const double * rate_weights,
                                    const unsigned int * pattern_weights,
                                    const unsigned int * freqs_indices,
                                    double * persite_lnl,
                                    unsigned int attrib);
void pll_core_root_likelihood_vector2(unsigned int states,
                                      unsigned int sites,
                                      unsigned int rate_cats,
                                      const double * clv,
                                      const unsigned int * scaler,
                                      double * const * frequencies,
                                      const double * rate_weights,
                                      const unsigned int * pattern_weights,
                                      const unsigned int * freqs_indices,
                                      double * persite_lh,
                                      unsigned int attrib);
double pll_core_root_loglikelihood_second_form(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               const double * clv,
                                               const unsigned int * scaler,
                                               double * const * frequencies,
                                               const double * rate_weights,
                                               const unsigned int * pattern_weights,
                                               const unsigned int * freqs_indices,
                                               double * persite_lnl,
                                               unsigned int attrib);
int pll_core_update_pmatrix2(double ** pmatrix,
                             unsigned int states,
                             unsigned int rate_cats,
                             const double * rates,
                             const double * branch_lengths,
                             const unsigned int * matrix_indices,
                             const unsigned int * param_indices,
                             double * const * eigenvals,
                             double * const * eigenvecs,
                             double * const * inv_eigenvecs,
                             unsigned int count,
                             unsigned int attrib);
void pll_core_update_partial_ii_norm(unsigned int states,
                                 unsigned int sites,
                                 unsigned int rate_cats,
                                 double * parent_clv,
                                 unsigned int * parent_scaler,
                                 const double * left_clv,
                                 const double * right_clv,
                                 const double * left_matrix,
                                 const double * right_matrix,
                                 const unsigned int * left_scaler,
                                 const unsigned int * right_scaler,
                                 unsigned int attrib,
                                 bool normalize);
int pll_core_update_pmatrix_16x16_jc69(double ** pmatrix,
                                       unsigned int states,
                                       unsigned int rate_cats,
                                       const double * rates,
                                       const double * branch_lengths,
                                       const unsigned int * matrix_indices,
                                       const unsigned int * param_indices,
                                       unsigned int count,
                                       unsigned int attrib);
int pll_core_update_pmatrix_16x16_jc69_matrix_second_form(double ** pmatrix,
                                                          unsigned int states,
                                                          unsigned int rate_cats,
                                                          const double * rates,
                                                          const double * branch_lengths,
                                                          const unsigned int * matrix_indices,
                                                          const unsigned int * param_indices,
                                                          unsigned int count,
                                                          unsigned int attrib);
void fill_parent_scaler2(unsigned int scaler_size,
                         unsigned int * parent_scaler,
                         const unsigned int * left_scaler,
                         const unsigned int * right_scaler);
void pll_core_update_parent_clv(unsigned int states,
                                unsigned int sites,
                                unsigned int rate_cats,
                                double * parent_clv,
                                unsigned int * parent_scaler,
                                const double * left_clv,
                                const double * right_clv,
                                const double * branch_lengths,
                                const double * left_matrix,
                                const double * right_matrix,
                                const unsigned int * left_scaler,
                                const unsigned int * right_scaler,
                                unsigned int attrib,
                                bool normalize);
#endif /* core_likelihood_hpp */
