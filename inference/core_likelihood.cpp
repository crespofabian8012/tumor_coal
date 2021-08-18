//
//  core_likelihood.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 11/08/2021.
//
#include "core_likelihood.hpp"
#include <iostream>
extern "C"
    {
#include "libpll/pll.h"
#include "libpll/pll_optimize.h"
#include "libpll/pllmod_common.h"
#include <libpll/pll_msa.h> //for pllmod_msa_empirical_frequencies
#include <libpll/pll_tree.h>
#include <libpll/pllmod_util.h>
#include "libpll/pll_optimize.h"
#include <libpll/pllmod_common.h>
#include <libpll/pllmod_algorithm.h>
    }
using namespace std;
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
                                    unsigned int attrib)
{
    unsigned int i,j,k;
    double logl = 0;
    const double * freqs = NULL;
    
    double term, term_r;
    double site_lk;
    
    unsigned int states_padded = states;
    
#ifdef HAVE_SSE3
    if (attrib & PLL_ATTRIB_ARCH_SSE)
    {
        if (states == 4)
        {
            return pll_core_root_loglikelihood_4x4_sse(sites,
                                                       rate_cats,
                                                       clv,
                                                       scaler,
                                                       frequencies,
                                                       rate_weights,
                                                       pattern_weights,
                                                       freqs_indices,
                                                       persite_lnl);
        }
        else
        {
            return pll_core_root_loglikelihood_sse(states,
                                                   sites,
                                                   rate_cats,
                                                   clv,
                                                   scaler,
                                                   frequencies,
                                                   rate_weights,
                                                   pattern_weights,
                                                   freqs_indices,
                                                   persite_lnl);
        }
        /* this line is never called, but should we disable the else case above,
         then states_padded must be set to this value */
        states_padded = (states+1) & 0xFFFFFFFE;
    }
#endif
#ifdef HAVE_AVX
    if (attrib & PLL_ATTRIB_ARCH_AVX)
    {
        if (states == 4)
        {
            return pll_core_root_loglikelihood_4x4_avx(sites,
                                                       rate_cats,
                                                       clv,
                                                       scaler,
                                                       frequencies,
                                                       rate_weights,
                                                       pattern_weights,
                                                       freqs_indices,
                                                       persite_lnl);
        }
        else
        {
            return pll_core_root_loglikelihood_avx(states,
                                                   sites,
                                                   rate_cats,
                                                   clv,
                                                   scaler,
                                                   frequencies,
                                                   rate_weights,
                                                   pattern_weights,
                                                   freqs_indices,
                                                   persite_lnl);
        }
        /* this line is never called, but should we disable the else case above,
         then states_padded must be set to this value */
        states_padded = (states+3) & 0xFFFFFFFC;
    }
#endif
#ifdef HAVE_AVX2
    if (attrib & PLL_ATTRIB_ARCH_AVX2)
    {
        if (states == 4)
        {
            return pll_core_root_loglikelihood_4x4_avx(sites,
                                                       rate_cats,
                                                       clv,
                                                       scaler,
                                                       frequencies,
                                                       rate_weights,
                                                       pattern_weights,
                                                       freqs_indices,
                                                       persite_lnl);
        }
        else
        {
            return pll_core_root_loglikelihood_avx2(states,
                                                    sites,
                                                    rate_cats,
                                                    clv,
                                                    scaler,
                                                    frequencies,
                                                    rate_weights,
                                                    pattern_weights,
                                                    freqs_indices,
                                                    persite_lnl);
        }
        /* this line is never called, but should we disable the else case above,
         then states_padded must be set to this value */
        states_padded = (states+3) & 0xFFFFFFFC;
    }
#endif
    
    /* iterate through sites */
    for (i = 0; i < sites; ++i)
    {
        term = 0;
        for (j = 0; j < rate_cats; ++j)
        {
            freqs = frequencies[freqs_indices[j]];
            term_r = 0;
            for (k = 0; k < states; ++k)
            {
                term_r += clv[k] * freqs[k];
            }
            
            term += term_r * rate_weights[j];
            
            clv += states_padded;
        }
        
        site_lk = term;
        
        assert(site_lk>0);
        /* compute site log-likelihood and scale if necessary */
        site_lk = log(site_lk);
        if (scaler && scaler[i]){
            site_lk += scaler[i] * log(PLL_SCALE_THRESHOLD);
           // std::cout<< "scaling" << std::endl;
        }
        
        site_lk *= pattern_weights[i];
        
        /* store per-site log-likelihood */
        if (persite_lnl)
            persite_lnl[i] = site_lk;
        
        logl += site_lk;
    }
    return logl;
}

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
                                      unsigned int attrib)
{
    unsigned int i,j,k;
    //double logl = 0;
    const double * freqs = NULL;
    
    double term, term_r;
    //double site_lk;
    
    unsigned int states_padded = states;
    
#ifdef HAVE_SSE3
    if (attrib & PLL_ATTRIB_ARCH_SSE)
    {
        if (states == 4)
        {
            pll_core_root_likelihood_vec_4x4_sse(sites,
                                                 rate_cats,
                                                 clv,
                                                 scaler,
                                                 frequencies,
                                                 rate_weights,
                                                 pattern_weights,
                                                 freqs_indices,
                                                 persite_lh);
        }
        else
        {
            pll_core_root_likelihood_vec_sse(states,
                                             sites,
                                             rate_cats,
                                             clv,
                                             scaler,
                                             frequencies,
                                             rate_weights,
                                             pattern_weights,
                                             freqs_indices,
                                             persite_lh);
        }
        return;
        /* this line is never called, but should we disable the else case above,
         then states_padded must be set to this value */
        states_padded = (states+1) & 0xFFFFFFFE;
    }
#endif
#ifdef HAVE_AVX
    if (attrib & PLL_ATTRIB_ARCH_AVX)
    {
        if (states == 4)
        {
            pll_core_root_likelihood_vec_4x4_avx(sites,
                                                 rate_cats,
                                                 clv,
                                                 scaler,
                                                 frequencies,
                                                 rate_weights,
                                                 pattern_weights,
                                                 freqs_indices,
                                                 persite_lh);
        }
        else
        {
            pll_core_root_likelihood_vec_avx(states,
                                             sites,
                                             rate_cats,
                                             clv,
                                             scaler,
                                             frequencies,
                                             rate_weights,
                                             pattern_weights,
                                             freqs_indices,
                                             persite_lh);
        }
        return;
        /* this line is never called, but should we disable the else case above,
         then states_padded must be set to this value */
        states_padded = (states+3) & 0xFFFFFFFC;
    }
#endif
#ifdef HAVE_AVX2
    if (attrib & PLL_ATTRIB_ARCH_AVX2)
    {
        if (states == 4)
        {
            pll_core_root_likelihood_vec_4x4_avx(sites,
                                                 rate_cats,
                                                 clv,
                                                 scaler,
                                                 frequencies,
                                                 rate_weights,
                                                 pattern_weights,
                                                 freqs_indices,
                                                 persite_lh);
        }
        else
        {
            pll_core_root_likelihood_vec_avx2(states,
                                              sites,
                                              rate_cats,
                                              clv,
                                              scaler,
                                              frequencies,
                                              rate_weights,
                                              pattern_weights,
                                              freqs_indices,
                                              persite_lh);
        }
        return;
        /* this line is never called, but should we disable the else case above,
         then states_padded must be set to this value */
        states_padded = (states+3) & 0xFFFFFFFC;
    }
#endif
    
    /* iterate through sites */
    for (i = 0; i < sites; ++i)
    {
        term = 0;
        for (j = 0; j < rate_cats; ++j)
        {
            freqs = frequencies[freqs_indices[j]];
            term_r = 0;
            for (k = 0; k < states; ++k)
            {
                term_r += clv[k] * freqs[k];
            }
            
            term += term_r * rate_weights[j];
            
            clv += states_padded;
        }
        
        persite_lh[i] = term;
#if 0
        site_lk = term;
        
        /* compute site log-likelihood and scale if necessary */
        site_lk = log(site_lk);
        if (scaler && scaler[i])
            site_lk += scaler[i] * log(PLL_SCALE_THRESHOLD);
        
        site_lk *= pattern_weights[i];
        
        /* store per-site log-likelihood */
        if (persite_lnl)
            persite_lnl[i] = site_lk;
        
        logl += site_lk;
#endif
    }
}
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
                             unsigned int attrib)
{
    unsigned int i,n,j,k,m;
    unsigned int states_padded = states;
    double * expd;
    double * temp;
    
    double * evecs;
    double * inv_evecs;
    double * evals;
    double * pmat;
    
    
    expd = (double *)malloc(states * sizeof(double));
    temp = (double *)malloc(states*states*sizeof(double));
    
    for (i = 0; i < count; ++i)
    {
        assert(branch_lengths[i] >= 0);
        
        /* compute effective pmatrix location */
        for (n = 0; n < rate_cats; ++n)
        {
            pmat = pmatrix[matrix_indices[i]] + n*states*states_padded;
            
            evecs = eigenvecs[param_indices[n]];
            inv_evecs = inv_eigenvecs[param_indices[n]];
            evals = eigenvals[param_indices[n]];
            
            /* if branch length is zero then set the p-matrix to identity matrix */
            if (!branch_lengths[i])
            {
                for (j = 0; j < states; ++j)
                    for (k = 0; k < states; ++k)
                        pmat[j*states_padded + k] = (j == k) ? 1 : 0;
            }
            else
            {
                /* NOTE: in order to deal with numerical issues in cases when Qt -> 0, we
                 * use a trick suggested by Ben Redelings and explained here:
                 * https://github.com/xflouris/libpll/issues/129#issuecomment-304004005
                 * In short, we use expm1() to compute (exp(Qt) - I), and then correct
                 * for this by adding an identity matrix I in the very end */
                
                /* exponentiate eigenvalues */
                for (j = 0; j < states; ++j)
                    expd[j] = expm1(evals[j] * rates[n] * branch_lengths[i]);
                
                for (j = 0; j < states; ++j)
                    for (k = 0; k < states; ++k)
                        temp[j*states+k] = inv_evecs[j*states_padded+k] * expd[k];
                
                for (j = 0; j < states; ++j)
                {
                    for (k = 0; k < states; ++k)
                    {
                        pmat[j*states_padded+k] = (j==k) ? 1.0 : 0;
                        for (m = 0; m < states; ++m)
                        {
                            pmat[j*states_padded+k] +=
                            temp[j*states+m] * evecs[m*states_padded+k];
                        }
                    }
                }
            }
            
            for (j = 0; j < states; ++j)
                for (k = 0; k < states; ++k)
                    assert(pmat[j*states_padded+k] >= 0);
            
        }
    }
    
    free(expd);
    free(temp);
    return PLL_SUCCESS;
}
int pll_core_update_pmatrix_16x16_jc69(double ** pmatrix,
                                       unsigned int states,
                                       unsigned int rate_cats,
                                       const double * rates,
                                       const double * branch_lengths,
                                       const unsigned int * matrix_indices,
                                       const unsigned int * param_indices,
                                       unsigned int count,
                                       unsigned int attrib)
{
    unsigned int i,n;
    unsigned int states_padded = states;
    
    double * pmat;
    
    
    for (i = 0; i < count; ++i)
    {
        assert(branch_lengths[i] >= 0);
        
        /* compute effective pmatrix location */
        for (n = 0; n < rate_cats; ++n)
        {
            pmat = pmatrix[matrix_indices[i]] + n*states*states_padded;
            
            double t = branch_lengths[i] * rates[n];
            
            if (t < 1e-100)
            {
                for(unsigned int i=0; i< states*states;i++)
                {
                    if (i %states ==0)
                        pmat[i] = 1;
                    else
                        pmat[i] = 0;
                    
                }
                
            }
            else
            {
                //        #if 0
                //        double a =  (1 + 3*exp(-4*t/3) ) / 4;
                //        double b = (1 - a) / 3;
                //        #endif
                //
                //        double exptm1 = expm1(-4*t/3);
                //        double a = 1 + 3/4. * exptm1;
                //        double b = -exptm1/4;
                
#if 0
                double a =  (1 + 15*exp(-5*t/3) ) / 16;
                double b = (1 - a) / 15;
#endif
                
                double exptm1 = expm1(-5*t/3);
                double a = 1 + 15/16. * exptm1;
                double b = -exptm1/16;
                
                for(unsigned int i=0; i< states*states;i++)
                {
                    if (i %states ==0)
                        pmat[i] = a;
                    else
                        pmat[i] = b;
                    
                }
            }
            
        }
    }
    
    return PLL_SUCCESS;
}

void pll_core_update_partial_ii2(unsigned int states,
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
                                 unsigned int attrib)
{
    unsigned int i,j,k,n;
    
    unsigned int scale_mode;  /* 0 = none, 1 = per-site, 2 = per-rate */
    unsigned int site_scale;
    unsigned int init_mask;
    
    const double * lmat;
    const double * rmat;
    
    unsigned int span = states * rate_cats;
    
#ifdef HAVE_SSE3
    if (attrib & PLL_ATTRIB_ARCH_SSE)
    {
        pll_core_update_partial_ii_sse(states,
                                       sites,
                                       rate_cats,
                                       parent_clv,
                                       parent_scaler,
                                       left_clv,
                                       right_clv,
                                       left_matrix,
                                       right_matrix,
                                       left_scaler,
                                       right_scaler,
                                       attrib);
        return;
    }
#endif
#ifdef HAVE_AVX
    if (attrib & PLL_ATTRIB_ARCH_AVX)
    {
        pll_core_update_partial_ii_avx(states,
                                       sites,
                                       rate_cats,
                                       parent_clv,
                                       parent_scaler,
                                       left_clv,
                                       right_clv,
                                       left_matrix,
                                       right_matrix,
                                       left_scaler,
                                       right_scaler,
                                       attrib);
        return;
    }
#endif
#ifdef HAVE_AVX2
    if (attrib & PLL_ATTRIB_ARCH_AVX2)
    {
        pll_core_update_partial_ii_avx2(states,
                                        sites,
                                        rate_cats,
                                        parent_clv,
                                        parent_scaler,
                                        left_clv,
                                        right_clv,
                                        left_matrix,
                                        right_matrix,
                                        left_scaler,
                                        right_scaler,
                                        attrib);
        return;
    }
#endif
    
    /* init scaling-related stuff */
    if (parent_scaler)
    {
        /* determine the scaling mode and init the vars accordingly */
        scale_mode = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 2 : 1;
        init_mask = (scale_mode == 1) ? 1 : 0;
        const size_t scaler_size = (scale_mode == 2) ? sites * rate_cats : sites;
        
        assert(scaler_size == sites );
        /* add up the scale vectors of the two children if available */
        fill_parent_scaler2(scaler_size, parent_scaler, left_scaler, right_scaler);
    }
    else
    {
        /* scaling disabled / not required */
        scale_mode = init_mask = 0;
    }
    
    /* compute CLV */
    for (n = 0; n < sites; ++n)
    {
        lmat = left_matrix;
        rmat = right_matrix;
        site_scale = init_mask;
        
        for (k = 0; k < rate_cats; ++k)
        {
            unsigned int rate_scale = 1;
            for (i = 0; i < states; ++i)
            {
                double terma = 0;
                double termb = 0;
                for (j = 0; j < states; ++j)
                {
                    terma += lmat[j] * left_clv[j];
                    termb += rmat[j] * right_clv[j];
                }
                parent_clv[i] = terma*termb;
                
                if (parent_clv[i]==0){
                    printf("Ooops");
                    for (unsigned int k =0; k < states; k++){
                        
                        std::cout<< "lmat:" << lmat[k] << std::endl;
                        std::cout<< "left_clv:" << left_clv[k] << std::endl;
                        std::cout<< "right_clv:" << right_clv[k] << std::endl;
                        std::cout<< lmat[k] * left_clv[k]<< std::endl;
                        std::cout<< rmat[k] * right_clv[k]<< std::endl;
                    }
                }
                
                
                rate_scale &= (parent_clv[i] < PLL_SCALE_THRESHOLD);
                
                lmat += states;
                rmat += states;
            }
            
            /* check if scaling is needed for the current rate category */
            if (scale_mode == 2)
            {
                /* PER-RATE SCALING: if *all* entries of the *rate* CLV were below
                 * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
                if (rate_scale)
                {
                    for (i = 0; i < states; ++i)
                        parent_clv[i] *= PLL_SCALE_FACTOR;
                    parent_scaler[n*rate_cats + k] += 1;
                }
            }
            else
                site_scale = site_scale && rate_scale;
            
            parent_clv += states;
            left_clv   += states;
            right_clv  += states;
        }
        /* PER-SITE SCALING: if *all* entries of the *site* CLV were below
         * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
        if (site_scale)
        {
            parent_clv -= span;
            for (i = 0; i < span; ++i)
                parent_clv[i] *= PLL_SCALE_FACTOR;
            parent_clv += span;
            parent_scaler[n] += 1;
        }
    }
}
void fill_parent_scaler2(unsigned int scaler_size,
                         unsigned int * parent_scaler,
                         const unsigned int * left_scaler,
                         const unsigned int * right_scaler)
{
    unsigned int i;
    
    if (!left_scaler && !right_scaler)
        memset(parent_scaler, 0, sizeof(unsigned int) * scaler_size);
    else if (left_scaler && right_scaler)
    {
        memcpy(parent_scaler, left_scaler, sizeof(unsigned int) * scaler_size);
        for (i = 0; i < scaler_size; ++i)
            parent_scaler[i] += right_scaler[i];
    }
    else
    {
        if (left_scaler)
            memcpy(parent_scaler, left_scaler, sizeof(unsigned int) * scaler_size);
        else
            memcpy(parent_scaler, right_scaler, sizeof(unsigned int) * scaler_size);
    }
}
