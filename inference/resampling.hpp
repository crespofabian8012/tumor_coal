//
//  resampling.hpp
//  run
//
//  Created by  Seong-Hwan Jun
//

#ifndef resampling_hpp
#define resampling_hpp




#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>



void multinomialResampling(const gsl_rng *random, const std::vector<double> *probs, unsigned int N, unsigned int *indices);
void stratifiedResampling(const gsl_rng *random, const std::vector<double> *probs, unsigned int N, unsigned int *indices);
void systematicResampling(const gsl_rng *random, const std::vector<double> *probs, unsigned int N, unsigned int *indices);

void multinomialResamplingSortedUniform(const gsl_rng *random, unsigned int N, double *u);
void stratifiedResamplingSortedUniform(const gsl_rng *random, unsigned int N, double *u);
void systematicResamplingSortedUniform(const gsl_rng *random, unsigned int N, double *u);

void uniform(const gsl_rng *random, unsigned int N, double *ret);
void multinomialSampleIndices(const gsl_rng *random, unsigned int N, const std::vector<double> &normalized_probs, unsigned int *indices);
#endif /* resampling_hpp */
