//
//  resampling.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 15/03/2021.
//

#include "resampling.hpp"
void multinomialResampling(const gsl_rng *random, const std::vector<double> *probs, unsigned int N, unsigned int *indices)
{
    multinomialSampleIndices(random, N, *probs, indices);
}

void stratifiedResampling(const gsl_rng *random, const std::vector<double> *probs, unsigned int N, unsigned int *indices)
{
    // draw N uniform numbers from (n/N, (n+1)/N]
    double num_samples = (double)N;
    double *u = new double[N];
    for (int n = 0; n < N; n++) {
        double a = n/num_samples;
        double b = (n+1)/num_samples;
        u[n] = gsl_ran_flat(random, a, b);
    }
    
    // u's are already sorted -- determine the indices
    double sum = 0.0;
    int idx = 0;
    for (unsigned int n = 0; n < N; n++) {
        while (u[n] > (sum + (*probs)[idx])) {
            sum += (*probs)[idx];
            idx++;
        }
        indices[n] = idx;
    }
    
    delete [] u;
}

void systematicResampling(const gsl_rng *random, const std::vector<double> *probs, unsigned int N, unsigned int *indices)
{
    // draw a single uniform from (0, 1/N]
    double u = gsl_ran_flat(random, 0, 1./N);
    double u_i = u;
    double inc = 1./N;
    double sum = 0.0;
    int idx = 0;
    for (unsigned int n = 0; n < N; n++) {
        // find the index
        while (u_i > (sum + (*probs)[idx])) {
            sum += (*probs)[idx];
            idx++;
        }
        indices[n] = idx;
        u_i += inc;
    }
}

void multinomialResamplingSortedUniform(const gsl_rng *random, unsigned int N, double *uvec)
{
    uniform(random, N, uvec); // N
    std::sort(uvec, uvec + N); // N log N
}

void stratifiedResamplingSortedUniform(const gsl_rng *random, unsigned int N, double *u)
{
    // draw N uniform numbers from (n/N, (n+1)/N]
    double num_samples = (double)N;
    for (int n = 0; n < N; n++) {
        double a = n/num_samples;
        double b = (n+1)/num_samples;
        u[n] = gsl_ran_flat(random, a, b);
    }
}

void systematicResamplingSortedUniform(const gsl_rng *random, unsigned int N, double *u)
{
    double one_over_N = 1./N;
    double initial_unif = gsl_ran_flat(random, 0, one_over_N);
    double double_N = (double)N;
    
    for (int n = 0; n < N; n++) {
        u[n] = (n - 1)/double_N + initial_unif;
    }
}
void uniform(const gsl_rng *random, unsigned int N, double *ret)
{
    for (int n = 0; n < N; n++)
    {
        ret[n] = gsl_rng_uniform(random);
    }
}
void multinomialSampleIndices(const gsl_rng *random, unsigned int N, const std::vector<double> &normalized_probs, unsigned int *indices)
{
    // draw N uniform values
    double *uvec = new double[N];
    uniform(random, N, uvec); // N
    std::sort(uvec, uvec + N); // N log N
    double sum = 0.0;
    int idx = 0;
    for (int n = 0; n < N; n++)
    {
        while (uvec[n] > sum + normalized_probs[idx]) {
            sum += normalized_probs[idx];
            idx++;
        }
        if (idx > N) {
            sum = sum + 0;
        }
        indices[n] = idx;
    }
    delete [] uvec;
}
