//
//  random.h
//  TumorCellCoal
//
//  Created by Fausto Fabian Crespo Fernandez on 7/8/19.
//  Copyright © 2019 Fausto Fabian Crespo Fernandez. All rights reserved.
//


#ifndef random_h
#define random_h

#include <vector>



using namespace std;

long double    RandomUniform (long int *seed);
long double  RandomExponential (long double  mean, long int *seed);
long double     RandomGamma (long double  shape, long int *seed);
long double  RandomGamma1 (long double  s, long int *seed);
long double  RandomGamma2 (long double  s, long int *seed);
int RandomNegativeBinomial (long double  mean, long double  dispersion, long int *seed);
long double  RandomNormal(long double  mean, long double  stddev, long int *seed);
long double  RandomLogNormal( long double  mean, long double  dispersion, long int *seed);
void  RandomDirichlet2 (long double  *vectorGamma, int vectorSize, long double  *outputVector, long int *seed);
int RandomPoisson (long double  lambda, long int *seed);

int RandomBinomial (long double  prob, int numTrials, long int *seed);

int RandomUniformTo (int max, long int *seed);

long double  RandomLogUniform( long double  from, long double  to);

void  RandomDirichlet (long double  s, int vectorSize, vector<double> &outputVector, long int *seed);
long double  LogUniformDensity(long double  value, long double  from, long double  to);

int ChooseUniformState ( long double  *prob, long int *seed);

long double  randomUniformFromGsl();
void   randomDirichletFromGsl(int vectorSize,  double  alpha[], long double  *theta);

void  randomDirichletFromVector (vector<long double> alpha, vector<long double> &outputVector);
long double  randomGammaBoost( long double  shape, long double  scale );
long double  randomGammaBoost2( long double  mean , long double  variance );
long double  logMultinomialProbability(const  size_t size, const  double  *p, const unsigned int *n );
long double randomNormalBoost(long double mean, long double sigma);
long double  RandomExponentialStartingFrom (long double  lambda, long double from);
#endif /* random_h */
