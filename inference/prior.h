//
//  prior.h
//  TumorCellCoal
//
//  Created by Fausto Fabian Crespo Fernandez on 7/8/19.
//  Copyright © 2019 Fausto Fabian Crespo Fernandez. All rights reserved.
//

#ifndef prior_h
#define prior_h
void GenerateTimesFromPriors(int noisy, int numClones,Population **populations, long int seed);
void GenerateTimesFromPriors2(int noisy, int numClones,Population **populations,  long int seed);
void GenerateTimesFromPriors3(int noisy, int numClones,Population **populations,  long int seed);
void GenerateTimesFromPriors4(int noisy, int numClones,double *proportionsVector, Population **populations,  long int seed);
void GenerateTimesFromPriorsOriginal(int noisy, int numClones,double *proportionsVector, Population **populations, long int *seed);

double DirichletDensity(double *proportionsVector, double *concentrationVector, int sizeVector);
void GeneratePopSizesFromPriors(int noisy, int numClones,Population **populations,  long int *seed, double gammaParam, int totalPopSize);
//void GenerateEffectPopSizesFromPriors2(int noisy, int numClones,Population **populations,  long int *seed, double gammaParam, int oldestPopSize, double **proportionsVector, int doGenerateProportionsVector);
void GenerateEffectPopSizesFromPriors2(Chain *chain, int noisy, int numClones,Population **populations,  long int *seed, double gammaParam, int totalPopSize,  int doGenerateProportionsVector);

double generateTimeofOriginfromPrior( double mean, double dispersion, long int *seed);


#endif /* prior_h */
