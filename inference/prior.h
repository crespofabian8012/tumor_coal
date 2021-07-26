/*################################################################################
 ##
 ##   Copyright (C) 2018-2020 Fausto Fabian Crespo Fernandez
 ##
 ##   This file is part of the tumor_coal C++ library.
 ##
 ##   Licensed under the Apache License, Version 2.0 (the "License");
 ##   you may not use this file except in compliance with the License.
 ##   You may obtain a copy of the License at
 ##
 ##       http://www.apache.org/licenses/LICENSE-2.0
 ##
 ##   Unless required by applicable law or agreed to in writing, software
 ##   distributed under the License is distributed on an "AS IS" BASIS,
 ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 ##   See the License for the specific language governing permissions and
 ##   limitations under the License.
 ##
 ################################################################################*/

/*
 * prior functions
 */
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
