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

double   RandomUniform (long int *seed);
double RandomExponential (double mean, long int *seed);
double    RandomGamma (double shape, long int *seed);
double RandomGamma1 (double s, long int *seed);
double RandomGamma2 (double s, long int *seed);
int RandomNegativeBinomial (double mean, double dispersion, long int *seed);
double RandomNormal(double mean, double stddev, long int *seed);
double RandomLogNormal( double mean, double dispersion, long int *seed);
void  RandomDirichlet2 (double *vectorGamma, int vectorSize, double *outputVector, long int *seed);
int RandomPoisson (double lambda, long int *seed);

int RandomBinomial (double prob, int numTrials, long int *seed);

int RandomUniformTo (int max, long int *seed);

double RandomLogUniform( double from, double to);

void  RandomDirichlet (double s, int vectorSize, vector<double> &outputVector, long int *seed);
double LogUniformDensity(double value, double from, double to);

int ChooseUniformState (double *prob, long int *seed);

double randomUniformFromGsl();
void   randomDirichletFromGsl(int vectorSize, double alpha[], double *theta);

void  randomDirichletFromVector (vector<double> alpha, vector<double> &outputVector);

#endif /* random_h */
