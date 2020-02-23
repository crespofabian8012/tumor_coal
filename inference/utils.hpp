//
//  util.h
//  TumorCellCoal
//
//  Created by Fausto Fabian Crespo Fernandez on 7/8/19.
//  Copyright © 2019 Fausto Fabian Crespo Fernandez. All rights reserved.
//

#ifndef util_h
#define util_h

#include "data_types.hpp"
#include "population.hpp"
#include "tree_node.hpp"

void ReadParametersFromFastaFile(char *fileName, ProgramOptions &programOptions);
void ReadFastaFile(char *fileName, vector<vector<int> > &ObservedData,  char **ObservedCellNames, ProgramOptions &programOptions);

int CheckMatrixSymmetry(double matrix[4][4]);
int WhichNucChar (char nucleotide);
int ChooseUniformState (double *prob, long int *seed);
char WhichIUPAC (int allele1, int allele2);
char *WhichGenotypeFromIUPAC (int  iupac);
int WhichGenotypeChar (char nucleotide);
char WhichNuc (int nucleotide);
char WhichConsensusBinary (int allele1, int allele2);
char WhichMut (int state);
void  normalizeVector(double *vector, int length);
int compareIntDescending(const void *a, const void *b);
int CompareGenotypes (int a1, int a2, int b1, int b2);

double * expand_uniq_rates(int states, const double * uniq_rates, const int * rate_sym);

static void unscale(double * prob, unsigned int times);

double DistanceBetweenSFS(int* SFS1, int*SFS, int numSNVs,  int numSites);


int computeTajimaD(TreeNode **treeTips,  int numSites, int numCells);

double ComputeESS(double *weights, int numberWeights);
void computeUnfoldedISMSFS(int numSites,SiteStr* allSites,int numSNVs, int* SNVsites, int* SFS, int *numberDifferences);

void Initialize( double (*Eij)[4], double (*Mij)[4], double *freq,  ProgramOptions &programOptions ) ;

void computeGenotypesFreq(double freqs[10], pll_msa_t * msa);

void InitNumberNodes(double *TotalBirthRate, double *TotalDeathRate, int *TotalN,  Population **populations, ProgramOptions &programOptions) ;


double parameterMultiplierMCMCmove (double lengthInterval);

#endif /* util_h */
