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
 * utils functions
 */
#ifndef utils_h
#define utils_h
#include <vector>


#include "data_types.hpp"
#include "population.hpp"
#include "tree_node.hpp"


//class ProgramOptions;
//class TreeNode;
//class Population;
//using namespace std;
struct Utils
{
static void ReadParametersFromFastaFile(char *fileName, int &numCells, int &TotalNumSequences, int &numSites);
static void ReadFastaFile(char *fileName, vector<vector<int> > &ObservedData,  char **ObservedCellNames, ProgramOptions &programOptions);

static int CheckMatrixSymmetry(double matrix[4][4]);
static int WhichNucChar (char nucleotide);
static vector<int> SequenceToIntegers(char ** sequence, int length);
static vector<int> GenotypesToIntegers(char ** sequence, int length);
static int ChooseUniformState ( double *prob, long int *seed);
static char WhichIUPAC (int allele1, int allele2);
static char *WhichGenotypeFromIUPAC (int  iupac);
static int WhichGenotypeChar (char nucleotide);
static char WhichNuc (int nucleotide);
static char WhichConsensusBinary (int allele1, int allele2);
static char WhichMut (int state);
static void  normalizeVector(double *vector, int length);
static int compareIntDescending(const void *a, const void *b);
static int CompareGenotypes (int a1, int a2, int b1, int b2);

static double * expand_uniq_rates(int states, const double * uniq_rates, const int * rate_sym);

static void unscale(double * prob, unsigned int times);

static double DistanceBetweenSFS(int* SFS1, int*SFS, int numSNVs,  int numSites);


static int computeTajimaD(TreeNode **treeTips,  int numSites, int numCells);

static double ComputeESS(double *weights, int numberWeights);
//void computeUnfoldedISMSFS(int numSites,SiteStr* allSites,int numSNVs, int* SNVsites, int* SFS, int *numberDifferences);

static void Initialize( double (*Eij)[4], double (*Mij)[4], double *freq,  ProgramOptions &programOptions ) ;

static void computeGenotypesFreq(double freqs[10], pll_msa_t * msa);

static void InitNumberNodes(double *TotalBirthRate, double *TotalDeathRate, int *TotalN,  Population **populations, ProgramOptions &programOptions) ;


static double parameterMultiplierMCMCmove (double lengthInterval);
    
static   char * cb_serialize(const pll_rnode_t * node);

static  char * rTreeToNewick(const pll_rnode_t * root,
                             char * (*cb_serialize)(const pll_rnode_t *));
static char* appendCharToCharArray(char* array, char a);
};
#endif /* utils_h */
