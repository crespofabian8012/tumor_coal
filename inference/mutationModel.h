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

//  mutationModel.h
//  TumorCellCoal
//
//  Created by Fausto Fabian Crespo Fernandez on 7/8/19.
//  Copyright © 2019 Fausto Fabian Crespo Fernandez. All rights reserved.
//

#ifndef mutationModel_h
#define mutationModel_h

#include "data_types.hpp"
#include "tree_node.hpp"

#include <gsl/gsl_rng.h>
#include <boost/random.hpp>

void SimulateISM (TreeNode *treeRooot, int genome, int doISMhaploid, long int *seed,  std::vector<int> &DefaultModelSites, int numDefaultModelSites, std::vector<int> &AltModelSites, int numAltModelSites, long  double   totalTreeLength , int &numISMmutations, int numFixedMutations, int numSNVmaternal, int doSimulateFixedNumMutations,  int alphabet, std::vector<SiteStr> &allSites, int  &numMU, long double  cumMij[4][4], long double  mutationRate, const gsl_rng *randomGsl, boost::mt19937* rngBoost);

void SimulateISMDNAforSite (TreeNode *p, int genome, int site, int doISMhaploid, long int *seed, long double  totalTreeLength, std::vector<SiteStr> &allSites, int &numMU, long double  cumMij[4][4], long double  mutationRate, long double  &uniform, long double  &cumBranchLength, long double  &ran, const gsl_rng *randomGsl, boost::mt19937* rngBoost);

void SimulateISMforSite (TreeNode *p, int genome, int site, int doISMhaploid, long int *seed, long double  totalTreeLength, std::vector<SiteStr> &allSites, int &numMU, long double  cumMij[4][4], long double  mutationRate,  long double  &cumBranchLength, long double  &uniform, int &mutationAdded,const  gsl_rng *randomGsl, boost::mt19937* rngBoost);


void SimulateISMGenotypeforSite (TreeNode *p, int genome, int site, int doISMhaploid, long int *seed, long double  totalTreeLength, std::vector<SiteStr> &allSites, int  &numMU, long double  cumMij[4][4],long double  mutationRate, long double  &uniform, long double  &cumBranchLength, long double  &ran, const gsl_rng *randomGsl, boost::mt19937* rngBoost, bool &added );

void GetTreeDepthFirstSearchOrder (TreeNode *p, long double  mutationRate, int &nextAvailableIndex, std::vector<TreeNode *> &dfs, std::vector<long double> &cumSumScaledBranchLengths,  std::vector<long double> &scaledBranchLengths,  double &cumScaledTreeLength, long double scaledTotalTreeLength);

void PlaceMutationsOnBranchOnSite(TreeNode *p, int site, int numMutations, int  &numMU, std::vector<SiteStr> &allSites, long int *seed, const gsl_rng *randomGsl, boost::mt19937* rngBoost);
void PropagateGenotypesForSiteAfterPositionInPreorder(int pos, int site, std::vector<TreeNode*> preorder);
//void SimulateMk2 (TreeNode *p, int genome, long int *seed, vector<int> &AltModelSites, int  numAltModelSites, int doUserTree, int rateVarAmongSites, long double  altModelMutationRate, vector<SiteStr> &allSites, int &numMU);

//void SimulateMk2forSite (TreeNode *p, int genome, int site, long int *seed, int doUserTree, int rateVarAmongSites, long double  altModelMutationRate, vector<SiteStr> &allSites, int &numMU );

void SimulateFiniteDNA (TreeNode *p, int genome, long int *seed, int doJC, int doHKY, int doGTR, int doGTnR, long double& freqR, long double  &freqY, long double  &freqAG, long double  &freqCT,  double  titv,  double  freq[4],  double  Mij[4][4], int numAltModelSites, std::vector<int> &AltModelSites, std::vector<SiteStr> &allSites,  int rateVarAmongSites, long double  altModelMutationRate, int &numMU,  double  Root[],  double  Cijk[], const gsl_rng *randomGsl, boost::mt19937* rngBoost);
void AllocateCellStructure(CellStr* cell, int numCells, int numSites);

void EvolveSitesOnTree (TreeNode *treeRoot, int genome, long int *seed, int rateVarAmongSites, int numSites, std::vector<SiteStr> &allSites, int doGeneticSignatures, int alphaSites, int  propAltModelSites , int numDefaultModelSites, int numAltModelSites, std::vector<int> &DefaultModelSites, std::vector<int> &AltModelSites,  long double  totalTreeLength , int &numISMmutations, int numFixedMutations, int numSNVmaternal, int doSimulateFixedNumMutations,  int alphabet, int  &numMU, long double  cumMij[4][4], int altModel, long double  altModelMutationRate, int doUserTree, int doJC, int doHKY, int doGTR, int doGTnR, long double  freqR, long double  freqY, long double  freqAG, long double  freqCT, double  titv, double  freq[4], double  Mij[4][4],   double  Root[],  double  Cijk[], const gsl_rng *randomGsl, boost::mt19937* rngBoost);

void EvolveGenotypesOnTree (TreeNode *treeRoot, int totalSampleSize,  std::vector<int> numberOfSitesWithKMutations, int numberVariableSites, long int *seed, int rateVarAmongSites, int numSites, std::vector<SiteStr> &allSites, int doGeneticSignatures, int alphaSites, int  propAltModelSites , int numDefaultModelSites, int numAltModelSites, std::vector<int> &DefaultModelSites, std::vector<int> &AltModelSites,  long double  totalTreeLength , int &numISMmutations, int numFixedMutations, int numSNVmaternal, int doSimulateFixedNumMutations,  int alphabet, int  &numMU, long double  cumMij[4][4], int altModel, long double  altModelMutationRate, int doUserTree, int doJC, int doHKY, int doGTR, int doGTnR, long double  freqR, long double  freqY, long double  freqAG, long double  freqCT, double  titv, double  freq[4], double  Mij[4][4],   double  Root[],  double  Cijk[],
                            const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost  );

//void SimulateFiniteDNAforSite (TreeNode *p, int genome, int site,vector<SiteStr> &allSites,  long int *seed, int rateVarAmongSites, long double  altModelMutationRate, int &numMU, int doJC, int doHKY, int doGTR, int doGTnR, long double  beta,  long double  kappa, long double  freqR, long double  freqY, long double  freq[4],  double  Root[],  double  Cijk[]);

void SimulateFiniteDNAforSite (TreeNode *p, int genome, int site,std::vector<SiteStr> &allSites,  long int *seed, int rateVarAmongSites, long double  altModelMutationRate, int &numMU, int doJC, int doHKY, int doGTR, int doGTnR, long double  beta,  long double  kappa,  long double&  freqR,  long double&  freqY,  double  freq[4],  double  Root[],  double  Cijk[], const gsl_rng *randomGsl, boost::mt19937* rngBoost);


void FillSubstitutionMatrix (long double  ch_prob[4][4], long double  branchLength, int doJC, int doHKY, int doGTR, int doGTnR, long double  beta, long double  kappa, long double  freqR, long double  freqY,  double  freq[4],  double  Root[],  double  Cijk[]);

void FillGenotypeGTJCSubstitutionMatrix (long double  ch_prob[16][16], long double  branchLength, long double  theta,  double  freq[16], double maxEntry, double  maxPerRow[16]);

void FillRowQGTJCSubstitutionMatrix ( long double  row[16], int indexCurrentGenotype);

int GetNewGenotypeIdx ( int indexCurrentGenotype, int selectedIdx);
void JCmodel (long double  Pij[4][4], long double  branchLength, long double  beta);
void HKYmodel (long double  Pij[4][4], long double  branchLength, long double  kappa, long double  freqR, long double  freqY, long double  beta,  double  freq[4]);
void GTRmodel (long double  Pij[4][4], long double  branchLength,  double  Root[],  double  Cijk[]);


void SimulateTriNucFreqGenome (int cell, long int *seed, TreeNode *p, int alphabet, int doUserGenome, int numSites, std::vector<SiteStr> &allSites, int doGeneticSignatures,  double  cumfreq[4],  long double  *triNucFreq, const gsl_rng *randomGsl, boost::mt19937* rngBoost);

void EvolveDeletionsOnTree (TreeNode *p, int genome, std::vector<SiteStr> &allSites, int &numISMdeletions,int numSites,  long double totalTreeLength, long double mutationRate, long double deletionRate,int &numDEL,  long int *seed, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost);

void SimulateDeletionforSite (TreeNode *p, int genome, int site, std::vector<SiteStr> &allSites, long double totalTreeLength, long double mutationRate,long double &cumDeletionBranchLength,  long double deletionRate, int &numDEL, long int *seed);

void EvolveCNLOHonTree (TreeNode *p, int genome, int &numISMCNLOH,std::vector<SiteStr> &allSites, int numSites, long int *seed, long double CNLOHrate, long double mutationRate, long double totalTreeLength,long double &cumCNLOHbranchLength,   const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost );

void SimulateCNLOHforSite (TreeNode *p, int genome, int site,long int *seed,std::vector<SiteStr> &allSites,long double CNLOHrate, long double  &cumCNLOHbranchLength, long double totalTreeLength, long double mutationRate, int  &numCNLOH, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost);

void AllelicDropout (int numCells,std::vector<SiteStr> &allSites, int doADOcell, int doADOsite, int numSites,long double fixedADOrate, long double meanADOcell, long double varADOcell, long double meanADOsite,long double varADOsite, std::vector<std::shared_ptr<TreeNode>> &nodes,    long int *seed, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost);

void addAllelicDropoutToTree( int numCells, std::vector<std::shared_ptr<TreeNode>> &nodes,std::vector<SiteStr> &allSites, int numSites, std::vector<std::vector<double> > &alleleADOrate, long int *seed, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost );

void FillErrorMatrix (long double  error_prob[4][4], double Eij[4][4]);


void GenotypeError (std::vector<std::shared_ptr<TreeNode>> &nodes,std::vector<SiteStr> &allSites, int alphabet, int numSites, int numCells,  double meanGenotypingError,  double varGenotypingError, double genotypingError, double Eij[4][4], long int *seed, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost);


void SequenceError (std::vector<std::shared_ptr<TreeNode>> &nodes,std::vector<SiteStr> &allSites, int alphabet, int numSites, int numCells,  double meanGenotypingError, int sampleSize, std::vector<int> numberOfSitesWithKSequencingErrors, double Eij[4][4], long int *seed, int &numberSeqErrorsAdded, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost);
//static void AddGermlineVariation (TreeNode *treeRoot,long int *seed, int numSites, long double  SNPrate, std::vector<SiteStr> &allSites, int alphabet, int HEALTHY_ROOT, long double  cumMij[4][4]);
#endif /* mutationModel_h */
