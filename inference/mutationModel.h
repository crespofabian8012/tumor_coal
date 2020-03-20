//
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

void SimulateISM (TreeNode *treeRooot, int genome, int doISMhaploid, long int *seed,  vector<int> &DefaultModelSites, int numDefaultModelSites, vector<int> &AltModelSites, int numAltModelSites, long  double   totalTreeLength , int &numISMmutations, int numFixedMutations, int numSNVmaternal, int doSimulateFixedNumMutations,  int alphabet, vector<SiteStr> &allSites, int  &numMU, long double  cumMij[4][4], long double  mutationRate);

void SimulateISMDNAforSite (TreeNode *p, int genome, int site, int doISMhaploid, long int *seed, long double  totalTreeLength, vector<SiteStr> &allSites, int &numMU, long double  cumMij[4][4], long double  mutationRate, long double  &uniform, long double  &cumBranchLength, long double  &ran);

void SimulateISMforSite (TreeNode *p, int genome, int site, int doISMhaploid, long int *seed, long double  totalTreeLength, vector<SiteStr> &allSites, int &numMU, long double  cumMij[4][4], long double  mutationRate,  long double  &cumBranchLength, long double  &uniform, int &mutationAdded);

//void SimulateMk2 (TreeNode *p, int genome, long int *seed, vector<int> &AltModelSites, int  numAltModelSites, int doUserTree, int rateVarAmongSites, long double  altModelMutationRate, vector<SiteStr> &allSites, int &numMU);

//void SimulateMk2forSite (TreeNode *p, int genome, int site, long int *seed, int doUserTree, int rateVarAmongSites, long double  altModelMutationRate, vector<SiteStr> &allSites, int &numMU );

void SimulateFiniteDNA (TreeNode *p, int genome, long int *seed, int doJC, int doHKY, int doGTR, int doGTnR, long double& freqR, long double  &freqY, long double  &freqAG, long double  &freqCT,  double  titv,  double  freq[4],  double  Mij[4][4], int numAltModelSites, vector<int> &AltModelSites, vector<SiteStr> &allSites,  int rateVarAmongSites, long double  altModelMutationRate, int &numMU,  double  Root[],  double  Cijk[]);
void AllocateCellStructure(CellStr* cell, int numCells, int numSites);

void EvolveSitesOnTree (TreeNode *treeRoot, int genome, long int *seed, int rateVarAmongSites, int numSites, vector<SiteStr> &allSites, int doGeneticSignatures, int alphaSites, int  propAltModelSites , int numDefaultModelSites, int numAltModelSites, vector<int> &DefaultModelSites, vector<int> &AltModelSites,  long double  totalTreeLength , int &numISMmutations, int numFixedMutations, int numSNVmaternal, int doSimulateFixedNumMutations,  int alphabet, int  &numMU, long double  cumMij[4][4], int altModel, long double  altModelMutationRate, int doUserTree, int doJC, int doHKY, int doGTR, int doGTnR, long double  freqR, long double  freqY, long double  freqAG, long double  freqCT, double  titv, double  freq[4], double  Mij[4][4],   double  Root[],  double  Cijk[]);


//void SimulateFiniteDNAforSite (TreeNode *p, int genome, int site,vector<SiteStr> &allSites,  long int *seed, int rateVarAmongSites, long double  altModelMutationRate, int &numMU, int doJC, int doHKY, int doGTR, int doGTnR, long double  beta,  long double  kappa, long double  freqR, long double  freqY, long double  freq[4],  double  Root[],  double  Cijk[]);

void SimulateFiniteDNAforSite (TreeNode *p, int genome, int site,vector<SiteStr> &allSites,  long int *seed, int rateVarAmongSites, long double  altModelMutationRate, int &numMU, int doJC, int doHKY, int doGTR, int doGTnR, long double  beta,  long double  kappa,  long double&  freqR,  long double&  freqY,  double  freq[4],  double  Root[],  double  Cijk[]);

void FillSubstitutionMatrix (long double  ch_prob[4][4], long double  branchLength, int doJC, int doHKY, int doGTR, int doGTnR, long double  beta, long double  kappa, long double  freqR, long double  freqY,  double  freq[4],  double  Root[],  double  Cijk[]);
void JCmodel (long double  Pij[4][4], long double  branchLength, long double  beta);
void HKYmodel (long double  Pij[4][4], long double  branchLength, long double  kappa, long double  freqR, long double  freqY, long double  beta,  double  freq[4]);
void GTRmodel (long double  Pij[4][4], long double  branchLength,  double  Root[],  double  Cijk[]);


void SimulateTriNucFreqGenome (int cell, long int *seed, TreeNode *p, int alphabet, int doUserGenome, int numSites, vector<SiteStr> &allSites, int doGeneticSignatures,  double  cumfreq[4],  long double  *triNucFreq );

static void AddGermlineVariation (TreeNode *treeRoot,long int *seed, int numSites, long double  SNPrate, vector<SiteStr> &allSites, int alphabet, int HEALTHY_ROOT, long double  cumMij[4][4]);
#endif /* mutationModel_h */
