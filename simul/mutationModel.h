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

void SimulateISM (TreeNode *treeRooot, int genome, int doISMhaploid, long int *seed,  vector<int> &DefaultModelSites, int numDefaultModelSites, vector<int> &AltModelSites, int numAltModelSites, double totalTreeLength , int &numISMmutations, int numFixedMutations, int numSNVmaternal, int doSimulateFixedNumMutations,  int alphabet, vector<SiteStr> &allSites, int  &numMU, double cumMij[4][4], double mutationRate);

void SimulateISMDNAforSite (TreeNode *p, int genome, int site, int doISMhaploid, long int *seed, double totalTreeLength, vector<SiteStr> &allSites, int &numMU, double cumMij[4][4], double mutationRate, double &uniform, double &cumBranchLength, double &ran);

void SimulateISMforSite (TreeNode *p, int genome, int site, int doISMhaploid, long int *seed, double totalTreeLength, vector<SiteStr> &allSites, int &numMU, double cumMij[4][4], double mutationRate,  double &cumBranchLength, double &uniform, int &mutationAdded);

//void SimulateMk2 (TreeNode *p, int genome, long int *seed, vector<int> &AltModelSites, int  numAltModelSites, int doUserTree, int rateVarAmongSites, double altModelMutationRate, vector<SiteStr> &allSites, int &numMU);

//void SimulateMk2forSite (TreeNode *p, int genome, int site, long int *seed, int doUserTree, int rateVarAmongSites, double altModelMutationRate, vector<SiteStr> &allSites, int &numMU );

void SimulateFiniteDNA (TreeNode *p, int genome, long int *seed, int doJC, int doHKY, int doGTR, int doGTnR, double freqR, double freqY, double freqAG, double freqCT, double titv, double freq[4], double Mij[4][4], int numAltModelSites, vector<int> &AltModelSites, vector<SiteStr> &allSites,  int rateVarAmongSites, double altModelMutationRate, int &numMU, double Root[], double Cijk[]);
void AllocateCellStructure(CellStr* cell, int numCells, int numSites);

void EvolveSitesOnTree (TreeNode *treeRoot, int genome, long int *seed, int rateVarAmongSites, int numSites, vector<SiteStr> &allSites, int doGeneticSignatures, int alphaSites, int  propAltModelSites , int numDefaultModelSites, int numAltModelSites, vector<int> &DefaultModelSites, vector<int> &AltModelSites,  double totalTreeLength , int &numISMmutations, int numFixedMutations, int numSNVmaternal, int doSimulateFixedNumMutations,  int alphabet,    int  &numMU, double cumMij[4][4], int altModel, double altModelMutationRate, int doUserTree,int doJC, int doHKY, int doGTR, int doGTnR, double freqR, double freqY, double freqAG, double freqCT, double titv, double freq[4], double Mij[4][4],  double Root[], double Cijk[] );

void SimulateFiniteDNAforSite (TreeNode *p, int genome, int site,vector<SiteStr> &allSites,  long int *seed, int rateVarAmongSites, double altModelMutationRate, int &numMU, int doJC, int doHKY, int doGTR, int doGTnR, double beta,  double kappa, double freqR, double freqY, double freq[4], double Root[], double Cijk[]);

void FillSubstitutionMatrix (double ch_prob[4][4], double branchLength, int doJC, int doHKY, int doGTR, int doGTnR, double beta, double kappa, double freqR, double freqY, double freq[4], double Root[], double Cijk[]);
void JCmodel (double Pij[4][4], double branchLength, double beta);
void HKYmodel (double Pij[4][4], double branchLength, double kappa, double freqR, double freqY, double beta, double freq[4]);
void GTRmodel (double Pij[4][4], double branchLength, double Root[], double Cijk[]);


void SimulateTriNucFreqGenome (int cell, long int *seed, TreeNode *p, int alphabet, int doUserGenome, int numSites, vector<SiteStr> &allSites, int doGeneticSignatures, double cumfreq[4], double *triNucFreq );

static void AddGermlineVariation (TreeNode *treeRoot,long int *seed, int numSites, double SNPrate, vector<SiteStr> &allSites, int alphabet, int HEALTHY_ROOT, double cumMij[4][4]);
#endif /* mutationModel_h */
