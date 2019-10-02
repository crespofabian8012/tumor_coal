//
//  likelihood.h
//  TumorCellCoal
//
//  Created by Fausto Fabian Crespo Fernandez on 7/8/19.
//  Copyright © 2019 Fausto Fabian Crespo Fernandez. All rights reserved.
//

#ifndef likelihood_h
#define likelihood_h

static double LogConditionalLikelihoodTree(TreeNode  *tree, TreeNode *nodes, Population **populations, int numClones);

static double  LogConditionalLikelihoodSequences(Chain *chain, pll_msa_t * msa, char* NewickString, ProgramOptions *programOptions);

void compute_state_probs(unsigned int state, double  **clvp,  unsigned int _states,  int from, int to,
                         double _seq_error_rate,
                         double _dropout_rate);

void set_partition_tips(Chain *chain, pll_partition_t * partition, pll_msa_t * msa, ProgramOptions *programOptions);
#endif /* likelihood_h */
