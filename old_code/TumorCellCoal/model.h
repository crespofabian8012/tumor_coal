//
//  model.h
//  TumorCellCoal
//
//  Created by Fausto Fabian Crespo Fernandez on 7/8/19.
//  Copyright © 2019 Fausto Fabian Crespo Fernandez. All rights reserved.
//

#ifndef model_h
#define model_h


static double FmodelTstandard (double t, double TOrigin, double delta);

static double GstandardTmodel (double V, double TOrigin, double delta);
static double CalculateH (double t, double TOrigin, double delta);
double DensityTime(double delta, double u);
static double LogDensityCoalescentTimesForPopulation(TreeNode  *tree, TreeNode *nodes,  Population **populations, int numClones);
double LogProbNoCoalescentEventBetweenTimes(Population *popI,double from, double to, int numberActiveInd);

#endif /* model_h */
