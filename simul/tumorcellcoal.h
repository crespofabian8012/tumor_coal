//
//  tumorcellcoal.h
//  TumorCellCoal
//
//  Created by Fausto Fabian Crespo Fernandez on 6/5/19.
//  Copyright © 2019 Fausto Fabian Crespo Fernandez. All rights reserved.
//

#ifndef tumorcellcoal_h
#define tumorcellcoal_h

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ftw.h>

#include "definitions.hpp"

extern double   Qij[16], mr;
extern double   ***selectedSignature;
extern double ****geneticSignature;

extern double **signatureProbs;
extern double *triNucFreq;

#endif /* tumorcellcoal_h */
