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

//  tree_node.cpp
//  simul
//
//  Created by Fausto Fabian Crespo Fernandez on 2019-10-08.
//
#include <vector>
#include "tree_node.hpp"


void init_to_empty_str(char str[MAX_NAME])
{
    for (size_t i = 0; i < MAX_NAME; i++) {
        str[i] = 0;
    }
}

TreeNode::TreeNode(int numSites)
{
    left = NULL;
    right = NULL;
    anc1 = NULL;
    outgroup = NULL;
    //    nodeRight=NULL;
    //    nodeLeft=NULL;
    //    nodeBack=NULL;
    time = 0;
    timePUnits = 0;
    scaledByThetaTimeInputTreeUnits=0;
    timeInputTreeUnits=0;
    length = 0;
    lengthModelUnits = 0;
    index = 0;
    label = 0;
    isOutgroup = NO;
    nodeClass = 0;
    indexOldClone = 0;
    indexCurrentClone = 0;
    orderCurrentClone = 0;
    effectPopSize=0;
    isLeaf = NO;
    init_to_empty_str(cellName);
    init_to_empty_str(observedCellName);
    numberOfTipsSubTree=0;
    initSequenceVectors(numSites);
    
}
void TreeNode::initNumberTipsVector(int numberClones)
{
    for (size_t i = 0; i < numberClones; i++) {
        numberTipsByPopulation.push_back(0);
    }
}
void TreeNode::resetNumberTipsVector(int numberClones)
{
    for (size_t i = 0; i < numberClones; i++) {
        numberTipsByPopulation.at(i)=0;
    }
}
void TreeNode::initSequenceVectors(int size){
    for (size_t i = 0; i < size; i++) {
        
        maternalSequence.push_back(0);
        paternalSequence.push_back(0);
        genotypeSequence.push_back(0);
        numbersMutationsUnderSubtreePerSite.push_back(0);
        numbersMaternalMutationsPerSite.push_back(0);
        numbersPaternalMutationsPerSite.push_back(0);        
    }
    
}

