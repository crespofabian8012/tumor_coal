//
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

TreeNode::TreeNode()
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
    
}
void TreeNode::initNumberTipsVector(int numberClones)
{
    for (size_t i = 0; i < numberClones; i++) {
        numberTipsByPopulation.push_back(0);
    }
}