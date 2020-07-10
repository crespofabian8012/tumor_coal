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
 * tree node
 */

#ifndef tree_node_hpp
#define tree_node_hpp

#include <vector>

#include  "definitions.hpp"


class TreeNode
{
public:
    
    TreeNode *left, *right, *anc1, *outgroup;
    int         index, label, isOutgroup;
    long double      length, time,lengthModelUnits, timePUnits;
    long double      oldlength, oldtime, oldlengthModelUnits, oldtimePUnits;
    long double      timeInputTreeUnits;
    int         nodeClass; // 0: leaf, 1: internal, 2: 3: 4: outgroup 5: healthy cells (TODO: check later)
    int         indexOldClone, indexCurrentClone,orderCurrentClone;
    //indexCoalClone;
    long double      effectPopSize;
    long double      oldeffectPopSize;
    char        cellName[200];
    char        observedCellName[200];
    std::vector<int>        maternalSequence;
    std::vector<int>        paternalSequence;
    std::vector<int>        numbersMutationsUnderSubtreePerSite;
    std::vector<int>        numbersMaternalMutationsPerSite;
    std::vector<int>        numbersPaternalMutationsPerSite;
    bool         isLeaf;
    int numberOfTipsSubTree;
    std::vector<int> numberTipsByPopulation;
    TreeNode(int numSites);
    void initNumberTipsVector(int numberClones);
    void resetNumberTipsVector(int numberClones);
    void initSequenceVectors(int size);
    ~TreeNode()
    {
        delete left;
        delete right;
        delete anc1;
        delete outgroup;
        
    }
};

typedef TreeNode* pTreeNode;


#endif /* tree_node_hpp */
