//
//  tree_node.hpp
//  simul
//
//  Created by Fausto Fabian Crespo Fernandez on 2019-10-08.
//

#ifndef tree_node_hpp
#define tree_node_hpp

#include <vector>

#include "definitions.hpp"

using namespace std;

class TreeNode
{
public:

    TreeNode *left, *right, *anc1, *outgroup;
    int         index, label, isOutgroup;
    double      length, time,lengthModelUnits, timePUnits;
    double      oldlength, oldtime, oldlengthModelUnits, oldtimePUnits;
    double      timeInputTreeUnits;
    int         nodeClass; // 0: leaf, 1: internal, 2: 3: 4: outgroup 5: healthy cells (TODO: check later)
    int         indexOldClone, indexCurrentClone,orderCurrentClone;
    //indexCoalClone;
    double      effectPopSize;
    double      oldeffectPopSize;
    char        cellName[MAX_NAME];
    char        observedCellName[MAX_NAME];
    vector<int>        maternalSequence;
    vector<int>        paternalSequence;
    vector<int>        numbersMutationsUnderSubtreePerSite;
    vector<int>        numbersMaternalMutationsPerSite;
    vector<int>        numbersPaternalMutationsPerSite;
    bool         isLeaf;
    int numberOfTipsSubTree;
    vector<int> numberTipsByPopulation;
    TreeNode();
    void initNumberTipsVector(int numberClones);
    void resetNumberTipsVector(int numberClones);
};

typedef TreeNode* pTreeNode;


#endif /* tree_node_hpp */
