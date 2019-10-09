//
//  tree_node.hpp
//  simul
//
//  Created by Seong-Hwan Jun on 2019-10-08.
//

#ifndef tree_node_hpp
#define tree_node_hpp

#include <vector>
#include <libpll/pll_tree.h>

#include "definitions.hpp"

using namespace std;

class TreeNode
{
public:
//    pll_unode_t     *nodeLeft;
//    pll_unode_t     *nodeRight;
//    pll_unode_t     *nodeBack;
//    pll_tree_edge_t  *edgeLeft;
//    pll_tree_edge_t  *edgeRight;
//    pll_tree_edge_t  *edgeBack;

    TreeNode *left, *right, *anc1, *outgroup;
    int         index, label, isOutgroup;
    double      length, time,lengthModelUnits, timePUnits;
    int         nodeClass; // 0: leaf, 1: internal, 2: 3: 4: outgroup 5: healthy cells (TODO: check later)
    int         indexOldClone, indexCurrentClone,orderCurrentClone;
    //indexCoalClone;
    double      effectPopSize;
    char        cellName[MAX_NAME];
    char        observedCellName[MAX_NAME];
    vector<int>        maternalSequence;
    vector<int>        paternalSequence;
    vector<int>        numbersMutationsUnderSubtreePerSite;
    vector<int>        numbersMaternalMutationsPerSite;
    vector<int>        numbersPaternalMutationsPerSite;
    bool         isLeaf;
    
    TreeNode();
};

typedef TreeNode* pTreeNode;


#endif /* tree_node_hpp */
