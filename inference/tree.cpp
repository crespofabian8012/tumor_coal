//
//  tree.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 26/08/2021.
//

#include "tree.hpp"
#include "parsers.hpp"
#include <iterator>
RootedTree::RootedTree(const std::string &str, bool isFile ):pll_rooted_tree(Parsers::readRooted(str,isFile)){
    
    if (!pll_rooted_tree){
        std::cout << "Error creating rooted tree from newick" << std:: endl;
        
    }
    num_tips = pll_rooted_tree->tip_count;
    // setMissingBranchLengths();
}
RootedTree::RootedTree (const RootedTree& other) : num_tips(other.num_tips),
    pll_rooted_tree(pll_utils::cloneRTree(other.getPtr())), partition_brlens(other.partition_brlens)
{
}

RootedTree::RootedTree (RootedTree&& other) :num_tips(other.num_tips), pll_rooted_tree(other.pll_rooted_tree.release())
{
  other.num_tips = 0;
  swap(pll_rtree_tips, other.pll_rtree_tips);
  swap(partition_brlens, other.partition_brlens);
}
pll_utree_t *RootedTree::getUnrootedPtr()const{
    //unsigned int tip_count = 0;
    pll_rtree_t * rtree_copy = pll_rtree_parse_newick_string(pll_rtree_export_newick( pll_rooted_tree->root, NULL));
    pll_utree_t *utree =  pll_rtree_unroot(rtree_copy);
    pll_unode_t * vroot = utree->nodes[utree->tip_count+utree->inner_count-1];
    pll_utree_reset_template_indices(vroot, utree->tip_count);
    return(utree);
};
//this will only look outgroup_label in the left and right child nodes of root!
// is not efficient method
RootedTree& RootedTree::getSubtreeWithoutOutgroup(std::string &outgroup_label) const{
    
    std::unordered_set<std::string> leaves_labels = getLabels(true);
    RootedTree *subtree = nullptr;
    if (leaves_labels.find(outgroup_label) != leaves_labels.end()){
        
        pll_rnode_t* left = pll_rooted_tree->root->left;
        pll_rnode_t* right = pll_rooted_tree->root->right;
        if (pll_utils::isLeaf(left) && std::string(left->label).compare(outgroup_label)==0){
            pll_rtree_t* right_tree_copy = pll_rtree_parse_newick_string( pll_rtree_export_newick( right, NULL));
            
            std::unique_ptr<pll_rtree_t> u_ptr= std::make_unique<pll_rtree_t>(*right_tree_copy);
            subtree = new RootedTree(u_ptr);
            
            return(*subtree);
        }
        else if (pll_utils::isLeaf(right) && std::string(right->label).compare(outgroup_label)==0){
            
            pll_rtree_t* left_tree_copy = pll_rtree_parse_newick_string( pll_rtree_export_newick( left, NULL));
            std::unique_ptr<pll_rtree_t> u_ptr= std::make_unique<pll_rtree_t>(*left_tree_copy);
            subtree = new RootedTree(u_ptr);
            return(*subtree);
        }
        else{
            std::cout << "Outgroup is not direct child of root!";
        }
    }
    else{
        std::cout << "Ourgroup label not found!";
    }
    return(*subtree);
}
//pll_utree_t * RootedTree::getUnrootedSubtreeWithoutOutgroupPtr(std::string &outgroup_label)const{
//    
//    RootedTree& subtree = getSubtreeWithoutOutgroup(outgroup_label);
//    
//    return(pll_rtree_unroot(subtree->pll_rooted_tree));
//    
//};
std::unordered_set<std::string> RootedTree::getLabels(bool leavesOnly) const{
    std::unordered_set<std::string> result;

     for (auto node: (leavesOnly ? getLeaves() : getNodes())) {
       if (node->label) {
         result.insert(node->label);
       }
     }
     
    assert(!result.empty());

    return result;
}
std::vector<std::string> RootedTree::getLabelList(bool leavesOnly) const{
    std::vector<std::string> result;
    
    for (auto const& node: (leavesOnly ? getLeaves() : getNodes()))
      result.emplace_back(std::string(node->label));

    assert(!result.empty());

    return result;
}
std::vector<pll_rnode_t*> RootedTree::getLeaves() const
{
    std::vector<pll_rnode_t*> leaves(num_tips);
    
    if (leaves.empty() && num_tips > 0)
    {
      assert(num_tips == pll_rooted_tree->tip_count);

      leaves.assign(pll_rooted_tree->nodes, pll_rooted_tree->nodes + pll_rooted_tree->tip_count);
    }

//    for (auto it = pll_rooted_tree->nodes; it != pll_rooted_tree->nodes+ num_tips; ++it) {
//        leaves.push_back(*it);
//    }
    return leaves;
}
std::vector<pll_rnode_t*> RootedTree::getNodes() const
{
    std::vector<pll_rnode_t*> nodes(numNodes());

    for (auto it = pll_rooted_tree->nodes; it != pll_rooted_tree->nodes+ numNodes(); ++it) {
        nodes.push_back(*it);
    }
    return nodes;
}
void RootedTree::fixMissingBranchLengths(double new_brlen)
{
  pll_utree_t * utree= unRootedTreeCopy();
  if (utree)
    pllmod_utree_set_length_recursive(utree, new_brlen, 1);
}
pll_utree_t * RootedTree::unRootedTreeCopy() const
{
  pll_utree_t *utree = pll_rtree_unroot(pll_rooted_tree.get());
  return utree ? pll_utree_clone(utree) : nullptr;
}

//pll_rtree_t * RootedTree::rootedTreeCopy() const
//{
//
//   pll_utree_t * utree= unRootedTreeCopy();
//    if (utree){
//
//        pllmod_utree_root_inplace(utree);
//    }
//
//}
RootedTree::~RootedTree(){
}
