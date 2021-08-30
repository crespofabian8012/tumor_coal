//
//  tree.hpp
//  run
//


#ifndef tree_hpp
#define tree_hpp

#include "pll_utils.hpp"
#include <stdio.h>
#include <vector>
#include <set>
#include <unordered_set>

extern "C"
    {
#include "libpll/pll.h"
    }
namespace std
    {
    template<>
    struct default_delete<pll_rtree_t> {
        void operator()(pll_rtree_t* ptr) { pll_rtree_destroy(ptr, nullptr); }
    };
    }
class RootedTree
{
protected:
    size_t num_tips;
    
    std::unique_ptr<pll_rtree_t> pll_rooted_tree;
    std::vector<std::vector<double>> partition_brlens;
    
    mutable std::vector<pll_rnode_t*>  pll_rtree_tips;
    
    std::vector<pll_rnode_t*>  const& tip_nodes() const;
    std::vector<pll_rnode_t*>  subnodes() const;
    std::string outgroup_label;
    bool has_outgroup;
public:
    RootedTree() : pll_rooted_tree(nullptr) {}
    RootedTree(unsigned int tip_count, const pll_rnode_t& root) :
    num_tips(tip_count),
    pll_rooted_tree(pll_rtree_wraptree(pll_utils::cloneRNode(&root), tip_count)) {}
    RootedTree(const pll_rtree_t& pll_rtree) :
    num_tips(pll_rtree.tip_count), pll_rooted_tree(pll_utils::cloneRTree(&pll_rtree)) {}
    
    RootedTree(const pll_rtree_t& pll_rtree, bool usingNewick ) :
    num_tips(pll_rtree.tip_count), pll_rooted_tree( pll_rtree_parse_newick_string( pll_rtree_export_newick( pll_rtree.root, NULL))) {}
    
    RootedTree(const std::string &str, bool isFile = true);
    
    RootedTree(std::unique_ptr<pll_rtree_t>&  pll_rtree) :
    num_tips(pll_rtree ? pll_rtree->tip_count : 0), pll_rooted_tree(pll_rtree.release()) {}
    RootedTree(std::unique_ptr<pll_rtree_t>&&  pll_rtree) :
    num_tips(pll_rtree ? pll_rtree->tip_count : 0), pll_rooted_tree(pll_rtree.release()) {}
    
    RootedTree (const RootedTree& other);
    RootedTree& operator=(const RootedTree& other);
    RootedTree (RootedTree&& other);
    RootedTree& operator=(RootedTree&& other);
    
    size_t numTips() const { return num_tips; };
    size_t numInner() const { return num_tips - 1; };
    size_t numNodes() const { return numTips() + numInner(); };
    size_t numSubnodes() const { return numBranches() * 2; };
    size_t numBranches() const { return num_tips ? num_tips + num_tips - 2 : 0; };
    size_t numSplits() const { return numBranches() - num_tips; };
    
    pll_rtree_t *getPtr() {return pll_rooted_tree.get();}
    
    const pll_rtree_t *getPtr() const {return pll_rooted_tree.get();}
    
    std::unordered_set<std::string> getLabels(bool leavesOnly) const;
    std::vector<std::string> getLabelList(bool leavesOnly) const;
    
    std::vector<pll_rnode_t*> getNodes() const;
    std::vector<pll_rnode_t*> getLeaves() const;
    
    
    RootedTree & getSubtreeWithoutOutgroup(std::string &outgroup_label)const;
    pll_utree_t *getUnrootedPtr()const;
    pll_utree_t *unRootedTreeCopy() const;
    
    void fixMissingBranchLengths(double new_brlen);
    virtual ~RootedTree();
    
    friend std::ostream& operator<<(std::ostream& os, const RootedTree &tree)
    {
        char *newick = pll_rtree_export_newick(tree.getPtr()->root, 0);
        std::string str(newick);
        os << str;
        free(newick);
        return os;
    }
};

class CoalescentTree:public RootedTree
{
    public:
        CoalescentTree();
        
};
#endif /* tree_hpp */
