//
//  state.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 16/03/2021.
//


#include "state.hpp"
#include "pll_utils.hpp"
#include <memory>

State::State(long double currentModelTime){
    this->currentModelTime = currentModelTime;
}

State::State(unsigned int num_sites, const pll_partition_t *partition):
partition(partition),num_sites(num_sites)
{
    currentModelTime = 0.0;
}

State::State(unsigned int num_sites, const pll_partition_t *partition,  int sampleSize, pll_msa_t *msa, std::vector<int> &positions, ProgramOptions &programOptions):
partition(partition),num_sites(num_sites)
{
    currentModelTime = 0.0;
    
    initForest(  sampleSize, msa, positions,  programOptions);
}


void State::setForest(std::vector<pll_rtree_t *> forestPar){
    forest= forestPar;
}
void State::initForest( int sampleSize, pll_msa_t *msa, std::vector<int> &positions, ProgramOptions &programOptions){
    
    assert(sampleSize==positions.size());
    
    const  pll_partition_t *reference_partition = partition;

    //pll_rtree_t * tree;
    std::shared_ptr<pll_rnode_t>  p;
    TreeNode *t;
    for(size_t i=0; i< positions.size();++i){
        
        unsigned int sites_alloc =
            reference_partition->asc_bias_alloc ? reference_partition->sites + reference_partition->states : reference_partition->sites;

        unsigned int clv_size =
            sites_alloc * reference_partition->states_padded * reference_partition->rate_cats * sizeof(double);
        unsigned int scaler_size = (reference_partition->attributes & PLL_ATTRIB_RATE_SCALERS)
                                       ? sites_alloc * reference_partition->rate_cats
                                       : sites_alloc;

        
        p=std::make_shared<pll_rnode_t>();
        p->left =NULL;
        p->right=NULL;
        p->parent=NULL;
        p->label =new char[strlen(msa->label[i]) + 1]{};
        std::copy(msa->label[i], msa->label[i] + strlen(msa->label[i]), p->label);
        p->length = 0.0;
        //p->scaler_index = i;
        //p->clv_index = i;
        
        t = new TreeNode(msa->length);
        std::copy(p->label, p->label + strlen(p->label), t->observedCellName);
        t->isLeaf=true;
        t->indexSequenceMSA= i;
        t->time = 0.0;//model time
        t->timePUnits = 0.0;// model time * theta * x
        //t->genotypeSequence = ObservedData[i];
        p->data = t;
        
        
//        delete (reference_partition->clv);
//        delete (reference_partition->scale_buffer);
//        reference_partition->clv = reference_partition->clv[i];
//        reference_partition->scale_buffer = nullptr;
//        reference_partition->ln_likelihood = compute_ln_likelihood(node->clv, nullptr, reference_partition);
//
//        tree = new pll_rtree_t();
//        tree->root = p;
//        roots.push_back(p);
//        forest.push_back(p);
    }
    currentModelTime = 0.0;
}
std::vector<pll_rtree_t *> State::getForest() const{
    return forest;
}
pll_rtree_t * State::getTreeAt(int i) const{
    return forest[i];
}
std::shared_ptr<State> State::makeDeepCopy(){
    
    std::shared_ptr<State> newState = std::make_shared<State>(currentModelTime);
    std::vector<pll_rtree_t *> forestCopy;
    pll_rtree_t *tree;
    pll_rtree_t *treeCopy;
    
    for(size_t i=0; i< forest.size();++i){
        tree = forest[i];
        char * newick=pll_rtree_export_newick(tree->root, NULL);
        //pll_utree_graph_clone
        treeCopy =pll_rtree_parse_newick_string(newick);
        treeCopy = pll_utils::cloneRTree(tree);
      
        forestCopy.push_back(treeCopy);
    }
    newState->setForest(forestCopy);
    return newState;
}

State &State::operator=(const State &original){
    
    if (this == &original)
       return *this;

     partition = original.partition;
     currentModelTime = original.currentModelTime;
     roots = original.roots;

     return *this;
    
}

std::vector<std::shared_ptr<pll_rtree_t*>> State::propose(int i, int j, double height){
    
    std::vector<std::shared_ptr<pll_rtree_t*>> result;
    
    return result;
}



double State::likelihood_factor(std::shared_ptr<pll_rtree_t*> root){
    
    double result=0.0;
    
    return result;
    
}

void State::remove_roots(int i, int j){}

State::~State(){
    
    
    
}
