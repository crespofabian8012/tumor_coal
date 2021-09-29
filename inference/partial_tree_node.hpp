//
//  partial_tree_node.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 26/07/2021.
//

#ifndef partial_tree_node_hpp
#define partial_tree_node_hpp

#include "pll_buffer_manager.hpp"
#include "genotype_error_model.hpp"
//#include <stdio.h>
#include <memory>
#include <string>
#include <vector>
//#include <cstring>

class PartialTreeNode;

class PartialTreeEdge {
    PLLBufferManager *manager;
public:
    double length;
    double *pmatrix;
    unsigned int matrix_size;//number of elements in matrix
    std::shared_ptr<PartialTreeNode> child;
    //constructors
    
    PartialTreeEdge(PLLBufferManager *manager,
                    std::shared_ptr<PartialTreeNode> const& child,
                    double length, unsigned int pmatrix_elements, size_t alignment);
    
    PartialTreeEdge(PLLBufferManager *manager,
                    std::shared_ptr<PartialTreeNode> const& child,
                    double length, unsigned int pmatrix_elements, double * pmatrix, size_t alignment);
    
    PartialTreeEdge(const PartialTreeEdge& original);
    PartialTreeEdge& operator= (PartialTreeEdge rhs);
    PartialTreeEdge(PartialTreeEdge&& rhs);
    
    std::shared_ptr<PartialTreeEdge> Clone();
    std::unique_ptr<PartialTreeEdge> CloneUnique();
    void showPMatrix( unsigned int states, unsigned int rate_cats, unsigned int states_padded,  unsigned int float_precision);
    void freeMatrix();
    virtual ~PartialTreeEdge(){};
    
};

class PartialTreeNode {
    
    PLLBufferManager *manager;
    
public:
    
    PartialTreeNode(PLLBufferManager *manager,
                    std::shared_ptr<PartialTreeNode> const& leftChild,
                    std::shared_ptr<PartialTreeNode> const& rightChild,
                    double leftLength,
                    double rightLength,
                    unsigned int pmatrix_elements,
                    std::string label,
                    double height,
                    unsigned int clv_size,
                    unsigned int scale_buffer_size,
                    size_t alignment,
                    unsigned  int index,
                    unsigned  int index_population,
                    double *left_pmatrix,
                    double *right_pmatrix);
    
//    PartialTreeNode(PLLBufferManager *manager,
//                    std::unique_ptr<PartialTreeEdge> edge_l,
//                    std::unique_ptr<PartialTreeEdge> edge_r,
//                    std::string label,
//                    double height, unsigned int clv_size,
//                    unsigned int scale_buffer_size, size_t alignment,
//                    double *pclv,
//                    unsigned int *pscale_buffer,
//                    unsigned  int index,
//                    unsigned  int index_population,
//                    double             ln_likelihood);
    
    PartialTreeNode( PartialTreeNode& original);
    PartialTreeNode& operator= (PartialTreeNode rhs) ;
    PartialTreeNode& operator=( const PartialTreeNode& rhs );
    PartialTreeNode(PartialTreeNode&& rhs);
    PartialTreeNode& operator=( PartialTreeNode&& rhs );
    
    void getNormalizeCLV(unsigned int sites, unsigned int states, std::vector<double> &norm_clv);
    void showClV(unsigned int states, unsigned int rate_cats, unsigned int states_padded,
                 unsigned int sites, unsigned int float_precision)const;
    
    void showpClV(unsigned int states, unsigned int rate_cats, unsigned int states_padded,unsigned int sites, unsigned int float_precision)const;
    static void  unscale(double * prob, unsigned int times);
    void init(int clv_elements);
    std::shared_ptr<PartialTreeNode> Clone();
    
    std::unique_ptr<PartialTreeEdge> edge_l;
    std::unique_ptr<PartialTreeEdge> edge_r;
    
    //std::vector<double>& getCLV(){return clv;};
    void buildCLV(int tip_id,int numberStates, pll_msa_t *msa, GenotypeErrorModel *gtErrorModel,  bool normalize);
    void freeVectors();
    
    // double * getCLVPointer(){return pclv.get();};
    // unsigned int * getScaleBufferPointer(){return pscale_buffer.get();};
    //std::vector<unsigned int>& getScaleBuffer(){return scale_buffer;};
    
    std::string label;
    double height;
    double ln_likelihood;
    double ln_coal_likelihood;
  
    unsigned int index_population;
    unsigned int index;
    unsigned int number_nodes_cluster;
    unsigned int number_leaves_cluster;
    
    unsigned int clv_size;//number of elements in clv and scale_buffer
    
    double *pclv;
    unsigned int *pscale_buffer;
    
    int getIndexLeftChild() const{
         if (edge_l) return edge_l->child->index;
         else return -1;
    };
    int getIndexRightChild() const{
        if (edge_r) return edge_r->child->index;
        else return -1;
        
    };
    double likelihood_factor()const;
    
    PLLBufferManager* getManager()const{return manager;};
    // std::unique_ptr<unsigned int[]> pscale_buffer;
    // std::unique_ptr<double[]> pclv;
    // std::unique_ptr<double> pclv;
    //std::vector<double> clv;
    //std::vector<unsigned int> scale_buffer;
    
    virtual ~PartialTreeNode(){};

    };
template<typename... Ts>
auto delPartialTreeNode = [](PartialTreeNode* pTreeNode)
{
    PLLBufferManager *manager = pTreeNode->getManager();
    manager->clv_buffer.push(pTreeNode->pclv);
    manager->scale_buffer_buffer.push(pTreeNode->pscale_buffer);
    
    pTreeNode->pclv = nullptr;
    pTreeNode->pscale_buffer = nullptr;
    
    manager->pmatrix_buffer.push(pTreeNode->edge_l->pmatrix);
    manager->pmatrix_buffer.push(pTreeNode->edge_r->pmatrix);
    // pll_aligned_free(pmatrix);
    pTreeNode->edge_l->pmatrix = nullptr;
    pTreeNode->edge_r->pmatrix = nullptr;
    delete pTreeNode;
};


#endif /* partial_tree_node_hpp */
    
    
    
    
    
    
    
