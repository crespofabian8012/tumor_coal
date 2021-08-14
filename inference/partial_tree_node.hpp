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
 
  PartialTreeEdge(PLLBufferManager *manager, std::shared_ptr<PartialTreeNode> child,
                double length, unsigned int pmatrix_elements, size_t alignment);

  void showPMatrix( unsigned int states, unsigned int rate_cats, unsigned int states_padded,  unsigned int float_precision);
  ~PartialTreeEdge();

  std::shared_ptr<PartialTreeNode> child;

  double length;
  double *pmatrix;
  unsigned int matrix_size;//number of elements in matrix
};

class PartialTreeNode {

  PLLBufferManager *manager;

public:

  PartialTreeNode(PLLBufferManager *manager,
                std::shared_ptr<PartialTreeEdge> edge_l,
                std::shared_ptr<PartialTreeEdge> edge_r, std::string label,
                double height, unsigned int clv_size,
                unsigned int scale_buffer_size, size_t alignment,
                unsigned  int index);

  void getNormalizeCLV(unsigned int sites, unsigned int states, std::vector<double> &norm_clv);
  void showClV(unsigned int states, unsigned int rate_cats, unsigned int states_padded,
               unsigned int sites, unsigned int float_precision)const;

  void showpClV(unsigned int states, unsigned int rate_cats, unsigned int states_padded,unsigned int sites, unsigned int float_precision)const;
  static void  unscale(double * prob, unsigned int times);
  void init(int clv_elements);

  std::shared_ptr<PartialTreeEdge> edge_l;
  std::shared_ptr<PartialTreeEdge> edge_r;
    
  std::vector<double>& getCLV(){return clv;};
  void buildCLV(int tip_id,int numberStates, pll_msa_t *msa, GenotypeErrorModel *gtErrorModel,  bool normalize);
  double * getCLVPointer(){return clv.data();};
  unsigned int * getScaleBufferPointer(){return scale_buffer.data();};
  std::vector<unsigned int>& getScaleBuffer(){return scale_buffer;};

  std::string label;
  double height;
  double ln_likelihood;
  unsigned int index_population;
  unsigned int index;
    
  unsigned int clv_size;//number of elements in clv and scale_buffer
  
  double *pclv;
  unsigned int *pscale_buffer;
  std::vector<double> clv;
  std::vector<unsigned int> scale_buffer;
    
~PartialTreeNode();
 
  
};


#endif /* partial_tree_node_hpp */







