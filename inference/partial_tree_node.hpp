//
//  partial_tree_node.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 26/07/2021.
//

#ifndef partial_tree_node_hpp
#define partial_tree_node_hpp

#include "pll_buffer_manager.hpp"
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
                double length, unsigned int pmatrix_size, size_t alignment);

  ~PartialTreeEdge();

  std::shared_ptr<PartialTreeNode> child;

  double length;
  double *pmatrix;
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


  ~PartialTreeNode();

  std::shared_ptr<PartialTreeEdge> edge_l;
  std::shared_ptr<PartialTreeEdge> edge_r;

  std::string label;
  double height;
  double ln_likelihood;
  unsigned int index_population;
  unsigned int index;
    
  double *clv;
  unsigned int *scale_buffer;
};


#endif /* partial_tree_node_hpp */







