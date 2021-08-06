//
//  partial_tree_node.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 26/07/2021.
//

#include "partial_tree_node.hpp"

#include <iostream>
extern "C"
    {
#include "libpll/pll.h"
    }

PartialTreeEdge::PartialTreeEdge(PLLBufferManager *manager,
                             std::shared_ptr<PartialTreeNode> child,
                             double length, unsigned int pmatrix_size, size_t alignment)
    : manager(manager), child(child), length(length) {
  if (manager->pmatrix_buffer.empty()) {
    //pmatrix = (double *)std::malloc(pmatrix_size);
     pmatrix =  (double*)pll_aligned_alloc(pmatrix_size, alignment);
  } else {
      
      std::cout<< "reusing pmatrix from PLLBufferManager"<< std::endl;
    pmatrix = manager->pmatrix_buffer.top();
    manager->pmatrix_buffer.pop();

    std::memset(pmatrix, 0, pmatrix_size);
  }
}

PartialTreeEdge::~PartialTreeEdge() {
  manager->pmatrix_buffer.push(pmatrix);
    
  pll_aligned_free(pmatrix);
  pmatrix = nullptr;
}

PartialTreeNode::PartialTreeNode(PLLBufferManager *manager,
                             std::shared_ptr<PartialTreeEdge> edge_l,
                             std::shared_ptr<PartialTreeEdge> edge_r,
                             std::string label, double height,
                             unsigned int clv_size,
                             unsigned int scale_buffer_size
                             , size_t alignment, unsigned int index)
    : manager(manager), edge_l(edge_l), edge_r(edge_r), label(label),
      height(height), index(index) {
  if (manager->clv_buffer.empty()) {
    //clv = (double *)std::malloc(clv_size);
      clv =  (double*)pll_aligned_alloc(clv_size, alignment);
  } else {
    
    std::cout<< "reusing clv from PLLBufferManager"<< std::endl;
    clv = manager->clv_buffer.top();
    manager->clv_buffer.pop();

    std::memset(clv, 0, clv_size);
  }

  if (manager->scale_buffer_buffer.empty()) {
    scale_buffer = (unsigned int *)std::malloc(scale_buffer_size);
     
  } else {
      
       std::cout<< "reusing scale_buffer from PLLBufferManager"<< std::endl;
    scale_buffer = manager->scale_buffer_buffer.top();
    manager->scale_buffer_buffer.pop();

    std::memset(scale_buffer, 0, scale_buffer_size);
  }
          ln_likelihood=0.0;
}

PartialTreeNode::~PartialTreeNode() {
  manager->clv_buffer.push(clv);
  manager->scale_buffer_buffer.push(scale_buffer);

  pll_aligned_free(clv);
  clv = nullptr;
  scale_buffer = nullptr;
  

}
