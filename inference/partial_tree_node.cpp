//
//  partial_tree_node.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 26/07/2021.
//

#include "partial_tree_node.hpp"

PartialTreeEdge::PartialTreeEdge(PLLBufferManager *manager,
                             std::shared_ptr<PartialTreeNode> child,
                             double length, unsigned int pmatrix_size)
    : manager(manager), child(child), length(length) {
  if (manager->pmatrix_buffer.empty()) {
    pmatrix = (double *)std::malloc(pmatrix_size);
  } else {
    pmatrix = manager->pmatrix_buffer.top();
    manager->pmatrix_buffer.pop();

    std::memset(pmatrix, 0, pmatrix_size);
  }
}

PartialTreeEdge::~PartialTreeEdge() {
  manager->pmatrix_buffer.push(pmatrix);
  pmatrix = nullptr;
}

PartialTreeNode::PartialTreeNode(PLLBufferManager *manager,
                             std::shared_ptr<PartialTreeEdge> edge_l,
                             std::shared_ptr<PartialTreeEdge> edge_r,
                             std::string label, double height,
                             unsigned int clv_size,
                             unsigned int scale_buffer_size)
    : manager(manager), edge_l(edge_l), edge_r(edge_r), label(label),
      height(height) {
  if (manager->clv_buffer.empty()) {
    clv = (double *)std::malloc(clv_size);
  } else {
    clv = manager->clv_buffer.top();
    manager->clv_buffer.pop();

    std::memset(clv, 0, clv_size);
  }

  if (manager->scale_buffer_buffer.empty()) {
    scale_buffer = (unsigned int *)std::malloc(scale_buffer_size);
  } else {
    scale_buffer = manager->scale_buffer_buffer.top();
    manager->scale_buffer_buffer.pop();

    std::memset(scale_buffer, 0, scale_buffer_size);
  }
}

PartialTreeNode::~PartialTreeNode() {
  manager->clv_buffer.push(clv);
  manager->scale_buffer_buffer.push(scale_buffer);

  clv = nullptr;
  scale_buffer = nullptr;
}
