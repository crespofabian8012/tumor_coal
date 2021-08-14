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
                             double length, unsigned int pmatrix_elements, size_t alignment)
    : manager(manager), child(child), length(length), matrix_size(pmatrix_elements) {
  if (manager->pmatrix_buffer.empty()) {
    pmatrix = (double *)std::malloc(pmatrix_elements*sizeof(double));
    //pmatrix =  (double*)pll_aligned_alloc(pmatrix_size, PLL_ALIGNMENT_SSE);
  } else {
      
      std::cout<< "reusing pmatrix from PLLBufferManager"<< std::endl;
    pmatrix = manager->pmatrix_buffer.top();
    manager->pmatrix_buffer.pop();

    std::memset(pmatrix, 0, pmatrix_elements*sizeof(double));
  }
}
void PartialTreeEdge::showPMatrix( unsigned int states, unsigned int rate_cats, unsigned int states_padded,  unsigned int float_precision){
    
    unsigned int i,j,k;

    for (k = 0; k < rate_cats; ++k)
     {
       for (i = 0; i < states; ++i)
       {
         for (j = 0; j < states; ++j)
           printf("%+2.*f   ", float_precision, pmatrix[i*states_padded+j]);
         printf("\n");
       }
       printf("\n");
     }
}
PartialTreeEdge::~PartialTreeEdge() {
  manager->pmatrix_buffer.push(pmatrix);
    
 // pll_aligned_free(pmatrix);
  pmatrix = nullptr;
}

PartialTreeNode::PartialTreeNode(PLLBufferManager *manager,
                             std::shared_ptr<PartialTreeEdge> edge_l,
                             std::shared_ptr<PartialTreeEdge> edge_r,
                             std::string label, double height,
                             unsigned int clv_elements,
                             unsigned int scale_buffer_elements
                             , size_t alignment, unsigned int index)
    : manager(manager), edge_l(edge_l), edge_r(edge_r), label(label),
      height(height), index(index), clv_size(clv_elements)
 //, clv(clv_elements,0.0), scale_buffer(clv_elements,0.0)
{
//  if (manager->clv_buffer.empty()) {
    pclv = (double *)std::malloc(clv_elements*sizeof(double));
//    //clv =  (double*)pll_aligned_alloc(clv_size, PLL_ALIGNMENT_SSE);
//  } else {
//
//    std::cout<< "reusing clv from PLLBufferManager"<< std::endl;
//    clv = manager->clv_buffer.top();
//    manager->clv_buffer.pop();
//
//    std::memset(clv, 0, clv_size);
//  }
//
//  if (manager->scale_buffer_buffer.empty()) {
     pscale_buffer = (unsigned int *)std::malloc(scale_buffer_elements*sizeof(unsigned int));
//      //scale_buffer =  (unsigned int *)pll_aligned_alloc(clv_size, PLL_ALIGNMENT_SSE);
//
//  } else {
//
//       std::cout<< "reusing scale_buffer from PLLBufferManager"<< std::endl;
//    scale_buffer = manager->scale_buffer_buffer.top();
//    manager->scale_buffer_buffer.pop();
//
//    std::memset(scale_buffer, 0, scale_buffer_size);
//  }
         // init(clv_elements);
          ln_likelihood=0.0;
}
void PartialTreeNode::getNormalizeCLV(unsigned int sites, unsigned int states, std::vector<double> &norm_clv )
{

    norm_clv.resize(clv.size());
    assert(norm_clv.size()==clv.size());
    /*iterators*/
    std::vector<double>::iterator clvp = clv.begin();
    std::vector<double>::iterator norm_clvp = norm_clv.begin();
    size_t k;
    for (size_t j = 0; j < sites; ++j)
       {
         double tsum = 0;
         for (k = 0; k < states; ++k){
                tsum += clvp[k];
           }
         if (tsum == 0.0)
            std::cout<<" CLV for site " << j << " cannot be normalized! "<< std::endl;
        for ( k = 0; k < states; ++k)
           norm_clvp[k] = norm_clvp[k]/tsum;

         norm_clvp += states;
         clvp += states;

       }
    
    assert(clvp == clv.end());
    assert(norm_clvp == norm_clv.end());
    
}
void PartialTreeNode::showClV(unsigned int states, unsigned int rate_cats, unsigned int states_padded,
                              unsigned int sites, unsigned int float_precision)const{
    
    unsigned int i,j,k;

     double prob;

//     if ((clv_index < numberTips) &&
//         (partition->attributes & PLL_ATTRIB_PATTERN_TIP))
//       return;

     printf ("[ ");
     for (i = 0; i < sites; ++i)
     {
       printf("{");
       for (j = 0; j < rate_cats; ++j)
       {
         printf("(");
         for (k = 0; k < states-1; ++k)
         {
           prob = clv[i*rate_cats*states_padded + j*states_padded + k];
          //if (scale_buffer) PartialTreeNode::unscale(&prob, scale_buffer[i]);
           printf("%.*f,", float_precision, prob);
         }
         prob = clv[i*rate_cats*states_padded + j*states_padded + k];
         //if (scale_buffer) PartialTreeNode::unscale(&prob, scale_buffer[i]);
         printf("%.*f)", float_precision, prob);
         if (j < rate_cats - 1) printf(",");
       }
       printf("} ");
     }
     printf ("]\n");
}
void PartialTreeNode::showpClV(unsigned int states, unsigned int rate_cats, unsigned int states_padded,unsigned int sites, unsigned int float_precision)const {
    
    unsigned int i,j,k;
    double prob;
    const double * pclv2 = pclv;
    
     printf ("[ ");
     for (i = 0; i < sites; ++i)
     {
       printf("{");
       for (j = 0; j < rate_cats; ++j)
       {
         printf("(");
         for (k = 0; k < states-1; ++k)
         {
           prob = pclv2[i*rate_cats*states_padded + j*states_padded + k];
          //if (scale_buffer) PartialTreeNode::unscale(&prob, scale_buffer[i]);
           printf("%.*f,", float_precision, prob);
         }
         prob = pclv2[i*rate_cats*states_padded + j*states_padded + k];
         //if (scale_buffer) PartialTreeNode::unscale(&prob, scale_buffer[i]);
         printf("%.*f)", float_precision, prob);
         if (j < rate_cats - 1) printf(",");
       }
       printf("} ");
     }
     printf ("]\n");
   
    
}
void PartialTreeNode::buildCLV(int tip_id,int numberStates, pll_msa_t *msa, GenotypeErrorModel *gtErrorModel,  bool normalize) {
    
    auto clv_size = msa->length * numberStates;
    //auto clv_size = msa->length * numberStates*numberRateCats;
    //here we assume numberRateCats=1
    //auto clvp = clv.begin();

    auto seq = msa->sequence[tip_id];
    //auto charmap = _model.charmap();
    auto charmap = pll_map_gt10;
    
    for (size_t j = 0; j < msa->length; ++j)
    {
        auto charstate = (pll_state_t) seq[j];
        pll_state_t state = charmap ? charmap[(int) charstate] : charstate;
        
        gtErrorModel->computeStateErrorProbPT17(state, pclv);
        
//        if (true)
//        {
//              printf("state: %llu ", state);
//              printf("char: %c ", seq[j]);
//              for (size_t k = 0; k < numberStates; ++k)
//                printf("%lf ", pclv[k]);
//              printf("\n");
//            for (size_t k = 0; k < numberStates; ++k)
//              printf("%lf ", pclv[k]);
//            printf("\n");
//         }
        
        pclv += numberStates;
       //  clvp += numberStates*numberRateCats;
    }
    
    //restore the pointer to the beginning
    pclv-= clv_size;
}
 void  PartialTreeNode::unscale(double * prob, unsigned int times)
{
  unsigned int i;

  for (i = 0; i < times; ++i)
    *prob *= PLL_SCALE_THRESHOLD;
}
void PartialTreeNode::init(int clv_elements){
    for (size_t i = 0; i < clv_elements; i++) {
        
        clv.push_back(0.0);
        scale_buffer.push_back(0);
  
    }
}


PartialTreeNode::~PartialTreeNode() {
 // manager->clv_buffer.push(clv);
  //manager->scale_buffer_buffer.push(scale_buffer);

  //pll_aligned_free(clv);
  //pll_aligned_free(scale_buffer);
    
    pclv = nullptr;
    pscale_buffer = nullptr;
   // delete clv ;
  // delete[] scale_buffer ;

}
