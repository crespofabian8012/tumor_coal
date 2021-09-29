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
                                 std::shared_ptr<PartialTreeNode> const& child,
                                 double length, unsigned int pmatrix_elements, size_t alignment)
: manager(manager),  length(length),  matrix_size(pmatrix_elements), child(child) {
    if (manager->pmatrix_buffer.empty()) {
        pmatrix = (double *)std::malloc(pmatrix_elements*sizeof(double));
        //pmatrix =  (double*)pll_aligned_alloc(pmatrix_size, PLL_ALIGNMENT_SSE);
    } else {
        
        //std::cout<< "reusing pmatrix from PLLBufferManager"<< std::endl;
        pmatrix = manager->pmatrix_buffer.top();
        manager->pmatrix_buffer.pop();
        
        std::memset(pmatrix, 0, pmatrix_elements*sizeof(double));
    }
}
PartialTreeEdge::PartialTreeEdge(PLLBufferManager *manager,
                                 std::shared_ptr<PartialTreeNode> const& child,
                                 double length, unsigned int pmatrix_elements, double * pmatrixP, size_t alignment)
: manager(manager), length(length), matrix_size(pmatrix_elements), child(child) {
    
    if (manager->pmatrix_buffer.empty()) {
        pmatrix = (double *)std::malloc(pmatrix_elements*sizeof(double));
        //pmatrix =  (double*)pll_aligned_alloc(pmatrix_size, PLL_ALIGNMENT_SSE);
    } else {
        
        //std::cout<< "reusing pmatrix from PLLBufferManager"<< std::endl;
        pmatrix = manager->pmatrix_buffer.top();
        manager->pmatrix_buffer.pop();
    }
    memcpy(pmatrix, pmatrixP, pmatrix_elements * sizeof(double));
}

PartialTreeEdge::PartialTreeEdge(const PartialTreeEdge& original)
{
    manager = original.manager;
    length = original.length;
    matrix_size = original.matrix_size;
    manager            = original.manager;
    pmatrix = new double(*original.pmatrix);
    //shallow copy of child
    child = original.child;
    
}
PartialTreeEdge& PartialTreeEdge::operator= (PartialTreeEdge rhs) {
    std::swap(*this, rhs);
    return *this;
}
PartialTreeEdge::PartialTreeEdge(PartialTreeEdge&& rhs) : child(std::move(rhs.child)) {
    
    manager = rhs.manager;
    length = rhs.length;
    matrix_size = rhs.matrix_size;
    manager            = rhs.manager;
    pmatrix = new double(*rhs.pmatrix);
}
std::shared_ptr<PartialTreeEdge> PartialTreeEdge::Clone(){
    
    std::shared_ptr<PartialTreeEdge> result;
    
    result = std::make_shared<PartialTreeEdge>(
                                               manager, child->Clone(), length,matrix_size , pmatrix,  0);
    
    return result;
}
std::unique_ptr<PartialTreeEdge> PartialTreeEdge::CloneUnique(){
    
    std::unique_ptr<PartialTreeEdge> result;
    
    result = std::make_unique<PartialTreeEdge>(
                                               manager, child->Clone(), length,matrix_size , pmatrix,  0);
    
    return result;
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
void PartialTreeEdge::freeMatrix(){
    
    delete []pmatrix;
    
}
//PartialTreeEdge::~PartialTreeEdge() {
//    manager->pmatrix_buffer.push(pmatrix);
//    
//    // pll_aligned_free(pmatrix);
//    pmatrix = nullptr;
//}

PartialTreeNode::PartialTreeNode(PLLBufferManager *manager,
                                 std::shared_ptr<PartialTreeNode> const& leftChild,
                                 std::shared_ptr<PartialTreeNode> const& rightChild,
                                 double leftLength,
                                 double rightLength,
                                 unsigned int pmatrix_elements,
                                 std::string label,
                                 double height,
                                 unsigned int clv_elements,
                                 unsigned int scale_buffer_elements,
                                 size_t alignment,
                                 unsigned  int index,
                                 unsigned  int index_population,
                                 double *left_pmatrix,
                                 double *right_pmatrix)
: manager(manager), label(label),
height(height),   index_population(index_population),index(index), clv_size(clv_elements)
//, clv(clv_elements,0.0), scale_buffer(clv_elements,0.0)
{
    
    
    if (manager->clv_buffer.empty()) {
        // pclv = std::make_unique<double[]>(clv_elements);
        
        pclv = (double *)std::malloc(clv_elements*sizeof(double));
        //    //clv =  (double*)pll_aligned_alloc(clv_size, PLL_ALIGNMENT_SSE);
    } else {
        
        //std::cout<< "reusing clv from PLLBufferManager"<< std::endl;
        pclv = manager->clv_buffer.top();
        manager->clv_buffer.pop();
        
        std::memset(pclv, 0, clv_elements*sizeof(double));
    }
    
    if (manager->scale_buffer_buffer.empty()) {
        // pscale_buffer = std::make_unique<unsigned int[]>(clv_elements);
        
        pscale_buffer = (unsigned int *)std::malloc(scale_buffer_elements*sizeof(unsigned int));
        //scale_buffer =  (unsigned int *)pll_aligned_alloc(clv_size, PLL_ALIGNMENT_SSE);
        
    } else {
        
        //std::cout<< "reusing scale_buffer from PLLBufferManager"<< std::endl;
        pscale_buffer = manager->scale_buffer_buffer.top();
        manager->scale_buffer_buffer.pop();
        
        std::memset(pscale_buffer, 0, scale_buffer_elements*sizeof(unsigned int));
    }
   
    if (leftChild!= nullptr && leftLength >0){
        if (left_pmatrix!=nullptr){
           edge_l = std::make_unique<PartialTreeEdge>(manager, leftChild, leftLength, pmatrix_elements,left_pmatrix,  alignment);
            
        }
        else
          edge_l = std::make_unique<PartialTreeEdge>(
                                                manager, leftChild, leftLength, pmatrix_elements, alignment);
        
        
    }
    else
        edge_l = nullptr;
    if (rightChild!= nullptr && rightLength >0){
      if (right_pmatrix!=nullptr){
                edge_r = std::make_unique<PartialTreeEdge>(manager, rightChild, rightLength, pmatrix_elements,  right_pmatrix,  alignment);
                 
             }
    else
      edge_r = std::make_unique<PartialTreeEdge>(manager, rightChild, rightLength, pmatrix_elements, alignment);
    }
    else{
        edge_r = nullptr;
        
    }
    ln_likelihood=0.0;
    ln_coal_likelihood = 0.0;
    number_leaves_cluster = 0;
    number_nodes_cluster = 0;
    
    if (edge_l!=nullptr && edge_r!=nullptr){
        if(edge_l->child!=nullptr){
            number_nodes_cluster += edge_l->child->number_nodes_cluster + 1;
            number_leaves_cluster += edge_l->child->number_leaves_cluster ;
        }
        if(edge_r->child!=nullptr){
            number_nodes_cluster += edge_r->child->number_nodes_cluster;
            number_leaves_cluster += edge_r->child->number_leaves_cluster;
        }
        //   std::cout<< " number leaves cluster " << number_leaves_cluster << std::endl;
    }
    else{
        number_leaves_cluster = 1;
        number_nodes_cluster = 1;
    }
    
}
//PartialTreeNode::PartialTreeNode(PLLBufferManager *manager,
//                                 std::unique_ptr<PartialTreeEdge> edge_l,
//                                 std::unique_ptr<PartialTreeEdge> edge_r,
//                                 std::string label,
//                                 double height, unsigned int clv_size,
//                                 unsigned int scale_buffer_size, size_t alignment,
//                                 double *pclv,
//                                 unsigned int *pscale_buffer,
//                                 unsigned  int index,
//                                 unsigned  int index_population,
//                                 double             ln_likelihood):
//manager(manager), edge_l(std::move(edge_l)), edge_r(std::move(edge_r)), label(label),
//height(height),   ln_likelihood(ln_likelihood),  index_population(index_population), index(index),  clv_size(clv_size),  pclv(pclv), pscale_buffer(pscale_buffer) {
//
//}
PartialTreeNode::PartialTreeNode( PartialTreeNode& original): edge_l( new PartialTreeEdge( *original.edge_l ) ),  edge_r( new PartialTreeEdge( *original.edge_r ) )
{
    manager = original.manager;
    height =  original.height;
    clv_size =  original.clv_size;
    index =  original.index;
    index_population =  original.index_population;
    label =  original.label;
    ln_likelihood =  original.ln_likelihood;
    number_nodes_cluster = original.number_nodes_cluster;
    number_leaves_cluster = original.number_leaves_cluster;
    
    
    //edge_l = std::move(original.edge_l);
    //pclv(std::make_unique<double[]>(*original.pclv))
    // pclv(std::move(original.pclv));
    // pclv = std::move( original.pclv );
    // pscale_buffer = std::move( original.pscale_buffer );
    pclv = original.pclv;
    pscale_buffer = original.pscale_buffer;
    //pclv( new double[](original.pclv.get()));
    //pscale_buffer( original.pscale_buffer );
    
}
PartialTreeNode& PartialTreeNode::operator=( const PartialTreeNode& rhs )
{
    manager = rhs.manager;
    height =  rhs.height;
    clv_size =  rhs.clv_size;
    index =  rhs.index;
    index_population =  rhs.index_population;
    label =  rhs.label;
    ln_likelihood =  rhs.ln_likelihood;
    ln_coal_likelihood = rhs.ln_coal_likelihood;
    number_nodes_cluster = rhs.number_nodes_cluster;
    number_leaves_cluster = rhs.number_leaves_cluster;
    
    edge_l.reset( new PartialTreeEdge( *rhs.edge_l ) );
    edge_r.reset( new PartialTreeEdge( *rhs.edge_r ) );
    
    pclv = rhs.pclv;
    pscale_buffer = rhs.pscale_buffer;
    //pclv( new double[](original.pclv.get()));
    //pscale_buffer( original.pscale_buffer );
    return *this;
    
}
PartialTreeNode& PartialTreeNode::operator=( PartialTreeNode&& rhs )
{
    manager = rhs.manager;
    pclv = std::move( rhs.pclv );
    pscale_buffer = std::move( rhs.pscale_buffer );
    height =  rhs.height;
    clv_size =  rhs.clv_size;
    index =  rhs.index;
    index_population =  rhs.index_population;
    label =  rhs.label;
    ln_likelihood =  rhs.ln_likelihood;
    ln_coal_likelihood = rhs.ln_coal_likelihood;
    number_nodes_cluster = rhs.number_nodes_cluster;
    number_leaves_cluster = rhs.number_leaves_cluster;
    
    edge_l = std::move( rhs.edge_l );
    edge_r = std::move( rhs.edge_r );
    return *this;
}
//PartialTreeNode& PartialTreeNode::operator= (PartialTreeNode& rhs) {
//    std::swap(*this, rhs);
//    return *this;
//}
PartialTreeNode::PartialTreeNode(PartialTreeNode&& rhs){
    
    manager = rhs.manager;
    height =  rhs.height;
    clv_size =  rhs.clv_size;
    index =  rhs.index;
    index_population =  rhs.index_population;
    label =  rhs.label;
    ln_likelihood =  rhs.ln_likelihood;
    ln_coal_likelihood = rhs.ln_coal_likelihood;
    number_nodes_cluster = rhs.number_nodes_cluster;
    number_leaves_cluster = rhs.number_leaves_cluster;
    pclv = std::move( rhs.pclv );
    pscale_buffer = std::move( rhs.pscale_buffer );
    
    
}
void PartialTreeNode::getNormalizeCLV(unsigned int sites, unsigned int states, std::vector<double> &norm_clv )
{
    
    //    norm_clv.resize(clv.size());
    //    assert(norm_clv.size()==clv.size());
    //    /*iterators*/
    //    std::vector<double>::iterator clvp = clv.begin();
    //    std::vector<double>::iterator norm_clvp = norm_clv.begin();
    //    size_t k;
    //    for (size_t j = 0; j < sites; ++j)
    //    {
    //        double tsum = 0;
    //        for (k = 0; k < states; ++k){
    //            tsum += clvp[k];
    //        }
    //        if (tsum == 0.0)
    //            std::cout<<" CLV for site " << j << " cannot be normalized! "<< std::endl;
    //        for ( k = 0; k < states; ++k)
    //            norm_clvp[k] = norm_clvp[k]/tsum;
    //
    //        norm_clvp += states;
    //        clvp += states;
    //
    //    }
    //
    //    assert(clvp == clv.end());
    //    assert(norm_clvp == norm_clv.end());
    
}
void PartialTreeNode::showClV(unsigned int states, unsigned int rate_cats, unsigned int states_padded,
                              unsigned int sites, unsigned int float_precision)const{
    
    //unsigned int i,j,k;
    
    // double prob;
    
    //     if ((clv_index < numberTips) &&
    //         (partition->attributes & PLL_ATTRIB_PATTERN_TIP))
    //       return;
    
    //    printf ("[ ");
    //    for (i = 0; i < sites; ++i)
    //    {
    //        printf("{");
    //        for (j = 0; j < rate_cats; ++j)
    //        {
    //            printf("(");
    //            for (k = 0; k < states-1; ++k)
    //            {
    //                prob = clv[i*rate_cats*states_padded + j*states_padded + k];
    //                //if (scale_buffer) PartialTreeNode::unscale(&prob, scale_buffer[i]);
    //                printf("%.*f,", float_precision, prob);
    //            }
    //            prob = clv[i*rate_cats*states_padded + j*states_padded + k];
    //            //if (scale_buffer) PartialTreeNode::unscale(&prob, scale_buffer[i]);
    //            printf("%.*f)", float_precision, prob);
    //            if (j < rate_cats - 1) printf(",");
    //        }
    //        printf("} ");
    //    }
    //    printf ("]\n");
}
void PartialTreeNode::showpClV(unsigned int states, unsigned int rate_cats, unsigned int states_padded,unsigned int sites, unsigned int float_precision)const {
    
    unsigned int i,j,k;
    double prob;
    
    double * pclv2 = pclv;
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
    double * pclv2 = pclv;
  
    
    auto seq = msa->sequence[tip_id];
    //auto charmap = _model.charmap();
    auto charmap = pll_map_gt10;
    
    for (size_t j = 0; j < msa->length; ++j)
    {
        auto charstate = (pll_state_t) seq[j];
        pll_state_t state = charmap ? charmap[(int) charstate] : charstate;
        
        double sum_lh = gtErrorModel->computeStateErrorProbPT20(state, pclv2);
        
        assert(sum_lh >0);
        if (normalize){
            
          
            for (size_t k = 0; k < numberStates; ++k)
                      pclv2[k] *= 1.0* numberStates /sum_lh;
            //with this , \sum_{i=1}^{numberStates} stat_prob * clv[i] = 1
         
            
        }
        
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
        
        pclv2 += numberStates;
        //  clvp += numberStates*numberRateCats;
    }
    
    //restore the pointer to the beginning
    pclv2 -= clv_size;
}
void  PartialTreeNode::unscale(double * prob, unsigned int times)
{
    unsigned int i;
    
    for (i = 0; i < times; ++i)
        *prob *= PLL_SCALE_THRESHOLD;
}
void PartialTreeNode::init(int clv_elements){
    for (size_t i = 0; i < clv_elements; i++) {
        
        //clv.push_back(0.0);
        //scale_buffer.push_back(0);
        
    }
}
void PartialTreeNode::freeVectors(){
    //if (pclv!=nullptr) delete []pclv;
    //if (pscale_buffer!=nullptr) delete []pscale_buffer;
    pclv = nullptr;
    pscale_buffer = nullptr;
}
std::shared_ptr<PartialTreeNode> PartialTreeNode::Clone(){
    
    std::shared_ptr<PartialTreeNode> result;
    
    if (edge_l==nullptr && edge_r==nullptr){
        
        result=std::make_shared<PartialTreeNode>(*this);
        
        
        return result;
    }
    else{
        std::unique_ptr<PartialTreeEdge> leftEdge = edge_l->CloneUnique();
        std::unique_ptr<PartialTreeEdge> rightEdge = edge_r->CloneUnique();
        
        //        result=std::make_shared<PartialTreeNode>(manager, leftEdge, rightEdge, label, height,                                        clv_size,clv_size, 0, pclv, pscale_buffer,                                          index, index_population, ln_likelihood);
        result=std::make_shared<PartialTreeNode>(*this);
        
        return result;
    }
    return result;
}
double PartialTreeNode::likelihood_factor()const{
    
    double result=0.0;
    assert(edge_l && edge_r && "Root cannot be a leaf");
    
    std::shared_ptr<PartialTreeNode> left = edge_l->child;
    std::shared_ptr<PartialTreeNode> right = edge_r->child;
    
    double ln_m = ln_likelihood;
    double ln_l = left->ln_likelihood;
    double ln_r = right->ln_likelihood;
    
    
   // assert(ln_m <= 0 && ln_l <= 0 && ln_r <= 0 &&
   //       "Likelihood can't be more than 100%");
    
    result= ln_m - (ln_l + ln_r);
    
    return result;
    
}
//PartialTreeNode::~PartialTreeNode() {
//    
//    manager->clv_buffer.push(pclv);
//    manager->scale_buffer_buffer.push(pscale_buffer);
//    
//    //pll_aligned_free(clv);
//    //pll_aligned_free(scale_buffer);
//    
//    //if (pclv!=0) delete []pclv;
//    //if (pscale_buffer!=0) delete []pscale_buffer;
//    pclv = nullptr;
//    pscale_buffer = nullptr;
//    // delete clv ;
//    // delete[] scale_buffer ;
//   
//    
//}
