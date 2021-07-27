//
//  pll_buffer_manager.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 27/07/2021.
//

#ifndef pll_buffer_manager_hpp
#define pll_buffer_manager_hpp

#include <stack>
struct PLLBufferManager {
  std::stack<double *> clv_buffer;
  std::stack<double *> pmatrix_buffer;
  std::stack<unsigned int *> scale_buffer_buffer;
};
#endif /* pll_buffer_manager_hpp */
