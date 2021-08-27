//
//  parsers.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 26/08/2021.
//

#ifndef parsers_hpp
#define parsers_hpp

#include "pll_utils.hpp"
#include <string>

namespace Parsers {
 pll_rtree_t *readRooted(const std::string &newick, bool isFile);
}

#endif /* parsers_hpp */
