//
//  parsers.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 26/08/2021.
//

#include "parsers.hpp"


using namespace std;

pll_rtree_t *Parsers::readRooted(const std::string &newick, bool isFile)
{
  pll_rtree_t* tree;
  if (isFile)
     tree = pll_rtree_parse_newick(newick.c_str());
  else {
     tree =  pll_rtree_parse_newick_string(newick.c_str());
      
  }
  if (!tree) {
    std::string errorMessage;
    if (isFile) {
      errorMessage = "Error while reading rooted tree from file " + newick + ".\n";
    } else {
      errorMessage = "Error while reading rooted tree from string " + newick + ".\n";
    }
      std::cout << errorMessage << std::endl;
  }
    return tree;
}
