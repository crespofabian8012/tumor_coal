//
//  RF_distance_calculator.hpp
//  run
//
//  Form RaxML Robinson Foulds Calculator
//

#ifndef RF_distance_calculator_hpp
#define RF_distance_calculator_hpp

#include "tree.hpp"
#include <vector>
extern "C"
    {
    #include "libpll/pll.h"
    #include "libpll/pll_optimize.h"
    #include "libpll/pllmod_common.h"
    #include <libpll/pll_msa.h> //for pllmod_msa_empirical_frequencies
    #include <libpll/pll_tree.h>
    #include <libpll/pllmod_util.h>
    #include "libpll/pll_optimize.h"
    #include <libpll/pllmod_common.h>
    #include <libpll/pllmod_algorithm.h>
}
using namespace std;

class RFDistanceCalculator
{
public:
  RFDistanceCalculator(const std::vector<RootedTree>& trees, bool lowmem = false);
  virtual
  ~RFDistanceCalculator ();

public:
  size_t numTrees() const;
  size_t numTips() const;
  size_t RF(size_t i, size_t j) const;
  double RRF(size_t i, size_t j) const;

  double avgRF() const;
  double avgRRF() const;
  size_t numUniqueTrees() const;

private:
  size_t num_trees;
  size_t num_tips;
  std::vector<unsigned int> rfdist_mat;
  std::vector<std::vector<bool>> split_occurence;

  double avg_rf;
  double avg_rrf;
  size_t num_uniq_trees;

  void calcRFdistance(const std::vector<RootedTree>& trees);
  void calcRFdistanceLowmem(const std::vector<RootedTree>& trees);
  void addTreeSplits(size_t tree_idx, const RootedTree& tree,
                       bitv_hashtable_t * splits_hash);
  double maxRF() const;
};


#endif /* RF_distance_calculator_hpp */
