//
//  RF_distance_calculator.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 25/08/2021.
//

#include "RF_distance_calculator.hpp"

#include <numeric>
#include <iostream>
#include <stdexcept>

using namespace std;
RFDistanceCalculator::RFDistanceCalculator (const std::vector<RootedTree>& trees, bool lowmem) :
    num_trees(0), num_tips(0), avg_rf(0.0), avg_rrf(0.0), num_uniq_trees(0)
{
  if (trees.size() > 1)
  {
    num_trees = trees.size();
    num_tips = trees.at(0).numTips();
    if (lowmem)
      calcRFdistanceLowmem(trees);
    else
      calcRFdistance(trees);
  }
  else
    throw std::invalid_argument("Need at least 2 trees to compute RF distances! "
                        "Given: " + to_string(trees.size()));
}

RFDistanceCalculator::~RFDistanceCalculator ()
{
}

void RFDistanceCalculator::addTreeSplits(size_t tree_idx, const RootedTree& tree,
                                       bitv_hashtable_t * splits_hash)
{
  assert(splits_hash);
  assert(tree.numTips() == num_tips);
  assert(tree_idx < num_trees);
    

  pll_utree_t * utree = (&tree)->getUnrootedPtr();
 
  //pll_unode_t * vroot =  utree->nodes[tree.numTips()];
  pll_unode_t * vroot = utree->nodes[utree->tip_count+utree->inner_count-1];
  pll_split_t * splits = pllmod_utree_split_create(vroot,
                                                 utree->tip_count,
                                                 nullptr);
  int num_splits = utree->edge_count - utree->tip_count;
  for (size_t i = 0; i < num_splits; ++i)
  {
    bitv_hash_entry_t * e = pllmod_utree_split_hashtable_insert_single(splits_hash,
                                                                       splits[i],
                                                                       1.0);
    if (!e)
      std::cout<< "Cannot add a split into hashtable: " << std::endl;

    assert(e->bip_number <= split_occurence.size());

    /* new split -> create bit vector with tree occurrence flags */
    if (e->bip_number ==split_occurence.size())
      split_occurence.emplace_back(std::vector<bool>(num_trees));

    split_occurence[e->bip_number][tree_idx] = true;
  }

  pllmod_utree_split_destroy(splits);
}

void RFDistanceCalculator::calcRFdistance(const std::vector<RootedTree>& trees)
{
  bitv_hashtable_t * splits_hash = pllmod_utree_split_hashtable_create(num_tips, 0);

  if (!splits_hash)
  {
    assert(pll_errno);
    std::cout<< "Cannot create split hashtable"<< std::endl;
  }

  /* add splits from all trees into a hashtable */
  for (size_t i = 0; i < num_trees; ++i)
  {
    addTreeSplits(i, trees.at(i), splits_hash);
  }

  assert(split_occurence.size() == splits_hash->entry_count);

  /* now compute all pairwise RF distances using hashtable */
  rfdist_mat.resize(num_trees * num_trees, 0.);

  /* iterate over all splits in the hashtable */
  for (size_t i = 0; i < splits_hash->table_size; ++i)
  {
    bitv_hash_entry_t * e =  splits_hash->table[i];
    while (e != NULL)
    {
      assert(e->bip_number < splits_hash->entry_count);

      const auto& occ = split_occurence[e->bip_number];
      for(size_t j = 0; j < num_trees; j++)
      {
        if (occ[j])
        {
          for(size_t k = 0; k < num_trees; k++)
          {
            if (j == k)
              continue;

            if (!occ[k])
            {
              rfdist_mat[j * num_trees + k]++;
            }
          }
        }
      }
      e = e->next;
    }
  }

  pllmod_utree_split_hashtable_destroy(splits_hash);

  /* compute averages */
  avg_rf = std::accumulate(rfdist_mat.begin(), rfdist_mat.end(), 0);
  num_uniq_trees = 1;
  for(size_t j = 0; j < num_trees-1; j++)
  {
    bool uniq = true;
    for(size_t k = j+1; k < num_trees; k++)
    {
      auto rf = rfdist_mat[j * num_trees + k] + rfdist_mat[k * num_trees + j];
      rfdist_mat[j * num_trees + k] = rfdist_mat[k * num_trees + j] = rf;
      uniq &= (rf > 0);
    }

    if (uniq)
      num_uniq_trees++;
  }

    auto num_pairs = num_trees * (num_trees - 1) / 2;
  avg_rf /= num_pairs;
  avg_rrf = avg_rf / maxRF();
}

void RFDistanceCalculator::calcRFdistanceLowmem(const std::vector<RootedTree>& trees)
{
  assert(trees.size() > 0);

  /* extract splits from all trees */
  std::vector<pll_split_t *> splits(num_trees);
  for (size_t i = 0; i < num_trees; ++i)
  {
    const auto& tree = trees.at(i);
      pll_utree_t * utree = (&tree)->getUnrootedPtr();
      //pll_unode_t * vroot =  utree->nodes[tree.numTips()];
      pll_unode_t * vroot = utree->nodes[utree->tip_count+utree->inner_count-1];
      splits[i] = pllmod_utree_split_create(vroot,
                                          utree->tip_count,
                                          nullptr);
  }

  avg_rf = 0.0;
  avg_rrf = 0.0;
  num_uniq_trees = 1;
  size_t num_pairs = 0;

  for (size_t i = 0; i < num_trees-1; ++i)
  {
    bool uniq = true;
    for (size_t j = i+1; j < num_trees; ++j)
    {
      auto rf = pllmod_utree_split_rf_distance(splits[i], splits[j], num_tips);
//      if (rf >0)
//          std::cout << "RF dist between tree " <<i << " and tree " << j << " " << rf << std::endl;
      avg_rf += rf;

      // TODO: maxrf will be different for multifurcating trees
      double rrf = ((double) rf) / maxRF();
      avg_rrf += rrf;

      uniq &= (rf > 0);
      num_pairs++;
    }

    if (uniq)
      num_uniq_trees++;
  }

  for (auto s: splits)
    pllmod_utree_split_destroy(s);

  avg_rf /= num_pairs;
  avg_rrf /= num_pairs;
}

double RFDistanceCalculator::avgRF() const
{
  return avg_rf;
}

double RFDistanceCalculator::avgRRF() const
{
  return avg_rrf;
}

size_t RFDistanceCalculator::numUniqueTrees() const
{
  return num_uniq_trees;
}

double RFDistanceCalculator::maxRF() const
{
  return (double) 2 * (num_tips - 3);
}

size_t RFDistanceCalculator::numTrees() const
{
  return num_trees;
}

size_t RFDistanceCalculator::numTips() const
{
  return num_tips;
}

size_t RFDistanceCalculator::RF(size_t i, size_t j) const
{
  return rfdist_mat.empty() ? 0 :rfdist_mat[i * num_trees + j];
}

double RFDistanceCalculator::RRF(size_t i, size_t j) const
{
 // TODO: maxrf will be different for multifurcating trees
 return RF(i, j) / maxRF();
}
