//
//  msa.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 26/08/2021.
//

#include <stdexcept>
#include "msa.hpp"
#include "pll_utils.hpp"

//this function is copied from RAXML
bool MSA::checkMSA(const MSA& msa)
{
  bool msa_valid = true;

  /* check taxa count */
  if (msa.getSize() < 4)
  {
    std::cout << "\nERROR: Your alignment contains less than 4 sequences! " << std::endl;
    msa_valid = false;
  }

  /* check for duplicate taxon names */
  unsigned long stats_mask = PLLMOD_MSA_STATS_DUP_TAXA;

  pllmod_msa_stats_t * stats = pllmod_msa_compute_stats(msa.pll_msa,
                                                        4,
                                                        pll_map_nt, // map is not used here
                                                        NULL,
                                                        stats_mask);


  assert(stats);

  if (stats->dup_taxa_pairs_count > 0)
  {

    for (unsigned long c = 0; c < stats->dup_taxa_pairs_count; ++c)
    {
      auto id1 = stats->dup_taxa_pairs[c*2];
      auto id2 = stats->dup_taxa_pairs[c*2+1];
      std::cout << "ERROR: Sequences " << id1+1 << " and "
                << id2+1 << " have identical name: "
                << msa.label(id1) <<  std::endl;
    }
    std::cout << "\nERROR: Duplicate sequence names found: "
              << stats->dup_taxa_pairs_count <<  std::endl;

    msa_valid = false;
  }

  pllmod_msa_destroy_stats(stats);

  return msa_valid;
}


MSA::MSA(const pll_msa_t *pll_msa) :
    length(0), num_sites(pll_msa->length), states(0), pll_msa(nullptr)
{
  for (auto i = 0; i < pll_msa->count; ++i)
  {
    append(string(pll_msa->sequence[i], pll_msa->length), pll_msa->label ? pll_msa->label[i] : "");
  }

  update_pll_msa();
}

MSA::MSA(MSA&& other) : length(other.length), num_sites(other.num_sites),
    sequences(move(other.sequences)), labels(move(other.labels)),
    label_id_map(move(other.label_id_map)), weights(move(other.weights)),
    probs(move(other.probs)),
    states(other.states), pll_msa(other.pll_msa), dirty(other.dirty)
{
  other.length = other.num_sites = 0;
  other.pll_msa = nullptr;
  other.dirty = false;
};

MSA::~MSA()
{
  free_pll_msa();
}

MSA& MSA::operator=(MSA&& other)
{
  if (this != &other)
  {
    // release the current object’s resources
    free_pll_msa();
    weights.clear();
    sequences.clear();
    labels.clear();
    label_id_map.clear();

    // steal other’s resource
    length = other.length;
    num_sites = other.num_sites;
    pll_msa = other.pll_msa;
    weights = std::move(other.weights);
    sequences = std::move(other.sequences);
    labels = std::move(other.labels);
    label_id_map = std::move(other.label_id_map);
    probs = std::move(other.probs);
    //local_seq_ranges = std::move(other.local_seq_ranges);
    states = other.states;
    dirty = other.dirty;

    // reset other
    other.length = other.num_sites = other.states = 0;
    other.pll_msa = nullptr;
    other.dirty = false;
  }
  return *this;
}
void MSA::append(const string& sequence, const string& header)
{
  if(length && sequence.length() != (size_t) length)
    throw runtime_error{string("Tried to insert sequence to MSA of unequal length: ") + sequence};

  sequences.push_back(sequence);

  if (!header.empty())
  {
    labels.push_back(header);
    label_id_map[header] = labels.size() - 1;
  }

  if (!length)
  {
    length = sequence.length();
    if (!num_sites)
      num_sites = length;
  }

  dirty = true;
}

void MSA::free_pll_msa() noexcept
{
  if (pll_msa)
  {
    free(pll_msa->sequence);
    if (pll_msa->label)
      free(pll_msa->label);
    free(pll_msa);
    pll_msa = nullptr;
  }
}
void MSA::update_pll_msa() const
{
  if (!pll_msa)
  {
    pll_msa = (pll_msa_t *) calloc(1, sizeof(pll_msa_t));
    dirty = true;
  }

  assert(labels.empty() || labels.size() == sequences.size());

  if (dirty)
  {
    pll_msa->count = getSize();
    pll_msa->length = getLength();

    size_t i = 0;
    pll_msa->sequence = (char **) calloc(pll_msa->count, sizeof(char *));
    for (const auto& entry : sequences)
    {
      pll_msa->sequence[i] = (char *) entry.c_str();
      ++i;
    }
    assert(i == getSize());

    if (!labels.empty())
    {
      i = 0;
      pll_msa->label = (char **) calloc(pll_msa->count, sizeof(char *));
      for (const auto& entry : labels)
      {
        pll_msa->label[i] = (char *) entry.c_str();
        ++i;
      }
      assert(i == getSize());
    }

    dirty = false;
  }
}
