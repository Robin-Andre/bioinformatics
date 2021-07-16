#pragma once
#include "PllSplits.hpp"
extern "C" {
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
}
#include <exception>
#include <string>
#include <utility>

/*
 *  Class representing a phylogenetic tree
 */
class PllTree {
public:
  /*
   * Creates a PllTree by parsing an input string in newick format
   * @param newick_string: string encoding the tree in newick format
   */
  explicit PllTree(const std::string &newick_string);

  /*
   * Creates a PllTree by parsing an input string in newick format
   * and aligns its node indices with a preexisting tree
   * @param newick_string: string encoding the tree in newick format
   * @param alignment_tree: the tree to align node indices with
   */
  PllTree(const std::string &newick_string, const PllTree& alignment_tree);
  PllTree()=delete;
  /* Rule of 5 constructors/destructors */
  ~PllTree();
  PllTree(const PllTree &other);
  PllTree(PllTree &&other) : _tree{std::exchange(other._tree, nullptr)} {}
  PllTree &operator=(const PllTree &other) { return *this = PllTree(other); }
  PllTree &operator=(PllTree &&other) {
    std::swap(_tree, other._tree);
    return *this;
  }

  /* Getters */

  const pll_utree_t *tree() const { return _tree; }

  /*
   * Creates the PllSplitList corresponding to the PllTree
   *
   * @return The corresponding PllSplitList
   */
  PllSplitList makeSplits() const;

  /*
   * Aligns this PllTree's node indices with those of the provided tree
   *
   * @param other: The tree to align node indices to
   */
  void         alignNodeIndices(const PllTree &other);

  /*
   *
   * @return The tip count, i.e. the number of taxa in the tree
   */
  size_t getTipCount() const {return _tree->tip_count;}

private:
  pll_utree_t *_tree;
};
