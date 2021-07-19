#pragma once
#include <vector>
#include <queue>
#include "PllTree.hpp"

/*
 * Builds a unique map of the splits occuring in the set of trees, filtering duplicate splits and 
 * creating PllSplitLists adapted to reference upon the map in accordance to the corresponding trees 
 */
class UniquePllMap {
  public:
  explicit UniquePllMap(const std::vector<PllTree>& trees);
  explicit UniquePllMap(std::vector<PllSplit>& splits);
  ~UniquePllMap() {}

  /** Accesses the split given by an offset 
   * @param index the position of the split in the unique map
   * @return the split at the corresponding position
   */
  const inline PllSplit& operator[] (size_t index) const {
    return all_splits_unique[index];
  }
  /**
   * @return the vector of PllSplitLists corresponding to the trees of the constructor
   */
  const std::vector<PllSplitList>& vectors() const {
      return split_lists;
  }
  /**
   * @return the size of the map / the number of unique splits of a tree set
   */
  size_t size() const {
    return all_splits_unique.size();
  }

  private:
  /*
    This structure encapsulates both a split and the tree from which the split originates.
  */
  struct SplitReference {
    PllSplit split_temp;
    size_t tree_number;
    /*Note that the comparision operator is inverted. This is intentional to ensure that the underlying
      priority queue operates on the smallest elements first instead of the largest elements 
      (see default implementation of PQ https://en.cppreference.com/w/cpp/container/priority_queue)
      */
    bool operator < (SplitReference other) const {
      return !(split_temp < other.split_temp);
    }
  };
  /* The following attributes are designed primarily to help build the unique map. They should not be
     utilized after successful construction anymore. 
  */
  std::priority_queue<SplitReference> queue; //The priority queue 
  std::vector<size_t> tree_iterator; //Saves the number of how many splits of tree number (i) have been processed.
  std::vector<pll_split_t*> trees_as_pointers; //The trees as processed by the external pll-lib
  size_t limit; //How many nontrivial splits are expected per tree : |Taxa| - 3
  size_t map_pos; //The position of the last valid entry into the map

  /* These attributes are the result of the construction 
  */
  std::vector<PllSplit> all_splits_unique; //contains all unique splits found in all trees sorted ascending
  std::vector<PllSplitList> split_lists; /*contains the split lists corresponding to each tree adapted to map
                                           accordingly to the unique map */

  /*
    Attempts to push the next element of tree #ID onto the queue
    returns true if an element could be pushed, false otherwise
  */
  bool push(size_t ID);
  /*
    Empties the queue
  */
  void processQueue();
  /* 
    updates the corresponding PllSplitList and inserts the split into the unique map if it 
    is indeed unique.
  */
  void insertIntoMap(const SplitReference& split_and_id);
  /*
    helper method for updating PllSplitLists
  */
  void updateSplitList(const PllSplit& current_split, size_t ID);
};
