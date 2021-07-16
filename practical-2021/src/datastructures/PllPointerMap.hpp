#pragma once
#include <vector>
#include <queue>
#include "PllTree.hpp"

/*
 * TODO: @Robin: Think about renaming and add docs
 */
class PllPointerMap {
  public:
  explicit PllPointerMap(const std::vector<PllTree> &trees) : tree_iterator(std::vector<size_t>(trees.size())) {
    trees_as_pointers = std::vector<pll_split_t*>(trees.size());
    limit = trees[0].getTipCount() - 3;
    split_lists = std::vector<PllSplitList>(trees.size());
    all_splits_unique = std::vector<PllSplit>(trees.size() * limit); // IN the worst case every split is unique so limit * trees
    map_pos = 0;
    //std::cout << "Size: " << all_splits_unique.size() <<" Limit: " << limit << "\n";
    for(unsigned i = 0; i < trees.size(); ++i) {
      trees_as_pointers[i] = pllmod_utree_split_create(trees[i].tree()->vroot, static_cast<unsigned int>(trees[0].getTipCount()), nullptr);
      push(i); //Push the first element of every tree into the queue
    }
    processQueue();
    all_splits_unique.resize(map_pos);
  }
  //This is for testing purposes only
  explicit PllPointerMap(std::vector<PllSplit>& splits) : all_splits_unique(splits) {
    //dummy initializations for Softwipe
    map_pos = 0;
    limit = 0;
  }
  ~PllPointerMap() {}
  const PllSplit& operator[] (size_t index) const {
    return all_splits_unique[index];
  }
  std::vector<PllSplitList>& vectors() {
      return split_lists;
  }
  size_t size() const {
    return all_splits_unique.size();
  }
  private:
  struct SplitReference {
    PllSplit split_temp;
    size_t tree_number;
    bool operator < (SplitReference other) const {
      return !(split_temp < other.split_temp);
    }
  };
  std::priority_queue<SplitReference> queue;
  std::vector<size_t> tree_iterator;
  std::vector<pll_split_t*> trees_as_pointers;
  std::vector<PllSplit> all_splits_unique;
  std::vector<PllSplitList> split_lists;
  size_t limit;
  size_t map_pos;
  /*
    Attempts to push the next element of tree #ID onto the queue
    returns true if an element could be pushed, false otherwise
  */
  bool push(size_t ID) {
    //Limit is hardcoded in the constructor to tipcount - 3
    if(tree_iterator[ID] < limit) {
      queue.push({PllSplit(trees_as_pointers[ID][tree_iterator[ID]]), static_cast<size_t>(ID)});
      ++tree_iterator[ID];
      return true;
    }
    return false;
  }
  void processQueue() {

    while(!queue.empty()) {
      SplitReference top = queue.top();
      insertIntoMap(top);
      push(top.tree_number);
      queue.pop();
    }
  }
  void insertIntoMap(const SplitReference& split_and_id) {
    assert (map_pos <= all_splits_unique.size());
    if(map_pos > 0 && all_splits_unique[map_pos - 1] == split_and_id.split_temp) {
      //TODO rewrite since nothing is executed here
      //This is extremely weird but this if block represents when a duplicate split has been found
      //since we do not
    }
    else {
      all_splits_unique[map_pos] = std::move(split_and_id.split_temp);
      ++map_pos;
    }
    updateSplitList(all_splits_unique[map_pos - 1], split_and_id.tree_number);
  }
  void updateSplitList(const PllSplit& current_split, size_t ID) {
    assert(ID < tree_iterator.size());
    //Since the push method always increases tree_iterator and points to the next element to be inserted into queue
    //We have to subtract 1 to find the location as where to insert the element.
    split_lists[ID].push(current_split, map_pos - 1);
  }
  /*void verify(unsigned ID) {
    for(unsigned i = 0; i < split_lists[ID].getSplitCount(); ++i) {
      std::cout << split_lists[ID][i]->toString();
    }
  }*/

};
