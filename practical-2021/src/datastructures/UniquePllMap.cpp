#include "UniquePllMap.hpp"

  UniquePllMap::UniquePllMap(const std::vector<PllTree> &trees) : tree_iterator(std::vector<size_t>(trees.size())) {
    trees_as_pointers = std::vector<pll_split_t*>(trees.size());
    limit = trees[0].getTipCount() - 3;
    split_lists = std::vector<PllSplitList>(trees.size());
    // IN the worst case every split is unique so limit * trees
    all_splits_unique = std::vector<PllSplit>(trees.size() * limit);
    map_pos = 0;

    for(unsigned i = 0; i < trees.size(); ++i) {
      trees_as_pointers[i] = pllmod_utree_split_create(trees[i].tree()->vroot,
        static_cast<unsigned int>(trees[0].getTipCount()), nullptr);
      push(i); //Push the first element of every tree into the queue
    }

    processQueue();
    //If duplicates exist the map won't be filled completely, therefore it is resized
    all_splits_unique.resize(map_pos);
  }
  //This constructor is required for testing purposes only
  UniquePllMap::UniquePllMap(std::vector<PllSplit>& splits) : all_splits_unique(splits) {
    map_pos = 0;
    limit = 0;
  }
  bool UniquePllMap::push(size_t ID) {
    //only push elements from a tree if there are any left
    if(tree_iterator[ID] < limit) {
      queue.push({PllSplit(trees_as_pointers[ID][tree_iterator[ID]]), static_cast<size_t>(ID)});
      ++tree_iterator[ID];
      return true;
    }
    return false;
  }
  void UniquePllMap::processQueue() {
    while(!queue.empty()) {
      SplitReference top = queue.top();
      insertIntoMap(top);
      push(top.tree_number);
      queue.pop();
    }
  }
  void UniquePllMap::insertIntoMap(const SplitReference& split_and_id) {
    assert (map_pos <= all_splits_unique.size());
    //if the split is unique or the first we must insert it into the map and increment the counter
    if(map_pos == 0 || !(all_splits_unique[map_pos - 1] == split_and_id.split_temp)) {
      all_splits_unique[map_pos] = std::move(split_and_id.split_temp);
      ++map_pos;
    }
    //regardless of whether it was inserted or not the PllSplitList must be updated
    //and will definitively reference the last element of the map. (either because it was inserted
    //or because it was a duplicate and the counter has not been increased)
    updateSplitList(all_splits_unique[map_pos - 1], split_and_id.tree_number);
  }
  void UniquePllMap::updateSplitList(const PllSplit& current_split, size_t ID) {
    assert(ID < tree_iterator.size());
    //Since the push method always increases tree_iterator and points to the next element to be inserted into queue
    //We have to subtract 1 to find the location as where to insert the element.
    split_lists[ID].push(current_split, map_pos - 1);
  }
