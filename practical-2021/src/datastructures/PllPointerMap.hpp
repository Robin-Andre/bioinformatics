#pragma once
#include <vector>
#include <queue>
#include "PllTree.hpp"
class PllPointerMap {
  public:
  explicit PllPointerMap(const std::vector<PllTree> &trees) {
    treeIter = std::vector<size_t>(trees.size());
    trees_as_pointers = std::vector<pll_split_t*>(trees.size());
    limit = trees[0].getTipCount() - 3;
    //TODO this preallocation needs trimming it reads like shit
    lol = std::vector<PllSplitList>(trees.size());
    /*for(unsigned i = 0; i < trees.size(); ++i) {
        lol[i] = PllSplitList(limit);
    }*/
    all_splits_unique = std::vector<PllSplit>(trees.size() * limit); // IN the worst case every split is unique so limit * trees
    map_pos = 0;
    //std::cout << "Size: " << all_splits_unique.size() <<" Limit: " << limit << "\n";
    for(unsigned i = 0; i < trees.size(); ++i) {
      trees_as_pointers[i] = pllmod_utree_split_create(trees[i].tree()->vroot, trees[0].getTipCount(), nullptr);
      push(i); //Push the first element of every tree into the queue
    }
    processQueue();
    //verify(0);
    //std::cout << "---------------------\n";
    //verify(1);
  }
  ~PllPointerMap() {
    for(unsigned i = 0; i < trees_as_pointers.size(); ++i) {
      free(trees_as_pointers[i]);
    }
  }
  PllSplit operator[] (size_t index) const {
    std::cout << "Accessing at: " <<index << " maxsize: " << all_splits_unique.size() << "\n";
    return all_splits_unique[index];
  }
  std::vector<PllSplit> getMap() {
      return all_splits_unique;
  }
  std::vector<PllSplitList>& vectors() {
      return lol;
  }
  private: 
  struct SplitReference {
    PllSplit split_temp;
    int tree_number;
    bool operator < (SplitReference other) const {
      return !(split_temp < other.split_temp);
    }
  };
  std::priority_queue<SplitReference> queue;
  std::vector<size_t> treeIter;
  std::vector<pll_split_t*> trees_as_pointers;
  std::vector<PllSplit> all_splits_unique;
  std::vector<PllSplitList> lol;
  size_t limit;
  size_t map_pos;
  /*
    Attempts to push the next element of tree #ID onto the queue
    returns true if an element could be pushed, false otherwise
  */
  bool push(size_t ID) {
    //Limit is hardcoded in the constructor to tipcount - 3 
    if(treeIter[ID] < limit) {
      //std::cout << "Inserting: " << treeIter[ID] <<"\n";
      queue.push({PllSplit(trees_as_pointers[ID][treeIter[ID]]), ID});
      ++treeIter[ID];
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
  void updateSplitList(PllSplit& pointer, size_t ID) {
    //Since the push method always increases treeIter and points to the next element to be inserted into queue
    //We have to subtract 1 to find the location as where to insert the element. 
    auto oldID = treeIter[ID] - 1; 
    lol[ID].push(&pointer, map_pos - 1);    
  }
  void verify(unsigned ID) {
    for(unsigned i = 0; i < lol[ID].getSplitCount(); ++i) {
      //std::cout << lol[ID][i]->toString();
    }
  }

};