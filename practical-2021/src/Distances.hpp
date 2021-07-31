#pragma once
extern "C" {
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
}
#include "datastructures/PllSplits.hpp"
#include "datastructures/PllTree.hpp"
#include "datastructures/UniquePllMap.hpp"
#include "datastructures/SimilarityCache.hpp"
#include "io/IOData.hpp"
#include "Metric.hpp"
#include <vector>
#include <iostream>
#include <queue>
#include <future>


class Distances {


  /**
   * Provides an interface to calculate simple and generalized n-to-n Robinson-
   * Foulds-Distances for sets of PllTrees
   *
   */
public:

  /**
   * Calculate n-to-n Robinson-Foulds distance as described in https://doi.org/10.1093/bioinformatics/btaa614
   *
   * This sum is the arithmetic sum, not some other kind of sum that only
   * mathematicians have heard of.
   *
   * @param trees The trees for which distances are to be calculated
   * @param metric The desired metric (MSI, SPI, MCI)
   * @param mode Gives whether (not normalized) similarities, absolute or relative  distances
   *    are to be calculated
   * @return IOData struct containing the results of the calculation
   */

  static io::IOData computeGeneralizedDistances(const std::vector<PllTree>& trees,
                                                const GeneralizedMetric& metric,
                                                Mode mode) {


    PllSplit::setTipCount(trees[0].getTipCount());
    size_t tree_count = trees.size();
    assert(tree_count > 0);
    size_t tip_count = trees[0].getTipCount();
    assert(tip_count > 3);

    //precalculations
    phylomath::initCache();
    UniquePllMap _map(trees);
    const std::vector<PllSplitList>& tree_splits = _map.vectors();
    SimilarityCacheMatrix pairwise_cache(_map, metric);

    //Init result
    io::IOData result;
    result.mode = ModeString[mode];
    result.metric = metric.name();
    result.mean_dst = 0;
    result.number_of_unique_trees = tree_count;
    result.pairwise_distance_mtx = std::vector<std::vector<double>>(tree_count, std::vector<double>());

    //Required for (parallel) computation
    size_t dist_count = 0;
    size_t split_count = tree_splits[0].getSplits().size();
    std::vector<std::vector<double>> similarities =
      std::vector<std::vector<double>>(split_count, std::vector<double>(split_count));
    std::vector<std::future<double>> futures;
    bool all_trees_equal_to_last = true;

    for(size_t i = 0; i < tree_count; ++i){
      futures.clear();
      for(size_t j = i; j < tree_count; ++j){
        similaritiesForSplits(tree_splits[i], tree_splits[j], pairwise_cache, &similarities);
        futures.push_back(std::async(MaximumMatcher::match, similarities));
      }
      bool is_unique = true;
      for(size_t j = i; j < tree_count; ++j){
        double dist = futures[j-i].get();
        assert(dist >= 0.0);
        dist = (mode == SIMILARITY) ? dist :
          distanceFromSimilarity(tree_splits[i], tree_splits[j], dist, mode, metric);
        //TODO: Check near 0 because of numerical issues
        if (i != j && dist <= 0 && is_unique){
          is_unique = false;
          --result.number_of_unique_trees;
        }
        if(j == tree_count - 1 && dist > 0){
          all_trees_equal_to_last = false;
        }
        result.pairwise_distance_mtx[j].emplace_back(dist);
        if(i != j){
          result.mean_dst += dist;
          ++dist_count;
        }
      }
    }
    if (all_trees_equal_to_last){
      --result.number_of_unique_trees;
      assert(result.number_of_unique_trees == 0);
    } else {
      assert(result.number_of_unique_trees > 0);
    }
    result.mean_dst = (dist_count == 0) ? 0 : result.mean_dst  / static_cast<double>(dist_count);
    return result;
  }

  /**
   * Calculate n-to-n Robinson-Foulds distance
   *
   * @param trees The trees for which distances are to be calculated
   * @param mode Gives whether absolute or relative distances are to be calculated
   * @return IOData struct containing the results of the calculation
   */


  static io::IOData computeRFDistances(const std::vector<PllTree>& trees, Mode mode) {
    PllSplit::setTipCount(trees[0].getTipCount());
    size_t tree_count = trees.size();
    assert(tree_count > 0);
    size_t tip_count = trees[0].getTipCount();
    assert(tip_count > 3);
    //precalculations
    phylomath::initCache();
    RFMetric metric;
    UniquePllMap map(trees);
    const std::vector<PllSplitList>& tree_splits = map.vectors();
    //Init result
    io::IOData result;
    result.mode = ModeString[mode];
    result.metric = metric.name();
    result.mean_dst = 0;
    result.number_of_unique_trees = tree_count;
    result.pairwise_distance_mtx = std::vector<std::vector<double>>(tree_count, std::vector<double>());


    //Required for computation
    size_t dist_count = 0;
    bool all_trees_equal_to_last = true;
    for(size_t i = 0; i < tree_count; ++i){
      bool is_unique = true;
      for(size_t j = i; j < tree_count; ++j){
        double dist = metric.distanceOf(tree_splits[i], tree_splits[j], mode);
        assert(dist >= 0.0);
        //TODO: Check near 0 because of numerical issues
        if (i != j && dist <= 0.0 && is_unique){
          is_unique = false;
          --result.number_of_unique_trees;
        }
        if(j == tree_count - 1 && dist > 0){
          all_trees_equal_to_last = false;
        }
        result.pairwise_distance_mtx[j].emplace_back(dist);
        if(i != j){
          result.mean_dst += dist;
          ++dist_count;
        }
      }
    }
    if (all_trees_equal_to_last){
      --result.number_of_unique_trees;
      assert(result.number_of_unique_trees == 0);
    } else {
      assert(result.number_of_unique_trees > 0);
    }
    result.mean_dst = (dist_count == 0) ? 0 : result.mean_dst  / static_cast<double>(dist_count);
    return result;
  }

  /**
   * Returns the absolute or relative distance of two trees represented by the provieded
   * split lists based on their similarity
   *
   * @param l1 First split list
   * @param l2 Second split list
   * @param similarity The trees similarity
   * @param mode Gives whether absolute or relative distances are to be calculated
   * @param metric The considered metric
   * @return Distance of the trees represented by the split lists
   */

  static double distanceFromSimilarity( const PllSplitList& l1,
                                        const PllSplitList& l2,
                                        double similarity,
                                        Mode mode,
                                        const GeneralizedMetric& metric) {
    double max_value = metric.maximum(l1, l2);
    assert(max_value > 0);
    assert((max_value - 2 * similarity) > -0.000001);
    double dist = std::max(0.0, max_value - 2 * similarity);
    return (mode == RELATIVE) ? dist : (dist / max_value);
  }

  /**
   * Writes the pairwise similarities for all splits in the provided split lists
   * into the provieded result vector
   *
   * @param l1 First split list
   * @param l2 Second split list
   * @param cache SimilarityCache holding the pairwise values for all splits
   * @param result Similarities are written in that vector
   */

  static void similaritiesForSplits(const PllSplitList& l1, const PllSplitList& l2, const SimilarityCache& cache,
       std::vector<std::vector<double>>* result){
    assert(l1.getSplits().size() == l1.getSplits().size());
    size_t n = l1.getSplits().size();
    assert(result->size() == n);
    for(size_t i = 0; i < n; ++i){
      assert((*result)[i].size() == n);
      for(size_t j = 0; j < n; ++j) {
        (*result)[i][j] = cache.access(l1[i], l2[j]);
      }
    }
  }



};
