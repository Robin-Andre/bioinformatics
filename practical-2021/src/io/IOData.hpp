#ifndef INFORF_IODATA_H
#define INFORF_IODATA_H

//#include "json/json.hpp"
#include "json.hpp"
#include <string>
#include <vector>
#include <cfloat>
#include <regex>
#include <fstream>
#include <sstream>
#include <iostream>
#include "../Metric.hpp"

/*
 * Struct encapsulating result data and meta data for an RF computation
 */

namespace io {
struct IOData {
	double mean_dst = 0.;
	std::string metric;
	std::vector<std::string> taxa_names;
	// Invariant: numEntries = numTaxa - rowID
	std::vector<std::vector<double>> pairwise_distance_mtx;
	std::string git_revision;
	std::string cpuInformation;
	size_t number_of_unique_trees = 0;
	std::string mode = "ABSOLUTE";

	bool operator==(const IOData &rhs) const {
		// compare floating points manually with relative distance measure

		bool is_eq = comparePairwiseDistances(rhs);
		is_eq &= nearly_eq_floating(mean_dst, rhs.mean_dst);
		is_eq &= metric == rhs.metric;
		is_eq &= taxa_names == rhs.taxa_names;
		is_eq &= git_revision == rhs.git_revision;
		is_eq &= cpuInformation == rhs.cpuInformation;
		is_eq &= number_of_unique_trees == rhs.number_of_unique_trees;
		is_eq &= mode == rhs.mode;
		return is_eq;
	}
	bool operator!=(const IOData &rhs) const {
		return !(rhs == *this);
	}

	bool comparePairwiseDistances(const IOData &other) const {
		bool equal_pairwise = pairwise_distance_mtx.size() == other.pairwise_distance_mtx.size();
		if (!equal_pairwise) {
			return false;
		}
		size_t i = 0;
		for (const auto &el : pairwise_distance_mtx) {
			const auto &othEl = other.pairwise_distance_mtx[i];
			equal_pairwise &= el.size() == othEl.size();
			if (!equal_pairwise) {
				return false;
			}
			for (size_t j = 0; j < el.size(); ++j) {
				equal_pairwise &= nearly_eq_floating(el[j], othEl[j]);
			}
			++i;
		}
		return equal_pairwise;
	}


	std::string toString() const{
		std::stringstream ss;
		ss << "mean_dst: " << mean_dst << std::endl;
		ss << "metric: " << metric << std::endl;
		ss << "number_of_unique_trees: " << number_of_unique_trees << std::endl;
		ss << "git_revision: " << git_revision << std::endl;
		//ss << "git_revision: " << git_revision << std::endl;
		ss << "distances: " << std::endl;
		for(size_t i = 0; i < pairwise_distance_mtx.size(); ++i){
			for(size_t j = 0; j < pairwise_distance_mtx[i].size(); ++j){
				ss << std::to_string(pairwise_distance_mtx[i][j]) << " ";
			}
			ss << std::endl;
		}
		return ss.str();
	}

  private:

	bool nearly_eq_floating(double a, double b) const {
		auto absA = std::abs(a);
		auto absB = std::abs(b);
		auto largest = (absA < absB) ? absB : absA;
		auto smallest = (absA < absB) ? absA : absB;
		return largest - smallest <= largest * static_cast<double>(FLT_EPSILON) * 1e5;
	}


};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(IOData,
                                   mean_dst,
                                   metric,
																	 mode,
                                   taxa_names,
                                   pairwise_distance_mtx,
                                   git_revision,
                                   cpuInformation,
                                   number_of_unique_trees);
} // namespace io
#endif // INFORF_IODATA_H
