//
// Created by tbwsl on 4/28/21.
//

#include "IOData.h"
#include <cfloat>
#include <fstream>
#include <regex>
#include <string>

namespace io {
bool IOData::handle_pairwise_dst(size_t major,
                                 size_t minor,
                                 std::ifstream &stream,
                                 std::vector<double> &out) {
	std::string regex_str = R"((\d+) (\d+): (\d+) (([0-9]*[.])?[0-9]+))";
	const std::regex pairwise_regex(regex_str);
	std::string curr_line;
	std::smatch match;
	// read line and match if successful
	if (std::getline(stream, curr_line) && std::regex_search(curr_line, match, pairwise_regex)) {
		if (match.size() <= 5) {
			return false;
		}
		if (std::stoul(match[1]) != major || std::stoul(match[2]) != minor) {
			return false;
		}
		out[major] = std::stod(match[4]);
		return true;
	}
	return false;
}

bool IOData::parse_pairwise_file(const std::string &distances_path,
                                 size_t num_taxa,
                                 std::vector<std::vector<double>> &res) {
	assert(res.empty());
	std::ifstream pairwise_file(distances_path);
	if (!pairwise_file.is_open()) {
		return false;
	}
	// preallocate LOWER triangle matrix
	for (size_t row = 0; row < num_taxa; ++row) {
		res.emplace_back(row + 1);
	}
	// assume UPPER triangle matrix input
	for (size_t expect_major = 0; expect_major < num_taxa - 1; ++expect_major) {
		// loop over every entry
		for (size_t expect_minor = expect_major + 1; expect_minor < num_taxa; ++expect_minor) {
			if (!handle_pairwise_dst(
			        expect_major, expect_minor, pairwise_file, res[expect_minor])) {
				return false;
			}
		}
	}
	return true;
}

bool IOData::operator==(const IOData &rhs) const {
	// compare floating points manually with relative distance measure

	bool is_eq = comparePairwiseDistances(rhs);
	is_eq |= nearly_eq_floating(mean_rf_dst, rhs.mean_rf_dst);
	is_eq |= split_score_calc == rhs.split_score_calc;
	is_eq |= nearly_eq_floating(mean_modified_rf_dst, rhs.mean_modified_rf_dst);
	is_eq |= taxa_names == rhs.taxa_names;
	is_eq |= git_revision == rhs.git_revision;
	is_eq |= cpuInformation == rhs.cpuInformation;
	return is_eq;
}
bool IOData::operator!=(const IOData &rhs) const {
	return !(rhs == *this);
}

bool IOData::parse_raxml(const std::string &overview_file_path,
                         const std::string &distances_path,
                         IOData &out) {
	// parse 2 values from overview file
	size_t num_trees = 0;
	{
		std::ifstream overview(overview_file_path);
		if (!overview.is_open()) {
			return false;
		}
		// read entire file to string
		std::string file_content((std::istreambuf_iterator<char>(overview)),
		                         (std::istreambuf_iterator<char>()));
		overview.close();

		const std::regex unique_t_patt("Number of unique trees in this tree set: (\\d+)");
		const std::regex num_t_patt("Found (\\d+) trees in File");
		const std::regex mean_rf_dst("Average relative RF in this set: (([0-9]*[.])?[0-9]+)");
		std::smatch unique_t_mat;
		// check if found num unique trees;
		if (!std::regex_search(file_content, unique_t_mat, unique_t_patt)) {
			return false;
		}
		out.number_of_unique_trees = std::stoul(unique_t_mat[1].str());

		std::smatch mean_dst_match;
		if (!std::regex_search(file_content, mean_dst_match, mean_rf_dst)) {
			return false;
		}
		out.mean_rf_dst = stod(mean_dst_match[1].str());

		std::smatch num_t_mat;
		// check if found num unique trees;
		if (!std::regex_search(file_content, num_t_mat, num_t_patt)) {
			return false;
		}
		num_trees = std::stoul(num_t_mat[1].str());
	}
	// parse pairwise distances
	bool success = parse_pairwise_file(distances_path, num_trees, out.pairwise_distance_mtx);

	return success;
}
bool IOData::comparePairwiseDistances(const IOData &other) const {
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
bool IOData::nearly_eq_floating(double a, double b) {
	auto absA = std::abs(a);
	auto absB = std::abs(b);
	auto largest = (absA < absB) ? absB : absA;
	auto smallest = (absA < absB) ? absA : absB;
	return largest - smallest <= largest * static_cast<double>(FLT_EPSILON) * 1e5;
}

} // namespace io