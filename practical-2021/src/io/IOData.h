#ifndef INFORF_IODATA_H
#define INFORF_IODATA_H
#include "json/json.hpp"
#include <string>
#include <vector>

namespace io {
struct IOData {
	double mean_rf_dst = 0.;
	enum Metric { MCI, MSI, SPI, RF } split_score_calc;
	double mean_modified_rf_dst = 0.;
	std::vector<std::string> taxa_names;
	// Invariant: numEntries = numTaxa - rowID
	std::vector<std::vector<double>> pairwise_distance_mtx;
	std::string git_revision;
	std::string cpuInformation;
	size_t number_of_unique_trees = 0;

	bool operator==(const IOData &rhs) const;
	bool operator!=(const IOData &rhs) const;

	bool comparePairwiseDistances(const IOData& other) const;

	static bool parse_raxml(const std::string &overview_file_path,
	                        const std::string &distances_path,
	                        IOData &out);

  private:
	// parsing helpers
	static bool handle_pairwise_dst(size_t major,
	                                size_t minor,
	                                std::ifstream &stream,
	                                std::vector<double> &out);
	static bool parse_pairwise_file(const std::string &distances_path,
	                                size_t num_taxa,
	                                std::vector<std::vector<double>> &res);
	static bool nearly_eq_floating(double a, double b);
};
NLOHMANN_JSON_SERIALIZE_ENUM(
    IOData::Metric,
    {{IOData::MCI, "MCI"}, {IOData::MSI, "MSI"}, {IOData::SPI, "SPI"}, {IOData::RF, "RF"}});
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(IOData,
                                   mean_rf_dst,
                                   split_score_calc,
                                   mean_modified_rf_dst,
                                   taxa_names,
                                   pairwise_distance_mtx,
                                   git_revision,
                                   cpuInformation,
                                   number_of_unique_trees);
} // namespace io
#endif // INFORF_IODATA_H
