#pragma once
#include <vector>
#include "PllPointerMap.hpp"
class IntersectionCache {
    public:
        struct IntersectionTableEntry {
        size_t cut_AA;
        size_t cut_AB;
        size_t cut_BA;
        size_t cut_BB;
        
    };
    IntersectionCache(const PllPointerMap& map) {
      std::cout << "Map Size: " << map.size() << "\n";
      cache = std::vector(map.size(), std::vector<IntersectionTableEntry>(map.size()));
      for(unsigned i = 0; i < cache.size(); ++i) {
          //Theoretically we could make this a triangular matrix and make the access well formed i.e. A < B
          for(unsigned j = 0; j < cache[i].size(); ++j) {
            cache[i][j] = calculateIntersectionSize(map, i, j);
          }
      }
    }
    IntersectionTableEntry access(size_t i, size_t j) const {
        return cache[i][j];
    }

    private:

    IntersectionTableEntry calculateIntersectionSize(const PllPointerMap& map, size_t i, size_t j) {
        size_t intersect_AA = map[i].intersectionSize(map[j], Block_A, Block_A);
        size_t intersect_AB = map[i].intersectionSize(map[j], Block_A, Block_B);
        size_t intersect_BA = map[i].intersectionSize(map[j], Block_B, Block_A);
        size_t intersect_BB = map[i].intersectionSize(map[j], Block_B, Block_B);
        return IntersectionTableEntry({intersect_AA, intersect_AB, intersect_BA, intersect_BB});
    }
 
    std::vector<std::vector<IntersectionTableEntry>> cache;
};