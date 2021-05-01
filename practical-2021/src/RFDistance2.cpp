#include "RFDistance2.hpp"
#include <algorithm>
namespace measure {
size_t rf_distance(const PllSplitList& p1, const PllSplitList& p2) {
    //assert(p1.sorted()) //TODO Enable proper assertion
    //assert(p2.sorted()) //TODO
    //vector union <- set::union
    //vector intersect <- set::intersect
    std::vector<PllSplit> set_A = p1.getSplits();
    std::vector<PllSplit> set_B = p2.getSplits();
    std::vector<PllSplit> unionA_B;
    std::vector<PllSplit> intersectA_B;
    std::sort(set_A.begin(), set_A.end());
    std::sort(set_A.begin(), set_A.end());
    std::set_union(set_A.begin(), set_A.end(), set_B.begin(), set_B.end(), std::back_inserter(unionA_B));
    std::set_intersection(set_A.begin(), set_A.end(), set_B.begin(), set_B.end(), std::back_inserter(intersectA_B));
    return unionA_B.size() - intersectA_B.size();
}
int lol() {
    return 1;
}
} //namespace measure