#include "RFDistance2.hpp"
#include "io/FileReader.hpp"
#include <algorithm>
namespace measure {
/* A quick and dirty implementation of the standard RF distance as defined in the lecture. 
   right now there is no persistent element that needs to be stored as an attribute so it is not
   implemented as a class. This might change. */
size_t rf_distance(const PllSplitList& p1, const PllSplitList& p2) {
    //assert(p1.sorted()) //TODO Enable proper assertion once sorting is removed
    //assert(p2.sorted()) 
    /* Also I have no clue how to properly assert a sorted input :sweat_smile: Uncommenting the above lines
    will result in noncompilation (We should treat it as pseudocode for now :P) */

    //p1.print();

    /* It might seem really silly to copy the Pll-vectors and I would wholeheartedly agree.
    but as long as the sorting happens in this method I felt unsure to remove this. Perhaps
    const modifiers might be disregarded or other unspecified behaviour might happen. 
    Maybe nothing would happen at all. But it is late at night and I cannot be bothered to test
    what happens if I sort the reference instead of the copy */ 
    std::vector<PllSplit> set_A = p1.getSplits();
    std::vector<PllSplit> set_B = p2.getSplits();
    /* The vector allocation will be redundant when we implement a simple counter 
    and be relevant again once we are interested in the actual elements of union and intersection.
    this might also be a good starting point for hashmaps etc. */
    std::vector<PllSplit> unionA_B;
    std::vector<PllSplit> intersectA_B;

    /*TODO it makes no sense to sort the splitlists every ...single...time.
    we should sort them beforehand and simply expect that behaviour to stick */ 
    //TODO spontaneous disabling of sorting
    //std::sort(set_A.begin(), set_A.end());
    //std::sort(set_A.begin(), set_A.end());
    /*TODO need to implement a simple counter instead of collecting the entire vectors, we are not 
    interested in the intersection/union elements ...yet.  Maybe even do both the intersection
    and union in one sweep instead of two */
    std::set_union(set_A.begin(), set_A.end(), set_B.begin(), set_B.end(), std::back_inserter(unionA_B));
    std::set_intersection(set_A.begin(), set_A.end(), set_B.begin(), set_B.end(), std::back_inserter(intersectA_B));
    return unionA_B.size() - intersectA_B.size();
}
/* TODO this should never be called and maybe even removed, it makes absolutely no sense to recalculate
   the splits every ...single...time ... but then again lazy me wants this done so I can learn*/
size_t rf_distance(const PllTree& t1, const PllTree& t2) {
    return rf_distance(t1.makeSplits(), t2.makeSplits());
}


/*Rudimentary implementation, a onedimensional vector, no structure, maybe a struct might be more suited to approach Alexis
  output format. But then again a problem for another time*/
std::vector<size_t> full_calculation(const std::string& file) {
    std::vector<PllTree> tree_vector = io::readTreeFile(file);
    std::vector<size_t> result; // TODO This should be preallocated to the exact amount of expected results
    
    
    
    //Align the trees, really weird stuff happens when they are misaligned

    for(unsigned i = 1/* This i=1 is intentional */; i < tree_vector.size(); ++i){
        tree_vector[i].alignNodeIndices(tree_vector[0]);
    }


    //Atleast this is semi-straightforward now. Though I am not sure about the indices, needs testing
    for(unsigned i = 0; i < tree_vector.size() - 1 /*The -1 "should" not cause issues*/; ++i) {
        for(unsigned j = i + 1; j < tree_vector.size(); ++j) {
            std::cout << "Running rf_distance on: " << i << " , " << j;
            size_t temp_rf_dist = measure::rf_distance(tree_vector[i], tree_vector[j]);
            std::cout << " result: " << temp_rf_dist <<"\n";
            result.push_back(temp_rf_dist);
        }
        
    }
    return result;
}
} //namespace measure