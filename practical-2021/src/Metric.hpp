#pragma once
#include "PllSplits.hpp"
class Metrics {
    public: 
    virtual double evaluate(const PllSplit& s1, const PllSplit& s2) const = 0;
    virtual double maximum(const PllSplitList& plist1, const PllSplitList& plist2) const = 0;
};