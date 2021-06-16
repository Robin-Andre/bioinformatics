#pragma once
#include <mpfr.h>
//TODO maybe rename to A, B, also A is equal to everything is 1 so they are
//additionally ordered B -> A to that the enum number aligns with the
//numeric value of the corrsponding bits in the partition of the split
enum Partition{Block_B, Block_A};

static const size_t PRECISION = 500;
static const mpfr_rnd_t RND = MPFR_RNDN;
static const double UNIQUE_EPSILON = 0.00001;
