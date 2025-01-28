#ifndef APPROXIMATE
#define APPROXIMATE
#include <vector>
#include <tuple>
#include "basepair.hpp"
#include "matrix.hpp"
/*
Greedy algorithm to find approximate lowest z-score global structure in O(nlogn)+O(1.5n) time
*/
namespace approximate
{
    std::vector<std::shared_ptr<basepair::BasePair>> greedy_approximation(matrix::BasePairMatrix& mat);
}
#endif