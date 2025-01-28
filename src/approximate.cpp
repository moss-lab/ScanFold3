#include "approximate.hpp"

using namespace approximate;

std::vector<std::shared_ptr<basepair::BasePair>> approximate::greedy_approximation(matrix::BasePairMatrix& mat)
{
    typedef std::shared_ptr<basepair::BasePair> pair_ptr;
    std::vector<pair_ptr> pairs; //stores base pairs sorted by Znorm
    std::vector<pair_ptr> chosen_pairs; //stores final set of pairs
    size_t unpaired_bases = mat.getSequenceLength();
    //create a vector to store whether each base is paired/unpaired
    std::vector<bool> is_available(unpaired_bases, true);
    //get all pairs from the base pair matrix
    for(auto& row : mat.Matrix)
    {
        for(auto& pair : row)
        {
            auto ptr = pair.getptr();
            pairs.push_back(ptr);
        }
    }
    //sort pairs via z-norm in ascending order (lambda needed since these are pointers)
    std::sort(pairs.begin(), pairs.end(), 
        [](const basepair::BasePair* a, const basepair::BasePair* b) 
            {return a->getZNorm() < b->getZNorm();}
    );
    pair_ptr current_pair;
    size_t i_coord;
    size_t j_coord;
    //choose best pair until all bases are paired
    //unpaired is represented as pairing to self
    for(std::vector<pair_ptr>::iterator it = pairs.begin(); unpaired_bases > 0; ++it) 
    {
        current_pair = *std::next(it);
        i_coord = current_pair->icoord;
        j_coord = current_pair->jcoord;
        if(is_available[i_coord] && is_available[j_coord])
        {
            chosen_pairs.push_back(current_pair);
            unpaired_bases -= 2;
            is_available[i_coord] = false;
            is_available[j_coord] = false;
        }
    }
    return chosen_pairs;
}