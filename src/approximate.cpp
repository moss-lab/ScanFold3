#include "approximate.hpp"

using namespace approximate;

std::vector<base_pair_pointer> approximate::greedy_approximation(matrix::BasePairMatrix& mat)
{
    std::vector<base_pair_pointer> all_possible_pairs;  //stores base pairs sorted by Znorm
    std::vector<base_pair_pointer> chosen_pairs;        //stores final set of pairs
    size_t unpaired_bases = mat.getSequenceLength();
    //create a vector to store whether each base is paired/unpaired
    std::vector<bool> is_available(unpaired_bases, true);
    //get all scanned pairs from the base pair matrix
    size_t index = 0;
    for(auto& row : mat.Matrix)
    {
        for(auto pair : row)
        {
            if(pair->icoord > -1 && pair->jcoord > -1)
            {
                all_possible_pairs.push_back(pair);
            }
        }
        ++index;
    }
    std::cout << "number of possible pairs: " << all_possible_pairs.size() << std::endl;
    //sort pairs via z-norm in ascending order (lambda needed since these are pointers)
    std::sort(all_possible_pairs.begin(), all_possible_pairs.end(), 
        [](const base_pair_pointer a, const base_pair_pointer b) 
            {return a->getZNorm() < b->getZNorm();}
    );
    base_pair_pointer current_pair;
    size_t i_coord;
    size_t j_coord;
    //choose best pair until all bases are paired
    //unpaired is represented as pairing to self
    for(auto it = all_possible_pairs.begin(); unpaired_bases > 0 && it != all_possible_pairs.end(); ++it) 
    {
        current_pair = *it;
        i_coord = current_pair->icoord;
        j_coord = current_pair->jcoord;
        if(is_available[i_coord] && is_available[j_coord])
        {
            chosen_pairs.push_back(current_pair);
            if(i_coord == j_coord)
            {
                unpaired_bases -= 1;
            }
            else
            {
                unpaired_bases -= 2;
            }
            is_available[i_coord] = false;
            is_available[j_coord] = false;
        }
        //std::cout << std::endl;
    }
    if(unpaired_bases > 0)
    {
        std::cout << unpaired_bases << " bases could not be paired!" << std::endl;
    }
    //sort pairs by position in sequence
    std::sort(chosen_pairs.begin(), chosen_pairs.end(), 
        [](const base_pair_pointer a, const base_pair_pointer b) 
            {return a->icoord < b->icoord;}
    );
    return chosen_pairs;
}
py::list approximate::py_greedy_approximation(matrix::BasePairMatrix& mat)
{
    //python wrapper
    std::vector<base_pair_pointer> pair_vector = greedy_approximation(mat);
    py::list pair_list;
    for(auto pair : pair_vector)
    {
        pair_list.append(pair);
    }
    return pair_list;
}