#include "approximate.hpp"

using namespace approximate;

std::vector<base_pair_pointer> approximate::greedy_approximation(matrix::BasePairMatrix& mat)
{
    std::vector<base_pair_pointer> all_possible_pairs; //stores base pairs sorted by Znorm
    std::vector<base_pair_pointer> chosen_pairs; //stores final set of pairs
    size_t unpaired_bases = mat.getSequenceLength();
    std::cout << "unpaired_bases: " << unpaired_bases << std::endl;
    //create a vector to store whether each base is paired/unpaired
    std::vector<bool> is_available(unpaired_bases, true);
    size_t num_available = is_available.size();
    std::cout << "num_available: " << num_available << std::endl;
    //get all scanned pairs from the base pair matrix
    size_t index = 0;
    std::vector<base_pair_pointer> must_be_paired;
    for(auto& row : mat.Matrix)
    {
        if(row[0]->icoord < 0)
        {
            //this base was always paired with something when scanned
            //so it gets priority when pairing since leaving it unpaired would be an error
            //and something else may outcompete it for the pairs it has
            //it's still possible for it to be left unpaired, since this is always a possibility
            //(if two bases each pair with exactly one other and it's the same one)
            //but this should minimize the times that happens
            size_t this_nucleotide = row[0]->icoord;
            size_t window_size = mat.getWinSize();
            size_t zero = 0;
            size_t start = std::min(zero, (index - window_size));
            size_t end = std::min((index + window_size), (mat.Matrix.size()-1));
            for(size_t i = start; i <= end; i++)
            {
                for(auto pair : mat.Matrix[i])
                {
                    if(pair->icoord == this_nucleotide || pair->jcoord == this_nucleotide)
                    {
                        must_be_paired.push_back(pair);
                    }
                }
            }
        }
        for(auto pair : row)
        {
            if(pair->icoord > -1 && pair->jcoord > -1)
            {
                all_possible_pairs.push_back(pair);
            }
        }
        ++index;
    }
    std::cout << "all_possible_pairs: " << all_possible_pairs.size() << std::endl;
    //sort pairs via z-norm in ascending order (lambda needed since these are pointers)
    std::sort(all_possible_pairs.begin(), all_possible_pairs.end(), 
        [](const base_pair_pointer a, const base_pair_pointer b) 
            {return a->getZNorm() < b->getZNorm();}
    );
    std::cout << "must_be_paired: " << must_be_paired.size() << std::endl;
    std::sort(must_be_paired.begin(), must_be_paired.end(), 
        [](const base_pair_pointer a, const base_pair_pointer b) 
            {return a->getZNorm() < b->getZNorm();}
    );
    base_pair_pointer current_pair;
    size_t i_coord;
    size_t j_coord;
    //choose best pair until all bases are paired
    //unpaired is represented as pairing to self
    //TODO: first run the algorithm on the bases that 
    //can't be left unpaired to minimize erroneously 
    //unmatched bases, then on the rest
    
    std::cout << "looping over must_be_paired..." << std::endl;
    for(auto it = must_be_paired.begin(); it != all_possible_pairs.end(); ++it) 
    {
        current_pair = *it;
        current_pair->print();
        i_coord = current_pair->icoord;
        j_coord = current_pair->jcoord;
        if(is_available[i_coord] && is_available[j_coord])
        {
            std::cout << "available!" << std::endl;
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
    }
    
    for(auto it = all_possible_pairs.begin(); unpaired_bases > 0 && it != all_possible_pairs.end(); ++it) 
    {
        //current_pair = *std::next(it);
        current_pair = *it;
        //current_pair->print();
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
    return chosen_pairs;
}
py::list approximate::py_greedy_approximation(matrix::BasePairMatrix& mat)
{
    std::vector<base_pair_pointer> pair_vector = greedy_approximation(mat);
    py::list pair_list;
    for(auto pair : pair_vector)
    {
        pair_list.append(pair);
    }
    return pair_list;
}