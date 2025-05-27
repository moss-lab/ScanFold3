#include "approximate.hpp"

using namespace approximate;

std::vector<base_pair_pointer> approximate::greedy_approximation(matrix::BasePairMatrix& mat)
{
    std::vector<base_pair_pointer> all_possible_pairs;  //stores base pairs sorted by Znorm
    std::vector<base_pair_pointer> chosen_pairs;        //stores final set of pairs
    std::vector<base_pair_pointer> must_be_paired;      //stores pairs involving bases that  
                                                        //were never seen to be unpaired 
                                                        //while scanning
    size_t unpaired_bases = mat.getSequenceLength();
    //create a vector to store whether each base is paired/unpaired
    std::vector<bool> is_available(unpaired_bases, true);
    //get all scanned pairs from the base pair matrix
    size_t index = 0;
    for(auto& row : mat.Matrix)
    {
        if(row[0]->icoord < 0)
        {
            //this base was always paired with something when scanned
            //so it gets priority when pairing since leaving it unpaired would be an error
            //and something else may outcompete it for the pairs it has
            //it's still possible for it to be left unpaired, since this is always a possibility
            //(ie if two bases each pair with exactly one other and it's the same one)
            //but this should minimize the times that happens
            size_t this_nucleotide = index; //index of this base in the sequence and row number are the
                                            //same, regardless of window size
            size_t window_size = mat.getWinSize();
            //start and end contain the indices of windows this base occurs within
            //assuming step size of 1 but naturally it works for any other step size
            size_t start = std::min(static_cast<size_t>(0), (index - window_size));//size_t wraps around on negative
            size_t end = std::min((index + window_size), (mat.Matrix.size()-1));
            std::cout << "start: " << start << std::endl;
            std::cout << "end: " << end << std::endl;
            //loop over a window length in either direction
            for(size_t i = start; i <= end; ++i)
            {
                for(auto pair : mat.Matrix[i])
                //loop over pairs in row i
                {
                    //store all pairs containing the nucleotide we're looking at
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
    std::cout << "number of possible pairs: " << all_possible_pairs.size() << std::endl;
    //sort pairs via z-norm in ascending order (lambda needed since these are pointers)
    std::sort(all_possible_pairs.begin(), all_possible_pairs.end(), 
        [](const base_pair_pointer a, const base_pair_pointer b) 
            {return a->getZNorm() < b->getZNorm();}
    );
    std::cout << "number of bases that can't be unpaired: " << must_be_paired.size() << std::endl;
    std::sort(must_be_paired.begin(), must_be_paired.end(), 
        [](const base_pair_pointer a, const base_pair_pointer b) 
            {return a->getZNorm() < b->getZNorm();}
    );
    base_pair_pointer current_pair;
    size_t i_coord;
    size_t j_coord;
    //choose best pair until all bases are paired
    //unpaired is represented as pairing to self
    //first run the algorithm on the bases that 
    //can't be left unpaired to minimize erroneously 
    //unmatched bases, then on the rest 
    std::cout << "looping over bases that can't be unpaired..." << std::endl;
    for(auto it = must_be_paired.begin(); it != must_be_paired.end(); ++it) 
    {
        //this one doesn't stop early when all bases are paired since that would
        //complicate things too much, it should be short enough that the extra
        //inefficiency is irrelevant anyway
        current_pair = *it;
        current_pair->print();
        i_coord = current_pair->icoord;
        j_coord = current_pair->jcoord;
        //check if bases haven't been paired
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
    }
    //looping over everything else    
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