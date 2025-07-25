#ifndef MATRIX
#define MATRIX
/*
BasePairMatrix class
used to store and access data for all base pairs
*/
#include <vector>
#include <string>
#include <memory>
#include <tuple>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "window.hpp"
#include "basepair.hpp"

#include <filesystem>

namespace matrix {
    namespace py = pybind11;
    typedef std::shared_ptr<basepair::BasePair> base_pair_pointer;
    //create a typedef for the graph and edge weights
    typedef boost::property<boost::edge_weight_t, int> EdgeProperty;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, EdgeProperty> Graph;
    //struct to store all base pairs (vector of vectors w/ znorm getter function)
    class BasePairMatrix 
    {
        private:
            size_t i_length{0};
            size_t j_length{0};
            size_t window_size{120};
            size_t step_size{1};
            double max_znorm{0.0};
            basepair::BasePair _invalid{-2,-2,'N','N',1}; //reference to this is returned from get() when you try to access something that wasn't scanned
                                                //reason for this is to avoid potential memory leaks from creating a new BasePair w/i that scope
                                                //DO NOT CHANGE THIS! It can't be set to const unfortunately due to being returned, otherwise it would be
                                                //default sets coordinates to -1, this one sets them to -2
            //goes through the matrix, finds any bases that were never unpaired across windows
            //and calculates a zscore for them by refolding those windows with the base constrained to be
            //unpaired and an MFE generated through shuffling/refolding those windows
            void add_zscore_for_unpaired(std::vector<window::ScanFoldWindow> &windows);
        public:
            std::vector<std::vector<base_pair_pointer>> Matrix;
            //constructor from vector of ScanFoldWindow
            BasePairMatrix(std::vector<window::ScanFoldWindow> &windows);
            //constructor from tsv file name
            BasePairMatrix(std::string tsv_name);
            //functions
            //get window size
            size_t getWinSize();
            size_t getStepSize();
            //update base pair at newData's position or add if none
            //return false if i,j are outside of scanned windows
            //this is an error and probably indicates a mistake somewhere, so check for it
            bool update(basepair::BasePair& newData);
            //return znorm, or null if invalid i,j
            double getZNorm(int i, int j);
            //return average zscore, or null if invalid i,j
            double getAvgZScore(int i, int j);
            //return a copy of the largest ZNorm found within the matrix
            double getMaxZNorm();
            //access base pair
            //returns default base pair (i,j = null) if i,j falls outside scanning window, so always check for that
            //get() used in other getter functions
            //can also be used to modify a BasePair in the matrix, so be careful w/ using it
            //when using pay attention to make sure you don't try to use a BasePair retrieved with get after the BasePairMatrix is deleted
            base_pair_pointer get(int i, int j);
            //create Graph from the contents of Matrix
            std::unique_ptr<Graph> toGraph();
            //return max matching as a vector of pairs, in-place (pairs is input and output)
            //O(n^3)
            void getBestPairing(std::vector<base_pair_pointer>& pairs);
            //python wrapper for getBestPairing
            void getBestPairing(py::list &pairs);
            py::list py_getAllPairs();
            //call and print to csv to save memory
            //void getBestPairing(std::string& ofile);
            //return length of the sequence that went into the matrix (j_length)
            size_t getSequenceLength();
            //return sequence from matrix
            std::string getSequence();
            //store matrix as a .csv of znorm values
            void toCSV(std::string fname);
            //print out pairs
            void print();
    };
}
#endif