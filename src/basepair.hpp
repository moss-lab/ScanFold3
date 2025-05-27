#ifndef BASEPAIR
#define BASEPAIR
/*
class to define a single base pair and its metrics
*/
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <limits>
#include "shared.hpp"
namespace basepair {
    namespace py = pybind11;
    //class to store per base pair data
    class BasePair 
        : public std::enable_shared_from_this<BasePair>
    {
        private:
            //cumulative metrics, use getter functions for averaged and normalized
            double zscore {0.0};  
            double mfe {0.0};
            double ed {0.0};
            double pvalue {0.0};

        public:
            //coordinates of -1 mean a BasePair wasn't initialized with any data
            int icoord {-1};
            int jcoord {-1};
            char inuc {'N'};
            char jnuc {'N'};
            size_t win_size {1};
            size_t pairs_read {0};
            size_t seq_length {1};
            //constructors
            //default values
            BasePair();
            //window size and length only
            BasePair(size_t win, size_t len);
            //coord/nucleotides only
            BasePair(int i, int j, char ival, char jval, size_t len);
            //all values
            BasePair(int i, int j, char ival, char jval, double z, double m, double e, double p, size_t win, size_t len, size_t num);
            //all but pairs read 
            BasePair(int i, int j, char ival, char jval, double z, double m, double e, double p, size_t win, size_t len); 
            //functions
            //update values for a new instance of a pair
            void update(BasePair& newData);
            //replace with another BasePair
            void swap(BasePair& newData);
            //getter functions
            //note that getter functions return by value, it's not possible to overwrite a metric except through update()
            double getZNorm();
            double getZNorm() const;
            double getCovgZNorm();
            double getAvgZScore(); 
            double getMFENorm();
            double getAvgMFE();
            double getPValNorm();
            double getAvgPVal();
            double getEDNorm();
            double getAvgED();
            //functions needed to be used with FlowGraph
            unsigned int getStart();
            unsigned int getEnd();
            double getWeight();
            //print data for debug
            void print();
            void py_print();
            //print BasePair as a line to a csv
            void print(std::ofstream& ofile);
            void printError(); 
            std::shared_ptr<BasePair> getptr();             //get pointer to this pair
            bool operator < (const BasePair& other);   //by default, based on znorm
            bool operator > (const BasePair& other);   //by default, based on znorm
            bool operator <= (const BasePair& other);   //by default, based on znorm
            bool operator >= (const BasePair& other);   //by default, based on znorm
    };
    //functions making use of BasePair class
    std::vector<BasePair> filterBasePairs(std::vector<BasePair>& pairs, 
        double min=std::numeric_limits<double>::min(), 
        double max=std::numeric_limits<double>::max());
    py::list py_filterBasePairs(py::list& pairs, 
        double min=std::numeric_limits<double>::min(), 
        double max=std::numeric_limits<double>::max());
}
#endif