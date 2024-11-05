//compiled w/ g++ (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0
//requires boost in include path, use 'sudo apt-get install libboost-all-dev'
#ifndef FOLD
#define FOLD
#include <vector>
#include <fstream>
#include <ios>
#include <limits>
#include <string>
#include <sstream>
#include <utility>
#include <stack>
#include <tuple>
#include <stdexcept>
#include <iostream>
#include <memory>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>
#include <pybind11/pybind11.h>
namespace py = pybind11;

//create a typedef for the graph and edge weights
typedef boost::property<boost::edge_weight_t, int> EdgeProperty;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, EdgeProperty> Graph;
//struct to store per base pair data
class BasePair 
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
        //getter functions
        //note that getter functions return by value, it's not possible to overwrite a metric except through update()
        double getZNorm();
        double getAvgZScore(); 
        double getMFENorm();
        double getAvgMFE();
        double getPValNorm();
        double getAvgPVal();
        double getEDNorm();
        double getAvgED();
        //print data for debug
        void print();
        void printError(); 
};
//struct to store a window from ScanFold
struct ScanFoldWindow 
{
    //variable names same as in scanfold-scan results
    unsigned int Start;  //indexed to 1
    unsigned int End;    //indexed to 1
    double Temperature;
    double NativeMFE;
    double Zscore;
    double pvalue;
    double ED;
    std::string Sequence;
    std::string Structure;
    std::string centroid;
    size_t window_size;
    size_t sequence_length;
    //constructor
    ScanFoldWindow(std::string& line, size_t len);
    //create a vector of BasePair from this window and return it
    std::vector<BasePair> getPairs();
    //print window contents
    void print();
};
//read output file from scanfoldscan into a vector of ScanFoldWindows
//can be passed directly into constructor of BasePairMatrix
std::vector<ScanFoldWindow> readScanTSV(std::ifstream& infile);
//find pairs in a dot-bracket structure
std::vector<std::pair<int,int>> findPairsInDotBracket(const std::string& dbstructure);
//get window/step size from scanfold-scan output
size_t getWindowSize(std::ifstream& file);
size_t getStepSize(std::ifstream& file);
//get window from ScanFoldWindow
size_t getWindowSize(const ScanFoldWindow& window);
//get step size from vector of ScanFoldWindow
size_t getStepSize(const std::vector<ScanFoldWindow>& windows);
//get length of the overall sequence that was scanned
size_t getSequenceLength(std::ifstream& file);
//struct to store all base pairs (vector of vectors w/ znorm getter function)
class BasePairMatrix 
{
    private:
        std::vector<std::vector<BasePair>> Matrix;
        size_t i_length{0};
        size_t j_length{0};
        size_t window_size{120};
        double max_znorm{0.0};
        BasePair _invalid{-2,-2,'N','N',1}; //reference to this is returned from get() when you try to access something that wasn't scanned
                                            //reason for this is to avoid potential memory leaks from creating a new BasePair w/i that scope
                                            //DO NOT CHANGE THIS! It can't be set to const unfortunately due to being returned, otherwise it would be
                                            //default sets coordinates to -1, this one sets them to -2
    public:
        //constructor from vector of ScanFoldWindow
        BasePairMatrix(std::vector<ScanFoldWindow> &windows);
        //constructor from tsv file name
        BasePairMatrix(std::string tsv_name);
        //functions
        //update base pair at newData's position or add if none
        //return false if i,j are outside of scanned windows
        //this is an error and probably indicates a mistake somewhere, so check for it
        bool update(BasePair& newData);
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
        BasePair& get(int i, int j);
        //create Graph from the contents of Matrix
        std::unique_ptr<Graph> toGraph();
        //return max matching as a vector of pairs, in-place (pairs is input and output)
        //O(n^3)
        void getBestPairing(std::vector<BasePair>& pairs);
        //python wrapper for getBestPairing
        void getBestPairing(py::list &pairs);
        //store matrix as a .csv of znorm values
        void toCSV(std::string fname);
        //print out pairs
        void print();
};
//go to the last line of a file
std::streampos findLastLine(std::ifstream& file);
//swap integers
void swapInt(int* x, int* y);
//python bindings
PYBIND11_MODULE(fold, var)
{
    var.doc() = "Classes and functions needed for ScanFold-Fold";
    var.def("readScanTSV", &readScanTSV, "read the .tsv file produced by scanfold-scan into a C++ data structure");
    py::class_<BasePair>(var, "BasePair")
        .def(py::init<size_t, size_t>())
        .def(py::init<int, int, char, char, size_t>())
        .def(py::init<int, int, char, char, double, double, double, double, size_t, size_t, size_t>())
        .def(py::init<int, int, char, char, double, double, double, double, size_t, size_t>())
        .def("getZNorm", &BasePair::getZNorm)
        .def("getAvgZScore", &BasePair::getAvgZScore)
        .def("getMFENorm", &BasePair::getMFENorm)
        .def("getAvgMFE", &BasePair::getAvgMFE)  
        .def("getPValNorm", &BasePair::getPValNorm)
        .def("getAvgPVal", &BasePair::getAvgPVal)
        .def("getEDNorm", &BasePair::getEDNorm)
        .def("getAvgED", &BasePair::getAvgED)
        .def("print", &BasePair::print)
        .def_readwrite("i_coord", &BasePair::icoord)
        .def_readwrite("j_coord", &BasePair::jcoord)
        .def_readwrite("i_nucleotide", &BasePair::inuc)
        .def_readwrite("j_nucleotide", &BasePair::jnuc)
    ;
    py::class_<ScanFoldWindow>(var, "ScanFoldWindow")
        .def(py::init<std::string&, size_t>())
        .def("getPairs", &ScanFoldWindow::getPairs)
        .def("print", &ScanFoldWindow::print)
        .def_readonly("Start", &ScanFoldWindow::Start)
        .def_readonly("End", &ScanFoldWindow::End)
        .def_readonly("Temperature", &ScanFoldWindow::Temperature)
        .def_readonly("NativeMFE", &ScanFoldWindow::NativeMFE)
        .def_readonly("Zscore", &ScanFoldWindow::Zscore)
        .def_readonly("pvalue", &ScanFoldWindow::pvalue)
        .def_readonly("Sequence", &ScanFoldWindow::Sequence)
        .def_readonly("Structure", &ScanFoldWindow::Structure)
        .def_readonly("centroid", &ScanFoldWindow::centroid)
    ; 
    py::class_<BasePairMatrix>(var, "BasePairMatrix")
        .def(py::init<std::string>())
        .def("getZNorm", &BasePairMatrix::getZNorm)
        .def("getAvgZScore", &BasePairMatrix::getAvgZScore)
        .def("getBestPairing", py::overload_cast<py::list&>(&BasePairMatrix::getBestPairing))
        .def("toCSV", &BasePairMatrix::toCSV)
    ;
}
#endif
