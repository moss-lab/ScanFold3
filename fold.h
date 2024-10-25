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
//create a typedef for the graph and edge weights
typedef boost::property<boost::edge_weight_t, int> EdgeProperty;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, EdgeProperty> Graph;
//struct to store per base pair data
struct BasePair 
{
    int icoord {-1};
    int jcoord {-1};
    char inuc {'N'};
    char jnuc {'N'};
    double _zscore {0.0};  //cumulative z-score, use getter functions for averaged and normalized
    double _mfe {0.0};
    double _ed {0.0};
    double _pvalue {0.0};
    int win_size {0};
    int pairs_read {0};
    //constructors
    //default values
    BasePair();
    //coord/nucleotides only
    BasePair(int i, int j, char ival, char jval);
    //all values
    BasePair(int i, int j, char ival, char jval, double z, double m, double e, double p, int win, int num);
    //all but pairs read 
    BasePair(int i, int j, char ival, char jval, double z, double m, double e, double p, int win); 
    //functions
    //update values for a new instance of a pair
    void update(BasePair& newData);
    //deep copy from another BasePair
    void copy(BasePair& newData);
    //getter functions
    //note that getter functions return by value
    double getZNorm();
    double getAvgZScore(); 
    double getMFEnorm();
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
    int Start;  //indexed to 1
    int End;    //indexed to 1
    double Temperature;
    double NativeMFE;
    double Zscore;
    double pvalue;
    double ED;
    int window_size;
    std::string Sequence;
    std::string Structure;
    std::string centroid;
    //constructor
    ScanFoldWindow(std::string& line);
    //create a vector of BasePair from this window and return it
    std::vector<BasePair> getPairs();
    //print window contents
    void print();
};
//read output file from scanfoldscan into a vector of ScanFoldWindows
std::vector<ScanFoldWindow> readScanTSV(std::ifstream& infile);
//find pairs in a dot-bracket structure
std::vector<std::pair<int,int>> findPairsInDotBracket(const std::string& dbstructure);
//get window/step size from scanfold-scan output
int getWindowSize(std::ifstream& file);
int getStepSize(std::ifstream& file);
//get window from ScanFoldWindow
int getWindowSize(const ScanFoldWindow& window);
//fet step size from vector of ScanFoldWindow
int getStepSize(const std::vector<ScanFoldWindow>& windows);
//create log file
//struct to store all base pairs (vector of vectors w/ znorm getter function)
struct BasePairMatrix 
{
    std::vector<std::vector<BasePair>> Matrix;
    int i_length{0};
    int j_length{0};
    int window_size{120};
    double max_score{0.0};
    BasePair _invalid{-2,-2,'N','N'};   //reference to this is returned from get() when you try to access something that wasn't scanned
                                        //reason for this is to avoid potential memory leaks from creating a new BasePair w/i that scope
                                        //DO NOT CHANGE THIS! It can't be set to const unfortunately, otherwise it would be
                                        //default sets coordinates to -1, this one sets them to -2
    std::unique_ptr<Graph> graph_ptr;   //pointer to graph g, used when constructing it in toGraph
    //constructor
    BasePairMatrix(std::vector<ScanFoldWindow> &windows);
    //functions
    //update base pair at newData's position or add if none
    //return false if i,j are outside of scanned windows
    //this is an error and probably indicates a mistake somewhere, so check for it
    bool update(BasePair& newData);
    //return znorm, or null if invalid i,j
    double getZNorm(int i, int j);
    //return average zscore, or null if invalid i,j
    double getAvgZScore(int i, int j);
    //access base pair
    //returns default base pair (i,j = null) if i,j falls outside scanning window, so always check for that
    //get() used in other getter functions
    //can also be used to modify a BasePair in the matrix, so be careful w/ using it
    BasePair& get(int i, int j);
    //convert to graph
    void toGraph();
    //stores max matching w.r.t. normalized z-score of the matrix in BasePair vector 'pairs'
    //O(n^3)
    void matchPairs(std::vector<BasePair>& pairs); 
    //store matrix as a .csv of znorm values
    void toCSV();
    //print out pairs
    void print();
};
//read scanfoldscan output
//go to a specific line in a file
std::ifstream& goToLine(std::ifstream& file, unsigned int num);
//swap integers
void swap(int* x, int* y);
#endif