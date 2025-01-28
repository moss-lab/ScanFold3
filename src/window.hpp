#ifndef WINDOW
#define WINDOW
/*
ScanFoldWindow and associated functions
used to store each window produced by ScanFold-Scan
*/
#include <string>
#include <vector>
#include <cstddef>
#include <sstream>
#include "basepair.hpp"
#include "shared.hpp"
namespace window {
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
        std::vector<basepair::BasePair> getPairs();
        //print window contents
        void print();
    };
    //read output file from scanfoldscan into a vector of ScanFoldWindows
    //can be passed directly into constructor of BasePairMatrix
    std::vector<ScanFoldWindow> readScanTSV(std::ifstream& infile);
    //get window from ScanFoldWindow
    size_t getWindowSize(const ScanFoldWindow& window);
    //get step size from vector of ScanFoldWindow
    size_t getStepSize(const std::vector<ScanFoldWindow>& windows);
}
#endif