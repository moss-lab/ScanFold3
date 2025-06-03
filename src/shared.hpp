#ifndef SHARED
#define SHARED
/*
shared utility functions, should be lowest level in include hierarchy
*/

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stack>
#include <cctype>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace shared {
    namespace py = pybind11;
    //swap integers
    void swapInt(int* x, int* y);
    //find pairs in a dot-bracket structure
    std::vector<std::pair<size_t,size_t>> findPairsInDotBracket(const std::string& dbstructure);
    //python binding
    py::list py_findPairsInDotBracket(const std::string& dbstructure);
    //get dot-bracket structure from pairs
    void getDBFromPairs(std::vector<std::pair<size_t, size_t>> &pairs, std::string &dbstructure, char open_paren='(', bool iserr=false);
    void py_getDBFromPairs(py::list &pairs, std::string &dbstructure);
    //get window/step size from scanfold-scan output
    size_t getWindowSize(std::ifstream& file);
    size_t getStepSize(std::ifstream& file);
    //get length of the overall sequence that was scanned
    size_t getSequenceLength(std::ifstream& file);
    //go to the last line of a file
    std::streampos findLastLine(std::ifstream& file);
    // custom exceptions
    class Exception : public std::exception
    {
        private:
            std::string message;
        public:
            Exception(const char* msg) : message(msg) {}
            Exception(const std::string &msg) : message(msg) {}
            const char* what() const throw()
            {
                return message.c_str();
            }
    };
    bool check_whitespace(std::string&str);
}
#endif