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
namespace shared {
    //swap integers
    void swapInt(int* x, int* y);
    //find pairs in a dot-bracket structure
    std::vector<std::pair<int,int>> findPairsInDotBracket(const std::string& dbstructure);
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
            const char* what() const throw()
            {
                return message.c_str();
            }
    };
}
#endif 