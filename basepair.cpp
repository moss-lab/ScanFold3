#include "basepair.hpp"
using namespace basepair;

//BasePair constructor
//default values (see fold.h)
BasePair::BasePair() {}
//window size and sequence length only
BasePair::BasePair(size_t win, size_t len) : 
    win_size(win), seq_length(len)
    {
        if(seq_length < 1) {std::cerr << "warning: sequence length of 0 entered for " << icoord << ", " << jcoord << std::endl;}
        if(win_size < 1) {std::cerr << "warning: window size of 0 entered for " << icoord << ", " << jcoord << std::endl;}
    }
//set coord/nucleotides only
//pairs_read will be set to 0
BasePair::BasePair(int i, int j, char ival, char jval, size_t len) :
        icoord(i), jcoord(j), inuc(ival), jnuc(jval), seq_length(len) {}
//set all values (debug use only, pairs_read should be left to the constructor normally)
BasePair::BasePair(int i, int j, char ival, char jval, double z, double m, double e, double p, size_t win, size_t len, size_t num) :
    icoord(i),
    jcoord(j),
    inuc(ival),
    jnuc(jval),
    zscore(z),
    mfe(m),
    ed(e),
    pvalue(p),
    win_size(win),
    seq_length(len),
    pairs_read(num)
    {
        if(seq_length < 1) {std::cerr << "warning: sequence length of 0 entered for " << icoord << ", " << jcoord << std::endl;}
        if(win_size < 1) {std::cerr << "warning: window size of 0 entered for " << icoord << ", " << jcoord << std::endl;}
    }
//set all but pairs read (this is set to one then increased automatically whenever new data are added)
//this is the version that should normally be used
BasePair::BasePair(int i, int j, char ival, char jval, double z, double m, double e, double p, size_t win, size_t len) :
    icoord(i),          //int - first coordinate in a pair
    jcoord(j),          //int - second coordinate in a pair
    inuc(ival),         //character - nucleotide at the i-th position
    jnuc(jval),         //character - nucleotide at the j-th position
    zscore(z),         //double - z-score of all instances of this i,j pair added together
    mfe(m),            //double - mean free energy of all instances of this i,j pair added together
    ed(e),             //double - ensemble diversity of all instances of this i,j pair added together
    pvalue(p),         //double - p-value of all instances of this i,j pair added together
    win_size(win),      //size_t - size of the window that was scanned, ScanFold-Scan defaults to 120
    seq_length(len),    //size_t - length of sequence this base pair is found within, used for normalizing values
    pairs_read(1)     //all instances of this i,j pair that were found, used to get normalized values (do not edit directly!)
    {
        if(seq_length < 1) {std::cerr << "warning: sequence length of 0 entered for " << icoord << ", " << jcoord << std::endl;}
        if(win_size < 1) {std::cerr << "warning: window size of 0 entered for " << icoord << ", " << jcoord << std::endl;}
    }
std::shared_ptr<BasePair> BasePair::getptr()
{
    //get shared pointer to this edge, to add to associated vertices
    return shared_from_this();
}
bool BasePair::operator < (const BasePair& other)
{
    return (this->getZNorm() < other.getZNorm());
}   //by default, based on znorm
bool BasePair::operator > (const BasePair& other)
{
    return (this->getZNorm() > other.getZNorm());
}   //by default, based on znorm
bool BasePair::operator <= (const BasePair& other)
{
    return (this->getZNorm() <= other.getZNorm());
}   //by default, based on znorm
bool BasePair::operator >= (const BasePair& other)
{
    return (this->getZNorm() >= other.getZNorm());
}   //by default, based on znorm
void BasePair::update(std::shared_ptr<BasePair> newData) 
{
    /*
        add cumulative data for another instance of this i,j pair
        input: another BasePair object
        output: none
    */
    pairs_read += newData->pairs_read;
    zscore += newData->zscore;
    mfe += newData->mfe;
    ed += newData->ed;
    pvalue += newData->pvalue;
}
//metrics are stored as cumulative values of all pairs read, use getter functions to access
double BasePair::getZNorm()
{
    /*
    return coverage normalized z-score for this pair
    input: none
    output: normalized z-score
    */
    int distance_from_start = this->icoord + 1;
    int distance_from_end = this->seq_length - this->jcoord;
    //window_occurences tracks the number of windows this pair appeared in 
    auto window_occurences = std::min({distance_from_start, distance_from_end, (int)this->win_size});
    return zscore/window_occurences;
}
double BasePair::getZNorm() const
{
    return this->getZNorm();
}
double BasePair::getAvgZScore() {return zscore/win_size;}
/*
    return average z-score for this i,j pair. Inaccurate when less than one full window length was scanned, ie first/last 120
    for window size of 120
    input: none
    output: average z-score
*/
//average/normalized values for other metrics
double BasePair::getMFENorm() 
{
    /*
    return normalized mean free energy score for this i,j pair
    input: none
    output: normalized MFE score
    */ 
    int distance_from_start = this->icoord + 1;
    int distance_from_end = this->seq_length - this->jcoord;
    //window_occurences tracks the number of windows this pair appeared in 
    auto window_occurences = std::min({distance_from_start, distance_from_end, (int)this->win_size});
    return zscore/window_occurences;
}
double BasePair::getAvgMFE() {return mfe/win_size;}
/*
    return average mean free energy score for this i,j pair
    input: none
    output: average MFE score
*/
double BasePair::getPValNorm() 
{
    /*
    return normalized p-value for this i,j pair
    input: none
    output: normalized p-value
    */
    int distance_from_start = this->icoord + 1;
    int distance_from_end = this->seq_length - this->jcoord;
    //window_occurences tracks the number of windows this pair appeared in 
    auto window_occurences = std::min({distance_from_start, distance_from_end, (int)this->win_size});
    return zscore/window_occurences;
}
double BasePair::getAvgPVal() {return pvalue/win_size;}
/*
    return average p-value for this i,j pair
    input: none
    output: average p-value
*/
double BasePair::getEDNorm() 
{
    /*
    return normalized ensemble diversity score for this i,j pair
    input: none
    output: normalized ED score
    */
    int distance_from_start = this->icoord + 1;
    int distance_from_end = this->seq_length - this->jcoord;
    //window_occurences tracks the number of windows this pair appeared in 
    auto window_occurences = std::min({distance_from_start, distance_from_end, (int)this->win_size});
    return zscore/window_occurences;
}
double BasePair::getAvgED() {return ed/win_size;}
/*
    return average ensemble diversity score for this i,j pair
    input: none
    output: normalized ED score
*/
//functions needed to woth with FlowGraph
//need to return a start coordinate, end coordinate, and weight
unsigned int BasePair::getStart() {return this->icoord+1;}
unsigned int BasePair::getEnd() {return this->jcoord+1;}
//to create a structure based on a metric other than normalized z-score, change this function
double BasePair::getWeight() {return this->getZNorm();}

void BasePair::print() 
{
    //print data for debugging
    std::cout << "coords: " << this->icoord << ", " << this->jcoord << std::endl;
    std::cout << "nucleotides: " << this->inuc << ", " << this->jnuc << std::endl;
    std::cout << "pairs read: " << this->pairs_read << std::endl;
    std::cout << "cumulative values: " << std::endl;
    std::cout << "z-score: " << this->zscore << "\t";
    std::cout << "mfe: " << this->mfe << "\t";
    std::cout << "ED: " << this->ed << "\t";
    std::cout << "p-value: " << this->pvalue << "\t" << std::endl;
    std::cout << "average values: " << std::endl;
    std::cout << "z-score: " << this->getAvgZScore() << "\t";
    std::cout << "mfe: " << this->getAvgMFE() << "\t";
    std::cout << "ED: " << this->getAvgED() << "\t";
    std::cout << "p-value: " << this->getAvgPVal() << "\t" << std::endl;
    std::cout << "normalized values: " << std::endl;
    std::cout << "z-score: " << this->getZNorm() << "\t";
    std::cout << "mfe: " << this->getMFENorm() << "\t";
    std::cout << "ED: " << this->getEDNorm() << "\t";
    std::cout << "p-value: " << this->getPValNorm() << "\t" << std::endl;
    
}
void BasePair::print(std::ofstream& ofile)
{
    ofile << this->icoord << "," << this->jcoord << "," << this->getZNorm() << "," << this->getMFENorm() << std::endl;    
}
void BasePair::printError() 
{
    //print data for debugging
    std::cerr << "error encountered in BasePair!" << std::endl;
    std::cerr << "coords: " << this->icoord << ", " << this->jcoord << std::endl;
    std::cerr << "nucleotides: " << this->inuc << ", " << this->jnuc << std::endl;
    std::cerr << "pairs read: " << this->pairs_read << std::endl;
    std::cerr << "cumulative values: " << std::endl;
    std::cerr << "z-score: " << this->zscore << "\t";
    std::cerr << "mfe: " << this->mfe << "\t";
    std::cerr << "ED: " << this->ed << "\t";
    std::cerr << "p-value: " << this->pvalue << "\t" << std::endl;
    std::cerr << "average values: " << std::endl;
    std::cerr << "z-score: " << this->getAvgZScore() << "\t";
    std::cerr << "mfe: " << this->getAvgMFE() << "\t";
    std::cerr << "ED: " << this->getAvgED() << "\t";
    std::cerr << "p-value: " << this->getAvgPVal() << "\t" << std::endl;
    std::cerr << "normalized values: " << std::endl;
    std::cerr << "z-score: " << this->getZNorm() << "\t";
    std::cerr << "mfe: " << this->getMFENorm() << "\t";
    std::cerr << "ED: " << this->getEDNorm() << "\t";
    std::cerr << "p-value: " << this->getPValNorm() << "\t" << std::endl;
    std::cerr << "/error" << std::endl;
    
}