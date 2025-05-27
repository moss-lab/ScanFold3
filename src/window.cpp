#include "window.hpp"
using namespace window;

size_t window::getStepSize(const std::vector<ScanFoldWindow>& windows)
{ 
    /*
    get step size from vector of ScanFoldWindows
    requires vector has length >= 2
    input: vector of ScanFoldWindow
    output: step size used in ScanFold-Scan, default is 1
    */
    if(windows.size() < 2) {throw std::runtime_error("less than two windows in scanfold-scan output!");}
    size_t step_size = windows[1].Start - windows[0].Start;
    return step_size;
}
size_t window::getWindowSize(const ScanFoldWindow& window)
{
    /*
    get window size from a ScanFoldWindow
    input: ScanFoldWindow
    output: window size (default from ScanFold-Scan is 120)
    */
    size_t win_size = window.End - window.Start + 1;
    return win_size;
}
ScanFoldWindow::ScanFoldWindow(std::string& line, size_t len) :
    sequence_length(len) 
{
    /*
    constructor for ScanFoldWindow
    input: a non-header line from the .tsv produced by ScanFold-Scan
    */
    //read a line from scanfold-scan output into a string stream
    std::istringstream iss(line);
    //read everything stored in the string stream into ScanFoldWindow values
    iss >> Start >> End >> Temperature >> NativeMFE >> Zscore >> pvalue >> ED >> Sequence >> Structure >> centroid;
    window_size = Sequence.size();
}
std::vector<ScanFoldWindow> window::readScanTSV(std::ifstream& infile) 
{
    /*
    read the .tsv file from ScanFold-Scan into a vector of window objects
    input: .tsv from ScanFold-Scan
    output: vector of ScanFoldWindow
    */
    std::cout << "reading TSV..." << std::endl;
    std::string line;
    std::vector<ScanFoldWindow> windows; 
    //get length of overall sequence
    size_t sequence_length = shared::getSequenceLength(infile); 
    //skip header
    std::getline(infile, line);
    //read contents
    while (std::getline(infile, line)) 
    {
        //skip empty lines
        if(line.empty()) {continue;}
        //otherwise create a ScanFoldWindow for that line
        ScanFoldWindow Window(line, sequence_length);
        //and add to vector that will be returned by this function
        windows.push_back(Window);
    }
    //reset file stream
    infile.clear();
    infile.seekg(0);
    return windows;
}
void ScanFoldWindow::print()
{
    std::cout << this->Start << "\t" << this->End << "\t" << this->Temperature << "\t" << this->NativeMFE << "\t";
    std::cout << this->Zscore << "\t" << this->pvalue << "\t" << this->ED << "\t" << this->Sequence << "\t";
    std::cout << this->Structure << "\t" << this->centroid << std::endl;
}
std::vector<std::shared_ptr<basepair::BasePair>> ScanFoldWindow::getPairs() 
{
    /*
    read this window's MFE structure as a vector of BasePairs
    pairs adjusted based on the start coord of this window
    so the pairs it gives you are for the entire sequence that was scanned, 0-indexed 
    passing by value to avoid issues elsewhere
    input: none, but requires ScanFoldWindow was constructed w/ some data
    output: vector of BasePair
    */
    std::vector<std::pair<size_t,size_t>> pair_indices = shared::findPairsInDotBracket(this->Structure);  
    std::vector<std::shared_ptr<basepair::BasePair>> pairs;    //stores all pairs found in this window as BasePairs
    int zero_indexed_start_coord = this->Start - 1; //Start is indexed to 1 in the input file
    for (auto & pair: pair_indices) 
    {       
        //loop over each pair and create a BasePair for each
        int icoord = zero_indexed_start_coord + pair.first;
        int jcoord = zero_indexed_start_coord + pair.second;
        //get nucleotides corresponding to this pair
        char ival = this->Sequence[pair.first];
        char jval = this->Sequence[pair.second];
        //create a BasePair for this pair w/ data from the window
        std::shared_ptr<basepair::BasePair> base_pair = std::make_shared<basepair::BasePair>(icoord, jcoord, ival, jval, this->Zscore, this->NativeMFE, this->ED, this->pvalue, this->window_size, this->sequence_length);
        pairs.push_back(base_pair);
    }    
    return pairs;
}