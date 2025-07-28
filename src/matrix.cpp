#include "matrix.hpp"
using namespace matrix;
//BasePairMatrix constructor from scanfold-scan windows
BasePairMatrix::BasePairMatrix(std::vector<window::ScanFoldWindow> &windows)
{
    //will have one column for each nucleotide
    //and one row for each window, so the last window's worth
    //if a given pairing never appears it will have NaN for the z-score
    //TODO: make more memory efficient, based on other constructor
    std::cout << "constructing BasePairMatrix..." << std::endl;
    //find step sizes and window sizes
    this->window_size = window::getWindowSize(windows[0]);
    this->step_size = window::getStepSize(windows);
    std::cout << "window size: " << window_size << std::endl;
    //sequence length is where the last window ends (since it's indexed to 1)
    size_t sequence_length = 0;
    //make sure not to use windows.back() if for some reason it's empty because segfault
    if(!windows.empty()) {sequence_length = windows.back().End;}
    std::cout << "sequence length: " << sequence_length << std::endl;
    this->i_length = sequence_length - window_size + 1;
    this->j_length = sequence_length;
    //create empty BasePairMatrix
    this->Matrix.resize(j_length);
    std::cout << "initializing Matrix" << std::endl;
    for (int row = 0; row < j_length; row++) 
    {
        for(int col = 0; col < window_size; col++)
        {
            auto B = std::make_shared<basepair::BasePair>(window_size, sequence_length);
            Matrix[row].push_back(B);
        }
    }
    std::cout << "dimensions: " << Matrix.size() << "x" << Matrix.back().size() << std::endl;
    std::cout << "extracting pairs from windows..." << std::endl;
    std::vector<base_pair_pointer> pairs;    //stores pairs found across all windows
    //loop over all windows
    for(auto win : windows)
    {
        //get pairs from each window (shared pointers)
        auto this_window_pairs = win.getPairs();
        //add them to pairs, using iterators
        pairs.insert(pairs.end(), this_window_pairs.begin(), this_window_pairs.end());
    }
    std::cout << "num pairs: " << pairs.size() << std::endl;
    //add all BasePairs from pairs to BasePairMatrix
    std::cout << "adding pairs to matrix..." << std::endl;
    for (auto pair : pairs) 
    {
        bool errorcheck = this->update(*pair);
        if(!errorcheck) 
        {
            std::cerr << "error in updating matrix at " << pair->icoord << ", " << pair->jcoord << "!" << std::endl;
            pair->printError();
        }
        if(pair->getZNorm() > this->max_znorm)
        {
            this->max_znorm = pair->getZNorm();
        }
    }
}
/*
BasePairMatrix::BasePairMatrix(std::string tsv_name)
{
    std::cout << tsv_name << std::endl;
    //
    //this is the constructor called from w/i Python
    //
    std::ifstream infile(tsv_name);
    std::string line;
    //for memory usage reasons, only reading 1 line at a time

    //get length of overall sequence
    size_t sequence_length = shared::getSequenceLength(infile); 
    this->window_size = shared::getWindowSize(infile);
    this->i_length = sequence_length - window_size + 1;
    this->j_length = sequence_length;
    //create empty BasePairMatrix
    this->Matrix.resize(j_length);
    std::cout << "initializing Matrix" << std::endl;
    for (int row = 0; row < j_length; row++) 
    {
        for(int col = 0; col < window_size; col++)
        {
            basepair::BasePair B(window_size, sequence_length);
            Matrix[row].push_back(B);
        }
    }
    //skip header
    std::getline(infile, line);
    //read contents
    while (std::getline(infile, line)) 
    {
        //skip empty lines
        if(line.empty()) {continue;}
        //otherwise create a ScanFoldWindow for that line
        window::ScanFoldWindow Window(line, sequence_length);
        //turn to vector of BasePair
        std::vector<basepair::BasePair> pairs = Window.getPairs();
        //loop over pairs and add to matrix
        for (basepair::BasePair pair : pairs) 
        {
            bool errorcheck = this->update(pair);
            if(!errorcheck) 
            {
                std::cerr << "error in updating matrix at " << pair.icoord << ", " << pair.jcoord << "!" << std::endl;
                pair.printError();
            }
            if(pair.getZNorm() > this->max_znorm)
            {
                //this will be necessary later to ensure the graph has no negative weights
                this->max_znorm = pair.getZNorm();
            }
        }
    }
    unsigned int size = 0;
    for (auto col : this->Matrix)
    {
        for (auto pair : col) 
        {
            size += sizeof(pair);
        }
    }
    std::cout << "BasePair size: " << sizeof(this->Matrix[0][0]) << std::endl;
    std::cout << "matrix size: " << size << std::endl;
}
*/
BasePairMatrix::BasePairMatrix(std::string tsv_name)
{
    std::cout << "BP Matrix from tsv: ";
    std::cout << tsv_name << std::endl;
    std::cout << "c++ cwd: " << std::filesystem::current_path() << std::endl;
    std::ifstream ifile(tsv_name);
    if (!ifile.is_open())
    {
        std::cout << "couldn't open file " << tsv_name << std::endl;
    }
    ifile.open(tsv_name, std::ifstream::in);
    if (!ifile.is_open())
    {
        std::cout << "couldn't open file (fatal) " << tsv_name << std::endl;
    }
    auto windows = window::readScanTSV(ifile);
    //TODO: remove, test
    windows[0].print();
    /*
    BasePairMatrix new_matrix = BasePairMatrix(windows);    //this is fucked but I don't want to use an init method
    std::swap(*this, new_matrix);
    */
    std::cout << "constructing BasePairMatrix..." << std::endl;
    //find step sizes and window sizes
    this->window_size = window::getWindowSize(windows[0]);
    this->step_size = window::getStepSize(windows);
    std::cout << "window size: " << window_size << std::endl;
    //sequence length is where the last window ends (since it's indexed to 1)
    size_t sequence_length = 0;
    //make sure not to use windows.back() if for some reason it's empty because segfault
    if(!windows.empty()) {sequence_length = windows.back().End;}
    std::cout << "sequence length: " << sequence_length << std::endl;
    this->i_length = sequence_length - window_size + 1;
    this->j_length = sequence_length;
    //create empty BasePairMatrix
    this->Matrix.resize(j_length);
    std::cout << "initializing Matrix" << std::endl;
    bool exec_once = true;
    for (int row = 0; row < j_length; row++) 
    {
        for(int col = 0; col < window_size; col++)
        {
            auto B = std::make_shared<basepair::BasePair>(window_size, sequence_length);
            Matrix[row].push_back(B);
        }
    }
    std::cout << "dimensions: " << Matrix.size() << "x" << Matrix.back().size() << std::endl;
    std::cout << "extracting pairs from windows..." << std::endl;
    std::vector<base_pair_pointer> pairs;    //stores pairs found across all windows
    //loop over all windows
    for(auto win : windows)
    {
        //get pairs from each window (shared pointers)
        auto this_window_pairs = win.getPairs();
        //add them to pairs, using iterators
        pairs.insert(pairs.end(), this_window_pairs.begin(), this_window_pairs.end());
    }
    std::cout << "num pairs: " << pairs.size() << std::endl;
    //add all BasePairs from pairs to BasePairMatrix
    std::cout << "adding pairs to matrix..." << std::endl;
    for (auto pair : pairs) 
    {
        bool errorcheck = this->update(*pair);
        if(!errorcheck) 
        {
            std::cerr << "error in updating matrix at " << pair->icoord << ", " << pair->jcoord << "!" << std::endl;
            pair->printError();
        }
        if(pair->getZNorm() > this->max_znorm)
        {
            this->max_znorm = pair->getZNorm();
        }
    }
    //std::cout << "first pair in matrix: " << std::endl;
    //this->Matrix[0][0]->print();
}
size_t BasePairMatrix::getWinSize()
{
    return this->window_size;
}
size_t BasePairMatrix::getStepSize()
{
    return this->step_size;
}
base_pair_pointer BasePairMatrix::get(int i, int j) 
{
    /*
    i,j are coordinates for a NxN matrix, where N is sequence length
    what is stored is NxM, where M is window size, as anything outside that isn't scanned anyway
    need to make sure i is less than j if it isn't already (matrix would be mirrored so assume i,j = j,i)
    also need to convert j to the actual position, j is indexed to the diagonal of the NxN matrix
    */
    //ensure i<j, it doesn't matter if you swap them because i,j = j,i
    if(i > j) {shared::swapInt(&i, &j);} 
    //set jcoord to where it is in the data structure
    j = j - i;
    //return invalid pair (i,j = -1,-1) if jcoord is outside what was scanned
    if(j >= this->window_size) {std::cerr << "invalid pair" << std::endl; return _invalid.getptr();}
    return Matrix[i][j];
}
bool BasePairMatrix::update(basepair::BasePair& newData) 
{
    /*
    update this matrix with a BasePair object
    input: coordinates of the BasePair and the BasePair itself
    output: true if the coordinate was updated successfully, false for a failure
    a failure should only happen if you try to reference something that wasn't scanned
    */
    //returns false if i or j fall outside the size the matrix was given when constructed
    if(newData.icoord > j_length) {std::cerr << "i coord outside of sequence!" << std::endl;return false;}
    if(newData.jcoord > j_length) {std::cerr << "j coord outside of sequence!\n" << std::endl;return false;}
    auto pair = this->get(newData.icoord, newData.jcoord);
    //tried to reference something that wasn't scanned (coords set to -2)
    if(pair->icoord == -2 || pair->jcoord == -2) 
    {
        std::cerr << "coordinates " << newData.icoord << ", " << newData.jcoord << " outside of scanned window!\n";
        return false;
    } 
    //no existing base pair - default BasePair object
    if(pair->pairs_read == 0) 
    {   
        //in this case, set the pair to newData
        //std::swap(*pair, newData);
        pair->swap(newData);
        return true;
    }
    //existing base pair
    else 
    {
        pair->update(newData);
        return true;
    }
    //just in case
    std::cerr << "unknown error in BasePairMatrix.update()!" << std::endl;
    std::cerr << "coordinates: " << newData.icoord << ", " << newData.jcoord << std::endl;
    newData.printError();
    return false;
}
double BasePairMatrix::getZNorm(int i, int j) 
{
    /*
    get normalized z-score for a pair at a given coordinate
    input: i,j coordinate for the pair (in the overall sequence scanned, not per window)
    output: normalized z-score
    */
    auto pair = this->get(i,j);
    //tried to reference something that wasn't scanned
    if(pair->icoord == -2 || pair->jcoord == -2) {return NULL;} 
    return pair->getZNorm();    
}
double BasePairMatrix::getAvgZScore(int i, int j) 
{
    //get average z-score for a pair at a given coordinate
    auto pair = this->get(i,j);
    //tried to reference something that wasn't scanned
    if(pair->icoord == -2 || pair->jcoord == -2) 
    {
        std::cerr << "coordinates " << i << ", " << j << " outside of scanned window!\n";
        return false;
    } 
    return pair->getAvgZScore();
}
double BasePairMatrix::getMaxZNorm() {return max_znorm;}
size_t BasePairMatrix::getSequenceLength() {return this->j_length;}
void BasePairMatrix::toCSV(std::string fname)
{
    //print pairs that were scanned to a csv file
    //rows are nucleotides, columns are 1 window's length downstream of that nucleotide
    //first entry in a row is for the nucleotide being unpaired
    //open matrix.csv
    std::ofstream matrixFile;
    matrixFile.open(fname);
    std::cout << "writing csv..." << std::endl;
    for(auto line : this->Matrix)
    {
        for(auto pair : line)
        {
            if(pair->pairs_read > 0)
            {
                matrixFile << pair->getZNorm() << ",";
            }
            else
            {
                matrixFile << ",";
            }
        }
        matrixFile << "\n";
    }
    matrixFile.close();
}
void BasePairMatrix::print()
{
    std::ofstream ofile;
    ofile.open("matrix.txt");
    std::cout << "printing Matrix:" << std::endl;
    for(auto line : this->Matrix)
    {
        for(auto pair : line)
        {
            if(pair->pairs_read > 0)
            {
                ofile << pair->icoord << "\t" << pair->jcoord << std::endl;
            }
        }
    } 
}
py::list BasePairMatrix::py_getAllPairs()
{
    py::list all_pairs;
    for(auto row : this->Matrix)
    {
        for(auto pair : row)
        {
            if(pair->icoord >= 0)
            {
                all_pairs.append(pair);
            }
        }
    }
    return all_pairs;
}
void BasePairMatrix::add_zscore_for_unpaired(std::vector<window::ScanFoldWindow> &windows)
{
    //TODO:
    /*
    Go through the filled base pair matrix, finds any rows where the first entry is unfilled
    which means that base was never observed to be unpaired. Refold each window that base was 
    found within but with it constrained to be unpaired, then refold shuffled versions of those
    windows without that constraint to get a thermodynamic z-score for the base being left unpaired
    
    haven't implemented this yet because it's so rare that a base can't be paired (like 200 bases in the
    entire 170kb EBV genome) that it's definitely an edge case
    it's also an issue in scanfold2 and no one's complained yet, so, eh?
    */
}
/*
//call and print to csv to save memory
void BasePairMatrix::getBestPairing(std::string& ofile)
{
    std::ofstream outfile;
    outfile.open(ofile);
    std::vector<BasePair> pair_vec;
    this->getBestPairing(pair_vec); 
    for(auto pair : pair_vec)
    {
        pair.print(outfile);
    }
    outfile.close();
}
*/
std::string BasePairMatrix::getSequence()
{
    std::string sequence;
    for (auto &row : this->Matrix)
    {
        for (auto &pair : row)
        {
            if (pair->inuc == 'N')
            {
                continue;
            }
            sequence.push_back(pair->inuc);
        }
    }
    auto &last_line = this->Matrix[this->Matrix.size()-1];
    for (auto it = last_line.begin()+1;it != last_line.end();++it)
    {
        if ((*it)->inuc == 'N')
        {
            continue;
        }
        sequence.push_back((*it)->inuc);
    }
    return sequence;
}