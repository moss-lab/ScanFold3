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
std::unique_ptr<Graph> BasePairMatrix::toGraph()
{   
    /*
    convert the base pair matrix to a boost adjacency list, see typedef for Graph
    in fold.h for details
    input: none
    output: Graph
    */
    //TODO: benchmark this vs return by value to see if there's a difference
    std::cout << "converting BasePairMatrix to graph..." << std::endl;
    //get length of the sequence that was originally scanned, 
    //earlier we made one row for each nucleotide
    //each nucleotide gets two vertices
    int number_of_nucleotides = this->j_length; //one column for each nucleotide
    std::cout << "num nt: " << number_of_nucleotides << std::endl;
    if(number_of_nucleotides == 0) {throw std::runtime_error("Base pair matrix empty!");}
    /*
    to allow self-pairing (representing unpaired nucleotide), we need 2 vertices per nucleotide
    will mirror the graph in vertices above the number of nucleotides and add an edge b/t respective 
    vertices to represent an unpaired nucleotide at that position
    elsewhere unpaired nucleotides are represented as pairing to themselves.
    Adding number_of_nucleotides to a vertex gives you its matching vertex on the
    mirrored graph
    */
    int num_vertices = number_of_nucleotides*2;
    std::cout << "num_vertices: " << num_vertices << std::endl;
    std::unique_ptr<Graph> g = std::make_unique<Graph>(num_vertices);
    //iterate over the base pair matrix
    std::vector<std::vector<base_pair_pointer>>::iterator row;
    std::vector<base_pair_pointer>::iterator col;
    for(row = this->Matrix.begin(); row != this->Matrix.end(); row++){
        for(col = row->begin(); col != row->end(); col++)
        {
            basepair::BasePair current_pair = **col;
            //case that a certain pairing never appeared in scanning
            if(current_pair.pairs_read < 1) {continue;}
            //get zNorm for edge weight, make negative since algorithm is max weighted matching
            //TODO: replace this with actual residuals
            //subtract max_znorm so every edge has a weight <= 0, then multiple by -10000 to make all positive
            //and to not lose too much when converting to int
            //round
            //cast to int
            //side note: maximum_weighted_matching gets stuck in some kind of infinite loop if you give it weights as double
            //nothing abt that in the documentation, but that's why it's done this way here
            int weight = (int) std::round(-10000*((current_pair.getZNorm())-this->max_znorm));
            //store the coordinates of this pair
            int i_val = current_pair.icoord;
            int j_val = current_pair.jcoord;
            //case for unpaired nucleotides
            //if i==j, add edge for i, i+number_of_nucleotides so an edge exists b/t the two identical nucleotides in the identical graphs
            if(i_val == j_val)
            {
                boost::add_edge(i_val, i_val+number_of_nucleotides, EdgeProperty(weight), *g);
            }
            //case for paired nucleotides
            //else add edge for i,j and i+number_of_nucleotides,j+number_of_nucleotides to make graphs mirrored
            else
            {
                boost::add_edge(i_val, j_val, EdgeProperty(weight), *g);
                boost::add_edge(i_val+number_of_nucleotides, j_val+number_of_nucleotides, EdgeProperty(weight), *g);
            }
        }
    }
    return g;
}
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
void BasePairMatrix::getBestPairing(std::vector<base_pair_pointer>& pairs)
{
    /*
    find optimal pairing w.r.t normalized z-score
    usage of underlying function described in https://www.boost.org/doc/libs/1_86_0/libs/graph/doc/maximum_weighted_matching.html
    input: vector of BasePairs (will be modified by function)
    output: none
    */

    //get length of the sequence that was originally scanned, 
    //earlier we made one row for each window
    //each nucleotide gets two vertices
    int number_of_nucleotides = this->j_length;
    std::cout << "number of nucleotides: " << number_of_nucleotides << std::endl;
    if(number_of_nucleotides == 0) {throw std::runtime_error("Base pair matrix empty!");}
    //same as in toGraph()
    int n_vertices = number_of_nucleotides*2;
    //create vertex iterators
    boost::graph_traits<Graph>::vertex_iterator vi, vi_end;
    //create graph
    std::unique_ptr<Graph> g = this->toGraph(); 
    //create vectors to hold matching
    //for each nucleotide i, mate[i] is the vertex it matches to 
    //(itself if vertex >= number of nucleotides, since those should be the only
    //edges crossing over to the latter half of the graph)
    std::vector<boost::graph_traits<Graph>::vertex_descriptor> mate(n_vertices); 
    std::cout << "performing weighted matching..." << std::endl;
    maximum_weighted_matching(*g, &mate[0]);
    std::cout << "weighted matching complete!" << std::endl;
    std::cout << "num vertices: " << boost::num_vertices(*g) << std::endl;
    std::cout << "len of mate[]: " << mate.size() << std::endl;
    //iterate over vertices and find matches
    for (boost::tie(vi, vi_end) = vertices(*g); vi != vi_end; ++vi)
    {    
        //null_vertex() means nothing was matched there so ignore those
        //also check to make sure the vertex is less than what it's matched to so we don't count the same matching twice
        if(mate[*vi] != boost::graph_traits<Graph>::null_vertex() && *vi < mate[*vi])
        {
            //for each icoord, finds its match (jcoord)
            int icoord = *vi;
            int jcoord = mate[*vi];
            //skip over any coords on the mirrored portion of the graph
            if(icoord >= number_of_nucleotides)
            {
                continue;
            }
            //if jcoord falls outside of the sequence length that means that this nucleotide was unpaired
            //represent that as self pairing
            if(jcoord >= number_of_nucleotides) {jcoord = icoord;}
            //go back to the matrix and get the matched pair
            auto pair = this->get(icoord, jcoord);
            //this shouldn't be possible since this was checked for in toGraph(), but just in case
            if(pair->inuc == 'N' || pair->jnuc == 'N') {throw std::runtime_error("Invalid pair encountered during matching!");}
            //add pairs from matching to vector pairs
            //note that this doesn't check for redundant pairs
            pairs.push_back(this->get(icoord, jcoord));
        }
    }
}

void BasePairMatrix::getBestPairing(py::list &pairs)
{
    std::vector<base_pair_pointer> pair_vec;
    this->getBestPairing(pair_vec);
    for(auto pair : pair_vec)
    {
        pairs.append(*pair);
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