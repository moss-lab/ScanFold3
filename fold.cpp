//compiled w/ g++ (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0
//requires boost in include path, use 'sudo apt-get install libboost-all-dev'
#include "fold.h"

void BasePair::update(BasePair& newData) 
{
    /*
        add cumulative data for another instance of this i,j pair
        input: another BasePair object
        output: none
    */
    pairs_read += newData.pairs_read;
    _zscore += newData._zscore;
    _mfe += newData._mfe;
    _ed += newData._ed;
    _pvalue += newData._pvalue;
}
void BasePair::copy(BasePair& newData)
{
    /*
    copy data from another BasePair into this one, overwriting what's there
    input: another BasePair
    output: none
    */
    this->pairs_read = newData.pairs_read;
    this->_zscore = newData._zscore;
    this->_mfe = newData._mfe;
    this->_ed = newData._ed;
    this->_pvalue = newData._pvalue;
}
//metrics are stored as cumulative values of all pairs read, use getter functions to access
double BasePair::getZNorm() {return _zscore/pairs_read;}
/*
    return normalized z-score for this pair; undefined for default constructed BasePairs
    input: none
    output: normalized z-score
*/
double BasePair::getAvgZScore() {return _zscore/win_size;}
/*
    return average z-score for this i,j pair. Inaccurate when less than one full window length was scanned, ie first/last 120
    for window size of 120
    input: none
    output: average z-score
*/
//average/normalized values for other metrics
double BasePair::getMFEnorm() {return _mfe/pairs_read;}
/*
    return normalized mean free energy score for this i,j pair; undefined for default constructed BasePairs
    input: none
    output: normalized MFE score
*/
double BasePair::getAvgMFE() {return _mfe/win_size;}
/*
    return average mean free energy score for this i,j pair
    input: none
    output: average MFE score
*/
double BasePair::getPValNorm() {return _pvalue/pairs_read;}
/*
    return normalized p-value for this i,j pair; undefined for default constructed BasePairs
    input: none
    output: normalized p-value
*/
double BasePair::getAvgPVal(){return _pvalue/win_size;}
/*
    return average p-value for this i,j pair
    input: none
    output: average p-value
*/
double BasePair::getEDNorm() {return _ed/pairs_read;}
/*
    return normalized ensemble diversity score for this i,j pair; undefined for default constructed BasePairs
    input: none
    output: normalized ED score
*/
double BasePair::getAvgED() {return _ed/win_size;}
/*
    return average ensemble diversity score for this i,j pair
    input: none
    output: normalized ED score
*/
void BasePair::print() 
{
    //print data for debugging
    std::cout << "coords: " << this->icoord << ", " << this->jcoord << std::endl;
    std::cout << "nucleotides: " << this->inuc << ", " << this->jnuc << std::endl;
    std::cout << "pairs read: " << this->pairs_read << std::endl;
    std::cout << "cumulative values: " << std::endl;
    std::cout << "z-score: " << this->_zscore << "\t";
    std::cout << "mfe: " << this->_mfe << "\t";
    std::cout << "ED: " << this->_ed << "\t";
    std::cout << "p-value: " << this->_pvalue << "\t" << std::endl;
    std::cout << "average values: " << std::endl;
    std::cout << "z-score: " << this->getAvgZScore() << "\t";
    std::cout << "mfe: " << this->getAvgMFE() << "\t";
    std::cout << "ED: " << this->getAvgED() << "\t";
    std::cout << "p-value: " << this->getAvgPVal() << "\t" << std::endl;
    std::cout << "normalized values: " << std::endl;
    std::cout << "z-score: " << this->getZNorm() << "\t";
    std::cout << "mfe: " << this->getMFEnorm() << "\t";
    std::cout << "ED: " << this->getEDNorm() << "\t";
    std::cout << "p-value: " << this->getPValNorm() << "\t" << std::endl;
    
}
void BasePair::printError() 
{
    //print data for debugging
    std::cerr << "error encountered in BasePair!" << std::endl;
    std::cerr << "coords: " << this->icoord << ", " << this->jcoord << std::endl;
    std::cerr << "nucleotides: " << this->inuc << ", " << this->jnuc << std::endl;
    std::cerr << "pairs read: " << this->pairs_read << std::endl;
    std::cerr << "cumulative values: " << std::endl;
    std::cerr << "z-score: " << this->_zscore << "\t";
    std::cerr << "mfe: " << this->_mfe << "\t";
    std::cerr << "ED: " << this->_ed << "\t";
    std::cerr << "p-value: " << this->_pvalue << "\t" << std::endl;
    std::cerr << "average values: " << std::endl;
    std::cerr << "z-score: " << this->getAvgZScore() << "\t";
    std::cerr << "mfe: " << this->getAvgMFE() << "\t";
    std::cerr << "ED: " << this->getAvgED() << "\t";
    std::cerr << "p-value: " << this->getAvgPVal() << "\t" << std::endl;
    std::cerr << "normalized values: " << std::endl;
    std::cerr << "z-score: " << this->getZNorm() << "\t";
    std::cerr << "mfe: " << this->getMFEnorm() << "\t";
    std::cerr << "ED: " << this->getEDNorm() << "\t";
    std::cerr << "p-value: " << this->getPValNorm() << "\t" << std::endl;
    std::cerr << "/error" << std::endl;
    
}
//BasePair constructor
//default values (see fold.h)
BasePair::BasePair() {}
//set coord/nucleotides only
//pairs_read will be set to 0
BasePair::BasePair(int i, int j, char ival, char jval) :
        icoord{i}, jcoord{j}, inuc{ival}, jnuc {jval} {}
//set all values (debug use only, pairs_read should be left to the constructor normally)
BasePair::BasePair(int i, int j, char ival, char jval, double z, double m, double e, double p, int win, int num) :
        icoord{i},
        jcoord{j},
        inuc{ival},
        jnuc{jval},
        _zscore{z},
        _mfe{m},
        _ed{e},
        _pvalue{p},
        win_size{win},
        pairs_read{num} {}
//set all but pairs read (this is set to one then increased automatically whenever new data are added)
//this is the version that should normally be used
BasePair::BasePair(int i, int j, char ival, char jval, double z, double m, double e, double p, int win) :
        icoord{i},          //int - first coordinate in a pair
        jcoord{j},          //int - second coordinate in a pair
        inuc{ival},         //character - nucleotide at the i-th position
        jnuc{jval},         //character - nucleotide at the j-th position
        _zscore{z},         //double - z-score of all instances of this i,j pair added together
        _mfe{m},            //double - mean free energy of all instances of this i,j pair added together
        _ed{e},             //double - ensemble diversity of all instances of this i,j pair added together
        _pvalue{p},         //double - p-value of all instances of this i,j pair added together
        win_size{win},      //int - size of the window that was scanned, ScanFold-Scan defaults to 120
        pairs_read{1} {}    //all instances of this i,j pair that were found, used to get normalized values (do not edit directly!)
//BasePairMatrix constructor from scanfold-scan windows
BasePairMatrix::BasePairMatrix(std::vector<ScanFoldWindow> &windows)
{ 
    //will have one column for each nucleotide
    //and one row for each window, so the last window's worth
    //if a given pairing never appears it will have NaN for the z-score
    std::cout << "constructing BasePairMatrix..." << std::endl;
    //find step sizes and window sizes
    int stepSize = getStepSize(windows);
    std::cout << "step size: " << stepSize << std::endl;
    int winSize = getWindowSize(windows[0]);
    std::cout << "window size: " << winSize << std::endl;
    //sequence length is where the last window ends (since it's indexed to 1)
    int seqLength = 0;
    //make sure not to use windows.back() if for some reason it's empty because segfault
    if(!windows.empty()) {seqLength = windows.back().End;}
    std::cout << "sequence length: " << seqLength << std::endl;
    this->i_length = seqLength - winSize + 1;
    this->j_length = seqLength;
    //create empty BasePairMatrix
    this->Matrix.resize(j_length);
    std::cout << "initializing Matrix" << std::endl;
    for (int row = 0; row < j_length; row++) 
    {
        for(int col = 0; col < winSize; col++)
        {
            BasePair B;
            Matrix[row].push_back(B);
        }
    }
    std::cout << "dimensions: " << Matrix.size() << "x" << Matrix.back().size() << std::endl;
    std::cout << "extracting pairs from windows..." << std::endl;
    std::vector<BasePair> pairs;    //stores pairs found across all windows
    //loop over all windows
    for(auto win : windows)
    {
        //get pairs from each window
        std::vector<BasePair> this_window_pairs = win.getPairs();
        //add them to pairs, using iterators
        pairs.insert(pairs.end(), this_window_pairs.begin(), this_window_pairs.end());
    }
    //add all BasePairs from pairs to BasePairMatrix
    std::cout << "adding pairs to matrix..." << std::endl;
    for (auto pair : pairs) 
    {
        bool errorcheck = this->update(pair);
        if(!errorcheck) 
        {
            std::cerr << "error in updating matrix at " << pair.icoord << ", " << pair.jcoord << "!" << std::endl;
            pair.printError();
        }
        double pair_znorm = pair.getZNorm();
        if(max_znorm > pair_znorm)
        {
            max_znorm = pair_znorm;
        }
    }
}

BasePair& BasePairMatrix::get(int i, int j) 
{
    /*
    i,j are coordinates for a NxN matrix, where N is sequence length
    what is stored is NxM, where M is window size, as anything outside that isn't scanned anyway
    need to make sure i is less than j if it isn't already (matrix would be mirrored so assume i,j = j,i)
    also need to convert j to the actual position, j is indexed to the diagonal of the NxN matrix
    */
    //ensure i<j, it doesn't matter if you swap them because i,j = j,i
    if(i > j) {swap(&i, &j);} 
    //set jcoord to where it is in the data structure
    j = j - i;
    //return invalid pair (i,j = -1,-1) if jcoord is outside what was scanned
    if(j >= window_size) {std::cerr << "invalid pair" << std::endl; return _invalid;}
    return Matrix[i][j];
}
bool BasePairMatrix::update(BasePair& newData) 
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
    BasePair& pair = this->get(newData.icoord, newData.jcoord);
    //tried to reference something that wasn't scanned (coords set to -2)
    if(pair.icoord == -2 || pair.jcoord == -2) 
    {
        std::cerr << "coordinates " << newData.icoord << ", " << newData.jcoord << " outside of scanned window!\n";
        return false;
    } 
    //no existing base pair - default BasePair object
    if(pair.pairs_read == 0) 
    {   
        std::cout << "no pairs read for " << newData.icoord << ", " << newData.jcoord << ", creating new BasePair" << std::endl;
        //in this case, set the pair to newData
        pair = newData;
        return true;
    }
    //existing base pair
    else 
    {
        pair.update(newData);
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
    BasePair& pair = this->get(i,j);
    //tried to reference something that wasn't scanned
    if(pair.icoord == -2 || pair.jcoord == -2) {return NULL;} 
    return pair.getZNorm();    
}
double BasePairMatrix::getAvgZScore(int i, int j) 
{
    //get average z-score for a pair at a given coordinate
    BasePair& pair = this->get(i,j);
    //tried to reference something that wasn't scanned
    if(pair.icoord == -2 || pair.jcoord == -2) 
    {
        std::cerr << "coordinates " << i << ", " << j << " outside of scanned window!\n";
        return false;
    } 
    return pair.getAvgZScore();
}
void BasePairMatrix::toGraph()
{   
    /*
    convert the base pair matrix to a boost adjacency list, see typedef for Graph
    in fold.h for details
    input: none
    output: Graph
    */
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
    this->graph_ptr = std::make_unique<Graph>(num_vertices);
    //iterate over the base pair matrix
    std::vector<std::vector<BasePair>>::iterator row;
    std::vector<BasePair>::iterator col;
    for(row = this->Matrix.begin(); row != this->Matrix.end(); row++){
        for(col = row->begin(); col != row->end(); col++)
        {
            BasePair current_pair = *col;
            //case that a certain pairing never appeared in scanning
            if(current_pair.pairs_read < 1) {continue;}
            //get zNorm for edge weight, make negative since algorithm is max weighted matching
            //subtract max_score so every edge has a weight <= 0, then multiple by -10000 to make all positive
            //and to not lose too much when converting to int
            //round
            //cast to int
            //side note: maximum_weighted_matching gets stuck in some kind of infinite loop if you give it weights as double
            //nothing abt that in the documentation, but that's why it's done this way here
            int weight = (int) std::round(-10000*((current_pair.getZNorm())-this->max_znorm));
            std::cout << current_pair.icoord << ", " << current_pair.jcoord << " weight: " << weight << std::endl;
            //store the coordinates of this pair
            int i_val = current_pair.icoord;
            int j_val = current_pair.jcoord;
            //case for unpaired nucleotides
            //if i==j, add edge for i, i+number_of_nucleotides so an edge exists b/t the two identical nucleotides in the identical graphs
            if(i_val == j_val)
            {
                boost::add_edge(i_val, i_val+number_of_nucleotides, EdgeProperty(weight), *graph_ptr);
                std::cout << "adding self-edge: " << i_val << ", " << i_val+number_of_nucleotides << "\tweight: " << weight << std::endl; 
            }
            //case for paired nucleotides
            //else add edge for i,j and i+number_of_nucleotides,j+number_of_nucleotides to make graphs mirrored
            else
            {
                boost::add_edge(i_val, j_val, EdgeProperty(weight), *graph_ptr);
                boost::add_edge(i_val+number_of_nucleotides, j_val+number_of_nucleotides, EdgeProperty(weight), *graph_ptr);
                std::cout << "adding edge: " << i_val << ", " << j_val << "\tweight: " << weight << std::endl; 
                std::cout << "adding redundant edge: " << i_val+number_of_nucleotides << ", " << j_val+number_of_nucleotides << "\tweight: " << weight << std::endl; 
            }
        }
    }
}
void BasePairMatrix::toCSV()
{
    //print pairs that were scanned to a csv file
    //rows are nucleotides, columns are 1 window's length downstream of that nucleotide
    //first entry in a row is for the nucleotide being unpaired
    //open matrix.csv
    std::ofstream matrixFile;
    matrixFile.open("matrix.csv");
    std::cout << "writing csv..." << std::endl;
    for(auto line : this->Matrix)
    {
        for(auto pair : line)
        {
            if(pair.pairs_read > 0)
            {
                matrixFile << pair.getZNorm() << ",";
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
            if(pair.pairs_read > 0)
            {
                ofile << pair.icoord << "\t" << pair.jcoord << std::endl;
            }
        }
    } 
}
void BasePairMatrix::matchPairs(std::vector<BasePair>& pairs)
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
    //create graph if it hasn't been already
    if(!graph_ptr)
    {
        this->toGraph();
    } 
    //make sure 'pairs' is empty
    if(!pairs.empty()) {pairs.clear();}
    //create vectors to hold matching
    //for each nucleotide i, mate[i] is the vertex it matches to 
    //(itself if vertex >= number of nucleotides, since those should be the only
    //edges crossing over to the latter half of the graph)
    std::vector<boost::graph_traits<Graph>::vertex_descriptor> mate(n_vertices); 
    std::cout << "performing weighted matching..." << std::endl;
    maximum_weighted_matching(*graph_ptr, &mate[0]);
    std::cout << "weighted matching complete!" << std::endl;
    std::cout << "num vertices: " << boost::num_vertices(*graph_ptr) << std::endl;
    std::cout << "len of mate[]: " << mate.size() << std::endl;
    //iterate over vertices and find matches
    std::cout << "final pairs:" << std::endl;
    for (boost::tie(vi, vi_end) = vertices(*graph_ptr); vi != vi_end; ++vi)
    {    
        //null_vertex() means nothing was matched there so ignore those
        //also check to make sure the vertex is less than what it's matched to so we don't count the same matching twice
        if(mate[*vi] != boost::graph_traits<Graph>::null_vertex() && *vi < mate[*vi])
        {
            //for each icoord, finds its match (jcoord)
            int icoord = *vi;
            int jcoord = mate[*vi];
            std::cout << "coords: " << icoord << ", " << jcoord << std::endl;
            //skip over any coords on the mirrored portion of the graph
            if(icoord >= number_of_nucleotides)
            {
                std::cout << "skipped!" << std::endl;
                continue;
            }
            //if jcoord falls outside of the sequence length that means that this nucleotide was unpaired
            //represent that as self pairing
            if(jcoord >= number_of_nucleotides) {jcoord = icoord;}
            std::cout << "updated coords: " << icoord << ", " << jcoord << std::endl;
            //go back to the matrix and get the matched pair
            BasePair pair = this->get(icoord, jcoord);
            std::cout << "pair:" << std::endl;
            pair.print();
            //this shouldn't be possible since this was checked for in toGraph(), but just in case
            if(pair.inuc == 'N' || pair.jnuc == 'N') {throw std::runtime_error("Invalid pair encountered during matching!");}
            //add pairs from matching to vector pairs
            //note that this doesn't check for redundant pairs
            pairs.push_back(this->get(icoord, jcoord));
            std::cout << icoord << ", " << jcoord << std::endl;
            std::cout << std::endl;
        }
    }
}
//functions to read scanfold-scan file into base pair matrix
//go to a specific line
std::ifstream& goToLine(std::ifstream& file, unsigned int num) 
{
    // go to a specific line in a text file
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i) {file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');}
    return file;
}
void swap(int* x, int* y)
{
    /*
    swap two integers in-place w/ bitwise operations
    use different variables otherwise the integer becomes 0
    ok if x and y are same value
    */
    *x = *x ^ *y;
    *y = *x ^ *y;
    *x = *x ^ *y;
    return;
}

int getStepSize(const std::vector<ScanFoldWindow>& windows)
{ 
    /*
    get step size from vector of ScanFoldWindows
    requires vector has length >= 2
    input: vector of ScanFoldWindow
    output: step size used in ScanFold-Scan, default is 1
    */
    if(windows.size() < 2) {throw std::runtime_error("less than two windows in scanfold-scan output!");}
    int step_size = windows[1].Start - windows[0].Start;
    return step_size;
}

int getStepSize(std::ifstream& file) 
{
    /*
    get step size from the .tsv produced by ScanFold-Scan
    this will set the file stream back to the first line!
    input: ifstream
    output: step size used in ScanFold-Scan, default is 1
    */

    //ensure ifstream is on to the first (header) line
    file.clear();
    file.seekg(0);
    std::string line;           //store one line from the .tsv
    std::string line_1_start;   //'Start' column of first content line
    std::string line_2_start;   //'Start' column of second content line
    int step_size;
    //skip header
    std::getline(file, line);
    //read first line 
    std::getline(file, line);
    std::istringstream line_1(line);
    //tab-separated so this gets the first column, which is Start
    line_1 >> line_1_start;
    //read second line
    std::getline(file, line);
    std::istringstream line_2(line);
    line_2 >> line_2_start;

    step_size = std::stoi(line_2_start) - std::stoi(line_1_start);
    //reset file stream
    file.clear();
    //go back to header
    file.seekg(0);
    return step_size;
}

int getWindowSize(const ScanFoldWindow& window)
{
    std::cout << "calculating window size..." << std::endl;
    /*
    get window size from a ScanFoldWindow
    input: ScanFoldWindow
    output: window size (default from ScanFold-Scan is 120)
    */
    int win_size = window.End - window.Start + 1;
    return win_size;
}

int getWindowSize(std::ifstream& file) 
{
    /*
    get window size from an ifstream
    important: this will set the file stream back to the first line
    use it before you process the file
    input: .tsv from ScanFold-Scan
    output: window size (default from ScanFold-Scan is 120)
    */

    //ensure ifstream is on to the first (header) line
    file.clear();
    file.seekg(0);
    std::string line;
    std::string line_1_start;
    std::string line_1_end;
    int win_size;
    //skip header
    std::getline(file, line);
    //read first line start
    std::getline(file, line);
    std::istringstream iss(line);
    //read first line start
    iss >> line_1_start;
    //read first line end;
    iss >> line_1_end;
    win_size = std::stoi(line_1_end) - std::stoi(line_1_start) + 1; //add 1 since the .tsv is 1-indexed
    //reset file stream
    file.clear();
    //go back to header
    file.seekg(0);
    return win_size;
}

ScanFoldWindow::ScanFoldWindow(std::string& line) 
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

std::vector<ScanFoldWindow> readScanTSV(std::ifstream& infile) 
{
    /*
    read the .tsv file from ScanFold-Scan into a vector of window objects
    input: .tsv from ScanFold-Scan
    output: vector of ScanFoldWindow
    */
    std::cout << "reading TSV..." << std::endl;
    std::string line;
    std::vector<ScanFoldWindow> windows; 
    //skip header
    std::getline(infile, line);
    
    //read contents
    while (std::getline(infile, line)) 
    {
        //skip empty lines
        if(line.empty()) {continue;}
        //otherwise create a ScanFoldWindow for that line
        ScanFoldWindow Window(line);
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

std::vector<std::pair<int,int>> findPairsInDotBracket(const std::string& dbstructure) 
{
    //return a vector of std::pair containing i,j coordinates
    //unpaired nucleotides are returned as i,i
    std::stack<int> iloc;                   //indices of start parens
    std::stack<int> iloc_pk1;               //indices of pseudoknot pairs indicated by []
    std::stack<int> iloc_pk2;               //indices of pseudoknot pairs indicated by {}
    std::stack<int> iloc_pk3;               //indices of pseudoknot pairs indicated by <>
    std::vector<std::pair<int,int>> pairs;  //pairs found in dot-bracket structure
    int icoord;                             //index of first nucleotide in pair
    int jcoord;                             //index of last nucleotide in pair
    for (int idx = 0; idx < int(dbstructure.size()); idx++) 
    {   
        switch (dbstructure[idx])
        {
            /*
            in the case of an open paren, will add it to the respective stack and move on to the next iteration
            in the case of a close paren, it will pop the corresponding open paren from the respective stack and 
            write it to the pairs vector
            in the case of something unpaired, will return it as pairing to itself
            */
            case '(':
                //if an open paren is found, store it and move to next loop iteration
                iloc.push(idx);
                break;
            case ')': 
                //if close paren is found, set it to jcoord and set last found open paren to icoord
                icoord = iloc.top();
                iloc.pop();
                jcoord = idx;
                pairs.push_back(std::pair<int,int> (icoord, jcoord));
                break; 
            case '[':
                //do same for pseudoknotted pairs 
                iloc_pk1.push(idx);
                break;
            case ']': 
                icoord = iloc_pk1.top();
                iloc_pk1.pop();
                jcoord = idx;
                pairs.push_back(std::pair<int,int> (icoord, jcoord));
                break;
            case '{': 
                iloc_pk2.push(idx);
                break;
            case '}': 
                icoord = iloc_pk2.top();
                iloc_pk2.pop();
                jcoord = idx;
                pairs.push_back(std::pair<int,int> (icoord, jcoord));
                break;
            case '<': 
                iloc_pk3.push(idx);
                break;
            case '>': 
                icoord = iloc_pk3.top();
                iloc_pk3.pop();
                jcoord = idx;
                pairs.push_back(std::pair<int,int> (icoord, jcoord));
                break;
            default: 
                //if no paren was found, store it as an unpaired nucleotide (i,i)
                icoord = idx;
                jcoord = idx;
                pairs.push_back(std::pair<int,int> (icoord, jcoord));
        }
    }
    return pairs;
}
std::vector<BasePair> ScanFoldWindow::getPairs() 
{
    /*
    read this window's MFE structure as a vector of BasePairs
    pairs adjusted based on the start coord of this window
    so the pairs it gives you are for the entire sequence that was scanned, 0-indexed 
    passing by value to avoid issues elsewhere
    input: none, but requires ScanFoldWindow was constructed w/ some data
    output: vector of BasePair
    */
    std::vector<std::pair<int,int>> pair_indices = findPairsInDotBracket(this->Structure);  
    std::vector<BasePair> pairs;    //stores all pairs found in this window as BasePairs
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
        BasePair base_pair(icoord, jcoord, ival, jval, this->Zscore, this->NativeMFE, this->ED, this->pvalue, this->window_size);
        pairs.push_back(base_pair);
    }    
    return pairs;
}

//main: test cases, for dev usage
//TODO: remove in live release
int main (int argc, char *argv[])
{   
    std::cout << ":3" << std::endl;
    freopen("log.txt", "w", stdout);
    freopen("error.txt", "w", stderr);
    std::string scanfold_scan_fname = "./test/coronaframeshift/fs.1.win_120.stp_1.tsv";
    std::ifstream scanfold_scan(scanfold_scan_fname.c_str());    
    std::vector<ScanFoldWindow> scan_data = readScanTSV(scanfold_scan); 
    std::cout << "BEGIN: testing functions to read scanfold-scan input" << std::endl;
    //readScanTSV()
    std::cout << "first three lines of test tsv: " << std::endl;
    for(int i = 0; i < 3; i++)
    {
        std::cout << "window " << i << ":" << std::endl;
        std::cout << "start: " << scan_data[i].Start << "\nend: " << scan_data[i].End;
        std::cout << "\ntemperature: " << scan_data[i].Temperature << "\nmfe: " << scan_data[i].NativeMFE;
        std::cout << "\nz-score: " << scan_data[i].Zscore << "\npvalue: " << scan_data[i].pvalue;
        std::cout << "\nED: " << scan_data[i].ED << "\nwindow_size: " << scan_data[i].window_size;
        std::cout << "\nsequence: " << scan_data[i].Sequence << "\ncentroid: " << scan_data[i].centroid;
        std::cout << "\nstructure: " << scan_data[i].Structure << std::endl;
    }
    std::cout << "number of windows (should be 250): " << scan_data.size() << std::endl;
    //findPairsInDotBracket()
    std::cout << "testing dot-bracket to pair" << std::endl;
    std::vector<std::pair<int, int>> pairs_win1 = findPairsInDotBracket(scan_data[0].Structure);
    //no pseudoknot
    std::cout << scan_data[0].Structure << std::endl;
    std::string pair_back_to_db = "";
    for(int i = 0; i < 120; i++) 
    {
        pair_back_to_db += "-";
    }
    for(auto pair : pairs_win1)
    {
        if(pair.first == pair.second)
        {
            pair_back_to_db[pair.first] = '.';
        }
        else
        {
            pair_back_to_db[pair.first] = '(';
            pair_back_to_db[pair.second] = ')';
        }
    }
    std::cout << pair_back_to_db << std::endl;
    //getWindowSize()
    int winSize = getWindowSize(scanfold_scan);
    std::cout << "window size, expect 120:\n" << winSize << std::endl;
    //getStepSize()
    int stepSize = getStepSize(scanfold_scan);
    std::cout << "step size, expect 1:\n" << stepSize << std::endl;
    //BasePair
    BasePair bp1 = BasePair(0, 1, 'A', 'U', 1.0, 2.0, 3.0, 4.0, 120, 1);
    std::cout << bp1.getZNorm() << '\t' << bp1.getAvgZScore() << std::endl;
    std::cout << bp1.getMFEnorm() << '\t' << bp1.getAvgMFE() << std::endl;
    std::cout << bp1.getEDNorm() << '\t' << bp1.getAvgED() << std::endl;
    std::cout << bp1.getPValNorm() << '\t' << bp1.getAvgPVal() << std::endl;
    //BasePairMatrix
    //convert scan_data to BasePairMatrix:
    BasePairMatrix bpmatrix = BasePairMatrix(scan_data);
    bpmatrix.toCSV();
    bpmatrix.print();
    std::vector<BasePair> pairs;
    bpmatrix.matchPairs(pairs);
    std::ofstream ofile;
    ofile.open("matching.txt");
    std::vector<std::pair<int, int>> intpairs;
    for(auto pair : pairs)
    {
        std::pair<int, int> p = std::make_pair(pair.icoord, pair.jcoord);
        intpairs.push_back(p);
    }
    std::sort(intpairs.begin(), intpairs.end());
    ofile << "i coord\t" << "j coord\t"  << "z-norm" << std::endl;
    for(auto ip : intpairs)
    {
        double znorm = bpmatrix.getZNorm(ip.first, ip.second);
        ofile << ip.first << "\t" << ip.second << "\t" << znorm << std::endl;
    }
    ofile.close();
}
