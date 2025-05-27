#include "shared.hpp"
using namespace shared;

void shared::swapInt(int* x, int* y)
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

size_t shared::getWindowSize(std::ifstream& file) 
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
    size_t win_size;
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
    std::cout << win_size << std::endl;
    //reset file stream
    file.clear();
    //go back to header
    file.seekg(0);
    return win_size;
}
size_t shared::getStepSize(std::ifstream& file) 
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
    size_t step_size;
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
size_t shared::getSequenceLength(std::ifstream& file)
{
    /*
    get sequence length from an ifstream
    important: this will set the file stream back to the first line
    use it before you process the file
    input: .tsv from ScanFold-Scan
    output: length of overall sequence that was scanned (size_t)
    */
    std::string last_line_end;
    std::string line;
    size_t sequence_length;
    file.clear();
    //go to one line before EOF
    std::streampos last_line = shared::findLastLine(file);
    file.seekg(last_line, file.beg);
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> last_line_end;
    iss >> last_line_end;
    sequence_length = std::stoi(last_line_end);
    //reset file stream
    file.clear();
    //go back to header
    file.seekg(0);
    return sequence_length;  
}
std::streampos shared::findLastLine(std::ifstream& file)
{
    /*
    takes a ifstream and goes to the last line in it
    works by finding EOF character and going back until it reaches the next line break
    */
    //go last char before EOF
    file.seekg(-1, std::ios_base::end);
    int counter = 0;    //need to skip over the first instance of a newline
    while(true)
    {
        ++counter;
        //c stores the character we're on
        char c;
        file.get(c);
        //check to see if we're at the start of the file
        if( (int)file.tellg() <= 1)
        {
            //stop looping if so
            file.seekg(0);
            break;
        }
        //check to see if we're at a newline and the loop has run at least once
        else if(c == '\n' && counter > 1)
        {
            //stop if so
            break;
        }
        //continue looping in all other cases
        else
        {
            //move 2 characters back so that get will get the previous character
            file.seekg(-2, std::ios_base::cur);
        }
    }
    return file.tellg();
}
py::list shared::py_findPairsInDotBracket(const std::string& dbstructure)
{
    auto pairs = shared::findPairsInDotBracket(dbstructure);
    py::list py_pairs = py::cast(pairs);
    return py_pairs;
}
std::vector<std::pair<size_t,size_t>> shared::findPairsInDotBracket(const std::string& dbstructure) 
{
    //return a vector of std::pair containing i,j coordinates
    //unpaired nucleotides are returned as i,i
    std::stack<size_t> iloc;                   //indices of start parens
    std::stack<size_t> iloc_pk1;               //indices of pseudoknot pairs indicated by []
    std::stack<size_t> iloc_pk2;               //indices of pseudoknot pairs indicated by {}
    std::stack<size_t> iloc_pk3;               //indices of pseudoknot pairs indicated by <>
    std::vector<std::pair<size_t,size_t>> pairs;  //pairs found in dot-bracket structure
    int icoord;                             //index of first nucleotide in pair
    int jcoord;                             //index of last nucleotide in pair
    for (size_t idx = 0; idx < (dbstructure.size()); idx++) 
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
                pairs.push_back(std::pair<size_t,size_t> (icoord, jcoord));
                break; 
            case '[':
                //do same for pseudoknotted pairs 
                iloc_pk1.push(idx);
                break;
            case ']': 
                icoord = iloc_pk1.top();
                iloc_pk1.pop();
                jcoord = idx;
                pairs.push_back(std::pair<size_t,size_t> (icoord, jcoord));
                break;
            case '{': 
                iloc_pk2.push(idx);
                break;
            case '}': 
                icoord = iloc_pk2.top();
                iloc_pk2.pop();
                jcoord = idx;
                pairs.push_back(std::pair<size_t,size_t> (icoord, jcoord));
                break;
            case '<': 
                iloc_pk3.push(idx);
                break;
            case '>': 
                icoord = iloc_pk3.top();
                iloc_pk3.pop();
                jcoord = idx;
                pairs.push_back(std::pair<size_t,size_t> (icoord, jcoord));
                break;
            default: 
                //if no paren was found, store it as an unpaired nucleotide (i,i)
                icoord = idx;
                jcoord = idx;
                pairs.push_back(std::pair<size_t,size_t> (icoord, jcoord));
        }
    }
    return pairs;
}
/*
size_t shared::getDBFromPairs(std::vector<std::pair<size_t, size_t>> &pairs, std::string &dbstructure, size_t start, char open_paren)
{
    //C++ implementation of pairs_to_db() from ScanFold2
    
    takes a list of base pairs (tuple), a dot-bracket structure initialized
    to a list of '.', a start coordinate, and a type of parenthesis to use
    the last two shouldn't be used when calling this
    each iteration of the function will either:
        -do nothing but call a new iteration on the next pair (case of an
         unpaired base outside of a hairpin)
        -fill out a nested hairpin that has no non-nested pairs and call a 
         new iteration on the next pair after that hairpin
        -fill out a nested hairpin, call a new iteration on a non-nested
         hairpin within that nested hairpin, and if it is not a non-nested 
         hairpin itself call a new iteration on the next pair
        -fill out a non-nested hairpin and return without calling a new
         iteration (the previous iteration will handle that), calling a 
         new iteration on any non-nested hairpins it finds within that 
         non-nested hairpin
        -end if it reaches the end of the list of pairs
    in any case it returns index of the next pair in the list, unless it 
    has reached the end of the list, in which case it returns the index of 
    the last pair plus one
    
    //choose type of brackets for this helix
    char close_paren;
    char new_paren;
    switch (open_paren) 
    {
        
        order is (), [], {}, <>; then aA, bB, etc
        once it reaches zZ, further pairs are just set to unpaired
        
        case '(':
            close_paren = ')';
            new_paren = '[';
            break;
        case '[':
            close_paren = ']';
            new_paren = '{';
            break;
        case '{':
            close_paren = '}';
            new_paren = '<';
            break;
        case '<':
            close_paren = '>';
            new_paren = 'a';
            break;
        default:
            if (open_paren > 96 && open_paren < 122)
            {
                close_paren = open_paren - 32;
                new_paren = open_paren + 1;
            }
            else
            {
                open_paren = '.';
                close_paren = '.';
                new_paren = '.';
            }
    }
    size_t new_start = start + 1;
    //reached end of pairs
    if (new_start >= pairs.size())
    {
        return new_start;
    }
    auto outer_i = pairs[start].first;
    auto outer_j = pairs[start].second;
    if (outer_i == outer_j)
    {
        //unpaired base outside of a helix
        //call new iteration on next pair to advance
        auto next_index = getDBFromPairs(pairs, dbstructure, new_start, open_paren);
        return next_index;
    }
    //if i and j are pairs, set them so in dbstructure
    dbstructure[outer_i] = open_paren;
    dbstructure[outer_j] = close_paren;
    //fill out hairpin b/t i and j
    for (auto it = pairs.begin() + new_start; it != pairs.end(); ++it)
    {   
        //get current index
        auto current_pair_index = std::distance(pairs.begin(), it);
        auto current_i = (*it).first;
        auto current_j = (*it).second;
        if (current_i > outer_j)
        { 
            // reached the end of this helix
            // since each pair is represented once, the pair (outer_j, j') doesn't
            // exist so we don't need to worry about it, this pair will always
            // begin after the last helix ends
            // since this will be at the end of any non-nested helices, we also
            // reset 'parens'
            auto next_index = getDBFromPairs(pairs, dbstructure, current_pair_index, '(');
            return next_index;
        }
        else if (current_i == current_j)
        {
            // unpaired
            continue;
        }
        else if (current_i <= outer_i)
        { 
            // a previous pair that wasn't caught in a previous recursion
            // should be unnecessary since unlike previous versions each 
            // pair is represented only once as (i,j) where j > i, sorted by i
            // included as an error just in case
            std::ostringstream err_stream;
            err_stream << "error converting list of pairs to helix at ";
            err_stream << current_i << ", " << current_j;
            err_stream << "(0-indexed), pairs are likely out of order; try sorting them before calling getDBFromPairs()";
            std::string err_str = err_stream.str();
            throw shared::Exception(err_str);
        }
        else if (current_j < outer_j)
        {
            // nested pair
            dbstructure[current_i] = open_paren;
            dbstructure[current_j] = close_paren;
        }
        else if (current_j > outer_j)
        { 
            // non-nested pair, do recursion
            // function is called again with parens used to denote a pseudoknot within
            auto next_index = getDBFromPairs(pairs, dbstructure, current_pair_index, new_paren);
            //returning this may skip over any non-nested pairs after the nested pair
            //TODO: test to see if returning this breaks things
            //return next_index;
        }
        else if (current_j == outer_j || current_i == outer_j)
        {
            std::ostringstream err_stream;
            err_stream << "error converting list of pairs to helix at ";
            err_stream << current_i << ", " << current_j;
            err_stream << "(0-indexed), a pair is likely represented twice";
            std::string err_str = err_stream.str();
            throw shared::Exception(err_str);
        }
    }
    return new_start;
}

void shared::py_getDBFromPairs(py::list &pairs, std::string& dbstructure)
{
    //pairs should be BasePair objects
    std::vector<std::pair<size_t, size_t>> pair_vec;
    try 
    {
        for (const auto& pair : pairs)
        {
            size_t icoord;
            size_t jcoord;
            if (py::isinstance<py::tuple>(pair))
            //usage on python tuples
            {
                py::tuple pypair = pair.cast<py::tuple>();
                icoord = pypair[0].cast<size_t>();
                jcoord = pypair[1].cast<size_t>();
            }
            else
            {
                //usage on BasePair objects
                //chicanery necessary to use python class members in C++
                //better than a circular import to mess with BasePair directly at least
                py::object py_icoord = pair.attr("i_coord");
                py::object py_jcoord = pair.attr("j_coord");
                icoord = py_icoord.cast<size_t>();
                jcoord = py_jcoord.cast<size_t>();
            }
            auto p = std::make_pair(icoord, jcoord);
            pair_vec.push_back(p);
        }
    }
    catch (...)
    {
        throw shared::Exception("Error in getDBFromPairs! Make sure you're passing a list of tuples or list of BasePairs.");
    }
    auto end_idx = shared::getDBFromPairs(pair_vec, dbstructure);
}
*/
void shared::getDBFromPairs(std::vector<std::pair<size_t, size_t>> &pairs, std::string &dbstructure, char open_paren, bool iserr)
{
    /*
    turns a vector of pairs into a dot-bracket structure
    assumes vector is sorted from low to high via first coordinate in pair
    also assumes all nucleotides are paired only once, with the lower coordinate first
    bases paired to themselves are assumed to be unpaired
    */
    if (dbstructure.empty())
    {
        //no dbstructure given, need to create one
        //find largest index, then make a dbstructure of that length
        size_t max_j = 0;
        for (auto pair : pairs)
        {
            if (pair.first > max_j) {
                max_j = pair.first;
            }
            if (pair.second > max_j) {
                max_j = pair.second;
            }
        }
        max_j += 1;
        std::string new_dbstructure(max_j, '.');
        std::swap(dbstructure, new_dbstructure);
    }
    try
        {
        //vector to store any non-nested pairs encountered
        std::vector<std::pair<size_t, size_t>> non_nested_pairs;
        //choose type of brackets for this helix
        char close_paren;
        char new_paren;
        switch (open_paren) 
        {
            /*    
            order is (), [], {}, <>; then aA, bB, etc
            once it reaches zZ, further pairs are just set to unpaired
            */
            case '(':
                close_paren = ')';
                new_paren = '[';
                break;
            case '[':
                close_paren = ']';
                new_paren = '{';
                break;
            case '{':
                close_paren = '}';
                new_paren = '<';
                break;
            case '<':
                close_paren = '>';
                new_paren = 'a';
                break;
            default:
                if (open_paren >= 'a' && open_paren <= 'z')
                {
                    close_paren = open_paren - 32;
                    new_paren = open_paren + 1;
                } 
                else
                {
                    open_paren = '.';
                    close_paren = '.';
                    new_paren = '.';
                }
        }
        //look for nested helices
        size_t outer_i = pairs[0].first;
        size_t outer_j = pairs[0].second;
        for (auto pair : pairs)
        {
            auto current_i = pair.first;
            auto current_j = pair.second;
            if (current_i > current_j) {std::swap(current_i, current_j);}

            //unpaired nucleotides
            if (current_i == current_j)
            {
                continue;
            }
            else if (current_i < outer_i)
            {
                throw shared::Exception("Error in getDBFromPairs(): pairs out of order! Attempting to fix...");
            }
            //have passed a helix
            else if (current_i > outer_j)
            {
                outer_i = current_i;
                outer_j = current_j;
                dbstructure[current_i] = open_paren;
                dbstructure[current_j] = close_paren;
            }
            //nested base pair
            else if (current_i >= outer_i && current_j <= outer_j)
            {
                dbstructure[current_i] = open_paren;
                dbstructure[current_j] = close_paren;
            }
            //non-nested pair
            else if (current_i < outer_j && current_j > outer_j)
            {
                non_nested_pairs.push_back(pair);
            }
        }
        if (non_nested_pairs.size() > 0)
        //if any non-nested pairs were encountered, repeat with a new type of bracket on non-nested pairs
        {
            getDBFromPairs(non_nested_pairs, dbstructure, new_paren, iserr);
        }
    }
    catch (...)
    {
        py::print("Error encountered in getDBFromPairs! Attempting to fix input...");
        //validate input and try again, checks to make sure this section was only called once
        size_t max_jcoord = 0;
        std::vector<std::pair<size_t, size_t>> new_pairs;
        //make sure error handling only happens once
        if (iserr)
        {
            throw shared::Exception("Unknown Error in getDBFromPairs()!");
        }
        //sort pairs
        std::sort(pairs.begin(), pairs.end());
        //loop to remove redundant pairs
        size_t last_i = pairs[0].first;
        size_t last_j = pairs[0].second;
        new_pairs.push_back(pairs[0]);
        for (auto it = pairs.begin()+1; it != pairs.end(); ++it)
        {
            auto pair = *it;
            if (pair.second > max_jcoord)
            {
                //keep track of highest coordinate seen
                max_jcoord = pair.second;
            }
            if (pair.first == last_i && pair.second > last_j)
            {
                //found a coord for this i with a greater j
                new_pairs.pop_back();
                new_pairs.push_back(pair);
                last_j = pair.second;
            }
            else if (pair.first != last_i)
            {
                //new i coord
                last_i = pair.first;
                last_j = pair.second;
                new_pairs.push_back(pair);
            }
            //otherwise this pair gets ignored
        }
        //make sure dbstructure is the right size
        if (dbstructure.size() < max_jcoord)
        {            
            std::string new_dbstructure(max_jcoord, '.');
            py::print("Warning! Dot-bracket structure given to getDBFromPairs is too small! Generating new empty dot-bracket structure");
            py::print("Some pairs may have been dropped attempting to fix the issue, check output.");
            getDBFromPairs(new_pairs, new_dbstructure, '(', true);
            return;
        }
        //try again
        getDBFromPairs(new_pairs, dbstructure, open_paren, true);
    }
}
std::string shared::py_getDBFromPairs(py::list &pairs)
{ 
    std::string dbstructure = "";
    //pairs should be BasePair objects
    std::vector<std::pair<size_t, size_t>> pair_vec;
    try 
    {
        for (const auto& pair : pairs)
        {
            size_t icoord;
            size_t jcoord;
            if (py::isinstance<py::tuple>(pair))
            //usage on python tuples
            {
                py::tuple pypair = pair.cast<py::tuple>();
                icoord = pypair[0].cast<size_t>();
                jcoord = pypair[1].cast<size_t>();
            }
            else
            {
                //usage on BasePair objects
                //chicanery necessary to use python class members in C++
                //better than a circular import to mess with BasePair directly at least
                py::object py_icoord = pair.attr("i_coord");
                py::object py_jcoord = pair.attr("j_coord");
                icoord = py_icoord.cast<size_t>();
                jcoord = py_jcoord.cast<size_t>();
            }
            auto p = std::make_pair(icoord, jcoord);
            pair_vec.push_back(p);
        }
    }
    catch (...)
    {
        throw shared::Exception("Error in getDBFromPairs! Make sure you're passing a list of tuples or list of BasePairs.");
    }
    shared::getDBFromPairs(pair_vec, dbstructure);
    return dbstructure;
}