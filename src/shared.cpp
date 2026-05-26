#include "shared.hpp"
using namespace shared;

// keeps track of how deeply nested a pseudoknot matching a given symbol is
/*
std::unordered_map<char, int> pk_nesting = {
    {'(', 0}, {')', 0}, {'[', 1}, {']', 1},    
    {'{', 2}, {'}', 2}, {'<', 3}, {'>', 3},
    {'A', 4}, {'a', 4}, {'B', 5}, {'b', 5},
    {'C', 6}, {'c', 6}, {'D', 7}, {'d', 7},
    {'E', 8}, {'e', 8}, {'F', 9}, {'f', 9},    
    {'G', 10}, {'g', 10}, {'H', 11}, {'h', 11},
    {'I', 12}, {'i', 12}, {'J', 13}, {'j', 13},
    {'K', 14}, {'k', 14}, {'L', 15}, {'l', 15},
    {'M', 16}, {'m', 16}, {'N', 17}, {'n', 17},    
    {'O', 18}, {'o', 18}, {'P', 19}, {'p', 19},
    {'Q', 20}, {'q', 20}, {'R', 21}, {'r', 21},
    {'S', 22}, {'s', 22}, {'T', 23}, {'t', 23},
    {'U', 24}, {'u', 24}, {'V', 25}, {'v', 25},
    {'W', 26}, {'w', 26}, {'X', 27}, {'x', 27},
    {'Y', 28}, {'y', 28}, {'Z', 29}, {'z', 29}
};
*/
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
   if(file.eof())
   {
    std::cout << "empty" << std::endl;
   }
    std::string last_line_end;
    std::string line;
    size_t sequence_length;
    if (!file.is_open())
    {
        std::cout << "error opening file" << std::endl;
        throw shared::Exception("error opening file!"); 
    }
    file.clear();
    //go to one line before EOF
    std::streampos last_line = shared::findLastLine(file);
    file.seekg(last_line, file.beg);
    std::getline(file, line);
    /*
    while (check_whitespace(line))
    {
        std::cout << "getting new line..." << std::endl;
        std::cout << line << std::endl;
        line.clear();
        std::getline(file, line);
    }
    */
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
    if(pairs.size() < 1)
    {
        return;
    }
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
void shared::py_getDBFromPairs(py::list &pairs, std::string &dbstructure)
{ 
    std::cout << "py_getDBFromPairs" << std::endl;
    if (pairs.size() < 1)
    {
        std::cout << "pairs empty!" << std::endl;
        return;
    }
    //pairs should be BasePair objects or python tuples
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
}
bool shared::check_whitespace(std::string&str)
{
    for (char c : str)
    {
        if (!std::isspace(c))
        {
            return false;
        }
    }
    return true;
}
void shared::py_structureExtract(std::string& glob_sequence, std::string& glob_structure, 
                                 py::list& seq_extract, 
                                 py::list &struc_extract, 
                                 py::list &coords_extract)
{
    /*
        glob_sequence and glob_structure are read in from .dbn
        seq_extract, struc_extract, and coords_extract are lists that hold 
        sequences (str), structures (str), and coordinates (pair of ints) for each
        extracted structure

        loop through global structure, when encountering an open paren start building 
        a sequence and structure from what is encountered
        when all open parens that have been used have been matched w/ a close paren stop
        building that structure and append it to the extract lists
    */
    // keeps track of how many of each open parentheses have been encountered
    // when an open paren is found increase by one
    // when a matching close paren is found decrease by one
    // when all hit 0 the structure is finished
    std::unordered_map<char, int> paren_count = {
        {'(', 0}, {'[', 0}, {'{', 0}, 
        {'<', 0}, {'A', 0}, {'B', 0}, 
        {'C', 0}, {'D', 0}, {'E', 0}, 
        {'F', 0}, {'G', 0}, {'H', 0},
        {'I', 0}, {'J', 0}, {'K', 0}, 
        {'L', 0}, {'M', 0}, {'N', 0},     
        {'O', 0}, {'P', 0}, {'Q', 0}, 
        {'R', 0}, {'S', 0}, {'T', 0}, 
        {'U', 0}, {'V', 0}, {'W', 0}, 
        {'X', 0}, {'Y', 0}, {'Z', 0}
    };
    std::unordered_map<char, char> match_paren = {
        {'(', ')'}, {')', '('}, {'[', ']'}, {']', '['},    
        {'{', '}'}, {'}', '{'}, {'<', '>'}, {'>', '<'},
        {'A', 'a'}, {'a', 'A'}, {'B', 'b'}, {'b', 'B'},
        {'C', 'c'}, {'c', 'C'}, {'D', 'c'}, {'d', 'D'},
        {'E', 'e'}, {'e', 'E'}, {'F', 'f'}, {'f', 'F'},    
        {'G', 'g'}, {'g', 'G'}, {'H', 'h'}, {'h', 'H'},
        {'I', 'i'}, {'i', 'I'}, {'J', 'j'}, {'j', 'J'},
        {'K', 'k'}, {'k', 'K'}, {'L', 'l'}, {'l', 'L'},
        {'M', 'm'}, {'m', 'M'}, {'N', 'n'}, {'n', 'N'},    
        {'O', 'o'}, {'o', 'O'}, {'P', 'p'}, {'p', 'P'},
        {'Q', 'q'}, {'q', 'Q'}, {'R', 'r'}, {'r', 'R'},
        {'S', 's'}, {'s', 'S'}, {'T', 't'}, {'t', 'T'},
        {'U', 'u'}, {'u', 'U'}, {'V', 'v'}, {'v', 'V'},
        {'W', 'w'}, {'w', 'W'}, {'X', 'x'}, {'x', 'X'},
        {'Y', 'y'}, {'y', 'Y'}, {'Z', 'z'}, {'z', 'Z'}
    };
    // keeps track of which parentheses have been used so
    // the entire map doesn't need to be looped through every time
    std::unordered_set<char> used_parens;
    // records if a structure has been encountered and is being recorded
    bool is_structure = false;
    // records the parenthesis at surrent location while looping
    char curr_paren;
    //builds sequence/structure/coords to add to output
    std::vector<char> curr_seq;
    std::vector<char> curr_struc;
    std::tuple<int, int> curr_coords(0, 0); // 1 indexed
    for (int i = 0; i < glob_structure.length(); ++i)
    {
        curr_paren = glob_structure[i];
        if (curr_paren == '\n') {continue;}
        else if (curr_paren == ' ') {continue;}
        else
        {
            if (is_structure)
            {
                if (curr_paren == '.')
                {
                    curr_struc.push_back(curr_paren);
                    curr_seq.push_back(glob_sequence[i]);
                    continue;
                }
                else if (curr_paren == '(' 
                    or curr_paren == '[' 
                    or curr_paren == '{' 
                    or curr_paren == '<' 
                    or (curr_paren >= 'A' and curr_paren <= 'Z'))
                {
                    // record this paren has been encountered
                    paren_count.at(curr_paren) += 1;
                    curr_struc.push_back(curr_paren);
                    curr_seq.push_back(glob_sequence[i]);                    
                    // record if it's a new one
                    if (used_parens.find(curr_paren) == used_parens.end())
                    {
                        used_parens.insert(curr_paren);
                    }
                    continue;
                }
                // check for close paren
                else if (curr_paren == ')' 
                    or curr_paren == ']' 
                    or curr_paren == '}' 
                    or curr_paren == '>' 
                    or (curr_paren >= 'a' and curr_paren <= 'z'))
                {
                    char matched_paren = match_paren[curr_paren];
                    paren_count.at(matched_paren) -= 1;
                    curr_struc.push_back(curr_paren);
                    curr_seq.push_back(glob_sequence[i]);       
                    //check to see if the structure has ended
                    bool finished = true;
                    // loop over used parentheses and check if any aren't 0
                    for (const auto& paren : used_parens)
                    {
                        if (paren_count.at(paren) != 0) {finished = false;break;}
                    }
                    // if none aren't 0, record the structure and set is_structure to false
                    if (finished)
                    {
                        is_structure = false;
                        for (const auto& paren : used_parens) {paren_count.at(paren) = 0;}
                        used_parens.clear();
                        std::string final_seq(curr_seq.begin(), curr_seq.end());
                        std::string final_struc(curr_struc.begin(), curr_struc.end());
                        std::get<1>(curr_coords) = i+1;
                        seq_extract.append(final_seq);
                        struc_extract.append(final_struc);
                        coords_extract.append(curr_coords);
                        curr_seq.clear();
                        curr_struc.clear();
                    }
                    continue;
                }             
                else
                { 
                    // throw an error since this should never happen
                    std::stringstream ss;
                    ss << "Invalid character: " << curr_paren << " in .dbn file at position " << i;
                    std::string except = ss.str();
                    throw shared::Exception(except);
                }
            }
            else
            {
                // no structure being recorded, if an open paren is encountered start one
                // no structure encountered
                if (curr_paren == '.') {continue;}
                // check for close paren and raise error if found
                else if (curr_paren == ')' 
                    or curr_paren == ']' 
                    or curr_paren == '}' 
                    or curr_paren == '>' 
                    or (curr_paren >= 'a' and curr_paren <= 'z'))
                {
                    std::stringstream ss;
                    ss << "Close paren encountered with no open paren in .dbn file at position " << i;
                    std::string except = ss.str();
                    throw shared::Exception(except);
                }
                else if (curr_paren == '(' 
                    or curr_paren == '[' 
                    or curr_paren == '{' 
                    or curr_paren == '<' 
                    or (curr_paren >= 'A' and curr_paren <= 'Z'))
                {
                    is_structure = true;
                    // record this paren has been encountered
                    paren_count.at(curr_paren) += 1;
                    curr_struc.push_back(curr_paren);
                    curr_seq.push_back(glob_sequence[i]);                    
                    used_parens.insert(curr_paren);
                    std::get<0>(curr_coords) = i+1;
                    continue;
                }
                else
                { 
                    // throw an error since this should never happen
                    std::stringstream ss;
                    ss << "Invalid character: " << curr_paren << " in .dbn file at position " << i;
                    std::string except = ss.str();
                    throw shared::Exception(except);
                }
            }
        } 
    }
}
