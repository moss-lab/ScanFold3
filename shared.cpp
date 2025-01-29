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
std::vector<std::pair<int,int>> shared::findPairsInDotBracket(const std::string& dbstructure) 
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