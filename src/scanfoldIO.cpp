#include "scanfoldIO.hpp"
using namespace io;

void io::list_to_ct(std::vector<basepair::BasePair>& base_pair_list, std::string out_file, double filter,
        char strand, std::string name, size_t start_coordinate, size_t end_coordinate)
{
    /*
    takes a sorted vector of BasePair (low to high), writes it to a .ct file
    inputs:
        base_pair_list: sorted list of BasePair
        out_file: name of .ct file
        filter: number z-scores have to be lower than to be written
        strand: 1 for positive, -1 for negative
        start_coordinate: start coordinate of scan (0 indexed)
        end_coordinate: end coordinate of scan (0 indexed)
    outputs:
        writes a .ct file, starts at 1 regardless of genomic coordinates
    */
    std::ofstream ofile(out_file);
    ofile << base_pair_list.size() << '/t' << name << std::endl;
    if (strand == 1)
    {
        for (auto pair : base_pair_list)
        {
            if (pair.getZNorm() > filter)
            {
                continue;
            }
            auto icoord = pair.icoord+1 - start_coordinate;
            auto jcoord = pair.jcoord+1 - start_coordinate;
            auto inuc = pair.inuc;
            auto jnuc = pair.jnuc;
            if (icoord > jcoord)
            {
                std::swap(jcoord, jcoord);
                std::swap(inuc, jnuc);
            }
            else if (icoord == jcoord)
            {
                jcoord = 0;
            }
            //TODO: remove redundant base pairs
            ofile << icoord << ' ' << inuc << ' ' << icoord-1 << ' ' << icoord+1 << ' ' << jcoord << icoord;
        }
    }
    else if (strand == -1)
    {
        //iterate backwards
        //base coordinates on distance from end rather than distance from start
        for (auto it = base_pair_list.rbegin(); it != base_pair_list.rend(); ++it)
        {
            auto pair = *it;
            if (pair.getZNorm() > filter)
            {
                continue;
            }
            auto icoord = end_coordinate+1 - pair.icoord;
            auto jcoord = end_coordinate+1 - pair.jcoord;
            auto inuc = pair.inuc;
            auto jnuc = pair.jnuc;
            if (icoord > jcoord)
            {
                std::swap(icoord, jcoord);
                std::swap(inuc, jnuc);
            }
            else if (icoord == jcoord)
            {
                jcoord = 0;
            }
            //TODO: remove redundant base pairs
            ofile << icoord << ' ' << inuc << ' ' << icoord-1 << ' ' << icoord+1 << ' ' << jcoord << icoord;
        }
    }
}
void io::py_list_to_ct(py::list& base_pair_list, std::string out_file, double filter,
        char strand, std::string name, size_t start_coordinate, size_t end_coordinate)
{
    /*
    takes a sorted python list of BasePair (low to high), writes it to a .ct file
    inputs:
        base_pair_list: sorted list of BasePair
        out_file: name of .ct file
        filter: number z-scores have to be lower than to be written
        strand: 1 for positive, -1 for negative
        start_coordinate: start coordinate of scan (0 indexed)
        end_coordinate: end coordinate of scan (0 indexed)
    outputs:
        writes a .ct file, starts at 1 regardless of genomic coordinates
    */
    std::ofstream ofile(out_file);
    ofile << base_pair_list.size() << '/t' << name << std::endl;
    if (strand == 1)
    {
        for (auto pair : base_pair_list)
        {
            //py::object znormMethod = pair.attr("getZNorm");
            //py::object result = znormMethod();
            //auto znorm = result.cast<double>();
            py::object py_znorm = pair.attr("getZNorm")();
            auto znorm = py_znorm.cast<double>();
            if (znorm > filter)
            {
                continue;
            }
            auto icoord = pair.attr("i_coord").cast<size_t>()+1 - start_coordinate;
            auto jcoord = pair.attr("j_coord").cast<size_t>()+1 - start_coordinate;
            auto inuc = pair.attr("inuc").cast<char>();
            auto jnuc = pair.attr("jnuc").cast<char>();
            if (icoord > jcoord)
            {
                std::swap(jcoord, jcoord);
                std::swap(inuc, jnuc);
            }
            else if (icoord == jcoord)
            {
                jcoord = 0;
            }
            //TODO: remove redundant base pairs
            ofile << icoord << ' ' << inuc << ' ' << icoord-1 << ' ' << icoord+1 << ' ' << jcoord << icoord;
        }
    }
    else if (strand == -1)
    {
        //iterate backwards
        //base coordinates on distance from end rather than distance from start
        for (auto pair : py::reinterpret_borrow<py::iterable>(base_pair_list).attr("__reversed__")())
        {
            py::object py_znorm = pair.attr("getZNorm")();
            auto znorm = py_znorm.cast<double>();
            if (znorm > filter)
            {
                continue;
            }
            auto icoord = end_coordinate+1 - pair.attr("i_coord").cast<size_t>();
            auto jcoord = end_coordinate+1 - pair.attr("j_coord").cast<size_t>();
            auto inuc = pair.attr("inuc").cast<char>();
            auto jnuc = pair.attr("jnuc").cast<char>();
            if (icoord > jcoord)
            {
                std::swap(icoord, jcoord);
                std::swap(inuc, jnuc);
            }
            else if (icoord == jcoord)
            {
                jcoord = 0;
            }
            //TODO: remove redundant base pairs
            ofile << icoord << ' ' << inuc << ' ' << icoord-1 << ' ' << icoord+1 << ' ' << jcoord << icoord;
        }
    } 
}
void io::py_makedbn(py::list& pairs, std::string sequence, std::string header, std::string file_name)
{
    std::cout << "making dbn file: " << file_name << std::endl;
    std::cout << "num pairs: " << pairs.size() << std::endl;
    //create a dbn file from a list of BasePair, a sequence, a header, and a filename
    //only callable in Python
    size_t len_seq = sequence.size();
    std::string dbstructure(len_seq, '.');
    shared::py_getDBFromPairs(pairs, dbstructure);
    std::ofstream ofile(file_name);
    if (ofile.is_open())
    {
        ofile << header << std::endl;
        ofile << sequence << std::endl;
        ofile << dbstructure;
        ofile.close();
    }
    else
    {
        std::string err_string = "Error: could not open " + file_name;
        throw shared::Exception(err_string);
    }
}
void io::dbn_to_ct(std::string in_file, std::string out_file)
{
    /*
    takes a .dbn file with a structure to turn into a .ct file
    .dbn file header should be of format:
        >NAME_motif_NUM_coordinates:START-END
    where start/end coordinates are 1-indexed
    */
    std::ifstream dbnfile(in_file);
    std::ofstream ofile (out_file);
    std::string header;
    std::string sequence;
    std::string structure;
    if (!dbnfile.is_open()) 
    {
        std::string err_string = "Error: Unable to open motif dbn file: " + in_file;
        throw shared::Exception(err_string);
    }
    //read dbn file data
    std::getline(dbnfile, header);
    std::getline(dbnfile, sequence);
    std::getline(dbnfile, structure);
    dbnfile.close();
    //get list of pairs
    std::vector<std::pair<size_t,size_t>> pairs = shared::findPairsInDotBracket(structure);
    //get pairs not stored by findPairsInDotBracket and add them
    std::vector<std::pair<size_t,size_t>> j_i_pairs;
    for (const auto& pair : pairs)
    {
        if (pair.first < pair.second)
        {
            std::pair<size_t,size_t> new_pair = {pair.second, pair.first};
            j_i_pairs.push_back(new_pair);
        }
    }
    pairs.insert(
            pairs.end(),
            std::make_move_iterator(j_i_pairs.begin()),
            std::make_move_iterator(j_i_pairs.end())
            );
    //dont use j_i_pairs past this point
    //
    //sort pairs, defaults to first item in pairs
    std::sort(pairs.begin(), pairs.end());
    for (const auto& pair : pairs)
    {
        auto icoord = pair.first+1;     //1 indexed
        auto jcoord = pair.second+1;
        char inuc = sequence[pair.first];
        char jnuc = sequence[pair.second];
        if (icoord == jcoord)
        {
            jcoord = 0;
        }
        ofile << icoord << ' ' << inuc << ' ' << icoord-1 << ' ' << icoord+1 << ' ' << jcoord << ' ' << icoord << std::endl;
    }
    /*
    else if (strand == -1)
    {
       //find end coordinate from header
       std::size_t last_spacer = header.rfind('-');
       std::string end_coord_str = header.substr(last_spacer);
       int end_coord = stoi(end_coord_str);

       for(auto it = pairs.rbegin(); it != pairs.rend(); ++it)
       {
          auto pair = *it;
          auto icoord = end_coord+1 - pair.first;
          auto jcoord = end_coord+1 - pair.second;
          auto inuc = sequence[pair.first];
          auto jnuc = sequence[pair.second];
    */
    //get info from header
    /*
    int spacercount = 0;
    std::string name;
    std::string start_coord_str;
    std::string end_coord_str;
    name.reserve(100);
    start_coord_str.reserve(100);
    end_coord_str.reserve(100);
    for (char c : header)
    {
        if (c == '>' or c == '\n' or c == ' ') {continue;}
        if (c == '_')
        {
            ++spacercount;
            continue;
        }
        if (spacercount == 0)
        {
            name += c;
        }
        else if (spacercount == 3)
        {
            if (c == ':')
            {
                ++spacercount;
                continue;
            }
        }
        else if (spacercount == 4)
        {
            if (c == '-')
            {
                ++spacercount;
                continue;
            }
            start_coord_str += c;
        }
        else if (spacercount == 5)
        {
            end_coord_str += c
        }        
    }
    int start_coord = std::stoi(start_coord_str);
    int end_coord = std::stoi(end_coord_str);
    */
}

