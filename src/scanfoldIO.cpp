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
            auto icoord = pair.attr("i_coord")+1 - start_coordinate;
            auto jcoord = pair.attr("j_coord")+1 - start_coordinate;
            auto inuc = pair.attr("inuc");
            auto jnuc = pair.attr("jnuc");
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
            py::object py_znorm = pair.attr("getZNorm")();
            auto znorm = py_znorm.cast<double>();
            if (znorm > filter)
            {
                continue;
            }
            auto icoord = end_coordinate+1 - pair.attr("i_coord");
            auto jcoord = end_coordinate+1 - pair.attr("j_coord");
            auto inuc = pair.attr("inuc");
            auto jnuc = pair.attr("jnuc");
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
    //create a dbn file from a list of BasePair, a sequence, a header, and a filename
    //only callable in Python
    std::string dbstructure = shared::py_getDBFromPairs(pairs);
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