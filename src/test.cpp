#include "flow.hpp"
#include "matrix.hpp"
#include "shared.hpp"
#include "basepair.hpp"
#include "window.hpp"
#include <chrono>
#include <filesystem>

template <
    class result_t   = std::chrono::milliseconds,
    class clock_t    = std::chrono::steady_clock,
    class duration_t = std::chrono::milliseconds
>
auto since(std::chrono::time_point<clock_t, duration_t> const& start)
{
    return std::chrono::duration_cast<result_t>(clock_t::now() - start);
}

int main(int argc, char *argv[])
{
    namespace fs = std::filesystem;
    std::string dir_path(argv[1]);
    for (const auto & entry : fs::directory_iterator(dir_path))
    {
        std::ifstream scanfold_scan(entry.path()); 
        std::cout << "reading " << entry.path() << ":" << std::endl;
        auto start = std::chrono::steady_clock::now();
        std::vector<window::ScanFoldWindow> scan_data = window::readScanTSV(scanfold_scan); 
        std::cout << "done! Elapsed time (ms): " << since(start).count() << std::endl << std::endl;
        //BasePairMatrix
        //convert scan_data to BasePairMatrix:
        std::string argv1(argv[1]);
        start = std::chrono::steady_clock::now();
        std::cout << "creating BasePairMatrix..." << std::endl;
        auto bpmatrix = matrix::BasePairMatrix(scan_data);
        std::cout << "done! Elapsed time (ms): " << since(start).count() << std::endl << std::endl;

        std::ofstream ffile;
        std::cout << "creating FlowGraph..." << std::endl;
        start = std::chrono::steady_clock::now();
        flow::FlowGraph fgraph(bpmatrix);
        std::cout << "done! Elapsed time (ms): " << since(start).count() << std::endl << std::endl;

        std::cout << "performing weighted matching via FlowGraph..." << std::endl;
        start = std::chrono::steady_clock::now(); 
        auto flowpairs = fgraph.minimum_weighted_matching();
        std::cout << "done! Elapsed time (ms): " << since(start).count() << std::endl;
        
        ffile.open("flow_matching.txt", std::ios::app);
        ffile << entry.path() << std::endl;
        double znorm;
        for(auto pair : flowpairs)
        {
            znorm = bpmatrix.getZNorm(pair.first, pair.second);
            ffile << pair.first << "\t" << pair.second << "\t" << znorm << std::endl;
        }
        ffile << std::endl;
    }
    return 0;
}