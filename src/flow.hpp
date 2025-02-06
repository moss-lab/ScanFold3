#ifndef FLOW
#define FLOW
/*
functions/classes to define a flow graph and solve for maximum flow and minimum weight
*/
#include <vector>
#include <cstddef>
#include <memory>
#include <limits>
#include <queue>
#include <stdexcept>
#include <cstring>
#include <iostream>
#include <algorithm>
#include "basepair.hpp"
#include "matrix.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace flow 
{
    namespace py = pybind11;
    struct Edge;
    //a Vertex is a vector of pointers to Edges that originate at this vertex
    typedef std::vector<std::shared_ptr<Edge>> Vertex;
    struct Edge 
        : public std::enable_shared_from_this<Edge>
    {
        /*
        Edge stores directed edges with a weight, which can be negative
        they also store whether or not this edge can take flow (available), aka whether the capacity is
        1 or 0
        each Edge has an associated residual, which connects the same vertex but with opposite direction
        and negative weight
        residual_index stores the index of that residual in the vertex this edge ends at
        the residual itself is an Edge; in the residual, residual_index stores the index of this Edge
        */
        double weight;  //getting rid of this and just using the matrix is one optimization, didn't
                        //do it here for sake of simplicity
        size_t start_vertex;    //edge has index of its vertices; if anyone decides to change this
                                //to pointers to vertices make sure to use non-owning pointers
        size_t end_vertex;
        bool available{true};
        size_t residual_index{0};   //index of residual in the vertex it leads to
        //constructors
        Edge(size_t start_vertex, 
             size_t end_vertex, 
             double weight
            );
        Edge(size_t start_vertex, 
             size_t end_vertex, 
             double weight,
             bool available //necessary to create a residual, which starts unable to take flow
            );
        //methods
        std::shared_ptr<Edge> getptr(); //get owning pointer to this edge
        void swap_pairing_status();     //make this edge paired (if unpaired) or unpaired (vice versa)
        void set_residual_index(size_t index);
    };
    struct FlowGraph {
        std::vector<Vertex> vertices;
        long double net_weight;
        size_t source;
        size_t sink;
        size_t sequence_length;
        bool flow_has_been_called{false};

        FlowGraph(matrix::BasePairMatrix &mat);
        //add edge, vertices must already be in graph
        void add_edge(size_t start, size_t end, double weight);
        //minimum weights are stored in the capacities of edges after flow() is called
        //use this to return indices of matches as a vector of pairs
        //rerun will cause flow() to be called again if it already has been
        std::vector<std::pair<size_t, size_t>> minimum_weighted_matching(bool rerun = false);
        //python wrapper for minimum_weighted_matching
        void minimum_weighted_matching(py::list &py_pairs, bool rerun = false);
        //run algorithm to find min flow, retrieve w/ minimum_weight()
        long double shortest_path_faster_flow();
        //find the lowest weight path through the graph currently available
        void find_shortest_path(size_t graph_size, std::vector<double>& cumulative_weight,
                                std::vector<std::shared_ptr<Edge>>& predecessor, size_t start = 0);
        void push_to_search_queue(size_t vertex_index, std::queue<size_t>& search_queue, 
                                  std::vector<bool> on_queue);
        //used to ensure the matching found uses each base once and only once
        void roll_for_sanity_damage(std::vector<std::pair<size_t, size_t>>& matches);
    };
}
#endif