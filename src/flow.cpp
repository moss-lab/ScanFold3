#include "flow.hpp"
using namespace flow;

Edge::Edge(size_t start_vertex, size_t end_vertex, double weight) :
    start_vertex(start_vertex), end_vertex(end_vertex), weight(weight) {}
Edge::Edge(size_t start_vertex, size_t end_vertex, double weight, bool available) : 
    start_vertex(start_vertex), end_vertex(end_vertex), weight(weight), available(available) {}
//bool Edge::is_paired() {return this->paired}
std::shared_ptr<Edge> Edge::getptr()
{
    //get shared pointer to this edge, to add to associated vertices
    return shared_from_this();
}
void Edge::swap_pairing_status()
{
    //std::cout << "availability before pairing status swapped: " << this->available << std::endl;
    //swap pairing status, turning non-residual edge to residual (paired) edge or vice versa
    //aka swap availability (whether it can be used by the algorithm)
    this->available = !this->available;
    //std::cout << "availability after pairing status swapped: " << this->available << std::endl;
}
void Edge::set_residual_index(size_t index)
{
    this->residual_index = index;
}
FlowGraph::FlowGraph(matrix::BasePairMatrix &mat)
{
    /*
    Defines a bipartite flow graph created from a BasePairMatrix
    the graph has a source, a sink, and one vertex for each base on either side of the graph
    the first half and second half of the graph, not counting source/sink, represent the
    two sets of the bipartite graph
    each edge represents a base pair, with a boolean used to track if it's paired
    paired edges act as a residual: the direction is reversed and the weight is made negative
    of what it is
    
    a BasePairMatrix can be used to create a weighted bipartite graph where each vertex has
    at least one edge
    the goal is to find a minimum weighted matching on that graph, under the condition that
    for edge (i,j) vertices i and j should be matched to j and i respectively
    */
    this->sequence_length = mat.getSequenceLength(); //number of bases represented in mat
    this->vertices.resize(2*sequence_length + 2);   //need 2 vertices for each base and another
                                                    //1 each for source and sink
    //indices of source and sink (first/last vertex)
    this->source = 0;
    this->sink = 2*sequence_length + 1;
    //each base should exist twice, once in the 1st half and once in the second half of the
    //graph
    //the one in the 1st half should be connected to the source, the one in the 2nd half should
    //be connected to the sink
    //the vertex in the first half should also be connected to every base that it pairs with
    
    //create an edge to source and sink for each base except the ones in the last window
    //source/sink means the graph is 1-indexed rather than 0, so we need to add 1 to coords
    for(size_t i = 1; i <= sequence_length; ++i)
    {
        this->add_edge(source, i, 0);
        this->add_edge(i + sequence_length, sink, 0);
    }
    //fill in edges of vertices; all are guaranteed at least one due to how ScanFold-Scan works
    for(auto row : mat.Matrix)
    {
        //std::cout << "row read!" << std::endl;
        for(auto pair : row)
        {
            if(pair->icoord < 0)
            {
                //std::cout << "skipped " << pair.icoord << std::endl;
                continue;   //default value, so no data were read for this pair (didn't occur 
                            //across any window)
            }
            //get coords, 1-indexed
            size_t icoord = pair->icoord + 1;
            size_t jcoord = pair->jcoord + 1;
            //add i,j edge (all edges go to other half of graph)
            this->add_edge(icoord, jcoord + sequence_length, pair->getZNorm());
            //add parallel (j,i) edge if i and j are different
            if(icoord != jcoord)
            {
                this->add_edge(jcoord, icoord + sequence_length, pair->getZNorm());
            }
        }
    }
}
void FlowGraph::add_edge(size_t start, size_t end, double weight)
{
    std::ofstream ofile;
    ofile.open("edges.txt", std::ios::app);
    /*
    adds an edge and its residual (backwards edge with negative weight and no capacity)
    note that when edges and residuals are swapped in the algorithm, the way they are referred
    to also changes. So an edge that is used will be called a residual even though it technically isn't
    */
    //dynamically create Edge
    auto edgeptr {std::make_shared<Edge>(start, end, weight)};  
    ofile << "edge: " << start << " " << end << " " << weight << std::endl;
    //do the same w/ its residual
    auto resptr{std::make_shared<Edge>(end, start, -weight, false)};
    ofile << "res: " << end << " " << start << " " << -weight << std::endl;
    size_t edge_index = this->vertices[start].size();
    ofile << "edge_index: " << edge_index << std::endl;
    size_t residual_index = this->vertices[end].size();
    ofile << "residual_index: " << residual_index << std::endl;
    //make sure the edges can find each other
    edgeptr->set_residual_index(residual_index);
    resptr->set_residual_index(edge_index);
    //push copy of pointers to start and end
    this->vertices[start].push_back(edgeptr);
    this->vertices[end].push_back(resptr);
    ofile.close();
}
std::vector<std::pair<size_t, size_t>> FlowGraph::minimum_weighted_matching(bool rerun)
{
    //will return matches based as a pair of 0-indexed coordinates corresponding to the position of
    //the base in the sequence
    std::vector<std::pair<size_t, size_t>> matches;
    //run flow algorithm if haven't already or told to
    if(!this->flow_has_been_called || rerun) 
    {
        this->net_weight = this->shortest_path_faster_flow();
    }
    //loop over the first half of vertices
    for(size_t i = 1; i <= this->sequence_length; ++i)
    {
        //std::cout << std::endl << std::endl;
        //std::cout << "vertex: " << i << std::endl;
        //loop over edges for each vertex
        for(auto &edge : this->vertices[i])
        {
            //std::cout << std::endl;
            //check if edge is used
            if(!edge->available)
            {
                //std::cout << edge->start_vertex << "\t" << edge->end_vertex << std::endl;
                //since vertex 1 and 2 have been swapped by pairing the edge, end_vertex
                //is the vertex on the left half (connected to source) and start_vertex is
                //the one on the right half (connected to sink)
                //skip over connections to sink/source
                if(edge->start_vertex == this->sink
                   || edge->end_vertex == this->source
                  ) {
                    //std::cout << "skipped!" << std::endl;
                    continue;
                    }
                //vertex this edge points towards but on the same side of the graph
                size_t i_pairs_with = edge->end_vertex - this->sequence_length;
                //std::cout << "i: " << i << "\ti in edge: " << edge->start_vertex << std::endl;
                //std::cout << "\ti_pairs_with: " << i_pairs_with << "\tj in edge: " << edge->end_vertex << std::endl;
                //don't count edges twice
                if(i_pairs_with >= i)
                {
                    size_t i_index = i - 1;
                    size_t j_index = i_pairs_with - 1;
                    //std::cout << "i_index: " << i_index << "\tj_index: " << j_index << std::endl;
                    //pairs are 0-indexed like the rest of the program
                    std::pair<size_t, size_t> match(i_index, j_index);
                    matches.push_back(match);
                }
            }
        }
    }
    //this->roll_for_sanity_damage(matches);
    return matches;
}
//python wrapper
void FlowGraph::minimum_weighted_matching(py::list &py_pairs, bool rerun)
{
    std::vector<std::pair<size_t, size_t>> pairs = this->minimum_weighted_matching(rerun);
    for(auto pair : pairs)
    {
        py_pairs.append(pair);
    }
}
void FlowGraph::roll_for_sanity_damage(std::vector<std::pair<size_t, size_t>>& matches)
{
    /*
    check to make sure i,j = j,i for everything in match
    if either true, throw error
    we're also taking this change to make sure every base got matched
    */ 
    double num_vertices_matched = 0.0;
    for(auto pair : matches)
    {
        size_t i_index = pair.first;
        size_t j_index = pair.second;
        auto j_i_pair = matches[j_index];
        if(i_index != j_i_pair.second)
        {
            throw shared::Exception("Inconsistent matching in flow graph!");
        }
        if(i_index == j_index)
        {
            num_vertices_matched++;
        }
        else
        {
            num_vertices_matched += 0.5;
        }
    }
    if(num_vertices_matched < this->sequence_length)
    {
        throw shared::Exception("not all vertices matched!");
    }

}
long double FlowGraph::shortest_path_faster_flow()
{
    /*
    shortest path faster algorithm flow (modified Bellman-Ford)
    overall algorithm:
    1.  create vectors for cumulative weight (distance) initialized to infinity, lowest 
        weight path backwards (predecessor), and if a given vertex is on the queue 
        (on_queue); index in these vectors corresponds to the vertex's index in the graph
    2.  add source to queue
    3.  while something is in the queue: 
        a.  pop vertex v1 from queue and set on_queue[v1] to false
        b.  loop over edges in v1:
            A.  check if edge is outgoing, if not continue to next iteration
            B.  get the destination vertex, v2, from edge
            C.  if v2 has more cumulative weight than v1 + edge weight, 
                aka distance[v2] >  distance[v1] + edge->weight:
                    I.  set distance[v2] to distance[v1] + edge->weight
                    II. set predecessor[v2] to edge
                    III.if on_queue[v2] false:
                        i.  add v2 to queue
                        ii. set on_queue[v2] to true
    4.  for edge in path[sink] to path[source]:
        a.  remove forward edge and create residual edge for the vertex this edge originated from
            aka edge.swap_pairing_status
    5. repeat algorithm until there are no outgoing edges from source
    note that this may terminate early due to negative cycles
    to check for negative cycles, try to find the shortest path again and see if it's shorter
    1.  for edge in edges:
        a.  get v1, v2, and weight from edge
        b.  if distance[v1] isn't infinity, and distance[v1] + weight < distance[v2] then there's
            a negative cycle
    the presence of a negative cycle means we can't frame this as a flow problem and need to use
    a general maximum weighted matching algorithm, like Edmonds' Blossom algorithm

    once the algorithm is finished, final pairings can be found by searching for edges not 
    connected to source/sink where the edge has been used
    these are the residuals created when an edge is used by the algorithm
    */
    //set infinity to be max double for the purposes of setting initial cumulative paths
    const double INF = std::numeric_limits<double>::max();
    size_t graph_size = this->vertices.size();
    size_t source = 0;
    size_t sink = graph_size - 1;

    //cumulative weight along the shortest path to a vertex given by the index found so far
    std::vector<double> cumulative_weight(graph_size, INF);
    //previous edge the shortest path leading to a vertex given by the index
    std::vector<std::shared_ptr<Edge>> predecessor(graph_size); 
    //maximum_flow is used to know when all bases have been paired
    //it's equivalent to one for every 2 vertices, not counting source/sink, where each valid path
    //increases flow by 1
    //we want as close to a maximum matching as we can get here because we want to construct a model 
    //with all data from scanfold-scan, which has at least one predicted pairing per nucleotide 
    //(being unpaired is modeled as self-pairing). Allowing the model to leave any base unpaired would 
    //erroneously ignore any actual metrics attached to it being unpaired
    unsigned int maximum_flow = (graph_size - 2) / 2;
    long double total_weight = 0.0;
    unsigned int total_flow = 0;
    while(true)
    {
        //std::cout << std::endl;
        //find lowest weight path available
        this->find_shortest_path(graph_size, cumulative_weight, predecessor);
        //std::cout << "sink: " << sink << std::endl;
        //std::cout << "cumulative weight: " << cumulative_weight[sink] << std::endl;
        
        if(cumulative_weight[sink] == INF) break;   //no path to sink remains

        //max flow on this path, initialized to highest number possible (number of remaining pairs)
        unsigned int this_path_flow = maximum_flow - total_flow;
        //std::cout << "this_path_flow: " << this_path_flow << std::endl;
        //std::cout << "maximum_flow: " << maximum_flow << std::endl;
        //loop backwards over the found path and ensure it is valid; TODO: remove this after debugging
        for(size_t current_vertex = sink; current_vertex != source;)
        { 
            unsigned int cap = (predecessor[current_vertex]->available ? 1 : 0);
            this_path_flow = std::min(this_path_flow, cap); //should come out to 1
                                                            //unless something went 
                                                            //wrong 
            current_vertex = predecessor[current_vertex]->start_vertex;
        }
        //store that we found a valid path
        total_flow += this_path_flow;
        //std::cout << "total_flow: " << total_flow << std::endl;
        //track the weight from this path
        total_weight += cumulative_weight[sink] * this_path_flow;
        //std::cout << "path: " << std::endl;
        //go backwards along the lowest weight path
        for(size_t current_vertex = sink; current_vertex != source;)
        {
            //std::cout << "current vertex: " << current_vertex << "\t";
            //std::cout << "next vertex: " << predecessor[current_vertex]->start_vertex << std::endl;
            //std::cout << "availability: " << predecessor[current_vertex]->available << std::endl;
            //mark edges used in this path as unavailable
            predecessor[current_vertex]->swap_pairing_status();
            //mark backwards edges along this path as available
            size_t residual_index = predecessor[current_vertex]->residual_index;
            //std::cout << "residual:\nend vertex: " << vertices[current_vertex][residual_index]->end_vertex << std::endl;
            //std::cout << "start vertex: " << vertices[current_vertex][residual_index]->start_vertex << std::endl;
            //std::cout << "available: " << vertices[current_vertex][residual_index]->available << std::endl;
            vertices[current_vertex][residual_index]->swap_pairing_status();
            //move to the previous vertex along this path
            current_vertex = predecessor[current_vertex]->start_vertex;
        }
        //reset vectors; slightly more efficient to do this than make new vectors each iteration
        std::fill(cumulative_weight.begin(), cumulative_weight.end(), INF);
        std::for_each(predecessor.begin(), predecessor.end(), [](std::shared_ptr<Edge>& p){p.reset();});
        //std::cout << std::endl;
    }
    return total_weight;
}
void FlowGraph::find_shortest_path(size_t graph_size, std::vector<double>& cumulative_weight,
                        std::vector<std::shared_ptr<Edge>>& predecessor, size_t start)
{
    //find shortest path given the current state of the graph using shortest path faster alg.
    //cumulative_weight is the weight of the shortest path to the vertex at index
    //predecessor is the last edge on the shortest path leading to the vertex at index
    //will throw an error if it encounters a negative cycle

    //queue for search
    std::queue<size_t> search_queue;
    //whether the vertex at the index is currently on the queue
    std::vector<bool> on_queue(graph_size, false);
    //used to find negative cycles by tracking how many times each vertex has been relaxed
    //if it's still possible to relax a vertex more times then there are vertices in the graph
    //then there's a negative cycle
    std::vector<unsigned int> relaxations(graph_size, 0);
    //set weight at starting point to 0
    cumulative_weight[start] = 0.0;
    //start search
    this->push_to_search_queue(start, search_queue, on_queue);
    while(!search_queue.empty())
    {
        //pop a vertex index from the queue
        size_t vertex = search_queue.front();
        search_queue.pop();
        on_queue[vertex] = false;
        //search through its edges
        for(auto &edge : vertices[vertex])
        {
            size_t end_vertex = edge->end_vertex;
            double path_weight = edge->weight + cumulative_weight[vertex];
            if(edge->available == true && (path_weight < cumulative_weight[end_vertex]))
            {
                cumulative_weight[end_vertex] = path_weight;
                predecessor[end_vertex] = edge;
                //only add to search queue if this edge is part of a path the algorithm is building
                if(!on_queue[end_vertex])
                {
                    this->push_to_search_queue(end_vertex, search_queue, on_queue);
                    ++relaxations[end_vertex];
                    if(relaxations[end_vertex] >= graph_size)
                    {
                        //throw std::runtime_error("negative cycle in graph!");
                    }
                    
                }
            }
        }
    }
}
void FlowGraph::push_to_search_queue(size_t vertex_index, std::queue<size_t>& search_queue, 
                                     std::vector<bool> on_queue)
{
    search_queue.push(vertex_index);
    on_queue[vertex_index] = true;
}