
#include "p3dGraph.h"

#include <limits>
#include <boost/config.hpp>
#include <iostream>
#include <fstream>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
//#include <boost/property_map/property_map.hpp>

#define I2(i,j,N)       ( (j)*(N) + (i) ) 
#define _DIJKSTRA_IN INT_MAX

using namespace boost;
using namespace std;

// This is an example of an exported function.

float p3dShortestPath(float* cost, int dim, int src, int trgt) {
    typedef adjacency_list < listS, vecS, directedS, no_property, property < edge_weight_t, float> > graph_t;
    typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
    typedef std::pair<int, int> Edge;

    graph_t g;

    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++) {
            if (!(abs(cost[I2(i, j, dim)] - (float) 0.0) < 1E-3)) {
                boost::add_edge(i, j, cost[I2(i, j, dim)], g);
                boost::add_edge(j, i, cost[I2(i, j, dim)], g);
                //cout << cost[I2(i,j,dim)] << " ";
            }
        }
    property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
    std::vector<vertex_descriptor> p(num_vertices(g));
    std::vector<float> d(num_vertices(g));

    // Set src:
    vertex_descriptor s = vertex(src, g);

    // Call algorithm:
    dijkstra_shortest_paths(g, s,
            predecessor_map(boost::make_iterator_property_map(p.begin(), get(boost::vertex_index, g))).
            distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, g))));

    // Set destination:
    int dest = vertex(trgt, g);

    // Return the distance of the shortest path:
    return d[dest];
}
