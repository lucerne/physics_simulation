/**
 * @file shortest_path.cpp
 * Test script for using our templated Graph to determine shortest paths.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Simplex.hpp"
#include "Graph.hpp"

#include <vector>
#include <queue>
#include <fstream>
#include <cmath>


/** Calculate shortest path lengths in @a g from @a root.
 * @param[in,out] g Input graph
 * @param[in] root Root node
 * @return Maximum path length in the Graph
 *
 * Changes each node's value() to the length of the shortest path to that node
 * from @a root. @a root's value() is 0. Nodes unreachable from @a root should
 * have value() -1. */
int shortest_path_lengths(Graph<int>& g, Graph<int>::Node root) {
  // HW1 #4: YOUR CODE HERE
  
  int d = 0;
  std::queue<unsigned> q;
  // root is the oldest element in queue
  q.push(root.index());  
	
  /** assign distance to value() and return largest distance
    * LI : 
    * for 0 < i < j < queue.size(), node.value() of i <=  node.value() of j
    *  each node index is stored in queue at most once
    *
    * easier to set d using breath-first search
    */
  while (!q.empty()){
  	// exam the next in queue 
  	auto node = g.node( q.front() );  
  	q.pop();

	// maximum distance is found if queue is empty
	if (q.empty()) d = node.value();
  
  	// for each adjacent vertex of node 
  	for(auto it = node.edge_begin(); it != node.edge_end(); ++it){
  		auto node2 = (*it).node2();
  		
  		// if node is not marked
		if ( (node2.value() == 0) && (node2 != root)) {
  			node2.value() = node.value() + 1;
   			q.push(node2.index());  
  		}
	}
  }
  
  // set unreached node.value() to -1 
  for(auto it = g.node_begin(); it != g.node_end(); ++it){
  	auto node = *it;
  	if ( (node != root) && (node.value() == 0)) node.value() = -1;
  }
  return d;
}

float d;
CS207::Color CFunction(Graph<int>::Node n){
	auto nv = n.value();
//	return CS207::Color::make_heat( n.value()/d );	// woohoo! 
//    return CS207::Color( 1,0.5,0.2);
	return CS207::Color::make_hsv( (n.value()/d), 0.5, 1);
}


int main(int argc, char* argv[]) {
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  Graph<int> graph;

  // Read all Points and add them to the Graph
  std::ifstream nodes_file(argv[1]);
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    graph.add_node(p);

  // Read all Tetrahedra and add the edges to the Graph
  std::ifstream tets_file(argv[2]);
  Tetrahedron t;
  while (CS207::getline_parsed(tets_file, t)) {
    for (unsigned i = 1; i < Tetrahedron::size; ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(graph.node(t.n[i]), graph.node(t.n[j]));
  }

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;	
  
  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  // HW1 #4: YOUR CODE HERE
  // Use shortest_path_lengths to set the node values to the path lengths
  // Construct a Color functor and view with the SDLViewer
  d = (float) shortest_path_lengths(graph, graph.node(0));
  
  viewer.clear();
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), CFunction, node_map); 
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();
  
  return 0;
}
