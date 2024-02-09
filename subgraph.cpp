/**
 * @file subgraph.cpp
 * Test script for making a subgraph from our Graph
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
#include <cmath>
#include <fstream>
#include <queue>

/** Construct an induced subgraph of Graph @a g.
 * @param[in] g Initial graph
 * @param[in] node_pred Node predicate to determine nodes in @a subgraph
 * @param[out] subgraph The subgraph of @a g induced by @a node_pred
 *
 * If @a node_pred(n) returns true for a node n of Graph @a g, then a node
 * with that position is added to @a subgraph. All induced edges are edges of
 * @a subgraph. */
template <typename G, typename Predicate>
void induced_subgraph(const G& g, Predicate node_pred, G& subgraph) {
  // HW1 #4: YOUR CODE HERE
  for (auto it = g.node_begin(); it != g.node_end(); ++it){
  	auto node = *it;
  	if (node_pred(node) ) 
		subgraph.add_node(node.position(), node.value());
  }

  // Add edge to subgraph
  for (auto it = g.edge_begin(); it != g.edge_end(); ++it){
  	auto edge = *it;
  	auto node1 = edge.node1();
  	auto node2 = edge.node2();
  	if (node_pred(node1) && node_pred(node2) ){
  		subgraph.add_edge(node1, node2);
  		}
  }
}

/** Test predicate for HW1 #4 */
struct IndexPredicate {
  unsigned last_index;
  IndexPredicate(unsigned last_index)
    : last_index(last_index) {
  }
  template <typename NODE>
  bool operator()(const NODE& n) {
    return n.index() < last_index;
  }
};

// HW1 #4: YOUR CODE HERE
// Specify and write an interesting predicate, useful for induced_subgraph().
// Explain what your predicate is intended to do, and supply a function like
// test_induced_subgraph() to test it. If you'd like you may create new
// nodes and tets files.

/** Test predicate for HW1 #4 
  *
  * % is computationally impossible on large number of edges
  */
struct InterPredicate {
  InterPredicate(){
  }
  bool operator()(Graph<int>::Node n) {
	float ni = (float) n.index();
	return ( ( ni <= 15000 ) || ( 20000 < ni ) );
  }
};

CS207::Color CFunction(Graph<int>::Node n);
/** Test function for HW1 #4 */
template <typename G>
void inter_induced_subgraph(const G& graph, CS207::SDLViewer& viewer) {
  G subgraph;
  InterPredicate ip;
  induced_subgraph(graph, ip, subgraph);

  viewer.clear();
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(subgraph.node_begin(), subgraph.node_end(), CFunction, node_map);
  viewer.add_edges(subgraph.edge_begin(), subgraph.edge_end(), node_map);
  viewer.center_view();
}


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
//	return CS207::Color::make_heat( n.value()/d );	// woohoo! 
//  return CS207::Color( 1,0.5,0.2);
	return CS207::Color::make_hsv( pow( (1 - n.value()/d), 2), 0.5, 1);
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
 
  d = (float) shortest_path_lengths(graph, graph.node(0));
  
  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();
  
  // Test your induced subgraph function
  inter_induced_subgraph(graph, viewer);

  return 0;
}
