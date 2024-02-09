/**
 * @file viewer.cpp
 * Test script for the SDLViewer and Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point list
 *
 * Prints
 * A B
 * where A = number of nodes
 *       B = number of edges
 * and launches an SDLViewer to visualize the system
 */

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"

#include "Graph.hpp"

#include <fstream>

/** Quickie class to represent and read Tetrahedra
 */
struct Tetrahedron
{
  static constexpr unsigned num_points = 4;
  int n[num_points];

  /** Called to read in a Tetrahedron */
  friend std::istream& operator>>(std::istream& s, Tetrahedron& t) {
    for (unsigned i = 0; i < num_points; ++i)
      s >> t.n[i];
    return s;
  }
};

using namespace std;


int main(int argc, char* argv[])
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: viewer NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  Graph<int> graph;

  // Read all Points and add them to the Graph
  std::ifstream nodes_file(argv[1]);
  Point p;
  while (CS207::getline_parsed(nodes_file, p)){
    graph.add_node(p);
	}

  // Read all Tetrahedra and add the edges to the Graph
  std::ifstream tets_file(argv[2]);
  typedef Graph<int>::node_type Node;
  Tetrahedron t;
  while (CS207::getline_parsed(tets_file, t)) {
    for (unsigned i = 1; i < Tetrahedron::num_points; ++i)
      for (unsigned j = 0; j < i; ++j){
//  		std::cerr << "CHECK : add_edge p\n";  
        graph.add_edge(graph.node(t.n[i]), graph.node(t.n[j]));
    	}
  }

  // Print number of nodes and edges
//  std::cerr << "CHECK : num_edges \n"; 
  cout << graph.num_nodes() << " " << graph.num_edges() << endl;

  // Launch a viewer
  CS207::SDLViewer viewer;
  viewer.launch();
/*
  // Set the viewer
  viewer.draw_graph(graph);
  viewer.center_view();
*/  
  // Draw nodes
/*
  auto node_map = viewer.empty_node_map(graph);
  auto first = graph.node_begin();
  auto last = graph.node_end();
  auto node = *first; 
*/  
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map); 
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();

  return 0;
}


