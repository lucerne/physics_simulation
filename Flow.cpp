/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"
#include "Point.hpp"
#include "Simplex.hpp"
#include "Force.hpp"
#include "Constraint.hpp"
#include "Mesh.hpp"

#include <fstream>
#include <cmath>
#include <assert.h>

// HW2 #1 YOUR CODE HERE
// Define your Graph type; this is a placeholder.

template <typename M, typename V>
	struct node_pair {
		M mass;
		V velocity;
};

template <typename K, typename L>
	struct edge_pair {
		K K_ij;
		L L_ij;
};


typedef Graph<std::pair<double, Point>, std::vector<std::pair<double, double>>> GraphType;
typedef GraphType::Node GraphNode;
typedef GraphType::Edge GraphEdge;

typedef Mesh<std::pair<double, Point>, std::vector<std::pair<double, double>>> MeshType;
typedef MeshType::Node MeshNode;
typedef MeshType::Edge MeshEdge;

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g Graph
 * @param[in] t The current time (useful for time-dependent forces)
 * @param[in] dt The time step
 * @param[in] force Function object defining the force per node
 * @pre G::node_value_type supports ???????? YOU CHOOSE
 * @return the next time step (usually @a t + @a dt)
 *
 * @a force is called as @a force(n, @a t), where n is a node of the graph
 * and @a t is the current time parameter. @a force must return a Point
 * representing the node's force at time @a t.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F F_k, C constr) {
  
  // Update position with positions at t = n
  for(auto it = g.node_begin(); it != g.node_end(); ++it){
  	auto node = (*it);
  	auto pos = node.position();
  	
  	// Update Qn for each node at t = n;
  	double Tn = 0;
  	
  	for (auto it = g.edge_begin(); it != g.edge_end(); ++it){
  		auto edge = (*it).edge;
  		auto node2 = (*it).node2;
  		auto T_k = node2.area(); 
  		auto Qbar_k = node2.value();
  		
  		Tn = Tn + T_k;
  		Qn = T_k * Qbar_k;
  	}
  	Qn = Qn/Tn;
  
  // Impose constraint and reset node position and velocity at t = n+1
  for(auto it = g.node_begin(); it != g.node_end(); ++it){
  	auto node = (*it);
  	constr(node, t);
  }
  	
  // Update Qbar_k at t = n+1
  for(auto it = g.node_begin(); it != g.node_end(); ++it){	
    auto node = (*it);
    
  	// Retrieve flow
  	auto f_k = edge.flow();
  	auto T_k = node.area();
  	  
	  // Calculate and update the next Qn. 
	  Point Qbar_next = Qbar_ - vel * sum(f_k) / T_k.mag() ;
  }
  return t + dt;
}


int main(int argc, char* argv[]) {
  // check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  GraphType graph;

  // Read all Points and add them to the Graph
  std::ifstream nodes_file(argv[1]);
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    graph.add_node(p);

  // Read all mesh squares and add their edges to the Graph
  std::ifstream tets_file(argv[2]);
  Tetrahedron t;		// Reuse this type
  while (CS207::getline_parsed(tets_file, t)) {
    if (t.n[0] < graph.size() && t.n[1] < graph.size()
	&& t.n[2] < graph.size() && t.n[3] < graph.size()) {
      graph.add_edge(graph.node(t.n[0]), graph.node(t.n[1]));
      graph.add_edge(graph.node(t.n[0]), graph.node(t.n[2]));
//#if 0
      // Diagonal edges: include as of HW2 #2
      graph.add_edge(graph.node(t.n[0]), graph.node(t.n[3]));
      graph.add_edge(graph.node(t.n[1]), graph.node(t.n[2]));
//#endif
      graph.add_edge(graph.node(t.n[1]), graph.node(t.n[3]));
      graph.add_edge(graph.node(t.n[2]), graph.node(t.n[3]));
    }
  }

  // Mesh
  MeshType mesh;
  
  // Create nodes
  
  // Create edges of the mesh 
  for (eit = graph.edge_begin(); eit != graph.edge_end(); ++eit){
  	auto edge = (*eit);
  	auto node1 = edge.node1();
  	auto node2 = edge.node2();
  	
  	for (m_nit = mesh.node_begin(); m_nit != mesh.node_end(); ++m_nit){
  		auto m_node = (*m_nit);
  		if (m_node.has_triangle(node1, node2)) mesh.add_edge(
    
  // Create edges of the mesh 
  for (nit1 = mesh.node_begin(); nit1 != mesh.node_end(); ++nit1){
  	auto m_node1 = (*nit1);
  	
  	for (nit2 = mesh.node_begin(); nit2 != mesh.node_end(); ++nit2){
  		auto m_node2 = (*nit2);
  		// m_node1 is adjacent to m_node2 if they share a pair of graph node indices
  		if (m_node1.adjacent_to(m_node2)) {
  			mesh.add_edge(m_node1, m_node2);
  		}	
  	}
  }
  
  // Initialize node value, edge value
  
  
  // Set initial conditions for your nodes if necessary.
  m = 1.0/graph.num_nodes();
  c = 1.0/graph.num_nodes();
  
  // Set node mass and velocity 
 
  // Forces to use!							
  									
 
  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // Begin shallow water simulation
  double dt = 0.001;
  double t_start = 0;
  double t_end = 3.0;

  for (double t = t_start; t < t_end; t += dt) {
  
		symp_euler_step(graph, mesh, t, dt, f, c);
		
    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.set_label(t);
  }

  return 0;
}
