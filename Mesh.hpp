#ifndef CS207_MESH_HPP
#define CS207_MESH_HPP

/**
 * @file Mesh.hpp
 */

#include "CS207/Util.hpp"
#include "Point.hpp"
#include "Graph.hpp"
#include <utility>
#include <algorithm>
#include <vector>
#include <assert.h>
#include <math.h>

using namespace std::rel_ops;

template<typename ND, typename ED, typename TD>
class Mesh {

	public: 
	  /** Type of indexes and sizes. Return type of Node::index() and
      Graph::num_nodes(), argument type of Graph::node. */
  typedef unsigned size_type;
  
  typedef ND node_data_type;
  typedef ED edge_data_type;
  typedef TD tri_data_type;
  
  struct mesh_tri_data;
	struct mesh_node_data;
	struct mesh_edge_data;
  
	// Define Graph type
	typedef Graph<mesh_node_data, mesh_edge_data> MeshGraphType;
	typedef typename MeshGraphType::Node MeshNode;
	typedef typename MeshGraphType::Edge MeshEdge;
	typedef typename MeshGraphType::node_value_type mesh_node_value_type;
	typedef typename MeshGraphType::edge_value_type mesh_edge_value_type;
	typedef typename MeshGraphType::node_iterator mesh_node_iterator;
	typedef typename MeshGraphType::edge_iterator mesh_edge_iterator;
	typedef typename MeshGraphType::incident_iterator mesh_incident_iterator;
			
	// Define Triangle Graph type
	typedef Graph<mesh_tri_data, int> TriGraphType;
	typedef typename TriGraphType::Node TriNode;
	typedef typename TriGraphType::Edge TriEdge;
  typedef typename TriGraphType::node_value_type tri_node_value_type;
	typedef typename TriGraphType::edge_value_type tri_edge_value_type;
	typedef typename TriGraphType::node_iterator tri_node_iterator;
	typedef typename TriGraphType::edge_iterator tri_edge_iterator;
	typedef typename TriGraphType::incident_iterator tri_incident_iterator;

	typedef std::vector<MeshNode> mesh_node_type;
	typedef std::vector<TriNode> tri_node_type;

	struct mesh_tri_data{
		tri_data_type tri_data;
	  mesh_node_type three_mesh_node;  // 3 mesh graph nodes of a triangle
//  	tri_node_type tri_neigh;	// vector of 0-3 neighboring triangles	
		// vector of outward normals, assume that tri edges are initialized 0-3 neighboring triangles	
	};
/*
	class triangle{
		Mesh* m_;
		tri_data_type tri_data_;
		mesh_node_type three_mesh_node_;  // 3 mesh graph nodes of a triangle
		
		triangle(const Mesh* m) 
    : m_(const_cast<Mesh*>(m)){   
    }
	
		tri_data_type& value(){
			return m_->tri_graph
		}
	};
*/
	struct mesh_node_data{
  	node_data_type node_data;
  	tri_node_type many_tri_node;	// any number of triangles adjacent to a mesh graph node
	};
	
	struct mesh_edge_data {
  	edge_data_type edge_data;
  	tri_node_type two_tri_node;	 // 1-2 triangles adjacent to a mesh graph edge
	};

  /** Return the number of nodes in Mesh. */
  size_type num_nodes() const {
    return mesh_graph.num_nodes();
  }

  /** Return the number of edges in Mesh. */
  size_type num_edges() const {
    return mesh_graph.num_edges();
  }

  /** Return the number of edges in Mesh. */
  size_type num_triangles() const {
    return tri_graph.num_nodes();
  }

 /** Return a node from Mesh
   *
   * Complexity: O(1) 
   */	
	MeshNode node(size_type i) const{
		return mesh_graph.node(i);
	}

 /** Add a node to Mesh
   *
   * Complexity: O(1) 
   */	
	void add_node(Point& p){
		mesh_graph.add_node(p);
//		std::cerr << "q " << node.value().q.h << "\n";
	}			

 /** Return an edge proxy 
   *
   * Complexity: O(1) 
   */	
	MeshEdge add_edge(MeshNode& n1, MeshNode& n2){
		return mesh_graph.add_edge(n1, n2);
//		std::cerr << "q " << node.value().q.h << "\n";
	}			

		
 /** Add a node to Dual Graph
	 * @param[in] n1, n2, n3 are proxies of Mesh nodes, n1 != n2 != n3
   * @pre at least one n1, n2, n3 does not exist; triangle does not exist
   * 
   * @post new tri_graph.num_nodes() = old tri_graph.num_nodes() + 1;
   * @post every mesh_graph edge is part of a triangle
   * @post normal of a triangle is defined only if both triangles exist
   *
   * Add a node to Dual Graph. Then add each edge to Primal Graph. For each 
   * edge in Primal Graph, if edge does not exist, add edge. Otherwise, add edges
   * to Dual Graph.
   *
   * RI : either 1 or 2 triangles for each edge
   *
   * Complexity: O(degree of an node in Primal Graph) 
   */	
	TriNode add_triangle(MeshNode n1, MeshNode n2, MeshNode n3){
      
    // Initialize TriData to Tri graph       
      TriNode t1 = tri_graph.add_node(Point(0,0,0));

      t1.value().three_mesh_node.push_back(n1);
      t1.value().three_mesh_node.push_back(n2);
      t1.value().three_mesh_node.push_back(n3);  
      
//      std::cerr<< tri_node.value().Q_bar.h << "\n";
  
      //Update NodeData in Mesh graph
      n1.value().many_tri_node.push_back(t1);
			n2.value().many_tri_node.push_back(t1);
			n3.value().many_tri_node.push_back(t1);
/*
			
			//Update EdgeData in Mesh graph and Edges in Tri graph
			if (mesh_graph_.has_edge(n2, n3)) {
				auto e = mesh_graph_.add_edge(n2, n3);
				e.value().tri_node_2 = t1;
				
				auto t2 = e.value().tri_node_1;
				auto tri_e = tri_graph_.add_edge(t1, t2);

				
				TriEdgeData ted;
				ted.outward_normal = normal(n1,n2,n3);
				ted.tri_node = t2;
				t1.value().tri_edge_data.push_back(ted);
				
			} else{
				auto e = mesh_graph_.add_edge(n2, n3);
				e.value().tri_node_1 = t1;
			}	
*/			
			//Update EdgeData in Mesh graph and Edges in Tri graph
			// mesh edge e has a tri node
			// triangle node has 0-3 neighbors

			if (mesh_graph.has_edge(n2, n3)) {
				auto e = mesh_graph.add_edge(n2, n3);
				e.value().two_tri_node.push_back(t1);

				e = mesh_graph.add_edge(n3, n2);
				e.value().two_tri_node.push_back(t1);
												
				auto t2 = e.value().two_tri_node.at(0);
				tri_graph.add_edge(t1, t2);
				
				
			} else{
				auto e = mesh_graph.add_edge(n2, n3);
				e.value().two_tri_node.push_back(t1);

				e = mesh_graph.add_edge(n3, n2);
				e.value().two_tri_node.push_back(t1);				
			}	
			
			if (mesh_graph.has_edge(n3, n1)) {
				auto e = mesh_graph.add_edge(n3, n1);
				e.value().two_tri_node.push_back(t1);

				e = mesh_graph.add_edge(n1, n3);
				e.value().two_tri_node.push_back(t1);
												
				auto t2 = e.value().two_tri_node.at(0);
				tri_graph.add_edge(t1, t2);
				
				
			} else{
				auto e = mesh_graph.add_edge(n3, n1);
				e.value().two_tri_node.push_back(t1);

				e = mesh_graph.add_edge(n1, n3);
				e.value().two_tri_node.push_back(t1);				
			}	

			if (mesh_graph.has_edge(n1, n2)) {
				auto e = mesh_graph.add_edge(n1, n2);
				e.value().two_tri_node.push_back(t1);

				e = mesh_graph.add_edge(n2, n1);
				e.value().two_tri_node.push_back(t1);
												
				auto t2 = e.value().two_tri_node.at(0);
				tri_graph.add_edge(t1, t2);
				
				
			} else{
				auto e = mesh_graph.add_edge(n1, n2);
				e.value().two_tri_node.push_back(t1);

				e = mesh_graph.add_edge(n2, n1);
				e.value().two_tri_node.push_back(t1);				
			}	
			
      //Update EdgeData, NodeData, TriData for each edge of the triangle
      n1.value().many_tri_node.push_back(t1);
			n2.value().many_tri_node.push_back(t1);
			n3.value().many_tri_node.push_back(t1);
			
			return t1;
	}	
	
 /** Return an area from a tri Mesh: 2D
   *
   * Complexity: O(1) 
   */		
	double area(MeshNode& n1, MeshNode& n2, MeshNode& n3){
		Point p1 = n1.position();
		Point p2 = n2.position();
		Point p3 = n3.position();
		
		auto v = p1 - p2;
		auto w = p3 - p2;
		
		auto u = v.x * w.y - v.y * w.x;
		
		return sqrt(u*u)/2.0;
	}

 /** Return a 2D outward normal from a tri Mesh
   * @pre n1, n2, n3 are nodes of a triangle
   *
   * Complexity: O(1) 
   *
   * Return a normal pointing away from n1 
   */			
	Point normal(MeshNode& n1, MeshNode& n2, MeshNode& n3){
		auto p1 = n1.position();
		auto p2 = n2.position();
		auto p3 = n3.position();
		
		auto v = p2 - p1;
		auto w = p2 - p3;
		
		// mind only x, y components
		auto w_mag = sqrt ( w.x * w.x + w.y * w.y );
		auto v_dot_w = v.x * w.x + v.y * w.y;	
		
		auto v_1_mag = v_dot_w/w_mag;
		
	
		auto v_2 = v - w/w_mag * v_1_mag;
		auto v_2_mag = sqrt( v_2.x * v_2.x + v_2.y * v_2.y );
		
		auto e = n2.position() - n3.position(); 
		auto edge_length = sqrt( e.x * e.x + e.y * e.y );
		
		return Point( v_2.x/v_2_mag, v_2.y/v_2_mag, 0)*edge_length; 
	}

 /** Return an iterator
   *
   * Complexity: O(1) 
   */		
	tri_node_iterator triangle_begin() const{
		return tri_graph.node_begin();
	}

	tri_node_iterator triangle_end() const{
		return tri_graph.node_end();
	}

	mesh_node_iterator node_begin() const{
		return mesh_graph.node_begin();
	}

	mesh_node_iterator node_end() const{
		return mesh_graph.node_end();
	}	

	mesh_edge_iterator edge_begin() const{
		return mesh_graph.edge_begin();
	}

	mesh_edge_iterator edge_end() const{
		return mesh_graph.edge_end();
	}	

//	private: 
		//create a Mesh graph
		MeshGraphType mesh_graph;
		//create a Tri graph
		TriGraphType tri_graph;

}; // class Mesh
#endif

