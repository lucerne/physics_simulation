#ifndef CS207_CONSTRAINT_HPP
#define CS207_CONSTRAINT_HPP

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"
#include "Point.hpp"
#include "Simplex.hpp"
#include "Force.hpp"
#include "BoundingBox.hpp"

#include <fstream>
#include <cmath>
#include <assert.h>

/* Making a combined Constraint using C1, C2 constructors. Constraints are nested. */
template< typename C1, typename C2>
struct ConstraintPair{
	C1 constraint1;
	C2 constraint2;
	
	/* Making a combined Constraint using C1, C2 constructors */
	ConstraintPair(C1 c1, C2 c2) : constraint1(c1), constraint2(c2) {
	}
		
	template <typename G>
  	void operator()(G& g, double t) {
  		constraint1(g,t);
  		constraint2(g,t);
  	}	
};

/* Constrain a Node to a point */
struct NodeConstraint {
	typedef std::vector<std::pair<int, Point>> constr_node_type;
	constr_node_type constrained_nodes;

	NodeConstraint(constr_node_type node) : constrained_nodes(node) {
	}
	
	template <typename G>
	void operator()(G& g, double t) {
		(void) t;
		typename std::vector<std::pair<int, Point>>::iterator constr_it;
	
		for(auto nit = g.node_begin(); nit != g.node_end(); ++nit)
			for (constr_it = constrained_nodes.begin(); constr_it != constrained_nodes.end(); ++constr_it){
				auto node = (*nit);
				auto constr_node = (*constr_it);
	  		if (node.index() == constr_node.first) node.set_position( constr_node.second );
  		}	
  }
};


/** Plane Constraint: Global constraint on Graph */
struct Plane{

	template <typename G>
  	void operator()(G& g, double t) {
  		// Impose constraint and reset node position and velocity at t = n+1
  		for(auto nit = g.node_begin(); nit != g.node_end(); ++nit){
  			auto node = (*nit);
  		
  			if ( node.position().z < -0.75){
  		
  				// Set position to nearest point (same x, y) on plane z = -0.75
  				node.set_position( Point(node.position().x, node.position().y, -0.75) );
  				// Set z component of velocity to 0
  				node.value().velocity.z = 0;
  			}
  		}
  	}
};

/** Sphere Constraint: Local constraint on nodes */
struct Sphere{

	template <typename G>
  	void operator()(G& g, double t) {
  		(void) t;
  		
  		auto c = Point(0.5, 0.5, -0.5);
  		double r = 0.15; 
  		auto bb = BoundingBox(c, r); 	
  		
  		//std::cerr << "min, max " << bb.min() << " " << bb.max() << "\n";	
  		
  		CS207::Clock start;
  		start.start();
  		
  		int i = 0;
//  		for(auto nit = g.node_begin(); nit != g.node_end(); ++nit){
  		
  		for (auto nit = g.neighbor_begin(bb); nit != g.neighbor_end(bb); ++nit){
  			auto node = (*nit);  	
  			
  			++i;  
  			//std::cerr << "index " << node.index() << "\n";			  
  		
  			if ( (node.position()-c).mag() < r){
  				auto R = (node.position()-c) / (node.position()-c).mag();
  				auto v = node.value().velocity;
  		
  			// Set position to nearest point on sphere 
  				node.set_position(c + R * r);
  			
  			// Set the normal component of velocity to 0
  				node.value().velocity = v - R * R.dot(v);
  			
  			}  			
  		}
	 		std::cerr << "i " << i << "\n";
  	//	std::cerr << "computation time " << start.elapsed() << "\n";
		}
};

/** Neighbor Constraint 
	* @pre Two neighbors will not become indefinitely close
	*
	* If node j violates node i constraint, set node j position and value. There is no need 
	* to check if node i violates node j. 0 <= i < j < g.num_node()
 	*/
struct Neighbor{

	template <typename G>
  	void operator()(G& g, double t) {
  		(void) t;
  		for(auto nit = g.node_begin(); nit != g.node_end(); ++nit){
  			auto node_i = (*nit);
  			auto c = node_i.position();
  			
  			// Neighbor() constraint: initialize radius to a large value
  			node_i.value().radius = 10;
  	
  			// Calculate and update the next radius
  			for(auto eit = node_i.edge_begin(); eit != node_i.edge_end(); ++eit){
  				auto edge = (*eit);	
	  			if (node_i.value().radius > edge.length()) node_i.value().radius = edge.length();
  			}    
  			double r = node_i.value().radius;
  						
  			for(auto mit = nit; mit != g.node_end(); ++mit){
  				// do not want the same node
  				if (nit != mit){
	  				auto node_j = (*mit);
  				
	  				if ( (node_j.position() - c).mag() < r ) {
  						auto R = (node_j.position() - c) / (node_j.position() - c).mag();
  						auto v = node_j.value().velocity;
  					
  						// Set position to nearest point on sphere 
  						node_j.set_position(c + R * r );	

  						// Set the normal component of velocity to 0
	  					node_j.value().velocity = v - R * R.dot(v);
	  					node_j.value().velocity = Point(0,0,0);
  					}
  				}	
  			}
  		}
  }
};



  template< typename C1, typename C2>
  ConstraintPair<C1, C2> make_combined_constraint(C1 first, C2 second){ 
  	return ConstraintPair<C1, C2>(first, second);
  }

#endif
