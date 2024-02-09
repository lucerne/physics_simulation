/**
 * @file shallow_water.cpp
 * Implementation of a shallow water system using Mesh
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D point list (one per line) defined by three doubles
 * Second file: Triangles (one per line) defined by 3 indices into the point list
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

#include <cmath>
#include <fstream>

#include "Point.hpp"
#include "Simplex.hpp"
#include "Graph.hpp"
#include "Mesh.hpp"

// Standard gravity (average gravity at Earth's surface) in meters/sec^2
static constexpr double grav = 9.80665;

/** Water column characteristics */
struct QVar {
  double h;	  //< Height of fluid
  double hu;	//< Height times average x velocity of column
  double hv;	//< Height times average y velocity of column

  /** Default constructor.
   * A default water column is 1 unit high with no velocity. */
  QVar()
      : h(1), hu(0), hv(0) {
  }
  /** Construct the given water column. */
  QVar(double h_, double hu_, double hv_)
      : h(h_), hu(hu_), hv(hv_) {
  }
  
  /** Add @a QVar to this point. */
/*
  QVar& operator+(const QVar& q) {
    h += q.h;
    hu += q.hu;
    hv += q.hv;
    return *this;
  }  
  */
};

/*
inline QVar operator+(QVar a, const QVar& b) {
  return a += b;
}
*/

/** Function object for calculating shallow-water flux.
 *          |e
 *   T_k    |---> n = (nx,ny)   T_m
 *   qk     |                   qm
 *          |
 * @param[in] nx, ny Defines the 2D outward normal vector n = (@a nx, @a ny)
 *            from triangle T_k to triangle T_m. The length of n is equal to the
 *            the length of the edge, |n| = |e|.
 * @param[in] dt The time step taken by the simulation. Used to compute the
 *               Lax-Wendroff dissipation term.
 * @param[in] qk The values of the conserved variables on the left of the edge.
 * @param[in] qm The values of the conserved variables on the right of the edge.
 * @return The flux of the conserved values across the edge e
 */
struct EdgeFluxCalculator {
  QVar operator()(double nx, double ny, double dt,
                  const QVar& qk, const QVar& qm) {
    double e_length = std::sqrt(nx*nx + ny*ny);
    nx /= e_length;
    ny /= e_length;

    // The velocities normal to the edge
    double wm = (qm.hu * nx + qm.hv * ny) / qm.h;
    double wk = (qk.hu * nx + qk.hv * ny) / qk.h;

    // Lax-Wendroff local dissipation coefficient
    double vm = std::sqrt(grav*qm.h) + std::sqrt(qm.hu*qm.hu+qm.hv*qm.hv) / qm.h;
    double vk = std::sqrt(grav*qk.h) + std::sqrt(qk.hu*qk.hu+qk.hv*qk.hv) / qk.h;
    double a = dt * std::max(vm*vm, vk*vk);

    // Helper values
    double scale = 0.5 * e_length;
    double gh2 = 0.5 * grav * (qm.h*qm.h + qk.h*qk.h);

    // Simple flux with dissipation for stability
    return QVar(scale * (wm*qm.h  + wk*qk.h)           - a * (qm.h  - qk.h),
                scale * (wm*qm.hu + wk*qk.hu + gh2*nx) - a * (qm.hu - qk.hu),
                scale * (wm*qm.hv + wk*qk.hv + gh2*ny) - a * (qm.hv - qk.hv));
  }
};

/** Node position function object for use in the SDLViewer. */
struct NodePosition {
  template <typename NODE>
  Point operator()(const NODE& n) {
    // HW3B: You may change this to plot something other than the
    // positions of the nodes
    return n.position();
  }
};

/** Node position function object for use in the SDLViewer. */
struct NodeColor {
  template <typename NODE>
  CS207::Color operator()(const NODE& n) {
    // HW3B: You may change this to plot colors other than white
    (void) n;
    return CS207::Color(1,1,1);
  }
};


// HW3B: Placeholder for Mesh Type!
struct NodeData {
  QVar q;
};
  
struct EdgeData {
  QVar f_e; // flux magnitude, outward_normal
};

struct TriData  {
  QVar q_bar; 
  double area;	// area of a triangle
  std::vector<Point>	outward_normal;		// vector of outward normal
};

typedef Mesh<NodeData,EdgeData,TriData> MeshType;

/** Integrate a hyperbolic conservation law defined over the mesh m
 * with flux functor f by dt in time.
 */
template <typename FLUX>
double hyperbolic_step(MeshType& m, FLUX& f, double t, double dt) {
  // HW3B: YOUR CODE HERE
  // Step the finite volume model in time by dt.
  // Implement Equation 7 from your pseudocode here.
  
  // Calculate and update the stored flux
  for(auto it = m.triangle_begin(); it != m.triangle_end(); ++it){
  	auto tri1 = (*it);
  	auto qk = tri1.value().tri_data.q_bar;
	
  	// iterate 0-3 neighbors of triangle
  	// should iterate over 3 mesh nodes of triangle
  	auto node1 = tri1.value().three_mesh_node[0];
  	auto node2 = tri1.value().three_mesh_node[1];
  	auto node3 = tri1.value().three_mesh_node[2];
  	
  	// three mesh edges at the boundary, same order as normal
  	auto e1 = m.add_edge(node2, node3);
  	auto e2 = m.add_edge(node3, node1);
  	auto e3 = m.add_edge(node1, node2);
  	
  	// q value for triangle, used to make ghost triangles
  	auto hI = tri1.value().tri_data.q_bar.h;
  	// there aught to be a better way to do this! 
  	// q values for three edges of triangle
  	auto q1 = QVar();  	
  	// checking boundary 
  	if (e1.value().two_tri_node.size() == 1){
  		q1 = QVar(hI,0,0);
  	} else{
	  	auto tri2 = e1.value().two_tri_node[0];
  		if ( tri1 == tri2) tri2 = e1.value().two_tri_node[1];
  		q1 = tri2.value().tri_data.q_bar;
  	}

  	auto q2 = QVar();
  	if (e2.value().two_tri_node.size() == 1){
  		q2 = QVar(hI,0,0);
  	} else{
  		auto tri2 = e2.value().two_tri_node[0];
  		if ( tri1 == tri2) tri2 = e2.value().two_tri_node[1];
  		q2 = tri2.value().tri_data.q_bar;
  	}  	

  	auto q3 = QVar();
  	if (e3.value().two_tri_node.size() == 1){
  		q3 = QVar(hI,0,0);
  	} else{
			auto tri2 = e3.value().two_tri_node[0];  
  		if ( tri1 == tri2) tri2 = e3.value().two_tri_node[1];
  		q3 = tri2.value().tri_data.q_bar;
  	}  	
 	
  	// normals at the boundary, same order as mesh edges 
  	auto n1 = tri1.value().tri_data.outward_normal[0];
  	auto n2 = tri1.value().tri_data.outward_normal[1];
  	auto n3 = tri1.value().tri_data.outward_normal[2];

		// Update edge data			
		e1.value().edge_data.f_e = f(n1.x, n1.y, dt, qk, q1);
		e2.value().edge_data.f_e = f(n2.x, n2.y, dt, qk, q2);
		e3.value().edge_data.f_e = f(n3.x, n3.y, dt, qk, q3);
		
		
	/*		
		std::cerr << "Triangle " << tri1.index() << "@ " << t << "\n";		
		std::cerr << "Area " << tri1.value().tri_data.area << "\n";
		std::cerr << "Node positions " << " (" << node1.position() << ") (" << node2.position() << ") (" << node3.position() << ") \n";
		std::cerr << "Triangle QVar [Q_bar, water column characteristics] " << "h " << tri1.value().tri_data.q_bar.h << " hu " << tri1.value().tri_data.q_bar.hu << " hv " << tri1.value().tri_data.q_bar.hv << "\n";			

		
		// Edge 0
  		std::cerr << "Edge 0 " << " (" << e1.node1().position() << ") (" << e1.node2().position() << ") \n";
  		auto nor = tri1.value().tri_data.outward_normal[0];
  		std::cerr << "  Normal " << " (" << nor << ") \n";
  		
  		if (e1.value().two_tri_node.size() == 2){
  		std::cerr << "	Adj Triangles " << e1.value().two_tri_node[0].index() << " " << e1.value().two_tri_node[1].index() << "\n";
  		
  			if (tri1 == e1.value().two_tri_node[0]) {
  				std::cerr << "	Opposite Tri QVar " << "h= " << e1.value().two_tri_node[1].value().tri_data.q_bar.h << " hu= " << e1.value().two_tri_node[1].value().tri_data.q_bar.hu  << " hv= " << e1.value().two_tri_node[1].value().tri_data.q_bar.hv << "\n"; 
  				
  				auto f1 = f(nor.x, nor.y, dt, tri1.value().tri_data.q_bar, e1.value().two_tri_node[1].value().tri_data.q_bar);
  				std::cerr << "	Boundary flux " << "h= " << f1.h << " hu= " << f1.hu << " hv= " << f1.hv << "\n";
  			} else {
  				std::cerr << "	Opposite Tri QVar " << "h= " << e1.value().two_tri_node[0].value().tri_data.q_bar.h << " hu= " << e1.value().two_tri_node[0].value().tri_data.q_bar.hu  << " hv= " << e1.value().two_tri_node[0].value().tri_data.q_bar.hv << "\n"; 
  				
  				auto f1 = f(nor.x, nor.y, dt, tri1.value().tri_data.q_bar, e1.value().two_tri_node[0].value().tri_data.q_bar);
  				std::cerr << "	Boundary flux " << "h= " << f1.h << " hu= " << f1.hu << " hv= " << f1.hv << "\n";
  			}
  		} else {
  			std::cerr << "	Adj Triangles " << e1.value().two_tri_node[0].index() << "\n";
  			QVar qm = QVar(e1.value().two_tri_node[0].value().tri_data.q_bar.h,0,0);
  			std::cerr << "	Opposite Tri QVar " << "h= " << qm.h << " hu= " << 0 << " hv= " << 0 << "\n"; 
  			auto f1 = f(nor.x, nor.y, dt, tri1.value().tri_data.q_bar, qm);
  			std::cerr << "	Boundary flux" << "h= " << f1.h << " hu= " << f1.hu << " hv= " << f1.hv << "\n";
  		}
  		
  		// Edge 1
  		std::cerr << "Edge 1 " << " (" << e2.node1().position() << ") (" << e2.node2().position() << ") \n";
  		nor = tri1.value().tri_data.outward_normal[1];
  		std::cerr << "  Normal " << " (" << nor << ") \n";
  		if (e2.value().two_tri_node.size() == 2){
  		std::cerr << "	Adj Triangles " << e2.value().two_tri_node[0].index() << " " << e2.value().two_tri_node[1].index() << "\n";
  			if (tri1 == e2.value().two_tri_node[0]) {
  				std::cerr << "	Opposite Tri QVar " << "h= " << e2.value().two_tri_node[1].value().tri_data.q_bar.h << " hu= " << e2.value().two_tri_node[1].value().tri_data.q_bar.hu  << " hv= " << e2.value().two_tri_node[1].value().tri_data.q_bar.hv << "\n"; 
  				
  				auto f2 = f(nor.x, nor.y, dt, tri1.value().tri_data.q_bar, e2.value().two_tri_node[1].value().tri_data.q_bar);
  				std::cerr << "	Boundary flux " << "h= " << f2.h << " hu= " << f2.hu << " hv= " << f2.hv << "\n";
  			} else {
  				std::cerr << "	Opposite Tri QVar " << "h= " << e2.value().two_tri_node[0].value().tri_data.q_bar.h << " hu= " << e2.value().two_tri_node[0].value().tri_data.q_bar.hu  << " hv= " << e2.value().two_tri_node[0].value().tri_data.q_bar.hv << "\n"; 
  				
  				auto f2 = f(nor.x, nor.y, dt, tri1.value().tri_data.q_bar, e2.value().two_tri_node[0].value().tri_data.q_bar);
  				std::cerr << "	Boundary flux " << "h= " << f2.h << " hu= " << f2.hu << " hv= " << f2.hv << "\n";
  			}
  		} else {
  			std::cerr << "	Adj Triangles " << e2.value().two_tri_node[0].index() << "\n";
  			QVar qm = QVar(e2.value().two_tri_node[0].value().tri_data.q_bar.h,0,0);
  			std::cerr << "	Opposite Tri QVar " << "h= " << qm.h << " hu= " << 0 << " hv= " << 0 << "\n"; 
  			auto f2 = f(nor.x, nor.y, dt, tri1.value().tri_data.q_bar, qm);
  			std::cerr << "	Boundary flux" << "h= " << f2.h << " hu= " << f2.hu << " hv= " << f2.hv << "\n";
  		}
  	
  		// Edge 2
  		std::cerr << "Edge 2 " << " (" << e3.node1().position() << ") (" << e3.node2().position() << ") \n";
  		nor = tri1.value().tri_data.outward_normal[2];
  		std::cerr << "  Normal " << " (" << nor << ") \n";
  		
  		if (e3.value().two_tri_node.size() == 2){
  		std::cerr << "	Adj Triangles " << e3.value().two_tri_node[0].index() << " " << e3.value().two_tri_node[1].index() << "\n";
  			if (tri1 == e3.value().two_tri_node[0]) {
  				std::cerr << "	Opposite Tri QVar " << "h= " << e3.value().two_tri_node[1].value().tri_data.q_bar.h << " hu= " << e3.value().two_tri_node[1].value().tri_data.q_bar.hu  << " hv= " << e3.value().two_tri_node[1].value().tri_data.q_bar.hv << "\n"; 
  				
  				auto f3 = f(nor.x, nor.y, dt, tri1.value().tri_data.q_bar, e3.value().two_tri_node[1].value().tri_data.q_bar);
  				std::cerr << "	Boundary flux " << "h= " << f3.h << " hu= " << f3.hu << " hv= " << f3.hv << "\n";
  			} else {
  				std::cerr << "	Opposite Tri QVar " << "h= " << e3.value().two_tri_node[0].value().tri_data.q_bar.h << " hu= " << e3.value().two_tri_node[0].value().tri_data.q_bar.hu  << " hv= " << e3.value().two_tri_node[0].value().tri_data.q_bar.hv << "\n"; 
  				
  				auto f3 = f(nor.x, nor.y, dt, tri1.value().tri_data.q_bar, e3.value().two_tri_node[0].value().tri_data.q_bar);
  				std::cerr << "	Boundary flux " << "h= " << f3.h << " hu= " << f3.hu << " hv= " << f3.hv << "\n";
  			}
  		} else {
  			std::cerr << "	Adj Triangles " << e3.value().two_tri_node[0].index() << "\n";
  			QVar qm = QVar(e3.value().two_tri_node[0].value().tri_data.q_bar.h,0,0);
  			std::cerr << "	Opposite Tri QVar " << "h= " << qm.h << " hu= " << 0 << " hv= " << 0 << "\n"; 
  			auto f3 = f(nor.x, nor.y, dt, tri1.value().tri_data.q_bar, qm);
  			std::cerr << "	Boundary flux" << "h= " << f3.h << " hu= " << f3.hu << " hv= " << f3.hv << "\n";
  		}
		*/
  }
  
  
  
  
  // Update Qk
  for(auto it = m.triangle_begin(); it != m.triangle_end(); ++it){
  	auto tri1 = (*it);
  	
  	// should iterate over 3 mesh nodes of triangle
  	auto node1 = tri1.value().three_mesh_node[0];
  	auto node2 = tri1.value().three_mesh_node[1];
  	auto node3 = tri1.value().three_mesh_node[2];
  	
  	// three mesh edges at the boundary, same order as normal
  	auto e1 = m.add_edge(node2, node3);
  	auto e2 = m.add_edge(node3, node1);
  	auto e3 = m.add_edge(node1, node2);
  	
  	// flux
		auto f1 = e1.value().edge_data.f_e;
		auto f2 = e2.value().edge_data.f_e;
		auto f3 = e3.value().edge_data.f_e;

  	// area
  	auto area =  tri1.value().tri_data.area;	  	

  	// Summing all the q	
		QVar dq;	
			
		dq.h = f1.h + f2.h + f3.h;
		dq.hu = f1.hu + f2.hu + f3.hu;
		dq.hv = f1.hv + f2.hv + f3.hv;
		
  	// update Qk
  	tri1.value().tri_data.q_bar.h = tri1.value().tri_data.q_bar.h - dq.h * dt / area; 
  	tri1.value().tri_data.q_bar.hu = tri1.value().tri_data.q_bar.hu - dq.hu * dt / area;   
  	tri1.value().tri_data.q_bar.hv = tri1.value().tri_data.q_bar.hv - dq.hv * dt / area; 
	}
  
  return t + dt;
}

/** Convert the triangle-averaged values to node-averaged values for viewing. */
template <typename FLUX>
void post_process(MeshType& m, FLUX& f, double t, double dt){
  // HW3B: Post-processing step
  // Translate the triangle-averaged values to node-averaged values
  // Implement Equation 8 from your pseudocode here
  for(auto it = m.node_begin(); it != m.node_end(); ++it){
  	auto node = (*it);
  	unsigned k = 0;
		double qh_sum = 0;
  	double qhu_sum = 0;
  	double qhv_sum = 0;
  	double area_sum = 0;
  	
  	while(k != node.value().many_tri_node.size()){
  		auto tri = node.value().many_tri_node[k];
  		auto qk = tri.value().tri_data.q_bar;
  		auto tk = tri.value().tri_data.area;
  		qh_sum = qh_sum + qk.h * tk;
  		qhu_sum = qhu_sum + qk.hu * tk;
  		qhv_sum = qhv_sum + qk.hv * tk;
  		area_sum = area_sum + tk;
  		++k;
  	}  	
  	node.value().node_data.q.h = qh_sum/area_sum;
  	node.value().node_data.q.hu = qhu_sum/area_sum;
		node.value().node_data.q.hv = qhv_sum/area_sum;
		
		/*
	for(auto it = m.triangle_begin(); it != m.triangle_end(); ++it){
  		auto tri = (*it);
  		
   		auto n1 = tri.value().three_mesh_node[0];
  		auto n2 = tri.value().three_mesh_node[1];
  		auto n3 = tri.value().three_mesh_node[2]; 

  		auto e1 = m.add_edge(n2, n3);
  		auto e2 = m.add_edge(n3, n1);
  		auto e3 = m.add_edge(n1, n2);
		
		std::cerr << "Triangle " << tri.index() << "@ " << t << "\n";		
			std::cerr << "Area " << tri.value().tri_data.area << "\n";
			std::cerr << "Node positions " << " (" << n1.position() << ") (" << n2.position() << ") (" << n3.position() << ") \n";
			std::cerr << "Triangle QVar [Q_bar, water column characteristics] " << "h " << tri.value().tri_data.q_bar.h << " hu " << tri.value().tri_data.q_bar.hu << " hv " << tri.value().tri_data.q_bar.hv << "\n";			
			
			// Edge 0
  		std::cerr << "Edge 0 " << " (" << e1.node1().position() << ") (" << e1.node2().position() << ") \n";
  		auto nor = tri.value().tri_data.outward_normal[0];
  		std::cerr << "  Normal " << " (" << nor << ") \n";
  		
  		if (e1.value().two_tri_node.size() == 2){
  		std::cerr << "	Adj Triangles " << e1.value().two_tri_node[0].index() << " " << e1.value().two_tri_node[1].index() << "\n";
  		
  			if (tri == e1.value().two_tri_node[0]) {
  				std::cerr << "	Opposite Tri QVar " << "h= " << e1.value().two_tri_node[1].value().tri_data.q_bar.h << " hu= " << e1.value().two_tri_node[1].value().tri_data.q_bar.hu  << " hv= " << e1.value().two_tri_node[1].value().tri_data.q_bar.hv << "\n"; 
  				
  				auto f1 = f(nor.x, nor.y, dt, tri.value().tri_data.q_bar, e1.value().two_tri_node[1].value().tri_data.q_bar);
  				std::cerr << "	Boundary flux " << "h= " << f1.h << " hu= " << f1.hu << " hv= " << f1.hv << "\n";
  			} else {
  				std::cerr << "	Opposite Tri QVar " << "h= " << e1.value().two_tri_node[0].value().tri_data.q_bar.h << " hu= " << e1.value().two_tri_node[0].value().tri_data.q_bar.hu  << " hv= " << e1.value().two_tri_node[0].value().tri_data.q_bar.hv << "\n"; 
  				
  				auto f1 = f(nor.x, nor.y, dt, tri.value().tri_data.q_bar, e1.value().two_tri_node[0].value().tri_data.q_bar);
  				std::cerr << "	Boundary flux " << "h= " << f1.h << " hu= " << f1.hu << " hv= " << f1.hv << "\n";
  			}
  		} else {
  			std::cerr << "	Adj Triangles " << e1.value().two_tri_node[0].index() << "\n";
  			QVar qm = QVar(e1.value().two_tri_node[0].value().tri_data.q_bar.h,0,0);
  			std::cerr << "	Opposite Tri QVar " << "h= " << qm.h << " hu= " << 0 << " hv= " << 0 << "\n"; 
  			auto f1 = f(nor.x, nor.y, dt, tri.value().tri_data.q_bar, qm);
  			std::cerr << "	Boundary flux" << "h= " << f1.h << " hu= " << f1.hu << " hv= " << f1.hv << "\n";
  		}
  		
  		// Edge 1
  		std::cerr << "Edge 1 " << " (" << e2.node1().position() << ") (" << e2.node2().position() << ") \n";
  		nor = tri.value().tri_data.outward_normal[1];
  		std::cerr << "  Normal " << " (" << nor << ") \n";
  		if (e2.value().two_tri_node.size() == 2){
  		std::cerr << "	Adj Triangles " << e2.value().two_tri_node[0].index() << " " << e2.value().two_tri_node[1].index() << "\n";
  			if (tri == e2.value().two_tri_node[0]) {
  				std::cerr << "	Opposite Tri QVar " << "h= " << e2.value().two_tri_node[1].value().tri_data.q_bar.h << " hu= " << e2.value().two_tri_node[1].value().tri_data.q_bar.hu  << " hv= " << e2.value().two_tri_node[1].value().tri_data.q_bar.hv << "\n"; 
  				
  				auto f2 = f(nor.x, nor.y, dt, tri.value().tri_data.q_bar, e2.value().two_tri_node[1].value().tri_data.q_bar);
  				std::cerr << "	Boundary flux " << "h= " << f2.h << " hu= " << f2.hu << " hv= " << f2.hv << "\n";
  			} else {
  				std::cerr << "	Opposite Tri QVar " << "h= " << e2.value().two_tri_node[0].value().tri_data.q_bar.h << " hu= " << e2.value().two_tri_node[0].value().tri_data.q_bar.hu  << " hv= " << e2.value().two_tri_node[0].value().tri_data.q_bar.hv << "\n"; 
  				
  				auto f2 = f(nor.x, nor.y, dt, tri.value().tri_data.q_bar, e2.value().two_tri_node[0].value().tri_data.q_bar);
  				std::cerr << "	Boundary flux " << "h= " << f2.h << " hu= " << f2.hu << " hv= " << f2.hv << "\n";
  			}
  		} else {
  			std::cerr << "	Adj Triangles " << e2.value().two_tri_node[0].index() << "\n";
  			QVar qm = QVar(e2.value().two_tri_node[0].value().tri_data.q_bar.h,0,0);
  			std::cerr << "	Opposite Tri QVar " << "h= " << qm.h << " hu= " << 0 << " hv= " << 0 << "\n"; 
  			auto f2 = f(nor.x, nor.y, dt, tri.value().tri_data.q_bar, qm);
  			std::cerr << "	Boundary flux" << "h= " << f2.h << " hu= " << f2.hu << " hv= " << f2.hv << "\n";
  		}
  		
  		// Edge 2
  		std::cerr << "Edge 2 " << " (" << e3.node1().position() << ") (" << e3.node2().position() << ") \n";
  		nor = tri.value().tri_data.outward_normal[2];
  		std::cerr << "  Normal " << " (" << nor << ") \n";
  		if (e3.value().two_tri_node.size() == 2){
  		std::cerr << "	Adj Triangles " << e3.value().two_tri_node[0].index() << " " << e3.value().two_tri_node[1].index() << "\n";
  			if (tri == e3.value().two_tri_node[0]) {
  				std::cerr << "	Opposite Tri QVar " << "h= " << e3.value().two_tri_node[1].value().tri_data.q_bar.h << " hu= " << e3.value().two_tri_node[1].value().tri_data.q_bar.hu  << " hv= " << e3.value().two_tri_node[1].value().tri_data.q_bar.hv << "\n"; 
  				
  				auto f3 = f(nor.x, nor.y, dt, tri.value().tri_data.q_bar, e2.value().two_tri_node[1].value().tri_data.q_bar);
  				std::cerr << "	Boundary flux " << "h= " << f3.h << " hu= " << f3.hu << " hv= " << f3.hv << "\n";
  			} else {
  				std::cerr << "	Opposite Tri QVar " << "h= " << e3.value().two_tri_node[0].value().tri_data.q_bar.h << " hu= " << e3.value().two_tri_node[0].value().tri_data.q_bar.hu  << " hv= " << e3.value().two_tri_node[0].value().tri_data.q_bar.hv << "\n"; 
  				
  				auto f3 = f(nor.x, nor.y, dt, tri.value().tri_data.q_bar, e3.value().two_tri_node[0].value().tri_data.q_bar);
  				std::cerr << "	Boundary flux " << "h= " << f3.h << " hu= " << f3.hu << " hv= " << f3.hv << "\n";
  			}
  		} else {
  			std::cerr << "	Adj Triangles " << e3.value().two_tri_node[0].index() << "\n";
  			QVar qm = QVar(e3.value().two_tri_node[0].value().tri_data.q_bar.h,0,0);
  			std::cerr << "	Opposite Tri QVar " << "h= " << qm.h << " hu= " << 0 << " hv= " << 0 << "\n"; 
  			auto f3 = f(nor.x, nor.y, dt, tri.value().tri_data.q_bar, qm);
  			std::cerr << "	Boundary flux" << "h= " << f3.h << " hu= " << f3.hu << " hv= " << f3.hv << "\n";
  		}
		
		*/
		node.set_position(Point(node.position().x, node.position().y, node.value().node_data.q.h));
	}
}


int main(int argc, char* argv[]) {
  // Check arguments  
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TRIS_FILE\n";
    exit(1);
  }

	// Mesh class
  MeshType mesh;

  // Read all Points and add them to the Mesh
  std::ifstream nodes_file(argv[1]);
  Point p;
  while (CS207::getline_parsed(nodes_file, p)) {
    // HW3B: Need to implement add_node before this can be used!
    mesh.add_node(p);
  }

  // Read all mesh triangles and add their edges to the Mesh
  std::ifstream tris_file(argv[2]);
  Triangle t;
  while (CS207::getline_parsed(tris_file, t)) {
    // HW3B: Need to implement add_triangle before this can be used!
    mesh.add_triangle(mesh.node(t.n[0]), mesh.node(t.n[1]), mesh.node(t.n[2]));
  }

  std::cerr << "Happy \n"; 
  // Print out the stats
  std::cout << mesh.num_nodes() << " "
            << mesh.num_edges() << " "
            << mesh.num_triangles() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  // HW3B: Need to define Mesh::node_type before these can be used!
  auto node_map = viewer.empty_node_map(mesh.mesh_graph);
  viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                   NodeColor(), NodePosition(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
  viewer.center_view();


  // HW3B Initialization Functors

  /* Compute Gaussian depression initial conditions from Point */
  struct GaussianDepression {
    QVar operator()(const Point& p) {
      double x = p.x - 0.75;
      double y = p.y - 0.0;
      return QVar(1 - 0.75*std::exp(-80*(x*x+y*y)), 0, 0);
    }
  };

  /* Compute the column of water initial conditions from Point */
  struct WaterColumn {
    QVar operator()(const Point& p) {
      double x = p.x - 0.75;
      double y = p.y - 0.0;
      if (x*x + y*y < 0.15*0.15)
        return QVar(1.75, 0, 0);
      else
        return QVar(1, 0, 0);
    }
  };

  /* Compute the dam break initial conditions from Point */
  struct DamBreak {
    QVar operator()(const Point& p) {
      if (p.x < 0)
        return QVar(1.75, 0, 0);
      else
        return QVar(1, 0, 0);
    }
  };

  // HW3B: Initialization
  // Set the initial conditions
  // Perform any needed precomputation
  // Compute the minimum edge length and maximum water height for computing dt
  double max_h = 0;
  auto init_cond = GaussianDepression();
  for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it) {
    auto n = *it;
    n.value().node_data.q = init_cond(n.position());
    max_h = std::max(max_h, n.value().node_data.q.h);
  }

	double min_edge_length = 10;
  for (auto it = mesh.edge_begin(); it != mesh.edge_end(); ++it) {
    auto edge = *it;
    min_edge_length = std::min(min_edge_length, edge.length());
  }
  
	// Set the initial values of the triangles to the average of their nodes
  for (auto it = mesh.triangle_begin(); it != mesh.triangle_end(); ++it) {
    auto t = *it;  
    
    auto n1 = t.value().three_mesh_node[0];
    auto n2 = t.value().three_mesh_node[1];
    auto n3 = t.value().three_mesh_node[2];

    t.value().tri_data.area = mesh.area(n1,n2,n3);  
    
    t.value().tri_data.outward_normal.push_back(mesh.normal(n1,n2,n3));
    t.value().tri_data.outward_normal.push_back(mesh.normal(n2,n3,n1));
    t.value().tri_data.outward_normal.push_back(mesh.normal(n3,n1,n2));   
            
    t.value().tri_data.q_bar.h = (n1.value().node_data.q.h +
                       n2.value().node_data.q.h +
                       n3.value().node_data.q.h) / 3.0;
    t.value().tri_data.q_bar.hu = (n1.value().node_data.q.hu +
                       n2.value().node_data.q.hu +
                       n3.value().node_data.q.hu) / 3.0;
    t.value().tri_data.q_bar.hv = (n1.value().node_data.q.hv +
                       n2.value().node_data.q.hv +
                       n3.value().node_data.q.hv) / 3.0;
  }

  // HW3B: Timestep
  // CFL stability condition requires dt <= dx / max|velocity|
  // For the shallow water equations with u = v = 0 initial conditions
  //   we can compute the minimum edge length and maximum original water height
  //   to set the time-step

  // dambreak() : dt = 0.03
  double dt = 0.25 * min_edge_length / (std::sqrt(grav*max_h));
//  std::cerr << "min_edge_length " << min_edge_length << " max_h " << max_h << " dt " << dt << "\n";
  double t_start = 0.0;
  double t_end = dt*1600;
  
  // Begin the simulation
  EdgeFluxCalculator flux;
  for (double t = t_start; t < t_end; t += dt) {  

  	for(auto it = mesh.triangle_begin(); it != mesh.triangle_end(); ++it){
  		auto tri = (*it);
  		
   		auto n1 = tri.value().three_mesh_node[0];
  		auto n2 = tri.value().three_mesh_node[1];
  		auto n3 = tri.value().three_mesh_node[2]; 

  		auto e1 = mesh.add_edge(n2, n3);
  		auto e2 = mesh.add_edge(n3, n1);
  		auto e3 = mesh.add_edge(n1, n2);

/*  		
  		std::cerr << "Triangle " << tri.index() << "@ " << t << "\n";		
			std::cerr << "Area " << tri.value().tri_data.area << "\n";
			std::cerr << "Node positions " << " (" << n1.position() << ") (" << n2.position() << ") (" << n3.position() << ") \n";
			std::cerr << "Triangle QVar [Q_bar, water column characteristics] " << "h " << tri.value().tri_data.q_bar.h << " hu " << tri.value().tri_data.q_bar.hu << " hv " << tri.value().tri_data.q_bar.hv << "\n";			
			
			// Edge 0
  		std::cerr << "Edge 0 " << " (" << e1.node1().position() << ") (" << e1.node2().position() << ") \n";
  		auto nor = tri.value().tri_data.outward_normal[0];
  		std::cerr << "  Normal " << " (" << nor << ") \n";
  		
  		if (e1.value().two_tri_node.size() == 2){
  		std::cerr << "	Adj Triangles " << e1.value().two_tri_node[0].index() << " " << e1.value().two_tri_node[1].index() << "\n";
  		
  			if (tri == e1.value().two_tri_node[0]) {
  				std::cerr << "	Opposite Tri QVar " << "h= " << e1.value().two_tri_node[1].value().tri_data.q_bar.h << " hu= " << e1.value().two_tri_node[1].value().tri_data.q_bar.hu  << " hv= " << e1.value().two_tri_node[1].value().tri_data.q_bar.hv << "\n"; 
  				
  				auto f = flux(nor.x, nor.y, dt, tri.value().tri_data.q_bar, e1.value().two_tri_node[1].value().tri_data.q_bar);
  				std::cerr << "	Boundary flux " << "h= " << f.h << " hu= " << f.hu << " hv= " << f.hv << "\n";
  			} else {
  				std::cerr << "	Opposite Tri QVar " << "h= " << e1.value().two_tri_node[0].value().tri_data.q_bar.h << " hu= " << e1.value().two_tri_node[0].value().tri_data.q_bar.hu  << " hv= " << e1.value().two_tri_node[0].value().tri_data.q_bar.hv << "\n"; 
  				
  				auto f = flux(nor.x, nor.y, dt, tri.value().tri_data.q_bar, e1.value().two_tri_node[0].value().tri_data.q_bar);
  				std::cerr << "	Boundary flux " << "h= " << f.h << " hu= " << f.hu << " hv= " << f.hv << "\n";
  			}
  		} else {
  			std::cerr << "	Adj Triangles " << e1.value().two_tri_node[0].index() << "\n";
  			QVar qm = QVar(e1.value().two_tri_node[0].value().tri_data.q_bar.h,0,0);
  			std::cerr << "	Opposite Tri QVar " << "h= " << qm.h << " hu= " << 0 << " hv= " << 0 << "\n"; 
  			auto f = flux(nor.x, nor.y, dt, tri.value().tri_data.q_bar, qm);
  			std::cerr << "	Boundary flux" << "h= " << f.h << " hu= " << f.hu << " hv= " << f.hv << "\n";
  		}
  		
  		// Edge 1
  		std::cerr << "Edge 1 " << " (" << e2.node1().position() << ") (" << e2.node2().position() << ") \n";
  		nor = tri.value().tri_data.outward_normal[1];
  		std::cerr << "  Normal " << " (" << nor << ") \n";
  		if (e2.value().two_tri_node.size() == 2){
  		std::cerr << "	Adj Triangles " << e2.value().two_tri_node[0].index() << " " << e2.value().two_tri_node[1].index() << "\n";
  			if (tri == e2.value().two_tri_node[0]) {
  				std::cerr << "	Opposite Tri QVar " << "h= " << e2.value().two_tri_node[1].value().tri_data.q_bar.h << " hu= " << e2.value().two_tri_node[1].value().tri_data.q_bar.hu  << " hv= " << e2.value().two_tri_node[1].value().tri_data.q_bar.hv << "\n"; 
  				
  				auto f = flux(nor.x, nor.y, dt, tri.value().tri_data.q_bar, e2.value().two_tri_node[1].value().tri_data.q_bar);
  				std::cerr << "	Boundary flux " << "h= " << f.h << " hu= " << f.hu << " hv= " << f.hv << "\n";
  			} else {
  				std::cerr << "	Opposite Tri QVar " << "h= " << e2.value().two_tri_node[0].value().tri_data.q_bar.h << " hu= " << e2.value().two_tri_node[0].value().tri_data.q_bar.hu  << " hv= " << e2.value().two_tri_node[0].value().tri_data.q_bar.hv << "\n"; 
  				
  				auto f = flux(nor.x, nor.y, dt, tri.value().tri_data.q_bar, e2.value().two_tri_node[0].value().tri_data.q_bar);
  				std::cerr << "	Boundary flux " << "h= " << f.h << " hu= " << f.hu << " hv= " << f.hv << "\n";
  			}
  		} else {
  			std::cerr << "	Adj Triangles " << e2.value().two_tri_node[0].index() << "\n";
  			QVar qm = QVar(e2.value().two_tri_node[0].value().tri_data.q_bar.h,0,0);
  			std::cerr << "	Opposite Tri QVar " << "h= " << qm.h << " hu= " << 0 << " hv= " << 0 << "\n"; 
  			auto f = flux(nor.x, nor.y, dt, tri.value().tri_data.q_bar, qm);
  			std::cerr << "	Boundary flux" << "h= " << f.h << " hu= " << f.hu << " hv= " << f.hv << "\n";
  		}
  		
  		// Edge 2
  		std::cerr << "Edge 2 " << " (" << e3.node1().position() << ") (" << e3.node2().position() << ") \n";
  		nor = tri.value().tri_data.outward_normal[2];
  		std::cerr << "  Normal " << " (" << nor << ") \n";
  		if (e3.value().two_tri_node.size() == 2){
  		std::cerr << "	Adj Triangles " << e3.value().two_tri_node[0].index() << " " << e3.value().two_tri_node[1].index() << "\n";
  			if (tri == e3.value().two_tri_node[0]) {
  				std::cerr << "	Opposite Tri QVar " << "h= " << e3.value().two_tri_node[1].value().tri_data.q_bar.h << " hu= " << e3.value().two_tri_node[1].value().tri_data.q_bar.hu  << " hv= " << e3.value().two_tri_node[1].value().tri_data.q_bar.hv << "\n"; 
  				
  				auto f = flux(nor.x, nor.y, dt, tri.value().tri_data.q_bar, e2.value().two_tri_node[1].value().tri_data.q_bar);
  				std::cerr << "	Boundary flux " << "h= " << f.h << " hu= " << f.hu << " hv= " << f.hv << "\n";
  			} else {
  				std::cerr << "	Opposite Tri QVar " << "h= " << e3.value().two_tri_node[0].value().tri_data.q_bar.h << " hu= " << e3.value().two_tri_node[0].value().tri_data.q_bar.hu  << " hv= " << e3.value().two_tri_node[0].value().tri_data.q_bar.hv << "\n"; 
  				
  				auto f = flux(nor.x, nor.y, dt, tri.value().tri_data.q_bar, e3.value().two_tri_node[0].value().tri_data.q_bar);
  				std::cerr << "	Boundary flux " << "h= " << f.h << " hu= " << f.hu << " hv= " << f.hv << "\n";
  			}
  		} else {
  			std::cerr << "	Adj Triangles " << e3.value().two_tri_node[0].index() << "\n";
  			QVar qm = QVar(e3.value().two_tri_node[0].value().tri_data.q_bar.h,0,0);
  			std::cerr << "	Opposite Tri QVar " << "h= " << qm.h << " hu= " << 0 << " hv= " << 0 << "\n"; 
  			auto f = flux(nor.x, nor.y, dt, tri.value().tri_data.q_bar, qm);
  			std::cerr << "	Boundary flux" << "h= " << f.h << " hu= " << f.hu << " hv= " << f.hv << "\n";
  		}

*/			
  	}
  	
    // Increment the simulation from t to t + dt
    hyperbolic_step(mesh, flux, t, dt);
	  
    // Update viewer with nodes' new positions
    post_process(mesh, flux, t, dt);  
		
    // HW3B: Need to define Mesh::node_type before these can be used!   
    viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                     CS207::DefaultColor(), NodePosition(), node_map);
    viewer.set_label(t);

    // These lines slow down the animation for small meshes.
    // Feel free to remove them or tweak the constants.
    if (mesh.num_nodes() < 100)
      CS207::sleep(0.05);
  }

  return 0;
}