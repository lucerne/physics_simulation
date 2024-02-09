 short pseudocode for the shallow water time inte- gration, Equation (8), and post-processing, Equation (9). 


/** Change a graph and a dual graph's nodes values Qbar_k, Q_n according to shallow 
 * water equation. 
 * @param[in,out] pg Graph
 * @param[in,out] dg Dual Graph 
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

double shallow_water_integration(double t, double dt) {
  // HW2 #1 YOUR CODE HERE
  
  // Update position with positions at t = n
  for(auto it = g.node_begin(); it != g.node_end(); ++it){
  	auto node = (*it);
  	auto pos = node.position();
  
  	// Skip two boundary points.
  	//if (node.position() != Point(0,0,0) && node.position() != Point(1,0,0)){
  	 
  		// Retrieve the velocity
  		auto vel = node.value().second;
  
	  	// Calculate and update the next position. 
	  	Point pos_next = pos + vel * dt;
  		node.set_position(pos_next);
  	//}
  }
  
  // Impose constraint and reset node position and velocity at t = n+1
  for(auto it = g.node_begin(); it != g.node_end(); ++it){
  	auto node = (*it);
  	constr(node, t);
  }
  	
  // Update velocity with positions at t = n+1
  for(auto it = g.node_begin(); it != g.node_end(); ++it){	
    auto node = (*it);
  	auto vel = node.value().second;
  	
  	// Skip two boundary points
  	//if (node.position() != Point(0,0,0) && node.position() != Point(1,0,0))
  	{  
  	  	// Calculate force
  		Point f_next = force(node,t);
  	 	
  	 	// Calculate and update the next velocity. 
  		Point v_next = vel + f_next * dt / m;
  		
  		// Not Sure:
  		// Shouldn't this be abstracted and accessed with a function rather than direct access?
  		node.value().second = v_next;
  	}
  }
  return t + dt;
}
