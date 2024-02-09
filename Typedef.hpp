#ifndef CS207_TYPEDEF_HPP
#define CS207_TYPEDEF_HPP

#include "Graph.hpp"
#include "Point.hpp"
#include "Simplex.hpp"


/* Gravity in meters/sec */
static constexpr double grav = 9.81;
/* mass */
double m;
/* viscosity */
double c;

/* Node attributes */
template <typename M, typename V>
	struct node_value {
		M mass;
		V velocity;
};

/* Edge attributes */
template <typename K, typename L>
	struct edge_value {
		K spring_constant;
		L spring_length;
};

/* Graph Types definitions */
typedef Graph<std::pair<double, Point>, edge_value<double, double>> GraphType;
typedef GraphType::Node Node;
typedef GraphType::Edge Edge;
typedef GraphType::node_value_type node_value_type;
typedef GraphType::edge_value_type edge_value_type;


#endif