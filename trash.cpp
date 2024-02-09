/* Constrain a Node to a point */
struct NodeConstraint {
	unsigned idx;
	Point p;
	NodeConstraint(unsigned _idx, const Point& _p)
		: idx(_idx), p(_p) {
	}
	template <typename NODE>
	void operator()(NODE n, double t) {
		if (n.index() == idx) {
			n.position() = p;
			n.value().second = Point(0,0,0);
		}
	}
};