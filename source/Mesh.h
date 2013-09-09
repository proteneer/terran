#ifndef MESH_H_
#define MESH_H_

#include <vector>

using namespace std;

struct double2 {
    double x;
    double y;
};
 
struct double3 {
    double x;
    double y;
    double z;
};
 
struct int2 {
    int x;
    int y;
};

struct bool2 {
	bool x;
	bool y;
};

namespace Terran {

class Mesh {

public:
	Mesh(std::vector<double> X, std::vector<double> Y, std::vector<std::vector<double> > Z, double xp, double yp);

	// Used for defining edge locations (clockwise)
	enum Direction {UP, RIGHT, DOWN, LEFT};
		
	struct Edge {
		bool operator== (const Edge &e1) const;
		bool operator!= (const Edge &e1) const;
		Direction orientation() const;
		int2 coords;
		Direction extend;
	};

	// Edge orientation is as follows:
	// if e extends UP, then the edge faces the LEFT direction
	// if e extends RIGHT, then the edge faces in the UP direction
	// if e extends DOWN, then the edge faces in the RIGHT direction 
	// if e extends LEFT, then the edge faces in the DOWN direction
	// Flags:
	// reset the visits array
	std::vector<std::vector<double2> > contour(double height);

	unsigned int xSize() const;

	unsigned int ySize() const;

	bool xPeriodic() const;

	bool yPeriodic() const;

	void test();

private:

	std::vector<double2> path(Edge startEdge, double height);

	// find the next logical edge that crosses height h
	// stores the cut location in cut
	Edge findNextEdge(Edge start, double height, double2 &cut) const;

	// move a point in a given direction, taking into account of periodic boundary conditions
	int2 move(const int2 &p, Direction d) const;

	// mark the edge as having been visited in the visits_ array
	void markEdge(const Edge &u);

	Edge leftEdge(const Edge &u) const;

	Edge forwardEdge(const Edge &u) const;

	Edge rightEdge(const Edge &u) const;

	// returns true if the edge crosses height h
	bool edgeCross(const Edge &e, double h) const;
	// returns true if the edge crosses height h,
	// stores the resulting value in coords
	bool edgeCross(const Edge &e, double h, double2 &coords) const;


	/*  convention for labelling visits:
         _  _  _  _  _ .
        |_ |_ |_ |_ |_ |
	    |_ |_ |_ |_ |_ |
        |_ |_ |_ |_ |_ |
	    |_ |_ |_ |_ |_ |                
        |_ |_ |_ |_ |_ |
	    |_ |_ |_ |_ |_ |
		This array needs to be reset whenever a different height is queried
	*/
	double xp_;
	double yp_;		

	std::vector<double> X_;
	std::vector<double> Y_;
	std::vector<std::vector<double> > Z_; 
	std::vector<std::vector<bool2> > visits_;

};

} // namespace Terran

#endif