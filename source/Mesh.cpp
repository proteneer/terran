#include "Mesh.h"
#include <stdexcept>
#include <math.h>
#include <iostream>

using namespace std;

namespace Terran {

Mesh::Mesh(vector<double> X, vector<double> Y, vector<vector<double> > Z, double xp, double yp) : 
	xp_(xp),
	yp_(yp),
	X_(X),
	Y_(Y),
	Z_(Z),
	visits_(X.size(), vector<bool2>(Y.size())) {

	if(X_.size() != Z.size()) {
		throw(std::runtime_error("Mesh::Mesh() - X size and Z size do not match!"));
	}

	if(Z.size() > 0) {
		if(Y_.size() != Z[0].size()) 
			throw(std::runtime_error("Mesh::Mesh() - Y size and Z size do not match!"));
	}

};	

bool Mesh::Edge::operator== (const Edge &e1) const {
	if(e1.coords.x != this->coords.x)
		return false;
	if(e1.coords.y != this->coords.y)
		return false;
	if(e1.extend != this->extend)
		return false;
	return true;
}

bool Mesh::Edge::operator!= (const Edge &e1) const {
	if(e1.coords.x != this->coords.x)
		return true;
	if(e1.coords.y != this->coords.y)
		return true;
	if(e1.extend != this->extend)
		return true;
	return false;
}

Mesh::Direction Mesh::Edge::orientation() const {
	return static_cast<Direction>((extend + 3) % 4);
};

vector<vector<double2> > Mesh::contour(double height) {
	// initialize visits array
	for(int x=0; x < visits_.size(); x++) {
		for(int y=0; y < visits_[x].size(); y++) {
			visits_[x][y].x = 0;
			visits_[x][y].y = 0;
		}
	}

	vector<vector<double2> > paths;
	for(int x=0; x < visits_.size(); x++) {
		for(int y=0; y < visits_[x].size(); y++) {

			// if the domain is aperiodic then the
			// boundary points will never need to have their edges searched
			// since we assume it is not possible for an edge to be on an edge
			if(!xPeriodic() && (x == visits_.size()-1))
				break;
			if(!yPeriodic() && (y == visits_[x].size()-1))
				break;

			// start a path only if the edge has not been visited and it is a cutting edge
			if(!visits_[x][y].x) {
				Edge edge;
				edge.coords.x = x;
				edge.coords.y = y;
				edge.extend = RIGHT;
				if(edgeCross(edge, height))
					paths.push_back(path(edge, height));
			}

			if(!visits_[x][y].y) {
				Edge edge;
				edge.coords.x = x;
				edge.coords.y = y;
				edge.extend = UP;
				if(edgeCross(edge, height))
					paths.push_back(path(edge, height));
			}
		}
	}
	return paths;
}

unsigned int Mesh::xSize() const {
	return X_.size();
}

unsigned int Mesh::ySize() const {
	return Y_.size();
}

vector<double2> Mesh::path(Edge startEdge, double height) {
	// mark the start edge as having been visited
	markEdge(startEdge);
	// the cut point is not computed as this will be taken care of
	// at the end of the loop
	vector<double2> list;
	Edge nextEdge = startEdge;

	// todo: deal with boundary edges later
	do {
		double2 cut;
		nextEdge = findNextEdge(nextEdge, height, cut);
		list.push_back(cut);
		markEdge(nextEdge);
	} while(nextEdge != startEdge);

	return list;
}

bool Mesh::xPeriodic() const {
	return xp_ != 0;
};

bool Mesh::yPeriodic() const {
	return yp_ != 0;
}

// find the next logical edge that crosses height h
// stores the cut location in cut
Mesh::Edge Mesh::findNextEdge(Mesh::Edge start, double height, double2 &cut) const {

	// TODO: Modify for four cut case
	Edge next;
	if(edgeCross(leftEdge(start), height, cut)) {
		next = leftEdge(start);
	} else if (edgeCross(forwardEdge(start), height, cut)) {
		next = forwardEdge(start);
	} else if (edgeCross(rightEdge(start), height, cut)) {
		next = rightEdge(start);
	} else {
		throw std::runtime_error("Mesh::findNextEdge() - No valid edge found!");
	}
	return next;
}

// move a point in a given direction, taking into account of periodic boundary conditions
int2 Mesh::move(const int2 &p, Direction d) const {
	// these must be ints because later on we do -1 % size!

	int2 q = p;

	if(d == UP) {
		q.y += 1;
	} else if(d == RIGHT) {
		q.x += 1;
	} else if(d == DOWN) {
		q.y -= 1;
	} else if(d == LEFT) {
		q.x -= 1;
	} else {
		throw(std::runtime_error("bad point extension!"));
	}

	// take into account of periodic conditions here
	// xSize() is unsigned, so need to add to make sure a negative number
	// isn't left of the mod operator
	if(xPeriodic()) {
		q.x = (q.x+xSize()) % xSize();
	} else if(q.x > xSize() - 1 || q.x < 0) {
		throw(std::runtime_error("move() - x out of bounds"));
	};

	if(yPeriodic()) {
		q.y = (q.y+ySize()) % ySize();
	} else if(q.y > ySize() - 1 || q.y < 0) {
		throw(std::runtime_error("move() - y out of bounds"));
	};

	return q;
}

// mark the edge as having been visited in the visits_ array
void Mesh::markEdge(const Edge &u) {
	// need to adjust so we know where to mark
	int2 pos = u.coords;
	if(u.extend == LEFT || u.extend == DOWN)
		pos = move(pos, u.extend);

	if(u.extend == LEFT || u.extend == RIGHT)
		visits_[pos.x][pos.y].x = true;
	else 
		visits_[pos.x][pos.y].y = true;
}

Mesh::Edge Mesh::leftEdge(const Edge &u) const {
	Edge v;
	v.coords = u.coords;
	v.extend = static_cast<Direction>((u.extend + 3) % 4);
	return v;
};

Mesh::Edge Mesh::forwardEdge(const Edge &u) const {
	Edge v;
	v.coords = move(u.coords, u.orientation());
	v.extend = u.extend;
	return v;
};

Mesh::Edge Mesh::rightEdge(const Edge &u) const {
	Edge v;
	v.coords = move(u.coords, u.orientation());
	v.coords = move(v.coords, u.extend);
	v.extend = static_cast<Direction>((u.extend + 1) % 4);
	return v;
};

// returns true if the edge crosses height h
bool Mesh::edgeCross(const Edge &e, double h) const {
	int2 a = e.coords;
	int2 b = move(a, e.extend);
	double z1 = Z_[a.x][a.y];
	double z2 = Z_[b.x][b.y];

	if((z1 < h && z2 > h) || (z1 > h && z2 < h))
		return true;
	else
		return false;
};

// returns true if the edge crosses height h,
// stores the resulting value in coords
bool Mesh::edgeCross(const Edge &e, double h, double2 &coords) const {
	
	// calculate the edge values
	int2 first = e.coords;
	int2 second = move(first, e.extend);
	int2 a,b;
	// if extends left, we need to reverse
	if(e.extend == LEFT || e.extend == DOWN) {
		a = second;
		b = first;
	} else {
		a = first;
		b = second;
	}
	// which one is left, which one is right?
	double3 ap, bp;
	ap.x = X_[a.x];
	ap.y = Y_[a.y];
	ap.z = Z_[a.x][a.y];
	bp.x = X_[b.x];
	bp.y = Y_[b.y];
	bp.z = Z_[b.x][b.y];
	if((ap.z < h && bp.z > h) || (ap.z > h && bp.z < h)) {
		double dx = ap.x-bp.x;
		double dy = ap.y-bp.y;
		if(xPeriodic())
			dx -= floor(dx/xp_+0.5)*xp_;
		if(yPeriodic())
			dy -= floor(dy/yp_+0.5)*yp_;		
		const double dz = fabs(ap.z-bp.z);
		const double dh = fabs(ap.z-h);
		coords.x = ap.x + fabs(dx)*(dh/dz);
		coords.y = ap.y + fabs(dy)*(dh/dz);
		return true;
	} else {
		return false;
	}
}

}

