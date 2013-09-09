#include <vector>
#include <stdexcept>
#include <math.h>
#include <iostream>

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

using namespace std;

class Mesh {

public:
	Mesh(vector<double> X, vector<double> Y, vector<vector<double> > Z, double xp, double yp) : 
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

	// Used for defining edge locations (clockwise)
	enum Direction {UP, RIGHT, DOWN, LEFT};
		
	struct Edge {
		
		bool operator== (const Edge &e1) const {
			if(e1.coords.x != this->coords.x)
				return false;
			if(e1.coords.y != this->coords.y)
				return false;
			if(e1.extend != this->extend)
				return false;
			return true;
		}

		bool operator!= (const Edge &e1) const {
			if(e1.coords.x != this->coords.x)
				return true;
			if(e1.coords.y != this->coords.y)
				return true;
			if(e1.extend != this->extend)
				return true;
			return false;
		}

		Direction orientation() const {
			return static_cast<Direction>((extend + 3) % 4);
		};

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
	vector<vector<double2> > contour(double height) {
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
				// boundary points will never need to have their edges
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

	unsigned int xSize() const {
		return X_.size();
	}

	unsigned int ySize() const {
		return Y_.size();
	}

	void test() {
		{
			double height = 0.5;

			{
			Edge e;
			e.coords.x = 0;
			e.coords.y = 1;
			e.extend = RIGHT;


			cout << e.coords.x << " " << e.coords.y << " " << e.extend << endl;

			vector<double2> p = path(e, height);
			//cout << p.size() << endl;

			contour(height);
		}

		}
		
	}
	
private:

	vector<double2> path(Edge startEdge, double height) {
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

	bool xPeriodic() const {
		return xp_ != 0;
	};

	bool yPeriodic() const {
		return yp_ != 0;
	}

	// find the next logical edge that crosses height h
	// stores the cut location in cut
	Edge findNextEdge(Edge start, double height, double2 &cut) const {

		// TODO: Modify for four cut case
		Edge next;
		if(edgeCross(leftEdge(start), height, cut)) {
			cout << "left" << endl;
			next = leftEdge(start);
		} else if (edgeCross(forwardEdge(start), height, cut)) {
						cout << "forward" << endl;
			next = forwardEdge(start);
		} else if (edgeCross(rightEdge(start), height, cut)) {
						cout << "right" << endl;
			next = rightEdge(start);
		} else {
			throw std::runtime_error("Mesh::findNextEdge() - No valid edge found!");
		}
		return next;
	}

	// move a point in a given direction, taking into account of periodic boundary conditions
	int2 move(const int2 &p, Direction d) const {
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
	void markEdge(const Edge &u) {
		// need to adjust so we know where to mark
		int2 pos = u.coords;
		if(u.extend == LEFT || u.extend == DOWN)
			pos = move(pos, u.extend);

		if(u.extend == LEFT || u.extend == RIGHT)
			visits_[pos.x][pos.y].x = true;
		else 
			visits_[pos.x][pos.y].y = true;
	}

	Edge leftEdge(const Edge &u) const {
		Edge v;
		v.coords = u.coords;
		v.extend = static_cast<Direction>((u.extend + 3) % 4);
		return v;
	};

	Edge forwardEdge(const Edge &u) const {
		Edge v;
		v.coords = move(u.coords, u.orientation());
		v.extend = u.extend;
		return v;
	};

	Edge rightEdge(const Edge &u) const {
		Edge v;
		v.coords = move(u.coords, u.orientation());
		v.coords = move(v.coords, u.extend);
		v.extend = static_cast<Direction>((u.extend + 1) % 4);
		return v;
	};

	// returns true if the edge crosses height h
	bool edgeCross(const Edge &e, double h) const {
		int2 a = e.coords;

		cout << "edgeCross properties: " << endl;
		cout << a.x << " " << a.y << endl;

		int2 b = move(a, e.extend);

		cout << b.x << " " << b.y << endl;
		
		double z1 = Z_[a.x][a.y];
		double z2 = Z_[b.x][b.y];

		if((z1 < h && z2 > h) || (z1 > h && z2 < h))
			return true;
		else
			return false;
	};

	// returns true if the edge crosses height h,
	// stores the resulting value in coords
	bool edgeCross(const Edge &e, double h, double2 &coords) const {
		int2 a = e.coords;
		int2 b = move(a, e.extend);
		double3 ap, bp;
		ap.x = X_[a.x];
		ap.y = Y_[a.y];
		ap.z = Z_[a.x][a.y];

		bp.x = X_[b.x];
		bp.y = X_[b.y];
		bp.z = Z_[b.x][b.y];

		if((ap.z < h && bp.z > h) || (ap.z > h && bp.z < h)) {
			double dx = fabs(ap.x-bp.x);
			double dy = fabs(ap.y-bp.y);
			if(xPeriodic())
				dx -= floor(dx/xp_+0.5)*xp_;
			if(yPeriodic())
				dy -= floor(dy/yp_+0.5)*yp_;		
			const double dz = fabs(ap.z-bp.z);
			const double dh = fabs(ap.z-h);
			coords.x = ap.x + dx*(dh/dz);
			coords.y = ap.y + dy*(dh/dz);
			return true;
		} else {
			return false;
		}
	}


	/*  convention for labelling visitss:
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

	vector<double> X_;
	vector<double> Y_;
	vector<vector<double> > Z_; 
	vector<vector<bool2> > visits_;


};


int main() {

/*
	const int numX = 4;
	const int numY = 3;

	double xp = 0;
	double yp = 0;

	vector<double> X(numX);
	vector<double> Y(numY);
	vector<vector<double> > Z(numX, vector<double>(numY));

	X[0] = 0;
	X[1] = 1;
	X[2] = 2;
	X[3] = 3;
	Y[0] = 0;
	Y[1] = 1;
	Y[2] = 2;
	for(int x=0; x < numX; x++) {
		for(int y=0; y < numY; y++) {
			Z[x][y] = 0;
		}
	}
	Z[1][1] = 1;
	Z[2][1] = 1;

	Mesh m(X,Y,Z,xp,yp);

	try{
		m.test();
	} catch(const exception &e) {
		cout << e.what() << endl;
	}

*/

	const int numX = 3;
	const int numY = 2;

	double xp = 4;
	double yp = 3;

	vector<double> X(numX);
	vector<double> Y(numY);
	vector<vector<double> > Z(numX, vector<double>(numY));

	X[0] = 0;
	X[1] = 1;
	X[2] = 2;
	Y[0] = 0;
	Y[1] = 1;
	for(int x=0; x < numX; x++) {
		for(int y=0; y < numY; y++) {
			Z[x][y] = 0;
		}
	}
	Z[1][1] = 1;
	Z[2][1] = 1;

	Mesh m(X,Y,Z,xp,yp);

	try{
		m.test();
	} catch(const exception &e) {
		cout << e.what() << endl;
	}



}