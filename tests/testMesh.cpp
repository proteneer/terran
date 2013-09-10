#include <Mesh.h>
#include <vector>
#include <iostream>
#include <stdexcept>

using namespace std;
using namespace Terran;

void testPeriodicMesh() {

	const int numX = 3;
	const int numY = 2;

	double xp = 3;
	double yp = 2;

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

	double height = 0.5;
	vector<vector<double2> > foo = m.contour(0.5);

	if(foo.size() != 1) {
		throw(std::runtime_error("Wrong number of paths!"));
	}

	if(foo[0].size() != 6) {
		throw(std::runtime_error("Wrong number of points on path!"));
	}

};

void testMesh() {
	
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

	double height = 0.5;
	vector<vector<double2> > foo = m.contour(0.5);

	if(foo.size() != 1) {
		throw(std::runtime_error("Wrong number of paths!"));
	}

	if(foo[0].size() != 6) {
		throw(std::runtime_error("Wrong number of points on path!"));
	}
}

int main() {

	try{
		testPeriodicMesh();
		testMesh();
	} catch(const exception &e) {
		cout << e.what() << endl;
	}
} // namespace Terran