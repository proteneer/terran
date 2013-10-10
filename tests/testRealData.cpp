#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>

#include <Terran.h>

using namespace std;
using namespace Terran;

vector<vector<double> > readData(string filename) {
	ifstream in_file(filename.c_str());
	vector<vector<double> > data;
	string line;
	while(getline(in_file, line)) {
		vector<double> point;
		stringstream ss(line);
		//cout << ss.str() << endl;
		while(!ss.eof()) {
			double angle;
			// convert to radians
			ss >> angle;
			angle = angle*(PI/180.0);
			point.push_back(angle);
		}
		data.push_back(point);
	}
	return data;
}

void print_vector(const vector<double> &v, string filename) {
	ofstream data(filename);
	for(int i=0; i < v.size(); i++) {
		data << v[i] << endl;
	}
}

int main() {

	vector<vector<double> > data = readData("dihed.txt");

	if(data.size() > 0) {
		cout << data.size() << "x" << data[0].size() << endl;
	} else {
		cout << "data empty" << endl;
	}
	
	vector<int> period(data[0].size(), 1);

	Cluster cc(data, period);
	
	/*
	for(int d = 0; d < cc.getNumDimensions(); d++) {
        cc.partition(d);
        vector<double> partitions = cc.getPartition(d);
		cout << d << " " << partitions.size() << endl;
    }*/
	ClusterTree ct(data, period);

	ct.setCurrentCluster();
    for(int i=0 ; i < ct.getCurrentCluster().getNumDimensions(); i++) {
        ct.partitionCurrentCluster(i);
        vector<double> partitions = ct.getCurrentCluster().getPartition(i);
		cout << i << " " << partitions.size() << endl;
    }
	ct.divideCurrentCluster();

	ofstream log("assign.txt");
	vector<int> out = ct.getAssignment();
	for(int i=0; i < out.size(); i++) {
		log << out[i] << endl;
	}



		/*
	for(int i=1; i < data.size(); i++) {
		if(data[i].size() != data[i-1].size()) {
			throw(std::runtime_error("data is not rectangular!"));
		}
	}

	// marginalize each dimension
	int numDims = data[0].size();
	int numPoints = data.size();

	for(int d=0; d < numDims; d++) {

	    vector<double> subset;
		for(int j=0; j < numPoints; j++) {
			subset.push_back(data[j][d]);
		}

		random_shuffle(subset.begin(), subset.end());
		subset.resize(3000);

		stringstream filename;
		filename << "data" << d;
		cout << filename.str() << endl;
		print_vector(subset, filename.str());

	}
	*/

}