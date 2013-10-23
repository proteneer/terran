#include "KMeans.h"
#include <iostream>

using namespace std;
using namespace Terran;

// 2 clusters.
void testSimpleNonPeriodic2D() {
   
    const int N = 10;
    const int D = 2;
    
    vector<vector<double> > dataset(N, vector<double>(2, 0));
    
    dataset[0][0] = -1.0;
    dataset[0][1] = -1.0;
    dataset[1][0] = -1.1;
    dataset[1][1] = -0.9;
    dataset[2][0] = -0.8;
    dataset[2][1] = -1.4;
    dataset[3][0] = -1.2;
    dataset[3][1] = -1.3;
    dataset[4][0] = -1.1;
    dataset[4][1] = -0.7;
   
    dataset[5][0] = 1.2;
    dataset[5][1] = 1.2;
    dataset[6][0] = 1.4;
    dataset[6][1] = 1.6; 
    dataset[7][0] = 0.8;
    dataset[7][1] = 1.1;
    dataset[8][0] = 1.5;
    dataset[8][1] = 1.0;
    dataset[9][0] = 0.9;
    dataset[9][1] = 1.1;
    
    vector<bool> period(2, false);
    
    KMeans km(dataset, period);
    
    vector<vector<double> > centers(2, vector<double>(2, 0));
    centers[0][0] = -1.3;
    centers[0][1] = -1.2;
    centers[1][0] = 1.1;
    centers[1][1] = 1.4;
     
    vector<int> assign = km.run(centers);
    
	for(int i=0; i < 5; i++) {
		if(assign[i] != 0) {
			throw(std::runtime_error("testSimpleNonPeriodic2D() - Wrong assignment"));
		}
	}

	for(int i=5; i < 10; i++) {
		if(assign[i] != 1) {
			throw(std::runtime_error("testSimpleNonPeriodic2D() - Wrong assignment"));
		}
	}
}

void testSimplePeriodic2D() {
   
    const int N = 10;
    const int D = 2;
    
    vector<vector<double> > dataset(N, vector<double>(2, 0));
    
    dataset[0][0] = -3.1;
    dataset[0][1] = -3.1;
    dataset[1][0] = -2.8;
    dataset[1][1] =  3.1;
    dataset[2][0] =  -2.9;
    dataset[2][1] =  -3.0;
    dataset[3][0] =  3.2;
    dataset[3][1] = -3.1;
    dataset[4][0] =  3.1;
    dataset[4][1] =  3.05;
   
    dataset[5][0] = 1.2;
    dataset[5][1] = 1.2;
    dataset[6][0] = 1.4;
    dataset[6][1] = 1.6; 
    dataset[7][0] = 0.8;
    dataset[7][1] = 1.1;
    dataset[8][0] = 1.5;
    dataset[8][1] = 1.0;
    dataset[9][0] = 0.9;
    dataset[9][1] = 1.1;
    
    vector<bool> period(2, true);
    
    KMeans km(dataset, period);
    
    vector<vector<double> > centers(2, vector<double>(2, 0));
    centers[0][0] = 3.0;
    centers[0][1] = 2.9;
    centers[1][0] = 1.1;
    centers[1][1] = 1.4;
     
    vector<int> assign = km.run(centers);
    
	for(int i=0; i < 5; i++) {
		if(assign[i] != 0) {
			throw(std::runtime_error("testSimplePeriodic2D() - Wrong assignment"));
		}
	}

	for(int i=5; i < 10; i++) {
		if(assign[i] != 1) {
			throw(std::runtime_error("testSimplePeriodic2D() - Wrong assignment"));
		}
	}
}

int main() {
	try {
		testSimpleNonPeriodic2D();
		testSimplePeriodic2D();
	} catch(const std::exception &e) {
		cout << e.what() << endl;
	}
}