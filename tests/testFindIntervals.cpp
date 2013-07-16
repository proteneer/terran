#include <findIntervals.h>
#include <vector>
#include <iostream>    
#include <stdexcept>

using namespace std;


void testFindPeriodicMaxima() {

    double tol = 1e-3;
    {
        vector<Param> p;
        p.push_back(Param(0.6,  1.0, 0.5));
        p.push_back(Param(0.4, -2.3, 0.5));
        partitionPeriodicGaussian(p);
        vector<double> results = findPeriodicMaxima(p);
        if(fabs(results[0]-1.0) > tol || fabs(results[1]-(-2.3)) > tol)
            throw(std::runtime_error("Finding Maxima 0 Failed."));
    }
    {
        vector<Param> p;
        p.push_back(Param(0.5, -0.2, 0.5));
        p.push_back(Param(0.5,  1.6, 0.6));
        partitionPeriodicGaussian(p);
        vector<double> results = findPeriodicMaxima(p);
        if(fabs(results[0]-(-0.188)) > tol || fabs(results[1]-(1.595)) > tol)
            throw(std::runtime_error("Finding Maxima 1 Failed."));

    }


}

int main() {
    try {
        testFindPeriodicMaxima();
        cout << "done" << endl;
    } catch(const std::exception &e) {
        cout << e.what() << endl;
    }
}


