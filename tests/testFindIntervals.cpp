#include <findIntervals.h>
#include <vector>
#include <iostream>    
#include <stdexcept>

using namespace std;


void testFindPeriodicMaxima() {
    double tol = 1e-3;

    // Bimodal
    {
        vector<Param> p;
        p.push_back(Param(0.6,  1.0, 0.5));
        p.push_back(Param(0.4, -2.3, 0.5));
        vector<double> results = findPeriodicMaxima(p);
        if(fabs(results[0]-1.0) > tol || fabs(results[1]-(-2.3)) > tol)
            throw(std::runtime_error("Finding Maxima 0 Failed."));
    }

    // Bimodal
    {
        vector<Param> p;
        p.push_back(Param(0.5, -0.2, 0.5));
        p.push_back(Param(0.5,  1.6, 0.6));
        vector<double> results = findPeriodicMaxima(p);
        if(fabs(results[0]-(-0.188)) > tol || fabs(results[1]-(1.595)) > tol)
            throw(std::runtime_error("Finding Maxima 1 Failed."));

    }

    // 2 Fused into 1
    {
        vector<Param> p;
        p.push_back(Param(0.5, 1.2, 0.5));
        p.push_back(Param(0.5, 1.6, 0.4));
        vector<double> results = findPeriodicMaxima(p);
        if(results.size() > 1)
            throw(std::runtime_error("Finding Maxima 2 Failed: Too many components!"));
        if(fabs(results[0]-1.473)>tol)
            throw(std::runtime_error("Finding Maxima 2 Failed: Wrong maxima"));
    }

    // 3 Fused into 1
    {
        vector<Param> p;
        p.push_back(Param(1.0/3.0,-0.2, 0.5));
        p.push_back(Param(1.0/3.0, 0.6, 0.6));
        p.push_back(Param(1.0/3.0, 1.6, 0.7));
        vector<double> results = findPeriodicMaxima(p);
        if(results.size() > 1)
            throw(std::runtime_error("Finding Maxima 3 Failed: Too many components!"));
        if(fabs(results[0]-(0.117)) > tol)
            throw(std::runtime_error("Finding Maxima 3 Failed."));
    }

    // 4 Fused into 2
    {
        vector<Param> p;
        p.push_back(Param(1.0/4.0,-2.1, 0.73));
        p.push_back(Param(1.0/4.0,-0.8, 0.44));
        p.push_back(Param(1.0/4.0, 0.8, 0.61));
        p.push_back(Param(1.0/4.0, 1.8, 0.45));
        vector<double> results = findPeriodicMaxima(p);
        if(results.size() > 2)
            throw(std::runtime_error("Finding Maxima 4 Failed: Too many components!"));
        if(fabs(results[0]-(-0.847)) > tol)
            throw(std::runtime_error("Finding Maxima 4 Failed."));
        if(fabs(results[1]-(1.670)) > tol)
            throw(std::runtime_error("Finding Maxima 4 Failed."));
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


