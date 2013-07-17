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
        if(results.size() != 1)
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
        if(results.size() != 1)
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
        if(results.size() != 2)
            throw(std::runtime_error("Finding Maxima 4 Failed: Too many components!"));
        if(fabs(results[0]-(-0.847)) > tol)
            throw(std::runtime_error("Finding Maxima 4 Failed."));
        if(fabs(results[1]-(1.670)) > tol)
            throw(std::runtime_error("Finding Maxima 4 Failed."));
    }

    // 4 Fused into 3
    {
        vector<Param> p;
        p.push_back(Param(1.0/4.0,-2.2, 0.73));
        p.push_back(Param(1.0/4.0,-0.8, 0.44));
        p.push_back(Param(1.0/4.0, 0.8, 0.61));
        p.push_back(Param(1.0/4.0, 1.8, 0.45));
        vector<double> results = findPeriodicMaxima(p);
        if(results.size() != 3)
            throw(std::runtime_error("Finding Maxima 5 Failed: Too many components!"));
        if(fabs(results[0]-(-2.141)) > tol)
            throw(std::runtime_error("Finding Maxima 5 Failed."));
        if(fabs(results[1]-(-0.835)) > tol)
            throw(std::runtime_error("Finding Maxima 5 Failed."));
        if(fabs(results[2]-( 1.671)) > tol)
            throw(std::runtime_error("Finding Maxima 5 Failed."));
    }

    // 4 Fused into 3, Almost saddle-point like
    {
        vector<Param> p;
        p.push_back(Param(1.0/4.0,-2.129, 0.73));
        p.push_back(Param(1.0/4.0,-0.8, 0.44));
        p.push_back(Param(1.0/4.0, 0.8, 0.61));
        p.push_back(Param(1.0/4.0, 1.8, 0.45));
        vector<double> results = findPeriodicMaxima(p);
        if(results.size() != 3)
            throw(std::runtime_error("Finding Maxima 5 Failed: Too many components!"));
        if(fabs(results[0]-(-1.968)) > tol)
            throw(std::runtime_error("Finding Maxima 5 Failed."));
        if(fabs(results[1]-(-0.843)) > tol)
            throw(std::runtime_error("Finding Maxima 5 Failed."));
        if(fabs(results[2]-( 1.670)) > tol)
            throw(std::runtime_error("Finding Maxima 5 Failed."));
    }

    // 7 Fused into 2, 1 Saddle-point
    {
        vector<Param> p;
        p.push_back(Param(1.0/7.0, -2.9, 0.3));
        p.push_back(Param(1.0/7.0, -2.1, 0.6));
        p.push_back(Param(1.0/7.0, -1.2, 1.2));
        p.push_back(Param(1.0/7.0, -0.2, 0.7));
        p.push_back(Param(1.0/7.0,  1.0, 0.4));
        p.push_back(Param(1.0/7.0,  1.3, 0.6));
        p.push_back(Param(1.0/7.0,  1.9, 0.4));
        vector<double> results = findPeriodicMaxima(p);
        if(results.size() != 2)
            throw(std::runtime_error("Finding Maxima 6 Failed: Too many components!"));
        if(fabs(results[0]-(-2.847)) > tol)
            throw(std::runtime_error("Finding Maxima 6 Failed."));
        if(fabs(results[1]-( 1.117)) > tol)
            throw(std::runtime_error("Finding Maxima 6 Failed."));
    }

    // Bimodal near boundary
    {
        vector<Param> p;
        p.push_back(Param(1.0/7.0, 3.1, 0.5));
        p.push_back(Param(1.0/7.0, -2.0, 0.5));
        partitionPeriodicGaussian(p);
        vector<double> results = findPeriodicMaxima(p);
        if(results.size() != 2)
            throw(std::runtime_error("Finding Maxima 6 Failed: Too many components!"));
        if(fabs(results[0]-(-3.074)) > tol)
            throw(std::runtime_error("Finding Maxima 6 Failed."));
        if(fabs(results[1]-(-2.109)) > tol)
            throw(std::runtime_error("Finding Maxima 6 Failed."));
    }
}

void testFindPeriodicMinima() {

    // Unimodal
    double tol = 1e-3;
    vector<Param> p;
    p.push_back(Param(1.0, 0.0, 0.5));
    partitionPeriodicGaussian(p);
    vector<double> results = findPeriodicMinima(p);
    if(results.size() != 1)
        throw(std::runtime_error("Finding Minima 1 Failed: Wrong number of minima"));
    if(fabs(results[0]-3.14159) > tol && fabs(results[0]-0) > tol)
        throw(std::runtime_error("Finding Minima 1 Failed."));  

}
/*
        for(int i=0; i< results.size(); i++) {
            cout << results[i] << endl;
        }
*/
int main() {
    try {
        testFindPeriodicMaxima();
        testFindPeriodicMinima();
        cout << "done" << endl;
    } catch(const std::exception &e) {
        cout << e.what() << endl;
    }
}
