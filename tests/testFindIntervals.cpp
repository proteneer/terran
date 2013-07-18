#include <findIntervals.h>
#include <vector>
#include <iostream>    
#include <stdexcept>

#include "util.h"

using namespace std;

void testFindPeriodicMaxima() {
    double tol = 1e-3;
    int testCount = 0;

    // Bimodal
    {
        vector<Param> p;
        p.push_back(Param(0.6,  1.0, 0.5));
        p.push_back(Param(0.4, -2.3, 0.5));
        vector<double> results = findPeriodicMaxima(p);
        vector<double> truth;
        truth.push_back(-2.3);
        truth.push_back( 1.0);
        Util::matchPoints(results, truth, tol, testCount++);
    }

    // Bimodal
    {
        vector<Param> p;
        p.push_back(Param(0.5, -0.2, 0.5));
        p.push_back(Param(0.5,  1.6, 0.6));
        vector<double> results = findPeriodicMaxima(p);
        vector<double> truth;
        truth.push_back(-0.188);
        truth.push_back( 1.595);
        Util::matchPoints(results, truth, tol, testCount++);
    }

    // 2 Fused into 1
    {
        vector<Param> p;
        p.push_back(Param(0.5, 1.2, 0.5));
        p.push_back(Param(0.5, 1.6, 0.4));
        vector<double> results = findPeriodicMaxima(p);
        vector<double> truth;
        truth.push_back(1.473);
        Util::matchPoints(results, truth, tol, testCount++);
    }

    // 3 Fused into 1
    {
        vector<Param> p;
        p.push_back(Param(1.0/3.0,-0.2, 0.5));
        p.push_back(Param(1.0/3.0, 0.6, 0.6));
        p.push_back(Param(1.0/3.0, 1.6, 0.7));
        vector<double> results = findPeriodicMaxima(p);
        vector<double> truth;
        truth.push_back(0.117);
        Util::matchPoints(results, truth, tol, testCount++);
    }

    // 4 Fused into 2
    {
        vector<Param> p;
        p.push_back(Param(1.0/4.0,-2.1, 0.73));
        p.push_back(Param(1.0/4.0,-0.8, 0.44));
        p.push_back(Param(1.0/4.0, 0.8, 0.61));
        p.push_back(Param(1.0/4.0, 1.8, 0.45));
        vector<double> results = findPeriodicMaxima(p);
        vector<double> truth;
        truth.push_back(-0.847);
        truth.push_back( 1.670);
        Util::matchPoints(results, truth, tol, testCount++);
    }

    // 4 Fused into 3
    {
        vector<Param> p;
        p.push_back(Param(1.0/4.0,-2.2, 0.73));
        p.push_back(Param(1.0/4.0,-0.8, 0.44));
        p.push_back(Param(1.0/4.0, 0.8, 0.61));
        p.push_back(Param(1.0/4.0, 1.8, 0.45));
        vector<double> results = findPeriodicMaxima(p);
        vector<double> truth;
        truth.push_back(-2.141);
        truth.push_back(-0.835);
        truth.push_back( 1.671);
        Util::matchPoints(results, truth, tol, testCount++);
    }

    // 4 Fused into 3, Almost saddle-point like
    {
        vector<Param> p;
        p.push_back(Param(1.0/4.0,-2.129, 0.73));
        p.push_back(Param(1.0/4.0,-0.8, 0.44));
        p.push_back(Param(1.0/4.0, 0.8, 0.61));
        p.push_back(Param(1.0/4.0, 1.8, 0.45));
        vector<double> results = findPeriodicMaxima(p);
        vector<double> truth;
        truth.push_back(-1.968);
        truth.push_back(-0.843);
        truth.push_back( 1.670);
        Util::matchPoints(results, truth, tol, testCount++);
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
        vector<double> truth;
        truth.push_back(-2.847);
        truth.push_back( 1.117);
        Util::matchPoints(results, truth, tol, testCount++);
    }

    // Bimodal near boundary
    {
        vector<Param> p;
        p.push_back(Param(1.0/7.0, 3.1, 0.5));
        p.push_back(Param(1.0/7.0, -2.0, 0.5));
        partitionPeriodicGaussian(p);
        vector<double> results = findPeriodicMaxima(p);
        vector<double> truth;
        truth.push_back(-3.074);
        truth.push_back(-2.109);
        Util::matchPoints(results, truth, tol, testCount++);
    }
}

void testFindPeriodicMinima() {

    const double tol = 1e-3;
    int testCount = 0;
    const double period = 2*PI;

/*
    // Unimodal
    {
        vector<Param> p;
        p.push_back(Param(1.0, 0.0, 0.5));
        vector<double> results = findPeriodicMinimaBS(p);
        vector<double> truth;
        truth.push_back(-period/2);
        Util::matchPeriodicPoints(results, truth, period, tol, testCount++);
    }

    // Unimodal
    {
        vector<Param> p;
        p.push_back(Param(1.0, 0.5, 1.0));
        vector<double> results = findPeriodicMinimaBS(p);
        vector<double> truth;
        truth.push_back(-2.642);
        Util::matchPeriodicPoints(results, truth, period, tol, testCount++);
    }
*/
    // Bimodal
    {
        vector<Param> p;
        p.push_back(Param(0.3, 0.5, 1.0));
        p.push_back(Param(0.7, -1.5, 0.5));
        vector<double> results = findPeriodicMinimaBS(p);
        vector<double> truth;
        partitionPeriodicGaussian(p);
        truth.push_back(-0.091);
        truth.push_back( 3.071);
        Util::matchPeriodicPoints(results, truth, period, tol, testCount++);
    }

/*
    // Bimodal
    {
        double tol2 = 1e-2;
        vector<Param> p;
        p.push_back(Param(0.5, 1.1, 0.3));
        p.push_back(Param(0.5, -2.9, 0.4));
        partitionPeriodicGaussian(p);
        vector<double> results = findPeriodicMinima(p);
        cout << results[0] << " " << results[1] << endl;
        vector<double> truth;
        truth.push_back(-0.632);
        truth.push_back( 2.108);
        Util::matchPeriodicPoints(results, truth, period, tol2, testCount++);
    }
    */


}
/*
        for(int i=0; i< results.size(); i++) {
            cout << results[i] << endl;
        }
*/
int main() {
    try {
        //testFindPeriodicMaxima();
        testFindPeriodicMinima();
        cout << "done" << endl;
    } catch(const std::exception &e) {
        cout << e.what() << endl;
    }
}
