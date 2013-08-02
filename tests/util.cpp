#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <math.h>
#include <iostream>
#include <fstream>

#include "util.h"
#include "MathFunctions.h"

using namespace std;
using Terran::Param;
using Terran::fabsp;


void Util::plotPeriodicGaussian(const vector<Param> &params, double period, string filename) {
    ofstream f1(filename.c_str());
    for(double xn = -period/2; xn < period/2; xn += 0.01) {
        f1 << xn << " " << periodicGaussianMixture(params, xn, period, 50) << endl;
    }
    ofstream f2((filename + "Dx").c_str());
    for(double xn = -period/2; xn < period/2; xn += 0.01) {
        f2 << xn << " " << periodicGaussianMixtureDx(params, xn, period, 50) << endl;
    }
}

// Naive O(N^2) method
void Util::matchParameters(const vector<Param> &p1, const vector<Param> &p2, double threshold) {
    if(p1.size() != p2.size()) {
        throw(std::runtime_error("Unequal Parameter Size!"));
    }
    list<Param> list1;
    list<Param> list2;
    for(int i=0; i<p1.size(); i++) {
        list1.push_back(p1[i]);
        list2.push_back(p2[i]);
    }
    for(list<Param>::iterator it = list1.begin(); it != list1.end(); it++) {
        bool found = false;
        for(list<Param>::iterator jt = list2.begin(); jt != list2.end(); jt++) {
            double pdiff = fabs(it->p - jt->p);
            double udiff = fabs(it->u - jt->u);
            double sdiff = fabs(it->s - jt->s);
            if(pdiff < threshold && udiff < threshold && sdiff < threshold) {
                list2.erase(jt);
                found = true;
                break;
            }
        }
        if(!found) 
            throw(std::runtime_error("Unpaired parameter found!"));
    }
}

void Util::matchPoints(const vector<double> &p1, const vector<double> &p2, double threshold, int count) {
    if(p1.size() != p2.size()) {
        stringstream ss;
        ss << "Unequal vector size!";
        if(count != -1)
            ss << " Test: " << count;
        throw(std::runtime_error("Unequal Vector Size!"));
    }
    vector<double> list1(p1);
    vector<double> list2(p2);
    sort(list1.begin(), list1.end());
    sort(list2.begin(), list2.end());
    for(int i=0; i < list1.size(); i++) {
        if(fabs(list1[i]-list2[i]) > threshold) { 
            stringstream ss;
            ss << "Nonmatching vector element found! ";
            if(count != -1)
                ss << " Test: " << count;
            throw(std::runtime_error(ss.str()));
        }
    }
}

void Util::matchPeriodicPoints(const vector<double> &p1, const vector<double> &p2, double period, double threshold, int count) {
    if(p1.size() != p2.size()) {
        stringstream ss;
        ss << "Unequal vector size!";
        if(count != -1)
            ss << " Test: " << count;
        throw(std::runtime_error("Unequal Vector Size!"));
    }
    vector<double> list1(p1);
    vector<double> list2(p2);
    sort(list1.begin(), list1.end());
    sort(list2.begin(), list2.end());
    for(int i=0; i < list1.size(); i++) {
        if(fabsp(list1[i],list2[i],period) > threshold) { 
            stringstream ss;
            ss << "Nonmatching vector element found! ";
            if(count != -1)
                ss << " Test: " << count;
            throw(std::runtime_error(ss.str()));
        }
    }
}
