#include "util.h"
#include <stdexcept>
#include <math.h>
#include <iostream>


using namespace std;

// Naive O(N^2) method
void Util::matchParameters(const vector<Param> &p1, const vector<Param> &p2, double threshold) {
    if(p1.size() != p2.size()) {
        throw(std::runtime_error("Invalid Parameter Size!"));
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
