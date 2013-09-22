#ifndef PARAM_H_
#define PARAM_H_

namespace Terran {

struct Param {

    Param() {} ;
    Param(double p, double u, double s) : p(p), u(u), s(s) {};

    // probability weight of this mixture
    double p; 
    // mean
    double u;
    // standard deviation
    double s;
};

}
#endif
