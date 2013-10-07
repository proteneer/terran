#ifndef PARAM_H_
#define PARAM_H_

#include <ostream>
#include <vector>

namespace Terran {

struct Param {

    Param() : p(0), u(0), s(0) {} ;
    Param(double p, double u, double s) : p(p), u(u), s(s) {};

    // probability weight of this mixture
    double p; 
    // mean
    double u;
    // standard deviation
    double s;
};

inline std::ostream& operator<<(std::ostream& os, const std::vector<Param>& p) {
	for(int i=0; i < p.size(); i++) {
		os << p[i].p << " " << p[i].u << " " << p[i].s << std::endl;
	}
	return os;
}

}
#endif
