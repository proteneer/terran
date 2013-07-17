#ifndef UTIL_H_
#define UTIL_H_

#include <Param.h>
#include <vector>
#include <list>

namespace Util {

void matchParameters(const std::vector<Param> &p1, const std::vector<Param> &p2, double threshold);

void matchPoints(const std::vector<double> &p1, const std::vector<double> &p2, double threshold, int count = -1);

void matchPeriodicPoints(const std::vector<double> &p1, const std::vector<double> &p2, double period, double threshold, int count = -1);

}

#endif
