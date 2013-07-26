#ifndef UTIL_H_
#define UTIL_H_

#include <Param.h>
#include <vector>
#include <list>

namespace Util {

void plotPeriodicGaussian(const std::vector<Terran::Param> &params, double period, int images);

void matchParameters(const std::vector<Terran::Param> &p1, const std::vector<Terran::Param> &p2, double threshold);

void matchPoints(const std::vector<double> &p1, const std::vector<double> &p2, double threshold, int count = -1);

void matchPeriodicPoints(const std::vector<double> &p1, const std::vector<double> &p2, double period, double threshold, int count = -1);

}

#endif
