#ifndef UTIL_H_
#define UTIL_H_

#include <Param.h>
#include <vector>
#include <list>

namespace Util {

void matchParameters(const std::vector<Param> &p1, const std::vector<Param> &p2, double threshold = 0.1);

}

#endif
