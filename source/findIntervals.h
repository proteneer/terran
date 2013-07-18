#ifndef FIND_INTERVALS_H_
#define FIND_INTERVALS_H_

#include "MathFunctions.h"
#include "Param.h"
#include <vector>

void partitionGaussian(const std::vector<Param> &params);

void partitionPeriodicGaussian(const std::vector<Param> &params, double period = 2*PI, int images = 30);

std::vector<double> findPeriodicMaxima(const std::vector<Param> &params, double period = 2*PI, int images = 30);

std::vector<double> findPeriodicMinima(const std::vector<Param> &params, double period = 2*PI, int images = 30);

std::vector<double> findPeriodicMinimaGS(const std::vector<Param> &params, double period = 2*PI, int images = 30);

std::vector<double> findPeriodicMinimaBS(const std::vector<Param> &params, double period = 2*PI, int images = 30);

#endif 
