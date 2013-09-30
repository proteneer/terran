#ifndef METHODS_GAUSSIAN_H_
#define METHODS_GAUSSIAN_H_

#include "MathFunctions.h"
#include "Methods.h"

namespace Terran {

class MethodsGaussian : public Methods {

public:

    explicit MethodsGaussian(const std::vector<Param> &params);

    // Partition the domain into disjoint intervals
    // std::vector<double> partition(double threshold) const;

    // Find the maxima for a gaussian mixture model
    std::vector<double> findMaxima() const;

    // Find the minima for a gaussian mixture model
    std::vector<double> findMinima() const;

private:

	struct Bracket {
		Bracket(float l, float m, float r) : left(l), middle(m), right(r) {};
		double left;
		double middle;
		double right;
	};

	vector<Bracket> minBrackets_;
	vector<Bracket> maxBrackets_;

};

}

#endif