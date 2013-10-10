#ifndef METHODS_H_
#define METHODS_H_

#include <vector>
#include "Param.h"


// An abstract class used to find minima and maxima of a mixture model
// Refactor into Methods class

namespace Terran {

class TERRAN_EXPORT Methods {

public:

    Methods(const std::vector<Param> &params) :
        params_(params) {
        
    }; 

    // Find the maxima
    virtual std::vector<double> findMaxima() const = 0;

    // Find the minima
    virtual std::vector<double> findMinima() const = 0;

protected:

	
	struct Bracket {
		Bracket(float l, float m, float r) : left(l), middle(m), right(r) {};
		double left;
		double middle;
		double right;
	};

    const std::vector<Param> params_;

};

}

#endif