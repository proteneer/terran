#include <vector>
#include "Param.h"

// An abstract class used to find minima and maxima of a mixture model
// Refactor into Methods class

namespace Terran {

class Methods {

public:

    Methods(const std::vector<Param> &params) :
        params_(params) {
        
    }; 

    // Find the maxima
    virtual std::vector<double> findMaxima() const = 0;

    // Find the minima
    virtual std::vector<double> findMinima() const = 0;

protected:

    const std::vector<Param> params_;

};

}
