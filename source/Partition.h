#include <vector>
#include "Param.h"

// An abstract class used to find minima and maxima of a mixture 
// model
class Partition {

public:

    Partition(const std::vector<Param> &params) :
        params_(params) {
        
    }; 

    // Partition the domain into disjoint intervals
    virtual std::vector<double> partition(double threshold) const = 0;

    // Find the maxima
    virtual std::vector<double> findMaxima() const = 0;

    // Find the minima
    virtual std::vector<double> findMinima() const = 0;

protected:

    const std::vector<Param> params_;

};
