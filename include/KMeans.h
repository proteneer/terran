#include <iostream>
#include <vector>

namespace Terran {

class KMeans {

    public:
    
        /* KMeans does not manage explicit the dataset object, so the user must take care
           to ensure its existence. Size is NxD. For periodic dimensions, the periodic mean
           is computed. */
        KMeans(const std::vector<std::vector<double> > &dataset, const std::vector<bool> &isPeriodic);
        ~KMeans();
        
        /* Run the k-means algorithm given centers and return assignments*/
        vector<int> run(const vector<double> &centers);
       
    private:
    
        const std::vector<std::vector<double> > &dataset_;
        const std::vector<bool> isPeriodic_;

}

}