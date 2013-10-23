#ifndef _KMEANS_H_
#define _KMEANS_H_

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
        std::vector<int> run(const std::vector<std::vector<double> > &initialCenters) const;
       
    private:
    
        std::vector<double> computeMean(const std::vector<int> &points) const;
        double distance(const std::vector<double> &p1, const std::vector<double> &p2) const;
        std::vector<int> assignPointsToCenters(const std::vector<std::vector<double> > &centers) const;

        const std::vector<std::vector<double> > &dataset_;
        const std::vector<bool> isPeriodic_;

};

} // namespace Terran

#endif // ifndef _KMEANS_H_