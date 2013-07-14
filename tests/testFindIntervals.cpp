#include <findIntervals.h>

using namespace std;

int main() {

    vector<Param> p;
    p.push_back(Param(0.6, 1.0, 0.5));
    p.push_back(Param(0.4, -2.3, 0.5));

    partitionPeriodicGaussian(p);

}


