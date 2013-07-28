#include <MathFunctions.h>
#include <iostream>
#include <stdexcept>
#include <sstream>

using namespace std;
using namespace Terran;

void tunersumrg() {
    double period = 2*PI;
    int maxImages = 100;
    double tolerance = 1e-9;
    for(double sk = 0.01; sk < period; sk += 0.05) {
        for(double uk = -period; uk < period; uk += 0.05) {
            // number of extra images radiating from the center
            for(double xn = -period; xn < period; xn+= 0.05) {
                // use 100 images as absolute truth
                double truth = 0;
                for(int r = -maxImages; r <= maxImages; r++) {
                    truth += r*gaussian(uk+r*period, sk, xn);
                }
                for(int i=0; i < 15; i++) {
                    double test = 0;
                    for(int r = -i; r <= i; i++) {
                        test += r*gaussian(uk+r*period, sk, xn);
                    }
                    if(fabs(test-truth) < tolerance) {
                        //if(i > 3) {
                            //cout << sk << " " << uk << " " << xn << " " << i << " " << fabs(test-truth) << endl;
                            cout << sk << " " << i << " " << fabs(test-truth) << endl;
                      //}
                        break;
                    }
                }
            }
        }
    }

}

// test the periodic gaussian to make sure the number of images is tuned automatically
void testPeriodicGaussian() {
    double period = 2*PI;
    double tolerance = 1e-9;
    // vary std devation
    for(double sk = 0.01; sk < period; sk += 0.05) {
        // vary mean
        for(double uk = -period/2; uk < period/2; uk += 0.05) {
            // vary point sample
            for(double xn = -period/2; xn < period/2; xn+= 0.05) {
                // use 100 images as absolute truth
                double truth = periodicGaussian(uk,sk,xn,100,period);
                double test = periodicGaussian(uk,sk,xn,period); 
                if(fabs(test-truth) > tolerance) {
                    stringstream msg;
                    msg << "testPeriodicGaussian failed (sk, uk, xn, diff) " << sk << " " << uk << " " << xn << " " << fabs(test-truth) << endl;
                    throw(std::runtime_error(msg.str()));
                }
            }
        }
    }
}

// images n
void tunePeriodicGaussianImages() {

    double period = 2*PI;
    int maxImages = 100;
    
    double tolerance = 1e-9;

    for(double sk = 0.01; sk < period; sk += 0.05) {
        // we can potentially get rid of this uk by considering when the mean is roughly within 7 sigma of the boundary (where we then proceed by adding one extra image to the left or the right depending on the situation)
        for(double uk = -period; uk < period; uk += 0.05) {
            // number of extra images radiating from the center
            for(double xn = -period; xn < period; xn+= 0.05) {
                // use 100 images as absolute truth
                double truth = periodicGaussian(uk,sk,xn,maxImages,period);
                for(int i=0; i < 15; i++) {
                    double test = periodicGaussian(uk,sk,xn,i,period); 
                    if(fabs(test-truth) < tolerance) {
                        if(i > 3) {
                            //cout << sk << " " << uk << " " << xn << " " << i << " " << fabs(test-truth) << endl;
                            cout << sk << " " << i << " " << fabs(test-truth) << endl;
                        }
                        break;
                    }
                }
            }
        }
    }

}

int main() {
    try {
        testPeriodicGaussian();
        cout << "done";
    } catch(const exception &e) {
        cout << e.what() << endl;
    }
}
