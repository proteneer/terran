#include "EMPeriodicGaussian.h"
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

namespace Terran {

EMPeriodicGaussian::EMPeriodicGaussian(const vector<double> &data, const vector<Param> &params, double period) : 
    EM(data, params), 
    period_(period) {
    for(int i=0; i<params.size(); i++) {
        if(params_[i].s > period_)
            throw(std::runtime_error("Cannot have s > period in parameters"));
        if(params_[i].u < -period_/2)
            throw(std::runtime_error("Cannot have u < -period/2"));
        if(params_[i].u > period_/2)
            throw(std::runtime_error("Cannot have u > period/2"));
    }
}

EMPeriodicGaussian::~EMPeriodicGaussian() {

}

void EMPeriodicGaussian::MStep() {
#ifndef NDEBUG
    testIntegrity();
#endif

    // Compute new mean
    for(int k=0; k<params_.size(); k++) {
        // find two points ak and bk such that
        // ak < bk, and dldu(ak) > dldu(bk)
        // we just divide the domain into 8 equal parts
        // and guess that its bound to cross at one point
        int increment = 16;

        // look for two points that bracket the point
        vector<double> signs(increment);
        vector<double> xvals(increment);
        for(int i=0; i<increment; i++) {
            double x = -period_/2+i*period_/increment;
            signs[i] = dldu(x,k);
            xvals[i] = x;
        }

        bool found = false;
        double ak,bk;
        for(int i=0; i<signs.size(); i++) {
            if(signs[i]>0 && signs[(i+1)%increment]<0) {
                found = true;
                ak = xvals[i];
                bk = xvals[(i+1)%increment];
            }
        }

        // need to search account of periodic domains as well

        if(!found) {
            throw(std::runtime_error("Fatal: Bracket not found!"));
        }

        double mk;
        double my;

        int iteration = 0;
        // bisection with periodic boundaries
        do {


            if(iteration > 1e3) {
            
                ofstream derp("Asdfasdf.txt");
                for(double a = -PI; a < PI; a+= 0.01) {
                    derp << a << " " << dldu(a,k) << endl;
                }
            
                throw(std::runtime_error("Error: MStep() (mean) maximum number of iterations reached."));
            }
             else
                 iteration++;
                  
            //assert(dldu(ak, k) > 0 || isnan(dldu(ak, k)));
            /*
            cout << ak << " " << bk << endl;
            cout << dldu(ak,k) << " " << dldu(bk,k) << endl;
            */
            assert(dldu(bk, k) < 0);
            
            if(bk < ak)
                mk = normalize((ak+period_+bk)/2.0);
            else
                mk = normalize((ak+bk)/2.0);
            my = dldu(mk, k);
            if(my > 0)
                ak = mk;
            else
                bk = mk;

            std::cout.precision(10);
            //cout <<  iteration << " " << ak << " " << bk << " " << my << endl;
        } while(fabs(my) > 1e-5);
        params_[k].u = mk;
    }

    // Compute new standard deviation and probability
    for(int k=0; k<params_.size(); k++) {
        // bisection root-finding algorithm
        double ak = 0.01;
        double bk = period_;
        double mk;
        double my;
        do {
            assert(dlds(ak, k) > 0 || isnan(dlds(ak, k)));
            assert(dlds(bk, k) < 0);
            mk = (ak+bk)/2.0;
            my = dlds(mk, k);
            if(my > 0)
                ak = mk;
            else
                bk = mk;
        } while(fabs(my) > 1e-5);
        params_[k].s = mk;

        // Compute new probability
        double denominatorSum = 0;
        for(int n=0; n<data_.size(); n++) {
            denominatorSum += pikn_[n][k];
        }
        params_[k].p = denominatorSum / data_.size();
    }

}

double EMPeriodicGaussian::qkn(int k, int n) const {
    double xn = data_[n];
    double pk = params_[k].p;
    double uk = params_[k].u;
    double sk = params_[k].s;
    return pk*periodicGaussian(uk,sk,xn,period_);
}

double EMPeriodicGaussian::dldu(double uk, int k) const {
    double sk = params_[k].s;
    double s_sum = 0;
    for(int n=0; n<data_.size(); n++) {
        double xn = data_[n];
        double summand = 0;
        double r_sum_rg = 0;
        double r_sum_g = 0;

        double left = uk-7*sk;
        double right = uk+7*sk;
        int lower = floor(left/period_+0.5);
        int upper = floor(right/period_+0.5);

        lower = min(-1, lower);
        upper = max( 1, upper);

        for(int r=lower; r <= upper; r++) {
            r_sum_g  += gaussian(uk+r*period_, sk, xn); 
            r_sum_rg += r*gaussian(uk+r*period_, sk, xn);
        }
        /*
        int images_ = 50;
        for(int r=-images_; r<=images_; r++) {
            r_sum_g  += gaussian(uk+r*period_, sk, xn); 
            r_sum_rg += r*gaussian(uk+r*period_, sk, xn);
        }*/
        summand = r_sum_rg/r_sum_g;
        if(isnan(summand))
            summand = 0;
        s_sum += pikn_[n][k]*(xn - uk - period_*summand);
    }
    return s_sum;
}

double EMPeriodicGaussian::dlds(double sk, int k) const {
    double s_sum = 0;
    const double uk = params_[k].u;
    for(int n=0; n<data_.size(); n++) {
        const double xn = data_[n];
        double p1 = 0;
        double p2 = 0;
        double p3 = 0;
        p1 = -sk*sk + (xn-uk)*(xn-uk);

        double r_sum_g = 0;
        double r_sum_rg = 0;
        double r_sum_rrg = 0;

        const int D = period_;

        double left = uk-7*sk;
        double right = uk+7*sk;
        int lower = floor(left/period_+0.5);
        int upper = floor(right/period_+0.5);
        lower = min(-1, lower);
        upper = max( 1, upper);
        for(int r=lower; r<=upper; r++) {
            double gauss = gaussian(uk+r*D, sk, xn);
            r_sum_g   += gauss;
            r_sum_rg  += r*gauss;
            r_sum_rrg += r*r*gauss;
        }
        p2 = -2*(xn-uk)*D*r_sum_rg / r_sum_g;
        p3 = D*D*r_sum_rrg / r_sum_g;
        if(isnan(p2))
            p2 = 0;
        if(isnan(p3))
            p3 = 0;
        s_sum += pikn_[n][k]*(p1+p2+p3);
    }
    return s_sum;
}

}
