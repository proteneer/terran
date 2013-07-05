#include "EMPeriodicGaussian.h"
#include <assert.h>

using namespace std;

EMPeriodicGaussian::EMPeriodicGaussian(const vector<double> &data, const vector<Param> &params, double period, int images) : 
    EM(data, params), 
    period_(period),
    images_(images) {

}

EMPeriodicGaussian::~EMPeriodicGaussian() {

}


#include <iostream>

void EMPeriodicGaussian::MStep() {

/*
    for(double u=-period_/2; u<period_/2; u+=0.01) {
        cout << u << " " << dldu(u, 0) << endl;
    }
*/


/*
    for(double s=0; s<6; s+= 0.01) {
        cout << s << " " << dlds(s, 0) << endl;
    }
  */  

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
    }
}

double EMPeriodicGaussian::qkn(int k, int n) const {
    double xn = data_[n];
    double pk = params_[k].p;
    double uk = params_[k].u;
    double sk = params_[k].s;
    return pk*periodicGaussian(uk,sk,xn,images_,period_);
}

double EMPeriodicGaussian::dldu(double uk, int k) const {
    double sk = params_[k].s;
    double s_sum = 0;
    for(int n=0; n<data_.size(); n++) {
        double xn = data_[n];
        double summand = 0;
        double r_sum_rg = 0;
        double r_sum_g = 0;
        for(int r=-images_; r<=images_; r++) {
            r_sum_g  += gaussian(uk+r*period_, sk, xn); 
            r_sum_rg += r*gaussian(uk+r*period_, sk, xn);
        }
        summand = r_sum_rg/r_sum_g;
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

        // how many images we need to use to approximate
        // the periodic gaussian
        const int R = images_;

        // domain length
        const int D = period_;

        for(int r=-R; r<=R; r++) {
            r_sum_g   += gaussian(uk+r*D, sk, xn);
            r_sum_rg  += r*gaussian(uk+r*D, sk, xn);
            r_sum_rrg += r*r*gaussian(uk+r*D, sk, xn);
        }

        p2 = -2*(xn-uk)*D*r_sum_rg / r_sum_g;
        p3 = D*D*r_sum_rrg / r_sum_g;
        s_sum += pikn_[n][k]*(p1+p2+p3);
    }
    return s_sum;
}
