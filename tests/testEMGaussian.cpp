// tests expectation maximization of gaussians for correctness

#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <EMGaussian.h>
#include <EMPeriodicGaussian.h>
#include <MathFunctions.h>
#include <MethodsGaussian.h>

#include "util.h"

using namespace std;
using namespace Terran;

void testUnimodalGaussian() {
    vector<Param> trueParams(1);
    trueParams[0].p = 1;
    trueParams[0].u = 12.3;
    trueParams[0].s = 4.7;
    vector<double> data;
	for(int i=0; i < 20000; i++) {
		data.push_back(gaussianMixtureSample(trueParams));
	}
    vector<Param> params;
    Param p;
    p.p = 1;
    p.u = 6.2;
    p.s = 8.1;
    params.push_back(p);
    EMGaussian em(data, params);
    em.run();
    vector<Param> optimizedParams = em.getParams();
    Util::matchParameters(trueParams, optimizedParams, 0.1);
}

void testBimodalGaussian() {
    int numSamples = 20000;
    vector<double> data;
    vector<Param> trueParams(2);
    trueParams[0].p = 0.4;
    trueParams[0].u = -3.4;
    trueParams[0].s = 1.2;

    trueParams[1].p = 0.6;
    trueParams[1].u = 7.4;
    trueParams[1].s = 6.2;

	for(int i=0; i < 20000; i++) {
		data.push_back(gaussianMixtureSample(trueParams));
	}

    vector<Param> params;
    {
        Param p;
        p.p = 0.5;
        p.u = 6.2;
        p.s = 8.1;
        params.push_back(p);
    }
    {
        Param p;
        p.p = 0.5;
        p.u = 0;
        p.s = 2.1;
        params.push_back(p);
    }
    EMGaussian em(data, params);
    em.run();
    vector<Param> optimizedParams = em.getParams();

    for(int i=0; i < optimizedParams.size(); i++) {
        cout << optimizedParams[i].p << " " << optimizedParams[i].u << " " << optimizedParams[i].s << endl;
    
    }

    Util::matchParameters(trueParams, optimizedParams, 0.2);
}

void testSimpleRunBimodal() {
    int numSamples = 20000;
    vector<double> data;
    vector<Param> trueParams(2);
    trueParams[0].p = 0.4;
    trueParams[0].u = -3.4;
    trueParams[0].s = 1.2;

    trueParams[1].p = 0.6;
    trueParams[1].u = 7.4;
    trueParams[1].s = 6.2;

	for(int i=0; i < 20000; i++) {
		data.push_back(gaussianMixtureSample(trueParams));
	}

    vector<double> subset(data);
    random_shuffle(subset.begin(), subset.end());
    subset.resize(3500);
    EMGaussian em(subset);

    bool rc = em.simpleRun(100);

	cout << rc << endl;

    vector<Param> optimizedParams = em.getParams();
    MethodsGaussian method(optimizedParams);
    vector<double> maxima = method.findMaxima();
    vector<double> minima = method.findMinima();

    vector<double> trueMaxima;
    trueMaxima.push_back(-3.25);
    trueMaxima.push_back( 7.4);

    vector<double> errorsMaxima;
    errorsMaxima.push_back(0.5);
    errorsMaxima.push_back(1.5);

	ofstream data2("gaussianH");
	for(int i=0; i < subset.size(); i++) {
		data2 << subset[i] << endl;
	}

	Util::plotGaussian(optimizedParams, -10, 15, "gaussianO");
	Util::plotGaussian(trueParams, -10, 15, "gaussianT");

		
    rc = em.simpleRun(100);
	vector<Param> optimizedParams2 = em.getParams();
	Util::plotGaussian(optimizedParams2, -10, 15, "gaussianO2");
	

    Util::matchPoints(maxima, trueMaxima, errorsMaxima);
    
	cout << "foo" << endl;
	vector<double> trueMinima;
    trueMinima.push_back(0.5);
    Util::matchPoints(minima, trueMinima, 0.3);

    if(rc != 1) {
        throw std::runtime_error("EM - Failed to converge");
    }    
}

void testUniSpecial() {
    int numSamples = 20000;
    vector<double> data;
    vector<Param> trueParams(1);
    trueParams[0].p = 1;
    trueParams[0].u = -3.4;
    trueParams[0].s = 1.2;
    for(int i=0; i < 20000; i++) {
        data.push_back(gaussianMixtureSample(trueParams));
    }

    vector<double> subset(data);
    random_shuffle(subset.begin(), subset.end());
    subset.resize(3500);
    EMGaussian em(subset);

    ofstream data2("gaussianH");
    for(int i=0; i < subset.size(); i++) {
        data2 << subset[i] << endl;
    }

    Util::plotGaussian(trueParams, -10, 15, "gaussianT");

    cout << trueParams[0].p << endl;

    //bool rc = em.simpleRun(2);

    bool rc = em.simpleRun(10);
    vector<Param> optimizedParams = em.getParams();
    Util::plotGaussian(optimizedParams, -10, 15, "gaussianO");
    cout << optimizedParams << endl;

}

int main() {
    try	{
        testUniSpecial();
		//srand(1);
        //cout << "testSimpleRunBimodal()" << endl;
		//testSimpleRunBimodal();
        /*
		cout << "testUnimodalGaussian()" << endl;
		srand(1);
        testUnimodalGaussian();
        cout << "testBimodalGaussian()" << endl;
		srand(1);
        testBimodalGaussian();
        */
    } catch( const std::exception &e ) {
        cout << e.what() << endl;
    }
    cout << "done" << endl;
}
