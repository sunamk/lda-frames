/*
 * Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
 * Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
 *
 */

#ifndef _DISTRIBUTIONS
#define _DISTRIBUTIONS

#include <vector>
#include <map>
#include <gsl/gsl_rng.h>


using namespace std;

class Distributions_t {

public:

    Distributions_t(long int seed);

    double randf(void);
    
    unsigned int randi(unsigned int min, unsigned int max);

    void normalizeLogMult(map<unsigned short int, double> &v);

    unsigned int sampleMultinomial(map<unsigned short int, double> &v);

    unsigned int sampleMultinomial(vector<double> &v, double sum, unsigned int first);
    
    vector<unsigned int> sampleMultinomial(map<vector<unsigned int>, unsigned int> &v, unsigned int sum);

    unsigned int sampleBernoulli(double p);

    double sampleBeta(double a, double b);

    void sampleDirichlet(double *sample, double *alphas, int dim);

    double sampleGamma(double a, double b);

    double digamma(double x);

    vector<double> stirling(unsigned int n);

    //Sample number of components m that a DP(alpha, G0) has after n samples.
    //This was first published by Antoniak (1974).
    unsigned int sampleAntoniak(double alpha, unsigned int n);

    ~Distributions_t();
    
    static double MININF;

private:
    gsl_rng * RANDOM_NUMBER; //random number generator
    vector<vector<double> > allss; //stored values of the stirling number
};

#endif

