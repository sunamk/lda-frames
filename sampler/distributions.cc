/*
 * Copyright (C) 2014 Jiri Materna <xmaterna@fi.muni.cz>
 * Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
 *
 */

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>
#include "distributions.h"

Distributions_t::Distributions_t(long int seed) {
    if (seed == 0){
        time_t t;
        time(&t);
        seed = (long int) t;
    }

    RANDOM_NUMBER = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(RANDOM_NUMBER,  seed);

    allss.push_back(vector<double>(1,1));
}


inline double Distributions_t::randf(void) {
    return  gsl_rng_uniform_pos(RANDOM_NUMBER);
}

unsigned int Distributions_t::randi(unsigned int min, unsigned int max) {
    return gsl_rng_uniform_int(RANDOM_NUMBER, max-min+1) + min;
}

void Distributions_t::normalizeLogMult(map<unsigned short int, double> &v) {
    //if (v.size() == 0) {
    //    cerr << "Normalizing empty vector." << endl;
    //    exit(10);
    // }
    double s = MININF;
    for (map<unsigned short int, double>::const_iterator it=v.begin(); it != v.end(); ++it) {
        if(s == MININF) {
            s = it->second;
        } else if (it->second == MININF) {
            continue;
        } else if (s - it->second > 16) {
            continue;
        } else if (s > it->second) {
            s +=  log(1 + exp(it->second - s));
        } else if (it->second - s > 16) {
            s = it->second;
        } else {
            s = it->second + log(1 + exp(s - it->second));
        }
    }
    for (map<unsigned short int, double>::const_iterator it=v.begin(); it != v.end(); ++it) {
        v[it->first] = exp(it->second - s);
    }
}

unsigned int Distributions_t::sampleMultinomial(map<unsigned short int, double> &v) {
    //if (v.size() == 0) {
    //    cerr << "Sampling empty vector." << endl;
    //    exit(10);
    //}
    double s = randf();
    for (map<unsigned short int, double>::const_iterator it=v.begin(); it != v.end(); ++it) {
        s -= it->second;
        if (s<0) return it->first;
    }
    return v.begin()->first;
}

unsigned int Distributions_t::sampleMultinomial(vector<double> &v, double sum, unsigned int first) {
    //if (v.size() == 0) {
    //    cerr << "Sampling empty vector." << endl;
    //    exit(10);
    //}
    double s = randf()*sum;
    unsigned int pos = first;
    for (vector<double>::const_iterator it=v.begin(); it != v.end(); ++it) {
        s -= *it;
        if (s<0) return pos;
        ++pos;
    }
    return first;
}

vector<unsigned int> Distributions_t::sampleMultinomial(map<vector<unsigned int>, unsigned int> &v, unsigned int sum) {
    double s = randf();
    for (map<vector<unsigned int>, unsigned int>::const_iterator it=v.begin();
                it != v.end(); ++it) {
        s-= (double) it->second / sum;
        if (s<0) return it->first;
    }
    return v.begin()->first;
}

unsigned int Distributions_t::sampleBernoulli(double p) {
    return gsl_ran_bernoulli(RANDOM_NUMBER, p);
}

double Distributions_t::sampleBeta(double a, double b) {
    return gsl_ran_beta(RANDOM_NUMBER, a, b);
}


void Distributions_t::sampleDirichlet(double *sample, double *alphas, int dim) {
    gsl_ran_dirichlet (RANDOM_NUMBER, dim, alphas, sample);
}

double Distributions_t::sampleGamma(double a, double b) {
  return gsl_ran_gamma_mt(RANDOM_NUMBER, a, b);
}


double Distributions_t::digamma(double x) {
    return gsl_sf_psi(x);
}

vector<double> Distributions_t::stirling(unsigned int n) {
    if (n > allss.size()) {
        for (unsigned int m = allss.size(); m < n; ++m) {
            unsigned int len = allss[m-1].size() + 1;
            allss.push_back(vector<double>(len, 0));
            double max = 0;
            for (unsigned int j = 0; j < len; ++j) {
                allss[m][j] += (j < len -1) ? allss[m - 1][j]*m: 0;
                allss[m][j] += (j == 0) ? 0 : allss[m - 1][j - 1];
                if (max < allss[m][j]) max = allss[m][j];
            }
            for (unsigned int j = 0; j < len; ++j) {
                allss[m][j] /= max;
            }
        }
    }
    return allss.at(n-1);
}

unsigned int Distributions_t::sampleAntoniak(double alpha, unsigned int n) {
    vector<double> p = stirling(n);
    double a = 1;
    double sum = 0;
    for (unsigned int m = 0; m < p.size(); ++m) {
        p[m] *= a;
        a *= alpha;
        sum += p[m];
    }
    return sampleMultinomial(p, sum, 1);        
}

double Distributions_t::MININF = -0/1;

Distributions_t::~Distributions_t() {
    //free random number generator
    gsl_rng_free(RANDOM_NUMBER);
}
