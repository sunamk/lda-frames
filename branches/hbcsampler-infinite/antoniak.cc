/*
 * Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
 * Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
 *
 */
#include <algorithm>
#include <iostream>
#include "stats.h"
#include "antoniak.h"

Antoniak_t::Antoniak_t() {
    allss.push_back(vector<double>(1,1));
}

vector<double> Antoniak_t::stirling(unsigned int n) {
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

unsigned int Antoniak_t::randBernoulli(double a) {
    double r = ranf();
    if (r < a) {
        return 1;
    } else {
        return 0;
    }
}

unsigned int Antoniak_t::randAntoniak(double alpha, unsigned int n) {
    vector<double> p = stirling(n);
    double a = 1;
    double sum = 0;
    for (unsigned int m = 0; m < p.size(); ++m) {
        p[m] *= a;
        a *= alpha;
        sum += p[m];
    }
    p.push_back(sum);

    double* pp = &p[0];
    
    return sample_Mult(pp, 1, p.size() - 1);
   
}

unsigned int Antoniak_t::randCRP(double alpha, unsigned int n) {
    double ainv = 1/alpha;
    double nt = 0;
    unsigned int R = n;
    double* p = (double*) malloc(sizeof(double)* (n+1));
    for (unsigned int i=0; i<n; ++i) p[i] = 0;

    for (unsigned int r=0; r<R; ++r) {
        for (unsigned int m=0; m<n; ++m) {
            for (unsigned int t=0; t<n; t++, nt += ainv) {
                p[m] += randBernoulli(1.0/(nt+1));
            }
        }
    }
    double sum = 0;    
    for (unsigned int i=0; i<n; ++i) {
        p[i] /= R;
        sum += p[i];
    }
    p[n] = sum;
    return sample_Mult(p, 1, n);

    free(p);
}



