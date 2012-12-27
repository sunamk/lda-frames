/*
 * Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
 * Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
 *
 */

#include <iostream>
#include <cstdlib>
#include <math.h>
#include "distributions.h"

double MININF = -0/1;

inline double randf(void) {
    return (double) rand()/(double) RAND_MAX;
}

void normalizeLogMult(map<unsigned short int, double> &v) {
    if (v.size() == 0) {
        cerr << "Normalizing empty vector." << endl;
        exit(10);
    }
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

unsigned int sampleMultinomial(map<unsigned short int, double> &v) {
    if (v.size() == 0) {
        cerr << "Sampling empty vector." << endl;
        exit(10);
    }
    double s = randf();
    for (map<unsigned short int, double>::const_iterator it=v.begin(); it != v.end(); ++it) {
        s -= it->second;
        if (s<0) return it->first;
    }
    return v.begin()->first;
}

