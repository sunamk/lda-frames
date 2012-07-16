/*
 * Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
 * Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
 *
 */

#ifndef _ANTONIAK
#define _ANTONIAK

#include <vector>

using namespace std;

/*
 * Sample number of components m that a DP(alpha, G0) has after n samples.
 * This was first published by Antoniak (1974).
 */

class Antoniak_t {
public:
    Antoniak_t();
    vector<double> stirling(unsigned int n);
    unsigned int randBernoulli(double a);
    unsigned int randAntoniak(double alpha, unsigned int n);
    unsigned int randCRP(double alpha, unsigned int n);

private:
    vector<vector<double> > allss;
};

#endif
