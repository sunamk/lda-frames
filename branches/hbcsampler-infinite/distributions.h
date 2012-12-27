/*
 * Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
 * Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
 *
 */

#ifndef _DISTRIBUTIONS
#define _DISTRIBUTIONS

#include <vector>
#include <map>

using namespace std;

double randf(void);

void normalizeLogMult(map<unsigned short int, double> &v);

unsigned int sampleMultinomial(map<unsigned short int, double> &v);


#endif

