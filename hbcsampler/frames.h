/*
 * Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
 * Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
 *
 */

#ifndef _FRAMES
#define _FRAMES

#include <set>
#include <vector>
using namespace std;

typedef vector<unsigned int> FrameKey_t;

class Frames_t {
public:
    Frames_t(unsigned int _S, unsigned int _R): S(_S), R(_R) {};

    FrameKey_t makeKey(unsigned int *f);
    
    FrameKey_t makeKey(unsigned int *f, unsigned int s, unsigned int r);

    bool inside(FrameKey_t frameKey);

    void insert(FrameKey_t frameKey);

    bool remove(FrameKey_t frameKey);

    static void printKey(FrameKey_t k);

    
private:
    unsigned int S;
    unsigned int R;
    set<FrameKey_t> frameSet;
};

#endif

