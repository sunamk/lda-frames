/*
 * Copyright (C) 2014 Jiri Materna <xmaterna@fi.muni.cz>
 * Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
 *
 */

#include "frames.h"
#include <iostream>


void Frames_t::setS(unsigned int slots) {
    S = slots;
}

FrameKey_t Frames_t::makeKey(vector<unsigned int> &f) {
    FrameKey_t k = f;
    return k;
}

FrameKey_t Frames_t::makeKey(vector<unsigned int> &f, unsigned int s, unsigned int r) {
    FrameKey_t k;
    for (unsigned int i = 0; i < S; ++i) {
        if (i == s-1) {
            k.push_back(r); 
        } else {
            k.push_back(f[i]);
        } 
    }
    return k;
}

bool Frames_t::inside(FrameKey_t frameKey) {
    set<FrameKey_t>::iterator  it;
    it = frameSet.find(frameKey);
    if (it == frameSet.end()) return false; else return true;
}

void Frames_t::insert(FrameKey_t frameKey) {
    frameSet.insert(frameKey);
}

bool Frames_t::remove(FrameKey_t frameKey) {
    return 1 == frameSet.erase(frameKey);
}

void Frames_t::printKey(FrameKey_t k) {
    for (FrameKey_t::const_iterator it = k.begin(); it != k.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;
}

void Frames_t::printAll(void) {
    for(set<FrameKey_t>::const_iterator it = frameSet.begin(); it!=frameSet.end(); ++it) {
        printKey(*it);
    }
}

void Frames_t::clear(void) {
    frameSet.clear();
}
