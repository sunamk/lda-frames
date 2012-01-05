#include "frames.h"
#include <iostream>


FrameKey_t Frames_t::makeKey(unsigned int *f) {
    FrameKey_t k;
    for (unsigned int i = 0; i < S; ++i) {
        k.push_back(f[i]);
    }
    return k;
}

FrameKey_t Frames_t::makeKey(unsigned int *f, unsigned int s, unsigned int r) {
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
