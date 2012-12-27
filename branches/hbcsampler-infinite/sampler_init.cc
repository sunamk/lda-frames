/*
 * Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
 * Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
 *
 */

#include <iostream>
#include "sampler.h"
#include "stats.h"

bool Sampler_t::initialize(bool recovery) {
    srand((unsigned)time(0));   // initialize random number generator
    cout << "Allocating memory..." << endl;

    frames = (unsigned int**) malloc(sizeof(unsigned int*) * U);
    for (unsigned int u=1; u<=U; ++u) {
        frames[u-1] = (unsigned int*) malloc(sizeof(unsigned int) * w[u-1].size());
    }

    roles = (unsigned int**) malloc(sizeof(unsigned int*) * F);
    for (unsigned int f=1; f<=F; ++f) {
        roles[f-1] = (unsigned int*) malloc(sizeof(unsigned int) * S);
        fc_f.push_back(0);
        fc_fsw.push_back(vector<vector<unsigned int> >(S,vector<unsigned int>(V, 0)));
    }

    post_phi = (double**) malloc(sizeof(double*) * U);
    for (unsigned int u=1; u<=U; ++u) {
        post_phi[u-1] = (double*) malloc(sizeof(double) * (F + 1));
    }

    post_theta = (double**) malloc(sizeof(double*) * R);
    for (unsigned int r=1; r<=R; ++r) {
        post_theta[r-1] = (double*) malloc(sizeof(double) * (V + 1));
    }

    for (unsigned int r=0; r<R; ++r) {
        beta.push_back((double*) malloc(sizeof(double) * (V + 1)));
    }

    post_omega = (double*) malloc(sizeof(double) * (R +1));

    if (!recovery) {
        cout << "Initializing variables..." << endl;
        initialize_beta();
        initialize_frames();
        initialize_roles();
        initialize_infinite_vars();
        initialize_post_phi();
        initialize_post_theta();
        initialize_post_omega();
    }

    initialized = true;
    return true;

}


void Sampler_t::initialize_frames(void) {
    for (unsigned int u=1; u<=U; ++u) {
        for (unsigned int t=1; t<=w[u-1].size(); ++t) {
            if (!infinite_F) {
                frames[u-1][t-1] = sample_MultSym(1, F);
                fc_f[frames[u-1][t-1]-1]++;
                for (unsigned int s=1; s<=S; ++s) {
                    fc_fsw[frames[u-1][t-1]-1][s-1][w[u-1][t-1][s-1]-1]++;
                }
            } else {
                frames[u-1][t-1] = 0;
            }
        }
    }
}

void Sampler_t::initialize_roles(void) {
    for (unsigned int f=1; f<=F; ++f) {
        do {
            for (unsigned int s=1; s<=S; ++s) {
                roles[f-1][s-1] = sample_MultSym(1, R);
            }
        } while (frameSet->inside(frameSet->makeKey(roles[f-1])));
        FrameKey_t fk = frameSet->makeKey(roles[f-1]);
        frameSet->insert(fk);
    }
}


void Sampler_t::initialize_post_phi(void) {
    for (unsigned int u=1; u<=U; ++u) {
        for (unsigned int f=1; f<=F; ++f) {
            post_phi[u-1][f-1] = 0;
        }
        post_phi[u-1][F] = 0;
    }
    if (!infinite_F) {
        resample_post_phi();
    }
}

void Sampler_t::initialize_beta(void) {
    for (unsigned int r=0; r<R; ++r) {
        beta[r][0]=V*beta0;
        for (unsigned int v = 1; v<=V; ++v) {
            beta[r][v] = beta0;
        }
    }
}

void Sampler_t::initialize_post_theta(void) {
    for (unsigned int r=1; r<=R; ++r) {
        for (unsigned int v=1; v<=V; ++v) {
            post_theta[r-1][v-1] = 0;
        }
        post_theta[r-1][V] = 0;
    }
    if (!infinite_F) {
        resample_post_theta();
    }
}

void Sampler_t::initialize_post_omega(void) {
    for (unsigned int r=1; r<=R; ++r) {
        post_omega[r-1] = 0;
    }
    post_omega[R] = 0;
    resample_post_omega();
}

void Sampler_t::initialize_infinite_vars(void) {
    for (unsigned int f=1; f<=F; ++f) {
        used_frames.insert(f);
        tau.insert(make_pair<unsigned int, bool>(f, 1.0/F+1));
    }
    tau.insert(make_pair<unsigned int, bool>(0, 1.0/(F+1))); //new component
    for (unsigned int r=1; r<=R; ++r) {
        used_roles.insert(r);
    }
}


