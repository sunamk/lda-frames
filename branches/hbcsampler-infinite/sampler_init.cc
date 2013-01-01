/*
 * Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
 * Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
 *
 */

#include <iostream>
#include "sampler.h"

bool Sampler_t::initialize(bool recovery) {
    cout << "Allocating memory..." << endl;

    dist = new Distributions_t(seed);

    for (unsigned int u=1; u<=U; ++u) {
        frames.push_back(vector<unsigned int>(w[u-1].size(), 0));
    }

    for (unsigned int f=1; f<=F; ++f) {
        roles.push_back(vector<unsigned int>(S, 0));
        fc_f.push_back(0);
        fc_fsw.push_back(vector<vector<unsigned int> >(S,vector<unsigned int>(V, 0)));
        alpha.push_back(alpha0);
    }

    for (unsigned int u=1; u<=U; ++u) {
        post_phi.push_back(vector<double>(F + 1, 0));
    }

    for (unsigned int r=1; r<=R; ++r) {
        post_theta.push_back(vector<double>(V + 1, 0));
        gamma[r] = gamma0;
    }
    gamma[0] = gamma0;

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
    bestPerplexity = 0;
    initialized = true;
    return true;

}


void Sampler_t::initialize_frames(void) {
    for (unsigned int u=1; u<=U; ++u) {
        for (unsigned int t=1; t<=w[u-1].size(); ++t) {
            if (!infinite_F) {
                frames[u-1][t-1] = dist->randi(1, F);
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
                roles[f-1][s-1] = dist->randi(1, R);
            }
        } while (frameSet.inside(frameSet.makeKey(roles[f-1])));
        FrameKey_t fk = frameSet.makeKey(roles[f-1]);
        frameSet.insert(fk);
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
    for (unsigned int v = 1; v<=V; ++v) {
        beta.push_back(beta0);
    }
    beta.push_back(V*beta0);
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
    for (unsigned int r=1; r<=R+1; ++r) {
        post_omega.push_back(0);
    }
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


