/*
 * Copyright (C) 2013 Jiri Materna <xmaterna@fi.muni.cz>
 * Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
 *
 */

#include <iostream>
#include <math.h>
#include "sampler.h"

bool Sampler_t::initialize(bool recovery) {
    cout << "Allocating memory..." << endl;

    dist = new Distributions_t(seed);

    for (unsigned int u=1; u<=U; ++u) {
        frames.push_back(vector<unsigned int>(w[u-1].size(), 0));
        if (testPhase) test_frames.push_back(vector<unsigned int>(test_w[u-1].size(), 0));
    }

    for (unsigned int f=1; f<=F; ++f) {
        roles.push_back(vector<unsigned int>(S, 0));
        fc_f.push_back(0);
        fc_fsw.push_back(vector<vector<unsigned int> >(S,vector<unsigned int>(V, 0)));
        alpha.push_back(alpha0);
    }
    alpha.push_back(alpha0*F);

    for (unsigned int u=1; u<=U; ++u) {
        post_phi.push_back(vector<double>(F + 1, 0));
    }

    for (unsigned int r=1; r<=R; ++r) {
        post_theta.push_back(vector<double>(V + 1, 0));
        gamma[r] = gamma0;
    }
    gamma[0] = gamma0;

    if (cores != 0) {
        omp_set_num_threads(cores);
    }

    if (!recovery) {
        cout << "Initializing variables..." << endl;
        initialize_beta();
        initialize_roles();
        initialize_frames();
        initialize_infinite_vars();
        initialize_post_phi();
        initialize_post_theta();
        initialize_post_omega();
    }
    bestPerplexity = 0;
    initialized = true;
    return true;

}


void Sampler_t::initialize_roles(void) {
    if (!infinite_F) {
        double boundary = 0;
        unsigned int f = 1;
        for (map<vector<unsigned int>, unsigned int>::const_iterator it=framePatterns.begin();
                it != framePatterns.end(); ++it) {
            boundary += (((double) it->second) / positions)*(F - framePatterns.size()) + 1;

            while (f <= round(boundary)) {
                do {
                   for (unsigned int s=1; s<=S; ++s) {
                        if (it->first.at(s-1) == 0) {
                            roles[f-1][s-1] = 0;
                        } else {
                            roles[f-1][s-1] = dist->randi(1, R);
                        }
                    }
                } while (frameSet.inside(frameSet.makeKey(roles[f-1])));
                FrameKey_t fk = frameSet.makeKey(roles[f-1]);
                frameSet.insert(fk);
                f++;
            }
        }
        
    } else {
        //not necessarily
        if (F != 1) { 
            cout << "Internal error." << endl;
            exit(10);
        }
        //sample frame for the first realization
        for (unsigned int s=1; s<=S; ++s) {
            if (w[0][0][s-1] == 0) {
                roles[0][s-1] = 0;
             } else {
                roles[0][s-1] = dist->randi(1, R);
             }
        }
        FrameKey_t fk = frameSet.makeKey(roles[0]);
        frameSet.insert(fk);
    }
}

void Sampler_t::initialize_frames(void) {
    for (unsigned int u=1; u<=U; ++u) {
        for (unsigned int t=1; t<=w[u-1].size(); ++t) {
            if (!infinite_F) {
                //do it better;
                do {
                    frames[u-1][t-1] = dist->randi(1, F);
                } while (!checkPattern(w[u-1][t-1], roles[frames[u-1][t-1]-1]));

                fc_f[frames[u-1][t-1]-1]++;
                for (unsigned int s=1; s<=S; ++s) {
                    if (w[u-1][t-1][s-1]>0) { //non-empty slot
                        fc_fsw[frames[u-1][t-1]-1][s-1][w[u-1][t-1][s-1]-1]++;
                    }
                }
            } else {
                frames[u-1][t-1] = 0;
            }
        }
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
    //TODO all frames and roles may not be used
    for (unsigned int f=1; f<=F; ++f) {
        used_frames.insert(f);
        tau.insert(make_pair<unsigned int, bool>(f, 1.0/F+1));
    }
    tau.insert(make_pair<unsigned int, bool>(0, 1.0/(F+1))); //new component
    for (unsigned int r=1; r<=R; ++r) {
        used_roles.insert(r);
    }
}


