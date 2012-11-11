/*
 * Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
 * Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <limits>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/errno.h>
#include <dirent.h>
#include <boost/tokenizer.hpp>

#include "sampler.h"
#include "stats.h"

void Sampler_t::resample_post_phi(void) {
    for (unsigned int u = 1; u <= U; ++u) {
        for (unsigned int f = 1; f<=F; ++f) {
            post_phi[u-1][f-1] = 0;
        }
        post_phi[u-1][F] = 0;
        for (unsigned int t = 1; t <= w[u-1].size(); ++t) {
            post_phi[u-1][F]++;
            post_phi[u-1][frames[u-1][t-1]-1]++;
        }
    }
}

void Sampler_t::resample_post_theta(void) {
    for (unsigned int r = 1; r <= R; ++r) {
        for (unsigned int v = 1; v <= V; ++v) {
            post_theta[r-1][v-1] = 0;
        }
        post_theta[r-1][V] = 0;
        for (unsigned int u = 1; u <= U; ++u) {
            for (unsigned int t = 1; t <= w[u - 1].size(); ++t) {
                for (unsigned int s = 1; s <= S; ++s) {
                    post_theta[r-1][V] += (r == roles[frames[u-1][t-1]-1][s-1]) ? 1 : 0;
                    post_theta[r-1][w[u-1][t-1][s-1]-1] += 
                        (r == roles[frames[u-1][t-1]-1][s-1]) ? 1 : 0;
                }
            }
        }
    }
}

void Sampler_t::resample_post_omega(void) {
    for (unsigned int f = 1; f <= F; ++f) {
        for (unsigned int s = 1; s <= S; ++s) {
             post_omega[roles[f-1][s-1]-1] = 0;
        }
    }
    post_omega[R] = 0;
    for (unsigned int f = 1; f <= F; ++f) {
        for (unsigned int s = 1; s <= S; ++s) {
            post_omega[roles[f-1][s-1]-1]++;
            post_omega[R]++;
        }
    }
}

void Sampler_t::resample_tau(void) {
    vector<unsigned int> idmap;
    unsigned int id = 0;
    double* fc = (double*) malloc(sizeof(double) * (used_frames.size() + 1));
    double* dirs = (double*) malloc(sizeof(double) * (used_frames.size() + 2));
    tables = 0;
    for(set<unsigned int>::const_iterator it = used_frames.begin(); it != used_frames.end(); ++it) {
        idmap.push_back(*it);
        fc[id] = 0;
        for (unsigned int u=1; u<=U; ++u) {
            if (post_phi[u-1][*it-1] > 1) {
                fc[id] += antoniak.randAntoniak(alpha * tau[*it], post_phi[u-1][*it-1]);
                //fc[id] += antoniak.randCRP(alpha * tau[*it], post_phi[u-1][*it-1]);
            } else {
                fc[id] += post_phi[u-1][*it-1];
            }
        }
        tables += fc[id];
        id++;
    }    
    fc[used_frames.size()] = delta;
    sample_Dir(dirs, fc, used_frames.size() + 1);
    for (unsigned int i=0; i<used_frames.size(); ++i) {
        tau[idmap[i]] = dirs[i];
    }
    tau[0] = dirs[used_frames.size()]; 
    free(fc);
    free(dirs);
}


void Sampler_t::resample_frames(void) {
    for (int u = 1; u <= (signed int) U; ++u) {
        double* post_frames = (double*) malloc(sizeof(double) * (F + 1));
        for (unsigned int t = 1; t <= w[u-1].size(); ++t) {

            //remove old values
            for (unsigned int s = 1; s <= S; ++s) {
                post_theta[roles[frames[u-1][t-1]-1][s-1]-1][V]--;
                post_theta[roles[frames[u-1][t-1]-1][s-1]-1][w[u-1][t-1][s-1]-1]--;
                fc_fsw[frames[u-1][t-1]-1][s-1][w[u-1][t-1][s-1]-1]--;
            }
            post_phi[u-1][F]--;
            post_phi[u-1][frames[u-1][t-1]-1]--;
            fc_f[frames[u-1][t-1]-1]--;

            //compute frame distribution
            post_frames[F] = 0;
            for (unsigned int f = 1; f <= F; ++f) {
                double prod = 0;
                post_frames[f-1] = 0;
                for (unsigned int s = 1; s <= S; ++s) {
                    prod += ldf_Mult_smooth(1, beta[roles[f-1][s-1]-1], w[u-1][t-1][s-1],
                            post_theta[roles[f-1][s-1]-1], 1, V);
                }
                post_frames[f-1] = prod + ldf_Mult_smooth(0, alpha, f, post_phi[u-1], 1, F);
            }
            normalizeLog(post_frames, 1, F);

            //sample frame
            frames[u-1][t-1] = sample_Mult(post_frames, 1, F);

            //update new values
            for (unsigned int s=1; s<=S; ++s) {
                post_theta[roles[frames[u-1][t-1]-1][s-1]-1][V]++;
                post_theta[roles[frames[u-1][t-1]-1][s-1]-1][w[u-1][t-1][s-1]-1]++;
                fc_fsw[frames[u-1][t-1]-1][s-1][w[u-1][t-1][s-1]-1]++;
            }
            post_phi[u-1][F]++;
            post_phi[u-1][frames[u-1][t-1]-1]++;
            fc_f[frames[u-1][t-1]-1]++;

        }
        free(post_frames);
    }
}

void Sampler_t::resample_frames_inf(void) {
    for (int u = 1; u <= (signed int) U; ++u) {
        double* post_frames = (double*) malloc(sizeof(double) * (F + 2));
        for (unsigned int t = 1; t <= w[u-1].size(); ++t) {
            
            bool tau_needs_resampling = false;

            //remove old values
            if (frames[u-1][t-1] > 0) { //already sampled?
                for (unsigned int s = 1; s <= S; ++s) {
                    post_theta[roles[frames[u-1][t-1]-1][s-1]-1][V]--;
                    post_theta[roles[frames[u-1][t-1]-1][s-1]-1][w[u-1][t-1][s-1]-1]--;
                    fc_fsw[frames[u-1][t-1]-1][s-1][w[u-1][t-1][s-1]-1]--;

                }
                post_phi[u-1][F]--;
                post_phi[u-1][frames[u-1][t-1]-1]--;
                fc_f[frames[u-1][t-1]-1]--;
            }

            //compute frame distribution
            for (unsigned int f = 1; f <= F; ++f) {
                double prod = 0;
                post_frames[f-1] = -1 * numeric_limits<double>::max(); //zero probability
                set<unsigned int>::iterator it = used_frames.find(f);

                //probability of used frames
                if (it != used_frames.end()) {

                    for (unsigned int s = 1; s <= S; ++s) {
                        prod += ldf_Mult_smooth(1, beta[roles[f-1][s-1]-1], w[u-1][t-1][s-1],
                                post_theta[roles[f-1][s-1]-1], 1, V);
                    }
                    post_frames[f-1] = prod + ldf_Mult_smooth(0, alpha*tau[f], f, post_phi[u-1], 1, F, 
                                       used_frames.size());
                }
            }
            //sample new frame
            vector<unsigned int> frame(S, 0);
            if (!sample_new_frame(frame, w[u-1][t-1])) {
                post_frames[F] = -1 * numeric_limits<double>::max();
            } else {
                double prod = log(alpha*tau[0]);
                for (unsigned int s = 1; s <= S; ++s) {
                    //prod += ldf_Mult_smooth(1, beta, w[u-1][t-1][s-1], opravit betu
                    //    post_theta[frame[s-1]-1], 1, V);
                    prod -= log(V);
                }
                post_frames[F] = prod;
            }

            normalizeLog(post_frames, 1, F + 1);

            //sample frame id
            unsigned int newFrame = sample_Mult(post_frames, 1, F + 1);
            
            //free unused frames and roles
            if (frames[u-1][t-1] != newFrame && frames[u-1][t-1] != 0 && fc_f[frames[u-1][t-1]-1] == 0) {
                frameSet->remove(frameSet->makeKey(roles[frames[u-1][t-1]-1]));
                for (unsigned int s=1; s<=S; ++s) {
                    post_omega[roles[frames[u-1][t-1]-1][s-1]-1]--;
                    post_omega[R]--;
                    if (post_omega[roles[frames[u-1][t-1]-1][s-1]-1] == 0 && infinite_R) {
                        unused_roles.insert(roles[frames[u-1][t-1]-1][s-1]);
                        used_roles.erase(roles[frames[u-1][t-1]-1][s-1]);
                    }
                }
                unused_frames.insert(frames[u-1][t-1]);
                used_frames.erase(frames[u-1][t-1]);
                tau_needs_resampling = true;
                tau.erase(frames[u-1][t-1]);
            }

            //create a new frame if required
            if (newFrame == F + 1) {
                if (frames[u-1][t-1]==0 && u==1 && t==1) {
                    used_frames.erase(1);
                    unused_frames.insert(1);
                    post_omega[0] -= 2;
                }
                frames[u-1][t-1] = createNewFrame(frame);
                tau[frames[u-1][t-1]] = tau[0];
                tau_needs_resampling = true;
                post_frames = (double*) realloc(post_frames, sizeof(double) * (F + 2));
                //resample_roles_inf();
            } else {
                frames[u-1][t-1] = newFrame;
            
                //remove possibly new roles of the unused new frame
                for (unsigned int s=1; s<=S; ++s) {
                    if (frame[s-1] != 0 && post_omega[frame[s-1]-1] == 0) {
                        used_roles.erase(frame[s-1]);
                        unused_roles.insert(frame[s-1]);
                    }
                }
            }
            

            //update new values
            for (unsigned int s=1; s<=S; ++s) {
                post_theta[roles[frames[u-1][t-1]-1][s-1]-1][V]++;
                post_theta[roles[frames[u-1][t-1]-1][s-1]-1][w[u-1][t-1][s-1]-1]++;
                fc_fsw[frames[u-1][t-1]-1][s-1][w[u-1][t-1][s-1]-1]++;
            }
            post_phi[u-1][F]++;
            post_phi[u-1][frames[u-1][t-1]-1]++;
            fc_f[frames[u-1][t-1]-1]++;
                    
            if(tau_needs_resampling) resample_tau();

        }
        free(post_frames);
    }
}

void Sampler_t::resample_roles(void) {

      for (set<unsigned int>::const_iterator fit = used_frames.begin(); fit!=used_frames.end(); ++fit) {
        double* post_roles = (double*) malloc(sizeof(double) * (R + 1));
        
        for (unsigned int s = 1; s <= S; ++s) {
            post_theta[roles[*fit-1][s-1]-1][V] -= fc_f[*fit-1];
            for(unsigned int v=1; v<=V; ++v) {
                post_theta[roles[*fit-1][s-1]-1][v-1] -= fc_fsw[*fit-1][s-1][v-1];
            }

            post_omega[roles[*fit-1][s-1]-1]--;
            post_omega[R]--;
            post_roles[R] = 0;
           
            FrameKey_t oldFrame; 
            oldFrame = frameSet->makeKey(roles[*fit-1]);

            for (unsigned int r = 1; r <= R; ++r) {
                double prod = 0;
                post_roles[r-1] = -1 * numeric_limits<double>::max(); //zero probability
                FrameKey_t newFrame = frameSet->makeKey(roles[*fit-1], s, r);
                bool inside = frameSet->inside(newFrame);
                if (newFrame == oldFrame || !inside) {
                    for (unsigned int v = 1; v<=V; ++v) {
                        prod += fc_fsw[*fit-1][s-1][v-1]*ldf_Mult_smooth(1, beta[r-1], v, 
                            post_theta[r-1], 1, V);
                    }
                    post_roles[r-1] = prod + ldf_Mult_smooth(0, gamma, r, post_omega, 
                        1, R);
                }
            }
            
            normalizeLog(post_roles, 1, R);

            roles[*fit-1][s-1] = sample_Mult(post_roles, 1, R);
        
            //Not necessarily
            FrameKey_t newFrame = frameSet->makeKey(roles[*fit-1]);
            if (frameSet->inside(newFrame) && oldFrame != newFrame){
                cerr << endl << "Error:  existing frame sampled [" << roles[*fit-1][s-1] << "]." << endl;
                for (unsigned int ii=0; ii<R; ++ii) {
                    cerr << ii+1 << ":" << post_roles[ii] << " ";
                }
                cout << endl;
                exit(10);
            }

            frameSet->remove(oldFrame);
            frameSet->insert(frameSet->makeKey(roles[*fit-1]));

            post_omega[roles[*fit-1][s-1]-1]++;
            post_omega[R]++;
            post_theta[roles[*fit-1][s-1]-1][V] += fc_f[*fit-1];
            for(unsigned int v=1; v<=V; ++v) {
                post_theta[roles[*fit-1][s-1]-1][v-1] += fc_fsw[*fit-1][s-1][v-1];
            }
        }
        free(post_roles);
    }
}

void Sampler_t::resample_roles_inf(void) {

    for (set<unsigned int>::const_iterator fit = used_frames.begin(); fit!=used_frames.end(); ++fit) {
        double* post_roles = (double*) malloc(sizeof(double) * (R + 2));
        
        for (unsigned int s = 1; s <= S; ++s) {
            post_theta[roles[*fit-1][s-1]-1][V] -= fc_f[*fit-1];
            for(unsigned int v=1; v<=V; ++v) {
                post_theta[roles[*fit-1][s-1]-1][v-1] -= fc_fsw[*fit-1][s-1][v-1];
            }

            post_omega[roles[*fit-1][s-1]-1]--;
            post_omega[R]--;

            post_roles[R+1] = 0;
           
            FrameKey_t oldFrame; 
            oldFrame = frameSet->makeKey(roles[*fit-1]);
            
            for (unsigned int r = 1; r <= R; ++r) {
                double prod = 0;
                post_roles[r-1] = -1 * numeric_limits<double>::max(); //zero probability
                FrameKey_t newFrame;
                bool inside;
                newFrame = frameSet->makeKey(roles[*fit-1], s, r);
                inside = frameSet->inside(newFrame);
                set<unsigned int>::iterator it = used_roles.find(r);

                //probability of used roles
                if ((newFrame == oldFrame || !inside) && it != used_roles.end()) {
                    for (unsigned int v = 1; v<=V; ++v) {
                        prod += fc_fsw[*fit-1][s-1][v-1]*ldf_Mult_smooth(1, beta[r-1], v, 
                            post_theta[r-1], 1, V);
                    }
                    //post_roles[r-1] = prod + ldf_Mult_smooth(0, gamma, r, post_omega, 
                    //    1, R, used_roles.size());
                    post_roles[r-1] = prod + log(post_omega[r-1] + gamma);
                }
            }

            //probability of a new role
            double prod = 0;
            for (unsigned int v = 1; v<=V; ++v) {
                prod -= fc_fsw[*fit-1][s-1][v-1]*log(V);
            }
            post_roles[R] = log(gamma) + prod;
            
            normalizeLog(post_roles, 1, R + 1);

            unsigned int newRole = sample_Mult(post_roles, 1, R + 1);

            //not necessarily
            set<unsigned int>::iterator it = used_roles.find(newRole);
            if (it == used_roles.end() && newRole != (R + 1)) {
                cerr << "Sampled disallowed role." << endl;
                exit(10);
            }
            
            //free unused role numbers
            if (roles[*fit-1][s-1] != newRole && post_omega[roles[*fit-1][s-1]-1] == 0) {
                unused_roles.insert(roles[*fit-1][s-1]);
                used_roles.erase(roles[*fit-1][s-1]);
            }
            //create a new role if required
            if (newRole == R + 1) {
                roles[*fit-1][s-1] = createNewRole();
                post_roles = (double*) realloc(post_roles, sizeof(double) * (R + 2));
            } else {
                roles[*fit-1][s-1] = newRole;
            }

            frameSet->remove(oldFrame);
            frameSet->insert(frameSet->makeKey(roles[*fit-1]));

            post_omega[roles[*fit-1][s-1]-1]++;
            post_omega[R]++;
            post_theta[roles[*fit-1][s-1]-1][V] += fc_f[*fit-1];
            for(unsigned int v=1; v<=V; ++v) {
                post_theta[roles[*fit-1][s-1]-1][v-1] += fc_fsw[*fit-1][s-1][v-1];
            }
        }
        free(post_roles);
    }
}

void Sampler_t::resample_hypers(void) {
    cout << "Sampling hyperparameters..." << endl;
    double bdelta = 5;
    double adelta = 0.1;
    double aalpha = 5;
    double balpha = 5.1;

    for (unsigned int iters=0; iters<10; iters++) {
        //sample delta
        double eta = sample_Bet(delta+1, tables);
        double bloge = bdelta - log(eta);
        double pie = 1.0 / (1.0-(tables*bloge/(delta+used_frames.size()-1)));
        int u = sample_Bern(pie);
        delta = sample_Gam(adelta+used_frames.size()-1+u, 1.0/bloge);

        //sample alpha
        double qs = 0;
        double qw = 0;
        for (unsigned int u=1; u<=U; ++u) {
            qs += sample_Bern(w[u-1].size() / (w[u-1].size() + alpha));
            qw += log(sample_Bet(alpha+1, w[u-1].size()));
        }
        alpha = sample_Gam(aalpha + tables - qs, 1.0/(balpha - qw));
    }
    cout << "...alpha = " << alpha << ", delta = " << delta << endl;
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


bool Sampler_t::loadData(string inputFileName) {

    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

    inputFile = inputFileName;

    ifstream ifs(inputFileName.c_str());
    if (!ifs.is_open()) {
        cout << "Cannot open file " << inputFileName << ".\n";
        return false;
    }

    string line;
    unsigned long int progress = 0;
    while (getline(ifs, line)) {
        progress++;
        boost::char_separator<char> sep("\t");
        tokenizer tokens(line, sep);
        vector<vector<unsigned int> > unit;
        unsigned long int itemcounter = 0;
        for (tokenizer::iterator tok_iter = tokens.begin();
            tok_iter != tokens.end(); ++tok_iter) {
            itemcounter++;
            boost::char_separator<char> sep(" ");
            tokenizer slots(*tok_iter, sep);

            vector<unsigned int> s;
            for (tokenizer::iterator slot_iter = slots.begin();
                slot_iter != slots.end(); ++slot_iter) {
                const unsigned int w = atoi(slot_iter->c_str());
                if (w == 0) {
                    cout << "Invalid input data: " << *slot_iter << endl;
                    ifs.close();
                    return false;
                }
                if (V < w) V = w;
                s.push_back(w);
            }
            if (S==0) {
                S = s.size();
            } else{
                if (S != s.size()) {
                    cout << "Inconsitent number of slots at line no. " << 
                        progress <<", item no. " << itemcounter << " (" << 
                        *tok_iter << "). \n";
                    ifs.close();
                    return false;
                }
            }
            unit.push_back(s);
        }
        w.push_back(unit);
    }
    U = w.size();

    ifs.close();


    if (F == 0) {
        cout << "F = automatic" << endl;
        infinite_F = true;
        F = 1;
    } else {
        cout << "F = " << F << endl;
    }
    if (R == 0) {
        infinite_R = true;
        //R = floor(exp(log(F)/S)) + 1; //minimum number of roles
        R = ceil(exp(log(F)/S)); //minimum number of roles
        cout << "R = automatic (min. " << R << ")" << endl;
    } else {
        cout << "R = " << R << endl;
    }

    
    if (F > pow(R, S)) {
        cout << "Number of frames (F) must be lower than or equal to the number of all " <<
                "combinations of possible semantic roles (R^S)." << endl;
        return false;
    }
    
    frameSet = new Frames_t(S);

    cout << "alpha = " << alpha << endl;
    cout << "beta0 = " << beta0 << endl;
    cout << "gamma = " << gamma << endl;
    cout << "delta = " << delta << endl;
    cout << "Lexical units = " << U << endl;
    cout << "Slots = " << S << endl;
    cout << "Vocabulary size = " << V << endl;


    return true;
}

bool Sampler_t::initialize(bool recovery) {
    setall(time(0),time(0));   /* initialize random number generator */
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
        //resample_tau();
    
    }

    initialized = true;
    return true;    

}

Sampler_t::~Sampler_t() {
    
    if (initialized) {

        for (unsigned int f=1; f<=F; ++f) {
            free(roles[f-1]);
        }
        free(roles);

        for (unsigned int r=1; r<=R; ++r) {
            free(post_theta[r-1]);
        }
        free(post_theta);

        free(post_omega);

        for (unsigned int u=1; u<=U; ++u) {
            free(post_phi[u-1]);
            free(frames[u-1]);
        }
        free(post_phi);
        free(frames);

        delete frameSet;
        
        for (unsigned int r=0; r<R; ++r) {
            free(beta[r]);
        } 
    }

}


void Sampler_t::sample(void) {
    cout << " frames..." << flush;
    if (infinite_F) {
        resample_frames_inf();
    } else {
        resample_frames();
    }
    cout << "roles..." << flush;
    if (infinite_R) {    
        resample_roles_inf();
    } else {
        resample_roles();
    }
    if (infinite_F) {
        cout << "tau..." << flush;
        resample_tau();
    }

    cout << "beta..." << flush;
    resample_beta(20);
}

bool Sampler_t::sampleAll(string outputDir, unsigned int iters, bool allSamples) {

    if (outputDir.at(outputDir.size()-1) != '/') outputDir += "/";
    int status = mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    if(iters < startIter) {
        cout << "All iterations already performed." << endl;
        return true;
    }   

    //TODO check outputDir
    if (status != 0 && errno != EEXIST) {
        cout << "Cannot create directory '" << outputDir << "'\n";
        return false;
    }

    //remove old files
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(outputDir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << outputDir << endl;
        return false;
    }

    while ((dirp = readdir(dp)) != NULL) {
        string fname = string(dirp->d_name);
        if (fname.length() > 5 && 
            fname.substr(fname.length()-5, 5).compare(".smpl") == 0) {
            if (remove((outputDir + fname).c_str()) != 0) {
                cout << "Cannot remove file " << outputDir + fname << endl;
                return false;
            }
        }
    }
    closedir(dp);

    for (unsigned int i = startIter; i < iters+1; ++i) {
        cout << "Iteration no. " << i << ":";
        cout << flush;
        sample();
        cout << "perplexity..." << flush;
        double p = perplexity();
        cout << " (" << used_frames.size() << " frames, " 
             << used_roles.size() << " roles).";
        cout << " Perplexity: " << p;
        cout << endl;
        /*
        if (i>100 && infinite_F) {
            resample_hypers();
        }*/


        stringstream ss;
        if (allSamples) {
            ss << i << "-";
        }

        if(infinite_F || infinite_R) {
            cout << "...packing frames and roles.";
            pack_FR();
            cout << endl;
        }

        if (!dump(outputDir + ss.str())) {
            return false;
        }
        if (!writeLog(outputDir, i, iters)) {
            return false;
        }
    }
    
    return true;
}

bool Sampler_t::dump(string prefix) {

    string ffn = prefix + "frames.smpl";
    string rfn = prefix + "roles.smpl";

    ofstream ffile(ffn.c_str());
    if (!ffile.is_open()) {
        cout << "Cannot open file '" << ffn << "\n";
        return false;
    }
    ofstream rfile(rfn.c_str());
    if (!rfile.is_open()) {
        cout << "Cannot open file '" << rfn << "\n";
        ffile.close();
        return false;
    }

    for (unsigned int u=1; u<=U; ++u) {
        for (unsigned int t=1; t<=w[u-1].size(); ++t) {
            ffile << frames[u-1][t-1];
            if (t != w[u-1].size()) ffile << " ";
        }
        if (u != U) ffile << endl;
    }

    for (set<unsigned int>::const_iterator it = used_frames.begin(); it != used_frames.end(); ++it) {
        for (unsigned int s=1; s<=S; ++s) {
            rfile << roles[*it-1][s-1];
            if (s != S) rfile << " ";
        }
        if (it != used_frames.end()) rfile << endl;
    }

    ffile.close();
    rfile.close();
    return true;
}

bool Sampler_t::writeLog(string outputDir, unsigned int citer, unsigned int aiter) {
    string fn = outputDir + "lda-frames.log";
    ofstream lfile(fn.c_str());
    if (!lfile.is_open()) {
        cout << "Cannot open file '" << fn << "\n";
        return false;
    }

    lfile << "Input file name:\t" << inputFile << endl;
    lfile << "Number of frames:\t" << F << endl;
    lfile << "Number of roles:\t" << R << endl;
    lfile << "Number of slots:\t" << S << endl;
    lfile << "Alpha:\t" << alpha << endl;
    lfile << "Beta0:\t" << beta0 << endl;
    lfile << "Gamma:\t" << gamma << endl;
    lfile << "Delta:\t" << delta << endl;
    lfile << "Required number of iterations:\t" << aiter << endl;
    lfile << "Last iteration:\t" << citer << endl; 
    lfile.close();

    return true;
}

bool Sampler_t::dumpBest(string outputDir) {
    vector<ifstream *> fsamples;
    vector<vector<unsigned int> > frames, roles;
    
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(outputDir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << outputDir << endl;
        return false;
    }

    while ((dirp = readdir(dp)) != NULL) {
        string fname = string(dirp->d_name);
        if (fname.length() > 11 && 
            fname.substr(fname.length()-11, 11).compare("frames.smpl") == 0) {
            fsamples.push_back(new ifstream((outputDir + fname).c_str()));
            if (!fsamples.back()->is_open()) {
                cout << "Cannot open '" << outputDir + fname << "'." << endl;
                for (size_t i = 0; i < fsamples.size(); ++i) {
                    fsamples.at(i)->close();
                    delete fsamples.at(i);
                }
                return false;
            }
            cout << outputDir + fname << endl;
        }
        
    }
    closedir(dp);

    for (unsigned int u = 1; u<=U; ++u) {
        vector<unsigned int> uvect;
        for (unsigned int t = 1; t<=w[u-1].size(); ++t) {
            string tmp;
            vector<unsigned int> bestFrames(F, 0);
            for (size_t i=0; i<fsamples.size(); ++i) {
                *fsamples.at(i) >> tmp;
                unsigned int f = atoi(tmp.c_str());
                if (f <= 0 || f > F) {
                    cout << "Some frame(s) is not in the range <1, " << F << "'." << endl;
                    for (size_t i = 0; i < fsamples.size(); ++i) {
                        fsamples.at(i)->close();
                        delete fsamples.at(i);
                    }
                    return false;
                }
                bestFrames[f]++;
            }
            
            unsigned int freq = 0, bestFrame = 1;
            for (unsigned int f = 1; f <= F; ++f) {
                if (bestFrames[f-1]>freq) {
                    freq = bestFrames[f-1];
                    bestFrame = f;
                }
            }
            uvect.push_back(bestFrame);
        }
        frames.push_back(uvect);
    }

    for (unsigned int u=1; u<=U; ++u) {
        for (unsigned int t=1; t<=w[u-1].size(); ++t) {
            cout << frames[u-1][t-1];
            if (t != w[u-1].size()) cout << " ";
        }
        if (u != U) cout << endl;
    }
    

    for (size_t i = 0; i < fsamples.size(); ++i) {
        fsamples.at(i)->close();
        delete fsamples.at(i);
    }
    return true;
}

void Sampler_t::printFrames(void) {
    for (unsigned int u=1; u<=U; ++u) {
        for (unsigned int t=1; t<=w[u-1].size(); ++t) {
            cout << frames[u-1][t-1];
            if (t != w[u-1].size()) cout << " ";
        }
        if (u != U) cout << endl;
    }
}



void Sampler_t::printRoles(void) {
    for (set<unsigned int>::const_iterator it = used_frames.begin(); it != used_frames.end(); ++it) {
        for (unsigned int s=1; s<=S; ++s) {
            cout << roles[*it-1][s-1];
            if (s != S) cout << " ";
        }
        if (it != used_frames.end()) cout << endl;
    }
}

bool Sampler_t::recoverParameters(string logDir) {

    string fname = logDir + "lda-frames.log";
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

    ifstream lfile(fname.c_str());
    if (!lfile.is_open()) {
        cout << "Cannot open log file '" << lfile << "\n";
        return false;
    }
    //parsing the log file
    string line;
    unsigned long int progress = 0;
    while (getline(lfile, line)) {
        progress++;
        boost::char_separator<char> sep("\t");
        tokenizer tokens(line, sep);
        vector<string> lineItems;
        for (tokenizer::iterator tok_iter = tokens.begin();
                tok_iter != tokens.end(); ++tok_iter) {
            lineItems.push_back(*tok_iter);
        }
        if (lineItems.size() == 2) {
            if (lineItems.at(0) == "Number of frames:") {
                F = atoi(lineItems.at(1).c_str());
                if (F <= 0) {
                    cout << "Wrong number of frames in the log file (" << lineItems.at(1) << ")." << endl;
                    return false;
                }
            }
            if (lineItems.at(0) == "Number of roles:") {
                R = atoi(lineItems.at(1).c_str());
                if (R <= 0) {
                    cout << "Wrong number of roles in the log file (" << lineItems.at(1) << ")." << endl;
                    return false;
                }
            }
            if (lineItems.at(0) == "Alpha:") {
                alpha = atof(lineItems.at(1).c_str());
                if (alpha <= 0) {
                    cout << "Wrong alpha in the log file (" << lineItems.at(1) << ")." << endl;
                    return false;
                }
            }
            if (lineItems.at(0) == "Beta0:") {
                beta0 = atof(lineItems.at(1).c_str());
                if (beta0 <= 0) {
                    cout << "Wrong beta0 in the log file(" << lineItems.at(1) << ")." << endl;
                    return false;
                }
            }
            if (lineItems.at(0) == "Gamma:") {
                gamma = atof(lineItems.at(1).c_str());
                if (gamma <= 0) {
                    cout << "Wrong gamma in the log file(" << lineItems.at(1) << ")." << endl;
                    return false;
                }
            }
            if (lineItems.at(0) == "Delta:") {
                delta = atof(lineItems.at(1).c_str());
                if (delta <= 0) {
                    cout << "Wrong delta in the log file(" << lineItems.at(1) << ")." << endl;
                    return false;
                }
            }
            if (lineItems.at(0) == "Last iteration:") {
                startIter = atoi(lineItems.at(1).c_str());
                if (startIter <= 0) {
                    cout << "Wrong last iteration number in the log file (" << 
                            lineItems.at(1) << ")." << endl;
                    return false;
                }
                cout << "The last iteration was " << startIter << "." << endl;
                startIter++;
            }
        }
        
    }

    lfile.close();
    return true;
}

bool Sampler_t::recoverData(string dataDir) {
    string ffname = dataDir + "frames.smpl";
    string rfname = dataDir + "roles.smpl";
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

    ifstream ffile(ffname.c_str());
    if (!ffile.is_open()) {
        cout << "Cannot open file '" << ffname << "'. (TODO: find the file " <<
                "with latest sample automatically.)\n";
        return false;
    }
    
    ifstream rfile(rfname.c_str());
    if (!rfile.is_open()) {
        cout << "Cannot open file '" << rfname << "'. (TODO: find the file " <<
                "with latest sample automatically.)\n";
        return false;
    }
    //parsing frames.smpl
    cout << "...parsing '" << dataDir << "frames.smpl'." << endl;
    string line;
    unsigned long int u = 0;
    while (getline(ffile, line)) {
        u++;
        if (u > U) {
            cout << "Wrong number of lexical units." << endl;
            ffile.close();
            rfile.close();
            return false;
        }
        boost::char_separator<char> sep(" ");
        tokenizer tokens(line, sep);
        unsigned long int t = 0;
        for (tokenizer::iterator tok_iter = tokens.begin();
                tok_iter != tokens.end(); ++tok_iter) {
            t++;
            if (t > w[u-1].size()) {
                cout << "Wrong number of realizations (line #" << u <<")." << endl;
                ffile.close();
                rfile.close();
                return false;
            }
            unsigned int frame = atoi(tok_iter->c_str());
            if (frame <= 0 || frame > F) {
                cout << "Wrong frame number (" << *tok_iter << ")." << endl;
                ffile.close();
                rfile.close();
                return false;
            }
            //recover data
            frames[u-1][t-1] = frame;
            used_frames.insert(frame);
            //tau[frame] = 0;

            fc_f[frame-1]++;
            for (unsigned int s=1; s<=S; ++s) {
                fc_fsw[frame-1][s-1][w[u-1][t-1][s-1]-1]++;
            }
        }
        if (t != w[u-1].size()) {
            cout << "Wrong number of realizations (line #" << u <<")." << endl;
            ffile.close();
            rfile.close();
            return false;
        }

    }
    if (u != U) {
        cout << "Wrong number of lexical units." << endl;
        ffile.close();
        rfile.close();
        return false;
    }

    //resample_tau();
    
    //parsing roles.smpl
    cout << "...parsing '" << dataDir << "roles.smpl'." << endl;
    unsigned long int f = 0;
    while (getline(rfile, line)) {
        f++;
        if (f > F) {
            cout << "Wrong number of frames." << endl;
            ffile.close();
            rfile.close();
            return false;
        }
        boost::char_separator<char> sep(" ");
        tokenizer tokens(line, sep);
        unsigned long int s = 0;
        for (tokenizer::iterator tok_iter = tokens.begin();
                tok_iter != tokens.end(); ++tok_iter) {
            s++;
            if (s > S) {
                cout << "Wrong number of slots (line #" << f <<")." << endl;
                ffile.close();
                rfile.close();
                return false;
            }
            unsigned int role = atoi(tok_iter->c_str());
            if (role <= 0 || role > R) {
                cout << "Wrong role number (" << *tok_iter << ")." << endl;
                ffile.close();
                rfile.close();
                return false;
            }
            //recover data
            roles[f-1][s-1] = role;
            used_roles.insert(role);
        }
        if (s != S) {
            cout << "Wrong number of slots (line #" << f <<")." << endl;
            ffile.close();
            rfile.close();
            return false;
        }

    }
    //recover data
    FrameKey_t fk = frameSet->makeKey(roles[f-1]);
    frameSet->insert(fk);
    if (f != F) {
        cout << "Wrong number of frames." << endl;
        ffile.close();
        rfile.close();
        return false;
    }

    cout << "...resampling phi and theta." << endl;
   
    initialize_post_phi();
    initialize_post_theta();
    initialize_post_omega();
    //resample_tau();

    ffile.close();
    rfile.close();

    return true;
}

unsigned int Sampler_t::createNewRole(void) {
    unsigned int role;
    if (!unused_roles.empty()) {
        set<unsigned int>::iterator it = unused_roles.begin();
        role = *it;
        unused_roles.erase(it);
    } else {
        R++;
        role = R;
        beta.push_back((double *) malloc(sizeof(double) * (V + 1)));
        post_theta = (double**) realloc(post_theta, sizeof(double*) * R);
        post_theta[R-1] = (double*) malloc(sizeof(double) * (V + 1));
        for (unsigned int v=1; v<=V; ++v) {
            post_theta[R-1][v-1] = 0;
        }
        post_theta[R-1][V] = 0;

        post_omega = (double*) realloc(post_omega, sizeof(double) * (R + 1));
        post_omega[R] = post_omega[R-1];
        post_omega[R-1] = 0;
        
    }
    beta[role-1][0] = V*beta0;
    for (unsigned int v=1; v<=V; ++v) beta[role-1][v] = beta0;

    used_roles.insert(role);
    return role;
}


bool Sampler_t::sample_new_frame(vector<unsigned int> &frame, vector<unsigned int> &pos) {
    for (unsigned int s=1; s<=S; ++s) {
        double* post_roles = (double*) malloc(sizeof(double) * (R + 2));
        for (unsigned int r = 1; r <= R; ++r) {
            post_roles[r-1] = -1 * numeric_limits<double>::max();
            set<unsigned int>::iterator it = used_roles.find(r);
            if (it != used_roles.end()) {
                post_roles[r-1] = 
                    ldf_Mult_smooth(1, beta[r-1], pos[s-1], post_theta[r-1], 1, V) +
                    log(gamma + post_omega[r-1]);
                    //ldf_Mult_smooth(0, gamma, r, post_omega, 1, R, used_roles.size());
            }
        }
        post_roles[R] = log(gamma) - log(V);
        normalizeLog(post_roles, 1, R + 1);
        unsigned int newRole = sample_Mult(post_roles, 1, R + 1);
        if (newRole == R + 1) {
            frame[s-1] = createNewRole();
        } else {
            frame[s-1] = newRole;
        }

        free(post_roles);
    }

    if(frameSet->inside(frameSet->makeKey(&frame[0]))) 
        return false;
    else
        return true;
}


unsigned int Sampler_t::createNewFrame(vector<unsigned int> &frame) {
    unsigned int frameId;
    if (!unused_frames.empty()) {
        set<unsigned int>::iterator it = unused_frames.begin();
        frameId = *it;
        unused_frames.erase(it);
    } else {
        F++;
        frameId = F;
        for (unsigned int u = 1; u <= U; ++u) {
            post_phi[u-1] = (double*) realloc(post_phi[u-1], sizeof(double) * (F + 1));
            post_phi[u-1][F] = post_phi[u-1][F-1];
            post_phi[u-1][F-1] = 0;
        }
        roles = (unsigned int**) realloc(roles, sizeof(unsigned int*) * F);
        roles[F-1] = (unsigned int*) malloc(sizeof(unsigned int) * S);
        fc_f.push_back(0);
        fc_fsw.push_back(vector<vector<unsigned int> >(S,vector<unsigned int>(V, 0)));
    }

    used_frames.insert(frameId);

    for (unsigned int s = 1; s <= S; ++s ) {
        roles[frameId-1][s-1] = frame[s-1];
        post_omega[roles[frameId-1][s-1]-1]++;
        post_omega[R]++;
    }

    FrameKey_t fk = frameSet->makeKey(roles[frameId-1]);
    frameSet->insert(fk);
    
    return frameId;
}

double Sampler_t::perplexity(void) {
    double loglik = 0;
    int words = 0;
       
    for (unsigned int u=1; u<=U; ++u) {
        for (unsigned int t=1; t <= w[u-1].size(); ++t) {
           unsigned int f = frames[u-1][t-1];
           loglik += log(post_phi[u-1][f-1] + alpha) -
                     log(post_phi[u-1][F] + used_frames.size()*alpha);
           for (unsigned int s=1; s<=S; ++s) {
                unsigned int r = roles[f-1][s-1];
                words++;
                loglik += log(post_theta[r-1][w[u-1][t-1][s-1]-1] + beta[r-1][w[u-1][t-1][s-1]]) -
                          log(post_theta[r-1][V] + beta[r-1][0]);
            }
        }
    }
    /*
    for (set<unsigned int>::const_iterator fit=used_frames.begin(); fit!=used_frames.end(); ++fit) {
        for (set<unsigned int>::const_iterator rit=used_roles.begin(); rit!=used_roles.end(); ++rit) {
            loglik += log(post_omega[*rit-1]+gamma) -
                      log(post_omega[R]+used_roles.size()*gamma); 
        }
    }*/
    /* 
    #pragma omp parallel for 
    for (unsigned int u=1; u<=U; ++u) {
        for (unsigned int t=1; t <= w[u-1].size(); ++t) {
            double tmp = 0;
            for (set<unsigned int>::const_iterator f =used_frames.begin(); f!=used_frames.end(); ++f) {
                double tmp2 = (post_phi[u-1][*f-1] + alpha)/(post_phi[u-1][F] + used_frames.size()*alpha);
                for (unsigned int s=1; s<=S; ++s) {
                    unsigned int r = roles[*f-1][s-1];
                    #pragma omp atomic
                    words++;
                    tmp2 *= (post_theta[r-1][w[u-1][t-1][s-1]-1] + beta[r-1][w[u-1][t-1][s-1]])/(post_theta[r-1][V] + beta[r-1][0]);
                }
                tmp += tmp2;
            
            }
            #pragma omp atomic            
            loglik += log(tmp);
        }
    }*/
    return exp(-loglik/words);
}

void Sampler_t::resample_beta(unsigned int iters) {
    for (set<unsigned int>::const_iterator rit = used_roles.begin(); rit!=used_roles.end(); ++rit) {
        for (unsigned int v=1; v<=V; ++v) {
            for (unsigned int iter = 0; iter < iters; ++iter) {
                double oldBeta = beta[*rit-1][v];
                beta[*rit-1][v] = max(
                    beta[*rit-1][v]*
                    (digamma(post_theta[*rit-1][v-1]+beta[*rit-1][v])-digamma(beta[*rit-1][v]))/
                    (digamma(post_theta[*rit-1][V]+beta[*rit-1][0])-digamma(beta[*rit-1][0])),
                    beta0);
                if (post_theta[*rit-1][V] == 0) beta[*rit-1][v] = beta0;
                beta[*rit-1][0] += beta[*rit-1][v] - oldBeta; 
            }
        
        }
    }
    
}

void Sampler_t::pack_FR(void) {

    //pack frames and roles    
    map<unsigned int, unsigned int> tmp_F, tmp_R;
    unsigned int f=0, r=0;
    for (set<unsigned int>::const_iterator it=used_frames.begin();
            it!= used_frames.end(); ++it) {
        tmp_F[*it] = ++f;
    }
    for (set<unsigned int>::const_iterator it=used_roles.begin();
            it!= used_roles.end(); ++it) {
        tmp_R[*it] = ++r;
    }
    for(unsigned int u=1; u<=U; ++u) {
        for (unsigned int t = 1; t <= w[u-1].size(); ++t) {
            frames[u-1][t-1] = tmp_F[frames[u-1][t-1]];
        }
    }
    unsigned int** tmp_roles; 
    fc_f.clear();
    fc_fsw.clear();
    tmp_roles = (unsigned int**) malloc(sizeof(unsigned int*) * used_frames.size());
    for (unsigned int f=1; f<=used_frames.size(); ++f) {
        tmp_roles[f-1] = (unsigned int*) malloc(sizeof(unsigned int) * S);
        fc_f.push_back(0);
        fc_fsw.push_back(vector<vector<unsigned int> >(S,vector<unsigned int>(V, 0)));
    }

    map<unsigned int, double> tmp_tau;
    tmp_tau[0] = F;
    frameSet->clear();
    for (set<unsigned int>::const_iterator it = used_frames.begin(); it != used_frames.end(); ++it) {
        tmp_tau[tmp_F[*it-1]] = tau[*it];
        for (unsigned int s=1; s<=S; ++s) {
            tmp_roles[tmp_F[*it]-1][s-1] = tmp_R[roles[*it-1][s-1]];
        }
        frameSet->insert(frameSet->makeKey(tmp_roles[tmp_F[*it]-1]));
    }
    tau = tmp_tau;

    for (unsigned int f=1; f<=F; ++f) {
        free(roles[f-1]);
    }
    free(roles);
    roles = tmp_roles;
    
    for(unsigned int u=1; u<=U; ++u) {
        for (unsigned int t = 1; t <= w[u-1].size(); ++t) {
            fc_f[frames[u-1][t-1]-1]++;
            for (unsigned int s=1; s<=S; ++s) {
                fc_fsw[frames[u-1][t-1]-1][s-1][w[u-1][t-1][s-1]-1]++;
            }
        }
    }

    //reduce number of frames and roles
    F = used_frames.size();
    R = used_roles.size();
    unused_frames.clear();
    unused_roles.clear();
    used_frames.clear();
    used_roles.clear();
    for (unsigned int f=1; f<=F; ++f) used_frames.insert(f);
    for (unsigned int r=1; r<=R; ++r) used_roles.insert(r);

    resample_post_phi();
    resample_post_theta();
    resample_post_omega();

}
