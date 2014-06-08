/*
 * Copyright (C) 2014 Jiri Materna <xmaterna@fi.muni.cz>
 * Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <sys/stat.h>
#include <sys/errno.h>
#include <dirent.h>
#include "sampler.h"

#define BOUNDPROB(x) (((x)<-300)?(-300):((((x)>300)?(300):(x))))

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
                    if (w[u-1][t-1][s-1] > 0) {
                        post_theta[r-1][V] += (r == roles[frames[u-1][t-1]-1][s-1]) ? 1 : 0;
                        post_theta[r-1][w[u-1][t-1][s-1]-1] += 
                            (r == roles[frames[u-1][t-1]-1][s-1]) ? 1 : 0;
                    }
                }
            }
        }
    }
}

void Sampler_t::resample_post_omega(void) {
    for (unsigned int r = 1; r <= R; ++r) {
        post_omega[r-1] = 0;
    }
    post_omega[R] = 0;
    for (unsigned int f = 1; f <= F; ++f) {
        for (unsigned int s = 1; s <= S; ++s) {
            if (roles[f-1][s-1] > 0) {
                post_omega[roles[f-1][s-1]-1]++;
                post_omega[R]++;
            }
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
                fc[id] += dist->sampleAntoniak(alpha0 * tau[*it], post_phi[u-1][*it-1]);
            } else {
                fc[id] += post_phi[u-1][*it-1];
            }
        }
        tables += fc[id];
        id++;
    }    
    fc[used_frames.size()] = delta;
    dist->sampleDirichlet(dirs, fc, used_frames.size() + 1);
    for (unsigned int i=0; i<used_frames.size(); ++i) {
        tau[idmap[i]] = dirs[i];
    }
    tau[0] = dirs[used_frames.size()]; 
    free(fc);
    free(dirs);
}


void Sampler_t::resample_frames(void) {
    #pragma omp parallel for
    for (int u = 1; u <= (signed int) U; ++u) {
        for (unsigned int t = 1; t <= w[u-1].size(); ++t) {

            map<unsigned short int, double> post_frames;

            //remove old values
            for (unsigned int s = 1; s <= S; ++s) {
                if (w[u-1][t-1][s-1] != 0) {
                    #pragma omp atomic
                    post_theta[roles[frames[u-1][t-1]-1][s-1]-1][V]--;
                    #pragma omp atomic
                    post_theta[roles[frames[u-1][t-1]-1][s-1]-1][w[u-1][t-1][s-1]-1]--;
                    #pragma omp atomic
                    fc_fsw[frames[u-1][t-1]-1][s-1][w[u-1][t-1][s-1]-1]--;
                }
            }
            post_phi[u-1][F]--;
            post_phi[u-1][frames[u-1][t-1]-1]--;
            #pragma omp atomic
            fc_f[frames[u-1][t-1]-1]--;

            //compute frame distribution
            for (set<unsigned int>::const_iterator fit = used_frames.begin(); fit != used_frames.end(); 
                    ++fit) {

                //skip frames that don't satisfy given pattern
                if (emptyFrames && !checkPattern(roles[*fit-1], w[u-1][t-1])) continue;

                double prod = 0;
                for (unsigned int s = 1; s <= S; ++s) {
                    if (roles[*fit-1][s-1] != 0) {
                        prod += BOUNDPROB(
                                log(post_theta[roles[*fit-1][s-1]-1][w[u-1][t-1][s-1]-1] + 
                                beta[w[u-1][t-1][s-1]-1]) -
                                log(post_theta[roles[*fit-1][s-1]-1][V] + beta[V])
                                ); 
                    }
                }
                post_frames[*fit] = prod + BOUNDPROB(log(post_phi[u-1][*fit-1] + alpha[*fit-1]));
            }

            //sample frame
            dist->normalizeLogMult(post_frames);
            frames[u-1][t-1] = dist->sampleMultinomial(post_frames);

            //update new values
            for (unsigned int s=1; s<=S; ++s) {
                if (w[u-1][t-1][s-1] != 0) {
                    #pragma omp atomic
                    post_theta[roles[frames[u-1][t-1]-1][s-1]-1][V]++;
                    #pragma omp atomic
                    post_theta[roles[frames[u-1][t-1]-1][s-1]-1][w[u-1][t-1][s-1]-1]++;
                    #pragma omp atomic
                    fc_fsw[frames[u-1][t-1]-1][s-1][w[u-1][t-1][s-1]-1]++;
                }
            }
            post_phi[u-1][F]++;
            post_phi[u-1][frames[u-1][t-1]-1]++;
            #pragma omp atomic
            fc_f[frames[u-1][t-1]-1]++;
        }
    }
}


void Sampler_t::resample_frames_inf(void) {
    for (int u = 1; u <= (signed int) U; ++u) {
        for (unsigned int t = 1; t <= w[u-1].size(); ++t) {

            map<unsigned short int, double> post_frames;
            bool tau_needs_resampling = false;

            //remove old values
            if (frames[u-1][t-1] > 0) { //already sampled?
                for (unsigned int s = 1; s <= S; ++s) {
                    if(w[u-1][t-1][s-1] != 0) {
                        post_theta[roles[frames[u-1][t-1]-1][s-1]-1][V]--;
                        post_theta[roles[frames[u-1][t-1]-1][s-1]-1][w[u-1][t-1][s-1]-1]--;
                        fc_fsw[frames[u-1][t-1]-1][s-1][w[u-1][t-1][s-1]-1]--;
                    }

                }
                post_phi[u-1][F]--;
                post_phi[u-1][frames[u-1][t-1]-1]--;
                fc_f[frames[u-1][t-1]-1]--;
            }

            //compute frame distribution
            for (set<unsigned int>::const_iterator fit = used_frames.begin(); fit !=used_frames.end(); 
                    ++fit) {

                //skip frames that don't satisfy given pattern
                if (emptyFrames && !checkPattern(roles[*fit-1], w[u-1][t-1])) continue;

                //probability of used frames
                double prod = 0;
                for (unsigned int s = 1; s <= S; ++s) {
                    if (w[u-1][t-1][s-1] !=0 ) {
                        prod += BOUNDPROB(
                                log(post_theta[roles[*fit-1][s-1]-1][w[u-1][t-1][s-1]-1] +
                                beta[w[u-1][t-1][s-1]-1]) -
                                log(post_theta[roles[*fit-1][s-1]-1][V] + beta[V])
                                );
                    }
                }
                post_frames[*fit] = prod + BOUNDPROB(log(post_phi[u-1][*fit-1] + alpha0*tau[*fit]));
            }
            //sample new frame
            vector<unsigned int> frame(S, 0);
            if (sample_new_frame(frame, w[u-1][t-1])) {
                double prod = BOUNDPROB(log(alpha0*tau[0]));
                for (unsigned int s = 1; s <= S; ++s) {
                    prod -= BOUNDPROB(log(V));
                }
                post_frames[F+1] = prod;
            }

            //sample frame id
            dist->normalizeLogMult(post_frames);
            unsigned int newFrame = dist->sampleMultinomial(post_frames);
            
            //free unused frames and roles
            if (frames[u-1][t-1] != newFrame && frames[u-1][t-1] != 0 && fc_f[frames[u-1][t-1]-1] == 0) {
                frameSet.remove(frameSet.makeKey(roles[frames[u-1][t-1]-1]));
                for (unsigned int s=1; s<=S; ++s) {
                    if(w[u-1][t-1][s-1] != 0) {
                        post_omega[roles[frames[u-1][t-1]-1][s-1]-1]--;
                        post_omega[R]--;
                        if (post_omega[roles[frames[u-1][t-1]-1][s-1]-1] == 0 && infinite_R) {
                            unused_roles.insert(roles[frames[u-1][t-1]-1][s-1]);
                            used_roles.erase(roles[frames[u-1][t-1]-1][s-1]);
                            gamma.erase(roles[frames[u-1][t-1]-1][s-1]);
                        }
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
            } else {
                frames[u-1][t-1] = newFrame;
                
                //remove possibly new roles of the unused new frame
                for (unsigned int s=1; s<=S; ++s) {
                    if(w[u-1][t-1][s-1] != 0) {
                        if (frame[s-1] != 0 && post_omega[frame[s-1]-1] == 0 && infinite_R) {
                            used_roles.erase(frame[s-1]);
                            unused_roles.insert(frame[s-1]);
                            gamma.erase(frame[s-1]);
                        }
                    }
                }
            }
            
            //update new values
            for (unsigned int s=1; s<=S; ++s) {
                if(w[u-1][t-1][s-1] != 0) {
                    post_theta[roles[frames[u-1][t-1]-1][s-1]-1][V]++;
                    post_theta[roles[frames[u-1][t-1]-1][s-1]-1][w[u-1][t-1][s-1]-1]++;
                    fc_fsw[frames[u-1][t-1]-1][s-1][w[u-1][t-1][s-1]-1]++;
                }
            }
            post_phi[u-1][F]++;
            post_phi[u-1][frames[u-1][t-1]-1]++;
            fc_f[frames[u-1][t-1]-1]++;
                    
            if(tau_needs_resampling) resample_tau();

        }
    }
}

void Sampler_t::resample_roles(void) {
   
    #pragma omp parallel for
    for (unsigned int f=1; f<=used_frames.size(); ++f) {
        set<unsigned int>::const_iterator fit = used_frames.begin();
        advance(fit, f-1);
    //for (set<unsigned int>::const_iterator fit = used_frames.begin(); fit!=used_frames.end(); ++fit) {
        //select frame pattern
        vector<unsigned int> pattern;
        if( fc_f[*fit-1] == 0) {
            //sample pattern 
            pattern = dist->sampleMultinomial(framePatterns, positions);
        } else {
            pattern = roles[*fit-1];
        }

        for (unsigned int s = 1; s <= S; ++s) {
            FrameKey_t oldFrame; 
            oldFrame = frameSet.makeKey(roles[*fit-1]);

            //empty slots
            if (pattern[s-1] == 0) {
                if (roles[*fit-1][s-1] != 0) {
                    #pragma omp critical
                    {
                        FrameKey_t newFrame = frameSet.makeKey(roles[*fit-1], s, 0);
                        bool inside = frameSet.inside(newFrame);
                        if (!inside) {
                            post_omega[roles[*fit-1][s-1]-1]--;
                            post_omega[R]--;
                            frameSet.remove(oldFrame);
                            frameSet.insert(newFrame);
                        }
                    }
                        roles[*fit-1][s-1] = 0;
                        continue;
                } else {
                    continue;
                }
            }


            map<unsigned short int, double> post_roles;

            if (roles[*fit-1][s-1] != 0) {
                #pragma omp critical
                {
                    post_theta[roles[*fit-1][s-1]-1][V] -= fc_f[*fit-1];
                    for(unsigned int v=1; v<=V; ++v) {
                        post_theta[roles[*fit-1][s-1]-1][v-1] -= fc_fsw[*fit-1][s-1][v-1];
                    }
                    post_omega[roles[*fit-1][s-1]-1]--;
                    post_omega[R]--;
                }
            }
           
            unsigned int possibilities = 0;
            for (set<unsigned int>::const_iterator rit = used_roles.begin(); rit != used_roles.end();
                    ++rit) {
                FrameKey_t newFrame = frameSet.makeKey(roles[*fit-1], s, *rit);
                bool inside;
                #pragma omp critical
                {
                    inside = frameSet.inside(newFrame);
                }
                if (newFrame == oldFrame || !inside) {
                    double prod = 0;
                    for (unsigned int v = 1; v<=V; ++v) {
                        prod += fc_fsw[*fit-1][s-1][v-1]*BOUNDPROB(                                   
                            log(post_theta[*rit-1][v-1] +
                            beta[v-1]) -
                            log(post_theta[*rit-1][V] + beta[V])
                            );
                    }
                    
                    post_roles[*rit] = prod + BOUNDPROB(log(post_omega[*rit-1] + gamma[*rit]));
                    possibilities++;
                }
            }

            //sample role
            if (possibilities > 0) {
                dist->normalizeLogMult(post_roles);
                roles[*fit-1][s-1] = dist->sampleMultinomial(post_roles);
                #pragma omp critical
                {
                    frameSet.remove(oldFrame);
                    frameSet.insert(frameSet.makeKey(roles[*fit-1]));
                    post_omega[roles[*fit-1][s-1]-1]++;
                    post_omega[R]++;
                    post_theta[roles[*fit-1][s-1]-1][V] += fc_f[*fit-1];
                    for(unsigned int v=1; v<=V; ++v) {
                        post_theta[roles[*fit-1][s-1]-1][v-1] += fc_fsw[*fit-1][s-1][v-1];
                    }
                }
            }
        }
    }
}

void Sampler_t::resample_roles_inf(void) {

    for (set<unsigned int>::const_iterator fit = used_frames.begin(); fit!=used_frames.end(); ++fit) {
        for (unsigned int s = 1; s <= S; ++s) {

            //empty slot
            if (roles[*fit-1][s-1] == 0) continue;

            map<unsigned short int, double> post_roles;

            post_theta[roles[*fit-1][s-1]-1][V] -= fc_f[*fit-1];
            for(unsigned int v=1; v<=V; ++v) {
                post_theta[roles[*fit-1][s-1]-1][v-1] -= fc_fsw[*fit-1][s-1][v-1];
            }

            post_omega[roles[*fit-1][s-1]-1]--;
            post_omega[R]--;

            FrameKey_t oldFrame; 
            oldFrame = frameSet.makeKey(roles[*fit-1]);
            
            for (set<unsigned int>::const_iterator rit = used_roles.begin(); rit != used_roles.end(); ++rit) {
                FrameKey_t newFrame = frameSet.makeKey(roles[*fit-1], s, *rit);
                bool inside = frameSet.inside(newFrame);

                //probability of used roles
                if (newFrame == oldFrame || !inside) {
                    double prod = 0;
                    for (unsigned int v = 1; v<=V; ++v) {
                        prod += fc_fsw[*fit-1][s-1][v-1]*BOUNDPROB(
                            log(post_theta[*rit-1][v-1] +
                            beta[v-1]) -
                            log(post_theta[*rit-1][V] + beta[V])
                            );
                    }
                    post_roles[*rit] = prod + BOUNDPROB(log(post_omega[*rit-1])); // + gamma[*rit]));
                }
            }

            //probability of a new role
            double prod = 0;
            for (unsigned int v = 1; v<=V; ++v) {
                prod -= fc_fsw[*fit-1][s-1][v-1]*log(V);
            }
            post_roles[R + 1] = prod + BOUNDPROB(log(gamma[0]));
            

            //Sample role id
            dist->normalizeLogMult(post_roles);
            unsigned int newRole = dist->sampleMultinomial(post_roles);

            //free unused role numbers
            if (roles[*fit-1][s-1] != newRole && post_omega[roles[*fit-1][s-1]-1] == 0) {
                unused_roles.insert(roles[*fit-1][s-1]);
                used_roles.erase(roles[*fit-1][s-1]);
                gamma.erase(roles[*fit-1][s-1]);
            }
            //create a new role if required
            if (newRole == R + 1) {
                roles[*fit-1][s-1] = createNewRole();
            } else {
                roles[*fit-1][s-1] = newRole;
            }

            frameSet.remove(oldFrame);
            frameSet.insert(frameSet.makeKey(roles[*fit-1]));

            post_omega[roles[*fit-1][s-1]-1]++;
            post_omega[R]++;
            post_theta[roles[*fit-1][s-1]-1][V] += fc_f[*fit-1];
            for(unsigned int v=1; v<=V; ++v) {
                post_theta[roles[*fit-1][s-1]-1][v-1] += fc_fsw[*fit-1][s-1][v-1];
            }
        }
    }
}


void Sampler_t::predict_test(void) {
    #pragma omp parallel for
    for (int u = 1; u <= (signed int) U; ++u) {
        for (unsigned int t = 1; t <= test_w[u-1].size(); ++t) {

            double prob = 0;             

            //compute frame distribution
            for (set<unsigned int>::const_iterator fit = used_frames.begin(); fit != used_frames.end(); 
                    ++fit) {

                //skip frames that don't satisfy given pattern
                if (emptyFrames && !checkPattern(roles[*fit-1], test_w[u-1][t-1])) continue;

                double prod = 0;
                for (unsigned int s = 1; s <= S; ++s) {
                    if (roles[*fit-1][s-1] != 0) {
                        prod += BOUNDPROB(
                                log(post_theta[roles[*fit-1][s-1]-1][test_w[u-1][t-1][s-1]-1] + 
                                beta[test_w[u-1][t-1][s-1]-1]) -
                                log(post_theta[roles[*fit-1][s-1]-1][V] + beta[V])
                                ); 
                    }
                }
                prod += BOUNDPROB(log(post_phi[u-1][*fit-1] + alpha[*fit-1]));
                if (prod >= prob || fit == used_frames.begin()) {
                    prob = prod;
                    test_frames[u-1][t-1] = *fit;
                }
            }

        }
    }
}


void Sampler_t::resample_hypers(unsigned int iters) {

    //sample beta
    for (unsigned int v=1; v<=V; ++v) {
        for (unsigned int iter = 0; iter < iters; ++iter) {
            double oldBeta = beta[v-1];
            double nsum = 0;
            double dsum = 0;
            for (set<unsigned int>::const_iterator rit = used_roles.begin(); 
                        rit!=used_roles.end(); ++rit) {
                nsum += dist->digamma(post_theta[*rit-1][v-1]+beta[v-1]); 
                dsum += dist->digamma(post_theta[*rit-1][V]+beta[V]);
            }
            beta[v-1] *=
                (nsum - used_roles.size()*dist->digamma(beta[v-1]))/
                (dsum - used_roles.size()*dist->digamma(beta[V]));
            if (beta[v-1]<=0) beta[v-1] = oldBeta;
            beta[V] += beta[v-1] - oldBeta;
        }
    }

    if (infinite_R) {
        //sample gamma
        double bgamma = 1.0;
        double agamma = 1.0;
        double gamma0_sum = 0;
        double n = 0;
        for (set<unsigned int>::const_iterator fit = used_frames.begin(); fit!=used_frames.end(); ++fit) {
            for (unsigned int s = 1; s <= S; ++s) {
                if (roles[*fit-1][s-1] != 0) n++;
            }
        }
        for (unsigned int i = 1; i <= iters; ++i) {
            double eta = dist->sampleBeta(gamma0 + 1, n); 
            double pi = agamma + used_roles.size() - 1;
            //double rate = 1.0 / bgamma - log(eta);
            double rate = bgamma - log(eta);
            pi = pi / (pi + rate * n);
    
            unsigned int cc = dist->sampleBernoulli(pi);
            if (cc == 1) {
                gamma0 = dist->sampleGamma(agamma + used_roles.size(), 1.0 / rate);
            } else {
                gamma0 = dist->sampleGamma(agamma + used_roles.size() - 1, 1.0 / rate);
            }
            if (i > iters/2) gamma0_sum += gamma0;

        }
        //cout << endl;
        gamma0 = 2*gamma0_sum / iters;
        for (set<unsigned int>::const_iterator rit=used_roles.begin(); rit != used_roles.end(); ++rit) {
            gamma[*rit] = 2*gamma0_sum / iters;
            //cout << gamma[*rit] << " "; 
        }
        gamma[0] = 2*gamma0_sum / iters;
        //cout << gamma[0] << endl;
    }

    if (infinite_F) {
        double bdelta = 1.0;
        double adelta = 1.0;
        double aalpha = 1.0;
        double balpha = 1.0;
        double delta_sum = 0;
        double alpha0_sum = 0;

        for (unsigned int i=1; i<=iters; i++) {

            //sample delta
            double eta = dist->sampleBeta(delta + 1, tables);
            double pi = adelta + used_frames.size() - 1;
            //double rate = 1.0 / bdelta - log(eta);
            double rate = bdelta - log(eta);
            pi = pi / (pi + rate * tables);
    
            unsigned int cc = dist->sampleBernoulli(pi);
            if (cc == 1) {
                delta = dist->sampleGamma(adelta + used_frames.size(), 1.0 / rate);
            } else {
                delta = dist->sampleGamma(adelta + used_frames.size() - 1, 1.0 / rate);
            }
            if (i > iters/2) delta_sum += delta;
    
            //sample alpha
            double sum_log_w = 0.0;
            double sum_s = 0.0;
            for (unsigned int u=1; u<=U; ++u) {
                sum_log_w += log(dist->sampleBeta(alpha0 + 1, w[u-1].size()));
                //sum_s += (double)dist->sampleBernoulli(w[u-1].size() / (w[u-1].size() + alpha0));
                sum_s += (double)dist->sampleBernoulli(w[u-1].size() / alpha0);
            }
            //rate = 1.0 / balpha - sum_log_w;
            rate = balpha - sum_log_w;
            alpha0 = dist->sampleGamma(aalpha + tables - sum_s, 1.0 / rate);
            if (i > iters/2) alpha0_sum += alpha0;

        }
        delta = 2*delta_sum / iters;
        alpha0 = 2*alpha0_sum / iters;
        //cout << endl << delta << " " << alpha0 << endl;
        
    } else {
        //sample alpha
        for (set<unsigned int>::const_iterator fit=used_frames.begin();
                fit!= used_frames.end(); ++fit) {
            for (unsigned int iter = 0; iter < iters; ++iter) {
                double oldAlpha = alpha[*fit-1];
                double nsum = 0;
                double dsum = 0;
                for (unsigned int u=1; u<=U; ++u) {
                    nsum += dist->digamma(post_phi[u-1][*fit-1]+alpha[*fit-1]);
                    dsum += dist->digamma(post_phi[u-1][F]+alpha[*fit-1]);
                }
                alpha[*fit-1] *=
                    (nsum - U*dist->digamma(alpha[*fit-1]))/
                    (dsum - U*dist->digamma(alpha[F]));
                if (alpha[*fit-1]<=0) alpha[*fit-1] = oldAlpha;
                alpha[F] += alpha[*fit-1] - oldAlpha;
            }
        }
        /*
        for (set<unsigned int>::const_iterator fit=used_frames.begin();
                fit!= used_frames.end(); ++fit) {
            cout << alpha[*fit-1] << " ";
        }
        cout << alpha[F] << endl; 
        */
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
    
    //dumpHypers();
}

bool Sampler_t::sampleAll(string outputDir, unsigned int iters, unsigned int burn_in, bool allSamples,
        bool no_hypers, bool no_perplexity, bool rm) {


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
    if (rm) {
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
    }

    string histfn = outputDir + "history.txt";
    ofstream histfile;
    if (startIter == 1) {
        histfile.open(histfn.c_str(), std::ofstream::out);
    } else {
        histfile.open(histfn.c_str(), std::ofstream::out | std::ofstream::app);
    }
    if (!histfile.is_open()) {
        cout << "Cannot open file '" << histfn << "\n";
        return false;
    }
    
    if (startIter == 1) {
        histfile << "iteration\tframes\troles\ttrain_perplexity";
        if (testPhase) histfile << "\ttest_perplexity";
        histfile << endl;
    }


    for (unsigned int i = startIter; i < iters+1; ++i) {
        cout << "Iteration no. " << i << ":";
        cout << flush;
        if (i>burn_in && !no_hypers) {
            if(reestimate_F) infinite_F = true;
            if(reestimate_R) infinite_R = true;
        }
        sample();
        if (i>burn_in && !no_hypers) {
            cout << "hyperparameters..." << flush;
            resample_hypers(100);
        }
        double train_p=0;
        double test_p=0;
        if (!no_perplexity) {
            cout << "perplexity..." << flush;
            train_p = perplexity(false);
            if (testPhase) {
                test_p = perplexity(true);
            }
            if (train_p < bestPerplexity || i==1) bestPerplexity = train_p;
        }
        cout << " (" << used_frames.size() << " frames, " 
             << used_roles.size() << " roles).";
        if (!no_perplexity) {
            cout << " Perplexity: " << train_p;
            if (testPhase) {
                cout << ", test perplexity: " << test_p;
            }
        }

        if(infinite_F || infinite_R) {
            pack_FR();
        }

        if (train_p==bestPerplexity) {
            cout << " *";
            if (!dump(outputDir)) {
                return false;
            }
            if (!writeLog(outputDir, i, iters)) {
                return false;
            }
        }

        stringstream ss;
        if (allSamples) {
            ss << i << "-";
            if (!dump(outputDir + ss.str())) {
                return false;
            }
        }
        cout << endl;

        histfile << i << "\t" << F << "\t" << R << "\t" << train_p;
        if (testPhase) histfile << "\t" << test_p;
        histfile << endl << flush;
    }

    histfile.close();
    
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
        post_theta.push_back(vector<double>(V + 1, 0));
        post_omega.push_back(0);
        post_omega[R] = post_omega[R-1];
        post_omega[R-1] = 0;
    }
    used_roles.insert(role);
    gamma[role] = gamma0;
    return role;
}


bool Sampler_t::sample_new_frame(vector<unsigned int> &frame, vector<unsigned int> &pos) {
    for (unsigned int s=1; s<=S; ++s) {
        if (pos[s-1]==0) {
            frame[s-1] = 0;
            continue;
        }
        map<unsigned short int, double> post_roles;
        for (set<unsigned int>::const_iterator rit=used_roles.begin(); rit!=used_roles.end();
                ++rit) {
            post_roles[*rit] = BOUNDPROB(
                log(post_theta[*rit-1][pos[s-1]-1] + beta[pos[s-1]-1]) -
                log(post_theta[*rit-1][V] + beta[V]) +
                //log(gamma[*rit] + post_omega[*rit-1])
                log(post_omega[*rit-1])
                );
        }
        if (infinite_R) {
            post_roles[R + 1] = log(gamma[0]) - log(V);
        }

        //sample role id
        dist->normalizeLogMult(post_roles);
        unsigned int newRole = dist->sampleMultinomial(post_roles);
        if (newRole == R + 1) {
            frame[s-1] = createNewRole();
        } else {
            frame[s-1] = newRole;
        }
    }

    if(frameSet.inside(frameSet.makeKey(frame))) 
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
            post_phi[u-1].push_back(0);
            post_phi[u-1][F] = post_phi[u-1][F-1];
            post_phi[u-1][F-1] = 0;
        }
        roles.push_back(vector<unsigned int>(S, 0));
        fc_f.push_back(0);
        fc_fsw.push_back(vector<vector<unsigned int> >(S,vector<unsigned int>(V, 0)));
    }

    used_frames.insert(frameId);

    for (unsigned int s = 1; s <= S; ++s ) {
        roles[frameId-1][s-1] = frame[s-1];
        if (frame[s-1] != 0) {
            post_omega[roles[frameId-1][s-1]-1]++;
            post_omega[R]++;
        }
    }

    FrameKey_t fk = frameSet.makeKey(roles[frameId-1]);
    frameSet.insert(fk);
    
    return frameId;
}


bool Sampler_t::checkPattern(vector<unsigned int> &u, vector<unsigned int> &v) {
    if (u.size() != v.size()) return false;
    for (unsigned int i=0; i<u.size(); ++i) if ((bool) u[i] != (bool) v[i]) return false;
    return true;
}

/*
double Sampler_t::perplexity(bool test) {
    vector<vector<vector<unsigned int> > > *words = &w; 
    vector<vector<unsigned int> > *fr = &frames;
    if(test) {
        words = &test_w;
        fr = &test_frames;
    } 

    double loglik = 0;
    int sum = 0;
    #pragma omp parallel for
    for (unsigned int u=1; u<=U; ++u) {
        double loglik_tmp = 0;
        int sum_tmp = 0;    
        for (unsigned int t=1; t <= words->at(u-1).size(); ++t) {
           unsigned int f = fr->at(u-1)[t-1];
           loglik_tmp += BOUNDPROB(log(post_phi[u-1][f-1]) -
                     log(post_phi[u-1][F]));
           for (unsigned int s=1; s<=S; ++s) {
                if (roles[f-1][s-1] > 0) {
                    unsigned int r = roles[f-1][s-1];
                    sum_tmp++;
                    loglik_tmp += BOUNDPROB(log(post_theta[r-1][words->at(u-1)[t-1][s-1]-1]) -
                              log(post_theta[r-1][V]));
                }
            }
        }
        #pragma omp atomic
        loglik += loglik_tmp;
        #pragma omp atomic
        sum += sum_tmp;
    }
    return exp(-loglik/sum);
}*/

double Sampler_t::perplexity(bool test) {
    vector<vector<vector<unsigned int> > > *words = &w; 
    if(test) {
        words = &test_w;
    } 
    double loglik = 0;
    int sum = 0;
    
    double tau_sum = 0;
    double alpha_sum = 0;
    for (set<unsigned int>::const_iterator fit = used_frames.begin(); fit != used_frames.end();
                    ++fit) {
        tau_sum += alpha0*tau[*fit];
        alpha_sum += alpha[*fit];
    }

    #pragma omp parallel for
    for (unsigned int u=1; u<=U; ++u) {
        double loglik_local = 0;
        int sum_local = 0;
        for (unsigned int t=1; t <= words->at(u-1).size(); ++t) {
            for (unsigned int s=1; s<=S; ++s) {
                double loglik_tmp = 0;
                for (set<unsigned int>::const_iterator fit = used_frames.begin(); fit != used_frames.end();
                        ++fit) {
                    if (!checkPattern(roles[*fit-1], w[u-1][t-1]) || roles[*fit-1][s-1] == 0) {
                         continue;
                    }

                    //loglik_tmp += post_phi[u-1][*fit-1]*(post_theta[roles[*fit-1][s-1]-1][words->at(u-1)[t-1][s-1]-1] + beta[words->at(u-1)[t-1][s-1]-1])
                    //              /
                    //              (post_phi[u-1][F]*(post_theta[roles[*fit-1][s-1]-1][V] + beta[V]));
                    if (infinite_F) {
                        loglik_tmp += (alpha0*tau[*fit] + post_phi[u-1][*fit-1])*(post_theta[roles[*fit-1][s-1]-1][words->at(u-1)[t-1][s-1]-1] + beta[words->at(u-1)[t-1][s-1]-1])
                                      /
                                      ((tau_sum + post_phi[u-1][F])*(post_theta[roles[*fit-1][s-1]-1][V] + beta[V]));
                    } else {
                        loglik_tmp += (alpha[*fit] + post_phi[u-1][*fit-1])*(post_theta[roles[*fit-1][s-1]-1][words->at(u-1)[t-1][s-1]-1] + beta[words->at(u-1)[t-1][s-1]-1])
                                      /
                                      ((alpha_sum + post_phi[u-1][F])*(post_theta[roles[*fit-1][s-1]-1][V] + beta[V]));

                    }


                }
                if (loglik_tmp !=0) {
                    loglik_local += BOUNDPROB(log(loglik_tmp));
                    sum_local++;
                }
            }
        }
        #pragma omp critical
        {
            loglik += loglik_local;
            sum += sum_local;
        }
    }
    return exp(-loglik/sum);
}

void Sampler_t::dumpHypers(void) {
    cout << "Alpha0: " << alpha0 << endl; 
    cout << "Alpha:" << endl; 
    for (vector<double>::const_iterator it = alpha.begin(); it != alpha.end(); ++it) {
        cout << *it << " ";
    }
    cout << endl;
    
    cout << "Beta0: " << beta0 << endl; 
    cout << "Beta:" << endl; 
    for (vector<double>::const_iterator it = beta.begin(); it != beta.end(); ++it) {
        cout << *it << " ";
    }
    cout << endl;
    
    cout << "Gamma0: " << gamma0 << endl; 
    cout << "Gamma:" << endl; 
    for (map<unsigned int, double>::const_iterator it = gamma.begin(); it != gamma.end(); ++it) {
        cout << "(" << it->first << "," << it->second << ") ";
    }
    cout << endl;
    
    cout << "Delta: " << delta << endl; 
    
    cout << "Tau:" << endl; 
    for (map<unsigned int, double>::const_iterator it = tau.begin(); it != tau.end(); ++it) {
        cout << "(" << it->first << "," << it->second << ") ";
    }
    cout << endl;
}



void Sampler_t::pack_FR(void) {

    //pack frames and roles    
    map<unsigned int, unsigned int> tmp_F, tmp_R;
    map<unsigned int, double >tmp_gamma;
    unsigned int f=0, r=0;
    tmp_gamma[0] = gamma[0];
    for (set<unsigned int>::const_iterator it=used_frames.begin();
            it!= used_frames.end(); ++it) {
        tmp_F[*it] = ++f;
    }
    for (set<unsigned int>::const_iterator it=used_roles.begin();
            it!= used_roles.end(); ++it) {
        tmp_R[*it] = ++r;
        tmp_gamma[r] = gamma[*it];
    }
    for(unsigned int u=1; u<=U; ++u) {
        for (unsigned int t = 1; t <= w[u-1].size(); ++t) {
            frames[u-1][t-1] = tmp_F[frames[u-1][t-1]];
        }
    }
    vector<vector<unsigned int> > tmp_roles;
    fc_f.clear();
    fc_fsw.clear();
    for (unsigned int f=1; f<=used_frames.size(); ++f) {
        tmp_roles.push_back(vector<unsigned int>(S, 0));
        fc_f.push_back(0);
        fc_fsw.push_back(vector<vector<unsigned int> >(S,vector<unsigned int>(V, 0)));
    }

    map<unsigned int, double> tmp_tau;
    tmp_tau[0] = F;
    frameSet.clear();
    for (set<unsigned int>::const_iterator it = used_frames.begin(); it != used_frames.end(); ++it) {
        tmp_tau[tmp_F[*it-1]] = tau[*it];
        for (unsigned int s=1; s<=S; ++s) {
            tmp_roles[tmp_F[*it]-1][s-1] = tmp_R[roles[*it-1][s-1]];
        }
        frameSet.insert(frameSet.makeKey(tmp_roles[tmp_F[*it]-1]));
    }
    tau = tmp_tau;
    gamma = tmp_gamma;
    roles = tmp_roles;
    
    for(unsigned int u=1; u<=U; ++u) {
        for (unsigned int t = 1; t <= w[u-1].size(); ++t) {
            fc_f[frames[u-1][t-1]-1]++;
            for (unsigned int s=1; s<=S; ++s) {
                if (w[u-1][t-1][s-1] != 0) {
                    fc_fsw[frames[u-1][t-1]-1][s-1][w[u-1][t-1][s-1]-1]++;
                }
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

Sampler_t::~Sampler_t() {
    if (initialized) {
        delete dist;
    }
}
