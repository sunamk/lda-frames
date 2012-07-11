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
#include "samplib.h"
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
            post_omega[roles[f-1][s-1]-1]++;
            post_omega[R]++;
        }
    }
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
                double tmp = 0;
                post_frames[f-1] = 0;
                for (unsigned int s = 1; s <= S; ++s) {
                    tmp += ldf_Mult_smooth(1, beta, w[u-1][t-1][s-1],
                            post_theta[roles[f-1][s-1]-1], 1, V);
                }
                post_frames[f-1] = tmp + ldf_Mult_smooth(0, alpha, f, post_phi[u-1], 1, F);
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
                double tmp = 0;
                post_frames[f-1] = 0;
                for (unsigned int s = 1; s <= S; ++s) {
                    tmp += ldf_Mult_smooth(1, beta, w[u-1][t-1][s-1],
                            post_theta[roles[f-1][s-1]-1], 1, V);
                }
                post_frames[f-1] = tmp + ldf_Mult_smooth(0, alpha, f, post_phi[u-1], 1, F);
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

void Sampler_t::resample_roles(void) {

    for (int f = 1; f <= (signed int) F; ++f) {
        double* post_roles = (double*) malloc(sizeof(double) * (R + 1));
        

        for (unsigned int s = 1; s <= S; ++s) {
            post_theta[roles[f-1][s-1]-1][V] -= fc_f[f-1];
            for(unsigned int v=1; v<=V; ++v) {
                post_theta[roles[f-1][s-1]-1][v-1] -= fc_fsw[f-1][s-1][v-1];
            }

            post_omega[roles[f-1][s-1]-1]--;
            post_omega[R]--;

            post_roles[R] = 0.0;
           
            FrameKey_t oldFrame; 
            oldFrame = frameSet->makeKey(roles[f-1]);
            
            for (unsigned int r = 1; r <= R; ++r) {
                double prod = 0.0;
                post_roles[r-1] = -1 * numeric_limits<double>::max(); //zero probability
                FrameKey_t newFrame;
                bool inside;
                {
                    newFrame = frameSet->makeKey(roles[f-1], s, r);
                    inside = frameSet->inside(newFrame);
                }

                if (newFrame == oldFrame || !inside) {
                    for (unsigned int v = 1; v<=V; ++v) {
                        prod += fc_fsw[f-1][s-1][v-1]*ldf_Mult_smooth(1, beta, v, 
                            post_theta[r-1], 1, V);
                    }
                    post_roles[r-1] = prod + ldf_Mult_smooth(0, gamma, r, post_omega, 
                        1, R);
                }
            }
            
            normalizeLog(post_roles, 1, R);

            {
                roles[f-1][s-1] = sample_Mult(post_roles, 1, R);
                frameSet->remove(oldFrame);
                frameSet->insert(frameSet->makeKey(roles[f-1]));
            }

            post_omega[roles[f-1][s-1]-1]++;
            post_omega[R]++;
            post_theta[roles[f-1][s-1]-1][V] += fc_f[f-1];
            for(unsigned int v=1; v<=V; ++v) {
                post_theta[roles[f-1][s-1]-1][v-1] += fc_fsw[f-1][s-1][v-1];
            }
        }
        free(post_roles);
    }
}

void Sampler_t::resample_roles_inf(void) {

    for (int f = 1; f <= (signed int) F; ++f) {
        double* post_roles = (double*) malloc(sizeof(double) * (R + 2));
        
        for (unsigned int s = 1; s <= S; ++s) {
            post_theta[roles[f-1][s-1]-1][V] -= fc_f[f-1];
            for(unsigned int v=1; v<=V; ++v) {
                post_theta[roles[f-1][s-1]-1][v-1] -= fc_fsw[f-1][s-1][v-1];
            }

            post_omega[roles[f-1][s-1]-1]--;
            post_omega[R]--;

            post_roles[R+1] = 0;
           
            FrameKey_t oldFrame; 
            oldFrame = frameSet->makeKey(roles[f-1]);
            
            for (unsigned int r = 1; r <= R; ++r) {
                double prod = 0.0;
                post_roles[r-1] = -1 * numeric_limits<double>::max(); //zero probability
                FrameKey_t newFrame;
                bool inside;
                newFrame = frameSet->makeKey(roles[f-1], s, r);
                inside = frameSet->inside(newFrame);
                set<unsigned int>::iterator it = used_roles.find(r);

                //probability of used roles
                if ((newFrame == oldFrame || !inside) && it != used_roles.end()) {
                    for (unsigned int v = 1; v<=V; ++v) {
                        prod += fc_fsw[f-1][s-1][v-1]*ldf_Mult_smooth(1, beta, v, 
                            post_theta[r-1], 1, V);
                    }
                    post_roles[r-1] = prod + ldf_Mult_smooth(0, gamma, r, post_omega, 
                        1, R);
                }
            }

            //probability of a new role
            double prod = 0.0;
            for (unsigned int v = 1; v<=V; ++v) {
                prod -= fc_fsw[f-1][s-1][v-1]*log(V);
            }
            post_roles[R] = log(gamma) + prod;
            
            normalizeLog(post_roles, 1, R + 1);
            unsigned int newRole = sample_Mult(post_roles, 1, R + 1);
            //cout << "New role: " <<  newRole << endl;

            //not necessarily
            set<unsigned int>::iterator it = used_roles.find(newRole);
            if (it == used_roles.end() && newRole != (R + 1)) {
                cerr << "Sampled disallowed role." << endl;
                exit(10);
            }
            
            //free unused role numbers
            if (roles[f-1][s-1] != newRole && post_omega[roles[f-1][s-1]-1] == 0) {
                unused_roles.insert(roles[f-1][s-1]);
                used_roles.erase(roles[f-1][s-1]);
            }
            //create a new role if required
            if (newRole == R + 1) {
                if (!createNewRole(roles[f-1][s-1])) {
                    cerr << "Error in creating new role." << endl;
                    exit(11);
                }
                post_roles = (double*) realloc(post_roles, sizeof(double) * (R + 2));
            } else {
                roles[f-1][s-1] = newRole;
            }

            frameSet->remove(oldFrame);
            frameSet->insert(frameSet->makeKey(roles[f-1]));

            post_omega[roles[f-1][s-1]-1]++;
            post_omega[R]++;
            post_theta[roles[f-1][s-1]-1][V] += fc_f[f-1];
            for(unsigned int v=1; v<=V; ++v) {
                post_theta[roles[f-1][s-1]-1][v-1] += fc_fsw[f-1][s-1][v-1];
            }
        }
        free(post_roles);
    }
}



void Sampler_t::initialize_frames(void) {
    for (unsigned int u=1; u<=U; ++u) {
        for (unsigned int t=1; t<=w[u-1].size(); ++t) {
            frames[u-1][t-1] = sample_MultSym(1, F);
            
            fc_f[frames[u-1][t-1]-1]++;
            for (unsigned int s=1; s<=S; ++s) {
                fc_fsw[frames[u-1][t-1]-1][s-1][w[u-1][t-1][s-1]-1]++;
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
    resample_post_phi();
}



void Sampler_t::initialize_post_theta(void) {
    for (unsigned int r=1; r<=R; ++r) {
        for (unsigned int v=1; v<=V; ++v) {
            post_theta[r-1][v-1] = 0;
        }
        post_theta[r-1][V] = 0;
    }
    resample_post_theta();
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
    }
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
        cout << "R = automatic" << endl;
        infinite_R = true;
        R = ceil(log(F)/log(S)); //minimum number of roles
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
    cout << "beta = " << beta << endl;
    cout << "gamma = " << gamma << endl;
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
        //frames[u-1] = (unsigned int*) malloc(sizeof(unsigned int) * (w[u-1].size() + 1));
        frames[u-1] = (unsigned int*) malloc(sizeof(unsigned int) * w[u-1].size());
    }

    roles = (unsigned int**) malloc(sizeof(unsigned int*) * F);
    for (unsigned int f=1; f<=F; ++f) {
        //roles[f-1] = (unsigned int*) malloc(sizeof(unsigned int) * (S + 1));
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
    
    post_omega = (double*) malloc(sizeof(double) * (R +1));
    
    if (!recovery) {
        cout << "Initializing variables..." << endl;
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
    }

}


void Sampler_t::sample(void) {
    if (infinite_F)
        resample_frames_inf();
    else
        resample_frames();
    
    if (infinite_R)    
        resample_roles_inf();
    else
        resample_roles();

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
        cout << "Iteration no. " << i << " (" << used_frames.size() << " frames, " 
             << used_roles.size() << " roles).\n";
        sample();
        stringstream ss;
        if (allSamples) {
            ss << i << "-";
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

    for (unsigned int f=1; f<=F; ++f) {
        for (unsigned int s=1; s<=S; ++s) {
            rfile << roles[f-1][s-1];
            if (s != S) rfile << " ";
        }
        if (f != F) rfile << endl;
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
    lfile << "Beta:\t" << beta << endl;
    lfile << "Gamma:\t" << gamma << endl;
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
    for (unsigned int f=1; f<=F; ++f) {
        for (unsigned int s=1; s<=S; ++s) {
            cout << roles[f-1][s-1];
            if (s != S) cout << " ";
        }
        if (f != F) cout << endl;
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
            if (lineItems.at(0) == "Beta:") {
                beta = atof(lineItems.at(1).c_str());
                if (beta <= 0) {
                    cout << "Wrong beta in the log file(" << lineItems.at(1) << ")." << endl;
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

    ffile.close();
    rfile.close();

    return true;
}

bool Sampler_t::createNewRole(unsigned int &role) {
    if (!unused_roles.empty()) {
        set<unsigned int>::iterator it = unused_roles.begin();
        role = *it;
        unused_roles.erase(it);
    } else {
        R++;
        role = R;
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

    used_roles.insert(role);
    
    return true;
}


bool Sampler_t::createNewFrame(unsigned int &frame) {
    if (!infinite_R && (used_frames.size() + 1 > pow(R,S))) {
         return false; //too many frames
    }

    if (!unused_frames.empty()) {
        set<unsigned int>::iterator it = unused_frames.begin();
        frame = *it;
        unused_frames.erase(it);
    } else {
        //TODO
        return false;
    }

    used_frames.insert(frame);

    if (!infinite_R) {
        do {
            for (unsigned int s=1; s<=S; ++s) {
                roles[frame-1][s-1] = sample_MultSym(1, R);
            }
        } while (frameSet->inside(frameSet->makeKey(roles[frame-1])));
        FrameKey_t fk = frameSet->makeKey(roles[frame-1]);
        frameSet->insert(fk);
    } else {
        //TODO
        return false;
    }
    
    return true;
}
