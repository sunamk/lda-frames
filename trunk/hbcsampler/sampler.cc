#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/errno.h>
#include <dirent.h>
#include <boost/tokenizer.hpp>

#include "sampler.h"
#include "samplib.h"
#include "stats.h"

using namespace std;

void Sampler_t::resample_post_phi(void) {
    double* tmp = (double*) malloc(sizeof(double) * (F + 1));
    for (unsigned int u = 1; u <= U; ++u) {
        for (unsigned int i = 1; i<=F; ++i) {
            tmp[i-1] = 0.0;
        }
        tmp[F] = 0.0;
        for (unsigned int t = 1; t <= w[u-1].size(); ++t) {
            tmp[F] += 1.0;
            tmp[frames[u-1][t-1] - 1] += 1.0;
        }
        sample_Delta(post_phi[u-1], tmp, F);
    }
    free(tmp);
}


void Sampler_t::resample_post_theta(void) {
    double* tmp = (double*) malloc(sizeof(double) * (V + 1));
    for (unsigned int r = 1; r <= R; ++r) {
        for (unsigned int i = 1; i <= V; ++i) {
            tmp[i-1] = 0.0;
        }
        tmp[V] = 0.0;
        for (unsigned int u = 1; u <= U; ++u) {
            for (unsigned int t = 1; t <= w[u - 1].size(); ++t) {
                for (unsigned int s = 1; s <= S; ++s) {
                    tmp[V] += (r == roles[frames[u-1][t-1]-1][s-1]) ? 1 : 0;
                    tmp[w[u-1][t-1][s-1]-1] += 
                        (r == roles[frames[u-1][t-1]-1][s-1]) ? 1 : 0;
                }
            }
        }
        sample_Delta(post_theta[r-1], tmp, V);
    }
    free(tmp);
}



void Sampler_t::resample_frames(void) {
    

    #pragma omp parallel for schedule(dynamic)
    for (int u = 1; u <= (signed int) U; ++u) {
        double* post_frames = (double*) malloc(sizeof(double) * (F + 1));
        for (unsigned int t = 1; t <= w[u-1].size(); ++t) {

            for (unsigned int s = 1; s <= S; ++s) {
                #pragma omp critical
                post_theta[roles[frames[u-1][t-1]-1][s-1]-1][V]--;
                #pragma omp critical
                post_theta[roles[frames[u-1][t-1]-1][s-1]-1][w[u-1][t-1][s-1]-1]--;
            }
            post_phi[u-1][F] -= 1.0;
            post_phi[u-1][frames[u-1][t-1]-1] -= 1.0;
            post_frames[F] = 0;

            for (unsigned int f = 1; f <= F; ++f) {
                double tmp = 0.0;
                post_frames[f-1] = 0.0;
                for (unsigned int s = 1; s <= S; ++s) {
                    tmp += ldf_Mult_smooth(0, beta, w[u-1][t-1][s-1], post_theta[roles[f-1][s-1]-1], 1, V);
                }
                post_frames[f-1] = tmp + ldf_Mult_smooth(0, alpha, f, post_phi[u-1], 1, F);
            }

            normalizeLog(post_frames, 1, F);
            for (unsigned int s=1; s<=S; ++s) {
                #pragma omp critical
                fc_fsw[frames[u-1][t-1]-1][s-1][w[u-1][t-1][s-1]-1]--;
            }
            #pragma omp critical
            fc_f[frames[u-1][t-1]-1]--;

            frames[u-1][t-1] = sample_Mult(post_frames, 1, F);

            for (unsigned int s=1; s<=S; ++s) {
                #pragma omp critical
                fc_fsw[frames[u-1][t-1]-1][s-1][w[u-1][t-1][s-1]-1]++;
            }
            #pragma omp critical
            fc_f[frames[u-1][t-1]-1]++;


            post_phi[u-1][F] += 1.0;
            post_phi[u-1][frames[u-1][t-1]-1] += 1.0;
            for (unsigned int s = 1; s<=S; ++s) {
                #pragma omp critical
                post_theta[roles[frames[u-1][t-1]-1][s-1]-1][V]++;
                #pragma omp critical
                post_theta[roles[frames[u-1][t-1]-1][s-1]-1][w[u-1][t-1][s-1]-1]++;
            }
        }
        free(post_frames);
    }
}

void Sampler_t::resample_roles(void) {

    double* vec = (double*) malloc(sizeof(double) * (R + 1));

    for (unsigned int r = 1; r <= R; ++r) {
        vec[r-1] = 1.0 / R;
    }
    vec[R] = 1.0;


    #pragma omp parallel for schedule(dynamic)
    for (int f = 1; f <= (signed int) F; ++f) {
        double* post_roles = (double*) malloc(sizeof(double) * (R + 1));
        for (unsigned int s = 1; s <= S; ++s) {

            #pragma omp critical
            post_theta[roles[f-1][s-1]-1][V] -= fc_f[f-1];
            for(unsigned int v=1; v<=V; ++v) {
                #pragma omp critical
                post_theta[roles[f-1][s-1]-1][v-1] -= fc_fsw[f-1][s-1][v-1];
            }

            post_roles[R] = 0.0;

            for (unsigned int r = 1; r <= R; ++r) {
                double prod = 0.0;
                post_roles[r-1] = 0.0;
                for (unsigned int u = 1; u <= U; ++u) {
                    for (unsigned int t=1; t <= w[u-1].size(); ++t) {
                        if ((unsigned int) f == frames[u-1][t-1]) prod += ldf_Mult_smooth(0, beta, w[u-1][t-1][s-1], post_theta[r-1], 1, V);
                    }
                }
                
                for (unsigned int v = 1; v<=V; ++v) {
                    for (unsigned int i = 0; i < fc_fsw[f-1][s-1][v-1]; ++i) {
                        prod += ldf_Mult_smooth(0, beta, v, post_theta[r-1], 1, V);
                    }
                }
                post_roles[r-1] = prod + ldf_Mult(0, r, vec, 1, R);
            }

            normalizeLog(post_roles, 1, R);
            roles[f-1][s-1] = sample_Mult(post_roles, 1, R);

            #pragma omp critical
            post_theta[roles[f-1][s-1]-1][V] += fc_f[f-1];
            for(unsigned int v=1; v<=V; ++v) {
                #pragma omp critical
                post_theta[roles[f-1][s-1]-1][v-1] += fc_fsw[f-1][s-1][v-1];
            }


        }
        free(post_roles);
    }
    free(vec);
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
        frames[u-1][(w[u-1].size())] = 0;
    }
}

void Sampler_t::initialize_roles(void) {
    for (unsigned int f=1; f<=F; ++f) {
        for (unsigned int s=1; s<=S; ++s) {
            roles[f-1][s-1] = sample_MultSym(1, R);
        }
        roles[f-1][S] = 0;
    }
}

void Sampler_t::initialize_post_phi(void) {
    for (unsigned int u=1; u<=U; ++u) {
        for (unsigned int f=1; f<=F; ++f) {
            post_phi[u-1][f-1] = 0.0;
        }
        post_phi[u-1][F] = 0.0;
    }
    resample_post_phi();
}



void Sampler_t::initialize_post_theta(void) {
    for (unsigned int r=1; r<=R; ++r) {
        for (unsigned int v=1; v<=V; ++v) {
            post_theta[r-1][v-1] = 0.0;
        }
        post_theta[r-1][V] = 0.0;
    }
    resample_post_theta();
}



bool Sampler_t::loadData(string inputFileName) {

    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

    ifstream ifs(inputFileName.c_str());
    if (!ifs.is_open()) {
        cout << "Cannot open file " << inputFileName << ".\n";
        return false;
    }

    string line;
    while (getline(ifs, line)) {
        boost::char_separator<char> sep("\t");
        tokenizer tokens(line, sep);
        vector<vector<unsigned int> > unit;
        for (tokenizer::iterator tok_iter = tokens.begin();
            tok_iter != tokens.end(); ++tok_iter) {

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
                    cout << "Inconsitent number of slots.\n";
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

    if (F > pow(R, S)) {
        cout << "Number of frames (F) must be lower than or equal to the number of all " <<
                "combinations of possible semantic roles (R^S)." << endl;
        return false;
    }

    cout << "F = " << F << endl;
    cout << "R = " << R << endl;
    cout << "alpha = " << alpha << endl;
    cout << "beta = " << beta << endl;
    cout << "Lexical units = " << U << endl;
    cout << "Slots = " << S << endl;
    cout << "Vocabulary size = " << V << endl;

    return true;
}

bool Sampler_t::initialize(void) {
    setall(time(0),time(0));   /* initialize random number generator */
    cout << "Allocating memory..." << endl;

    frames = (unsigned int**) malloc(sizeof(unsigned int*) * (1+(U)-(1)));
    for (unsigned int malloc_dim_1=1; malloc_dim_1<=U; malloc_dim_1++) {
        frames[malloc_dim_1-1] = (unsigned int*) malloc(sizeof(unsigned int) * (1+((w[malloc_dim_1-1].size()) + (1))-(1)));
    }

    post_phi = (double**) malloc(sizeof(double*) * (1+(U)-(1)));
    for (unsigned int malloc_dim_1=1; malloc_dim_1<=U; malloc_dim_1++) {
        post_phi[malloc_dim_1-1] = (double*) malloc(sizeof(double) * (1+((F) + (1))-(1)));
    }

    post_theta = (double**) malloc(sizeof(double*) * (1+(R)-(1)));
    for (unsigned int malloc_dim_1=1; malloc_dim_1<=R; malloc_dim_1++) {
        post_theta[malloc_dim_1-1] = (double*) malloc(sizeof(double) * (1+((V) + (1))-(1)));
    }    
    
    roles = (unsigned int**) malloc(sizeof(unsigned int*) * (1+(F)-(1)));
    for (unsigned int malloc_dim_1=1; malloc_dim_1<=F; malloc_dim_1++) {
        roles[malloc_dim_1-1] = (unsigned int*) malloc(sizeof(unsigned int) * (1+((S) + (1))-(1)));
        fc_f.push_back(0);
        fc_fsw.push_back(vector<vector<unsigned int> >(S,vector<unsigned int>(V, 0)));
    }    

    cout << "Initializing variables..." << endl;
    initialize_frames();
    initialize_roles();
    initialize_post_phi();
    initialize_post_theta();
    
    initialized = true;
    return true;    

}

Sampler_t::~Sampler_t() {
    
    if (initialized) {

        for (unsigned int malloc_dim_1=1; malloc_dim_1<=F; malloc_dim_1++) {
            free(roles[malloc_dim_1-1]);
        }
        free(roles);

        for (unsigned int malloc_dim_1=1; malloc_dim_1<=R; malloc_dim_1++) {
            free(post_theta[malloc_dim_1-1]);
        }
        free(post_theta);

        for (unsigned int malloc_dim_1=1; malloc_dim_1<=U; malloc_dim_1++) {
            free(post_phi[malloc_dim_1-1]);
        }
        free(post_phi);

        for (unsigned int malloc_dim_1=1; malloc_dim_1<=U; malloc_dim_1++) {
            free(frames[malloc_dim_1-1]);
        }
        free(frames);
    }
}


void Sampler_t::sample(void) {
    resample_frames();
    resample_roles();
}

bool Sampler_t::sampleAll(string outputDir, unsigned int iters) {
    if (outputDir.at(outputDir.size()-1) != '/') outputDir += "/";
    int status = mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
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


    for (unsigned int i = 0; i < iters; ++i) {
        cout << "Iteration no. " << i + 1 << ".\n";
        sample();
        stringstream ss;
        ss << i + 1;
        if (!dump(outputDir + ss.str())) {
            return false;
        }
    }
    return true;
}

bool Sampler_t::dump(string prefix) {

    string ffn = prefix + "-frames.smpl";
    string rfn = prefix + "-roles.smpl";

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

bool Sampler_t::dumpBest(string outputDir) {
    if (outputDir.at(outputDir.size()-1) != '/') outputDir += "/";
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

