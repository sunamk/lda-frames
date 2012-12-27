/*
 * Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
 * Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
 *
 */

#ifndef _SAMPLER
#define _SAMPLER

#include <vector>
#include <set>
#include <map>
#include <string>
#include "antoniak.h"
#include "frames.h"


using namespace std;

class Sampler_t {
public:
    Sampler_t(unsigned int _F,
              unsigned int _R,
              float _alpha,
              float _beta,
              float _gamma,
              float _delta): F(_F), R(_R), S(0), U(0), V(0),
                            alpha(_alpha), beta0(_beta), gamma(_gamma),
                            delta(_delta),
                            infinite_F(false), infinite_R(false),
                            startIter(1), minHyperIter(50), initialized(false)
                            {};

    bool loadData(string inputFileName);

    bool initialize(bool recovery);
    
    void sample(void);

    bool dump(string prefix);

    //bool dumpBest(string outputDir);
    
    bool sampleAll(string outputDir, unsigned int iters, bool allSamples);

    void printFrames(void);

    void printRoles(void);

    bool writeLog(string outputDir, unsigned int citer, unsigned int aiter);

    bool recoverParameters(string logDir);

    bool recoverData(string dataDir);
    
    //infinite LDA-frames
    unsigned int createNewFrame(vector<unsigned int> &frame);
    unsigned int createNewRole(void);

    ~Sampler_t();

private:
    unsigned int F;
    unsigned int R;
    unsigned int S;
    unsigned int U;
    unsigned int V;
    double alpha;
    double beta0;
    vector<vector<double> > beta;
    double gamma;
    double delta;

    unsigned int** frames;
    unsigned int** roles;

    double** post_phi;
    vector<vector<double> > post_theta;
    double* post_omega;

    bool infinite_F;
    bool infinite_R;

    Frames_t *frameSet;

    string inputFile;
    unsigned int startIter;
    unsigned int minHyperIter;
  
    vector<vector<vector<unsigned int> > > w;//inputData;
    vector<unsigned int> fc_f;
    vector<vector<vector<unsigned int> > > fc_fsw;

    //infinite LDA-frames
    set<unsigned int> unused_frames, used_frames, unused_roles, used_roles;
    map<unsigned int, double> tau;
    Antoniak_t antoniak;
    unsigned int tables;
    void pack_FR(void);

    //sampling
    void resample_post_phi(void);
    void resample_post_theta(void);
    void resample_post_omega(void);
    void resample_frames(void);
    void resample_frames_inf(void);
    void resample_roles(void);
    void resample_roles_inf(void);
    void resample_tau(void);
    void resample_hypers(unsigned int iters);
    //void resample_beta(unsigned int iters);
    bool sample_new_frame(vector<unsigned int> &frame, vector<unsigned int> &pos);


    //inititalization
    bool initialized;
    void initialize_frames(void);
    void initialize_roles(void);
    void initialize_beta(void);
    void initialize_post_phi(void);
    void initialize_post_theta(void);
    void initialize_post_omega(void);
    void initialize_infinite_vars(void);

    double perplexity(void);
};


#endif
