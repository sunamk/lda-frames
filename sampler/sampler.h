/*
 * Copyright (C) 2014 Jiri Materna <xmaterna@fi.muni.cz>
 * Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
 *
 */

#ifndef _SAMPLER
#define _SAMPLER

#include <vector>
#include <set>
#include <map>
#include <string>
#ifdef _OPENMP
    #include <omp.h>
#endif
#include "distributions.h"
#include "frames.h"

using namespace std;

class Sampler_t {
public:
    Sampler_t(unsigned int _F,
              unsigned int _R,
              float _alpha,
              float _beta,
              float _gamma,
              float _delta,
              long int _seed,
              bool _reestimate_F,
              bool _reestimate_R,
              unsigned int _cores,
              bool _testPhase): F(_F), R(_R), S(0), U(0), V(0),
                            alpha0(_alpha), beta0(_beta), gamma0(_gamma),
                            delta(_delta),
                            seed(_seed),
                            reestimate_F(_reestimate_F), reestimate_R(_reestimate_R),
                            cores(_cores), testPhase(_testPhase), infinite_F(false), infinite_R(false),
                            startIter(1), initialized(false)
                            {};

    ~Sampler_t();

    bool loadData(string inputFileName);
    
    bool loadTestData(string inputFileName);

    bool initialize(bool recovery);
    
    void sample(void);

    bool dump(string prefix);

    bool sampleAll(string outputDir, unsigned int iters, unsigned int burn_in, bool allSamples, 
            bool no_hypers, bool no_perplexity, bool rm);

    void printFrames(void);

    void printRoles(void);
    
    void printTest(void);

    bool writeLog(string outputDir, unsigned int citer, unsigned int aiter);

    bool recoverParameters(string logDir);

    bool recoverData(string dataDir, unsigned int burn_in);
    
    double perplexity(bool test);
    
    double bestPerplexity;
    unsigned int requiredIters;
    

private:
    unsigned int F;
    unsigned int R;
    unsigned int S;
    unsigned int U;
    unsigned int V;
    double alpha0;
    vector<double> alpha;
    double beta0;
    vector<double> beta;
    double gamma0;
    map<unsigned int, double> gamma;
    double delta;
    long int seed;
    Distributions_t *dist;

    vector<vector<unsigned int> > frames, test_frames;
    vector<vector<unsigned int> > roles;

    vector<vector<double> > post_phi;
    vector<vector<double> > post_theta;
    vector<double> post_omega;

    bool reestimate_F;
    bool reestimate_R;
    unsigned int cores;
    bool testPhase;
    bool infinite_F;
    bool infinite_R;

    Frames_t frameSet;

    string inputFile;
    unsigned int startIter;
  
    vector<vector<vector<unsigned int> > > w, test_w;//inputData;
    vector<unsigned int> fc_f;
    vector<vector<vector<unsigned int> > > fc_fsw;

    set<unsigned int> unused_frames, used_frames, unused_roles, used_roles;
    map<unsigned int, unsigned int> frame_count;
    map<unsigned int, double> tau;
    unsigned int tables;

    bool emptyFrames;
    bool checkPattern(vector<unsigned int> &u, vector<unsigned int> &v);
    map<vector<unsigned int>, unsigned int> framePatterns;
    unsigned int positions;

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
    bool sample_new_frame(vector<unsigned int> &frame, vector<unsigned int> &pos);
    void predict_test(void);

    //inititalization
    bool initialized;
    void initialize_frames(void);
    void initialize_roles(void);
    void initialize_beta(void);
    void initialize_post_phi(void);
    void initialize_post_theta(void);
    void initialize_post_omega(void);
    void initialize_infinite_vars(void);
    
    //infinite LDA-frames
    void pack_FR(void);
    unsigned int createNewFrame(vector<unsigned int> &frame);
    unsigned int createNewRole(void);

};

#endif
