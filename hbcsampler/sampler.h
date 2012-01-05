#ifndef _SAMPLER
#define _SAMPLER

#include <vector>
#include <string>
#include "frames.h"


using namespace std;

class Sampler_t {
public:
    Sampler_t(unsigned int _F,
              unsigned int _R,
              float _alpha,
              float _beta): F(_F), R(_R), S(0), U(0), V(0),
                            alpha(_alpha), beta(_beta),
                            initialized(false)
                            {};

    bool loadData(string inputFileName);

    bool initialize(void);
    
    void sample(void);

    bool dump(string prefix);

    bool dumpBest(string outputDir);
    
    bool sampleAll(string outputDir, unsigned int iters);

    void printFrames(void);

    void printRoles(void);

    ~Sampler_t();

private:
    unsigned int F;
    unsigned int R;
    unsigned int S;
    unsigned int U;
    unsigned int V;
    double alpha;
    double beta;

    unsigned int** frames;
    unsigned int** roles;

    double** post_phi;
    double** post_theta;

    Frames_t *frameSet;
  
    vector<vector<vector<unsigned int> > > w;//inputData;
    vector<unsigned int> fc_f;
    vector<vector<vector<unsigned int> > > fc_fsw;

    //sampling
    void resample_post_phi(void);
    void resample_post_theta(void);
    void resample_frames(void);
    void resample_roles(void);


    //inititalization
    bool initialized;
    void initialize_frames(void);
    void initialize_roles(void);
    void initialize_post_phi(void);
    void initialize_post_theta(void);


};


#endif
