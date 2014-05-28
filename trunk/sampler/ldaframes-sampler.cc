/*
 * Copyright (C) 2014 Jiri Materna <xmaterna@fi.muni.cz>
 * Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
 *
 */

#include <iostream>
#include <boost/program_options.hpp>
#include "distributions.h"
#include "sampler.h"


namespace po = boost::program_options;

void outputUsage(const po::options_description& desc, char* prog )
{
    
   std::cout << "Usage:\n  " << prog << " input_file output_directory\n\n" << desc;
}


int main(int argc, char **argv) {

    
    string inputFileName, outputDirectoryName, testFileName;

    unsigned int frames = 0;
    unsigned int roles = 0;
    unsigned int iters = 1000;
    unsigned int burn_in = 50;
    float alpha = 0.1;
    float beta = 0.1;
    float gamma = 0.1;
    float delta = 1.0;
    bool printResult = false;
    bool allSamples = false;
    bool reestimate_F = false;
    bool reestimate_R = false;
    bool no_hypers = false;
    bool no_perplexity = false;
    bool remove_old_samples = false;
    unsigned int cores = 1;
    bool testPhase = false;
    
    long int seed = 0;
    
    po::options_description desc("Allowed options");

    desc.add_options()
        ("help,h", "Print this help message.")
        ("input-file", po::value<string>(), "Path to the input file.")
        ("test-file,T", po::value<string>(), "Path to a test file. If the file is not provided, testing is skipped.")
        ("output-directory", po::value<string>(),"Path to the output directory." )
        ("frames,F", po::value<unsigned int>(), 
            "Number of frames (if the value is not specified, it is chosen automatically, starting from 1).")
        ("roles,R", po::value<unsigned int>(),
            "Number of roles (if the value is not specified, it is chosen automatically, starting from 1).")
        ("iters,I", po::value<unsigned int>(), "Number of iterations (default is 1000).")
        ("burn-in", po::value<unsigned int>(), "Number of burn-in iterations before estimating hyperparameters (default is 50).")
        
        ("alpha", po::value<float>(), "Alpha hyperparameter.")
        ("beta", po::value<float>(), "Beta parameter.")
        ("gamma", po::value<float>(), "Gamma parameter.")
        ("delta", po::value<float>(), "Delta parameter.")
        ("seed", po::value<long int>(), "Random number generator seed (0 = current time).")
        #ifdef _OPENMP
            ("cores,C", po::value<unsigned int>(), "number of cores (default is 1, 0 = all available cores). This feature works only when a fixed number of frames and roles is given).")
        #endif
        ("reestimate_F","Reestimate number of frames automatically.")
        ("reestimate_R","Reestimate number of roles automatically.")
        ("no_hypers", "Do not estimate hyperparameters.")
        ("no_perplexity", "Do not compute perplexity.")
        ("rm", "Remove old samples from the output directory.")
        ("all-samples,A", "Save samples of all iterations.")
        ("recovery", "Try to recover data and continue sampling.")
        ("print,P", "Print sampled frames and roles to stdout.")
    ;

    po::positional_options_description p;
    p.add("input-file", 1);
    p.add("output-directory", 1);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).
            options(desc).positional(p).run(), vm);
        po::notify(vm);
    } catch (std::exception& e) {
        cout << e.what() << endl << endl;
        outputUsage(desc, argv[0]);
        return 2;
    }

    if (vm.count("help")) {
        outputUsage(desc, argv[0]);
        return 1;
    }

    if (vm.count("print")) {
        printResult = true;
    }
    
    if (vm.count("all-samples")) {
        allSamples = true;
    }
    
    if (vm.count("reestimate_F")) {
        reestimate_F = true;
    }
    
    if (vm.count("reestimate_R")) {
        reestimate_R = true;
    }
    
    if (vm.count("no_hypers")) {
        no_hypers = true;
    }
    
    if (vm.count("no_perplexity")) {
        no_perplexity = true;
    }
    
    if (vm.count("rm")) {
        remove_old_samples = true;
    }

    if (vm.count("input-file"))
    {
        inputFileName = vm["input-file"].as<string>();
    } else {
        outputUsage(desc, argv[0]);
        return 2;
    }

    if (vm.count("test-file"))
    {
        testFileName = vm["test-file"].as<string>();
        testPhase = true;
    } 
    
    if (vm.count("output-directory"))
    {
        outputDirectoryName = vm["output-directory"].as<string>();
        if (outputDirectoryName.at(outputDirectoryName.size()-1) != '/') 
            outputDirectoryName += "/";
    } else {
        outputUsage(desc, argv[0]);
        return 2;
    }
    
    if (vm.count("frames")) {
        frames = vm["frames"].as<unsigned int>();
    }
    if (frames == 0) alpha = 5.0;
    
    if (vm.count("roles")) {
        roles = vm["roles"].as<unsigned int>();
    }
    
    if (vm.count("iters")) {
        iters = vm["iters"].as<unsigned int>();
    }
    if (vm.count("burn-in")) {
        burn_in = vm["burn-in"].as<unsigned int>();
    }
    
    if (vm.count("alpha")) {
        alpha = vm["alpha"].as<float>();
    }
    
    if (vm.count("beta")) {
        beta = vm["beta"].as<float>();
    }
    if (vm.count("gamma")) {
        gamma = vm["gamma"].as<float>();
    }
    if (vm.count("delta")) {
        delta = vm["delta"].as<float>();
    }
    if (vm.count("seed")) {
        seed = vm["seed"].as<long int>();
    }
    #ifdef _OPENMP 
        if (vm.count("cores")) {
            cores = vm["cores"].as<unsigned int>();
        }
    #endif


    Sampler_t sampler(frames, roles, alpha, beta, gamma, delta, seed, reestimate_F, reestimate_R, cores, testPhase);

    if (!vm.count("recovery")) {
        cout << "Number of iterations is " << iters << "." << endl;
        cout << "Loading input data..." << endl;
        if (!sampler.loadData(inputFileName)) return 3;

        if (testPhase) {
            cout << "Loading test data..." << endl;
            if (!sampler.loadTestData(testFileName)) return 3;
        }
        cout << "Initializing..." << endl;
        sampler.initialize(false);

    } else {
        //recovery
        cout << "Recovering parameters from log..." << endl;
        if (!sampler.recoverParameters(outputDirectoryName)) return 3;
        if (!vm.count("iters") && sampler.requiredIters!=0) {
            iters = sampler.requiredIters;
        }
        cout << "Required number of iterations is " << iters << "." << endl;
        cout << "Loading input data..." << endl;
        if (!sampler.loadData(inputFileName)) return 3;
        if (testPhase) {
            cout << "Loading test data..." << endl;
            if (!sampler.loadTestData(testFileName)) return 3;
        }
        sampler.initialize(true);
        cout << "Recovering sampled data..." << endl;
        if (!sampler.recoverData(outputDirectoryName, burn_in)) return 3;
        if (!no_perplexity) {
            cout << "Computing perplexity..." << endl;
            sampler.bestPerplexity = sampler.perplexity(false);
        }
    }

    sampler.sampleAll(outputDirectoryName, iters, burn_in, allSamples, no_hypers, 
                      no_perplexity, remove_old_samples);
    
    if (printResult) {
        sampler.printFrames();
        cout << endl << endl;
        if (testPhase) {
            sampler.printTest();
            cout << endl << endl;
        }
        sampler.printRoles();
        cout << endl;
    }

    return 0;

}
