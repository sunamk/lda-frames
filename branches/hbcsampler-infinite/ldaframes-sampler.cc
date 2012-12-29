/*
 * Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
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

    
    string inputFileName, outputDirectoryName;

    unsigned int frames = 0;
    unsigned int roles = 0;
    unsigned int iters = 1000;
    float alpha = 0.1;
    float beta = 0.1;
    float gamma = 0.1;
    float delta = 1.0;
    bool printResult = false;
    bool allSamples = false;
    long int seed = 0;
    
    po::options_description desc("Allowed options");

    desc.add_options()
        ("help", "Print this help message.")
        ("input-file", po::value<string>(), "Path to the input file.")
        ("output-directory", po::value<string>(),"Path to the output directory." )
        ("frames,F", po::value<unsigned int>(), 
            "Number of frames (if the value is not specified, it is chosen automatically).")
        ("roles,R", po::value<unsigned int>(),
            "Number of roles (if the value is not specified, it is chosen automatically).")
        ("iters,I", po::value<unsigned int>(), "Number of iterations (default is 1000).")
        ("alpha", po::value<float>(), "Alpha hyperparameter.")
        ("beta", po::value<float>(), "Beta parameter.")
        ("gamma", po::value<float>(), "Gamma parameter.")
        ("delta", po::value<float>(), "Delta parameter.")
        ("seed", po::value<long int>(), 
        "Random number generator seed (0 = current time).")
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

    if (vm.count("input-file"))
    {
        inputFileName = vm["input-file"].as<string>();
    } else {
        outputUsage(desc, argv[0]);
        return 2;
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


    Sampler_t sampler(frames, roles, alpha, beta, gamma, delta, seed);
    

    
    if (!vm.count("recovery")) {
        cout << "Loading input data..." << endl;
        if (!sampler.loadData(inputFileName)) return 3;
        cout << "Initializing..." << endl;
        sampler.initialize(false);

    } else {
        //recovery
        cout << "Recovering parameters from log..." << endl;
        if (!sampler.recoverParameters(outputDirectoryName)) return 3;
        cout << "Loading input data..." << endl;
        if (!sampler.loadData(inputFileName)) return 3;
        sampler.initialize(true);
        cout << "Recovering sampled data..." << endl;
        if (!sampler.recoverData(outputDirectoryName)) return 3;
    }

    sampler.sampleAll(outputDirectoryName, iters, allSamples);
    
    if (printResult) {
        sampler.printFrames();
        cout << endl << endl;
        sampler.printRoles();
        cout << endl;
    }

    return 0;

}
