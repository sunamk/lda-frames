/*
 * Copyright (C) 2014 Jiri Materna <xmaterna@fi.muni.cz>
 * Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
 *
 */

#include <iostream>
#include <boost/program_options.hpp>
#include "sampler.h"

namespace po = boost::program_options;

void outputUsage(const po::options_description& desc, char* prog )
{

   std::cout << "Usage:\n  " << prog << "input_file working_directory\n\n" << desc;
}

int main(int argc, char **argv) {

    string inputFileName, outputDirectoryName, testFileName;

    unsigned int frames = 0;
    unsigned int roles = 0;
    unsigned int iters = 1000;
    float alpha = 0.1;
    float beta = 0.1;
    float gamma = 0.1;
    float delta = 1.0;
    bool testPhase = false;
    double train_perplexity, test_perplexity;

    po::options_description desc("Allowed options");

    desc.add_options()
        ("help,h", "Print this help message.")
        ("input-file", po::value<string>(), "Path to the input file.")
        ("test-file,T", po::value<string>(), "Path to a test file. If the file is not provided, testing is skipped.")
        ("working-directory", po::value<string>(),"Path to the working directory." )
    ;

    po::positional_options_description p;
    p.add("input-file", 1);
    p.add("working-directory", 1);

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

    if (vm.count("working-directory"))
    {
        outputDirectoryName = vm["working-directory"].as<string>();
        if (outputDirectoryName.at(outputDirectoryName.size()-1) != '/')
            outputDirectoryName += "/";
    } else {
        outputUsage(desc, argv[0]);
        return 2;
    }

    Sampler_t sampler(frames, roles, alpha, beta, gamma, delta, 0, false, false, 1, false);  

    cout << "Reading parameters from log..." << endl;
    if (!sampler.recoverParameters(outputDirectoryName)) return 3;
    iters = sampler.requiredIters;
    cout << "Required number of iterations is " << iters << "." << endl;
    cout << "Loading input data..." << endl;
    if (!sampler.loadData(inputFileName)) return 3;
    if (testPhase) {
        cout << "Loading test data..." << endl;
        if (!sampler.loadTestData(testFileName)) return 3;
    }
    sampler.initialize(true);
    cout << "Recovering sampled data..." << endl;
    if (!sampler.recoverData(outputDirectoryName, 0, true)) return 3;
    //cout << "Dumping hyperparameters: " << endl;
    //sampler.dumpHypers(); 
    cout << "Computing train perplexity... ";
    train_perplexity = sampler.perplexity(false);
    cout << train_perplexity << endl;
    if (testPhase) {
        cout << "Computing test perplexity... ";
        test_perplexity = sampler.perplexity(true);
        cout << test_perplexity << endl;
    }
    
    

    return 0;

}




