#include <iostream>
#include <boost/program_options.hpp>

#include "sampler.h"

namespace po = boost::program_options;

void outputUsage(const po::options_description& desc, char* prog )
{
    
   std::cout << "Usage:\n  " << prog << " input_file output_directory\n\n" << desc;
}


int main(int argc, char **argv) {

    
    string inputFileName, outputDirectoryName;
    string framesFileName = "frames.sampled";
    string rolesFileName = "roles.sampled";

    unsigned int frames = 7;
    unsigned int roles = 4;
    unsigned int iters = 100;
    float alpha = 0.1;
    float beta = 0.1;
    
    po::options_description desc("Allowed options");

    desc.add_options()
        ("help", "produce help message")
        ("input-file", po::value<string>(), "input file name")
        ("output-directory", po::value<string>(),"output directory" )
        ("frames-file", po::value<string>(), "frames output file name")
        ("roles-file", po::value<string>(), "roles output file name")
        ("frames,F", po::value<unsigned int>(), "sampled frames")
        ("roles,R", po::value<unsigned int>(), "sampled roles")
        ("iters,I", po::value<unsigned int>(), "iterations")
        ("alpha", po::value<float>(), "alpha")
        ("beta", po::value<float>(), "beta")
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
    } else {
        outputUsage(desc, argv[0]);
        return 2;
    }
    
    if (vm.count("frames-file")) {
        framesFileName = vm["frames-file"].as<string>();
    }

    if (vm.count("roles-file")) {
        rolesFileName = vm["roles-file"].as<string>();
    }
    
    if (vm.count("frames")) {
        frames = vm["frames"].as<unsigned int>();
    }
    
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
    

    Sampler_t sampler(frames, roles, alpha, beta);
    
    cout << "Loading data..." << endl;
    if (!sampler.loadData(inputFileName)) return 3;

    cout << "Initializing..." << endl;
    sampler.initialize();
    cout << "Sampling..." << endl;
    sampler.sampleAll(outputDirectoryName, iters);
    
    sampler.printFrames();
    cout << endl << endl;
    sampler.printRoles();
    cout << endl;
    
    //cout << "Selecting the best sample..." << endl;
    //sampler.dumpBest(outputDirectoryName);

}
