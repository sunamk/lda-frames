/*
 * Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
 * Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
 *
 */

#include <fstream>
#include <boost/tokenizer.hpp>
#include <cmath>
#include "sampler.h"

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
    emptyFrames = false;
    positions = 0;
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

            vector<unsigned int> s, framePattern;
            for (tokenizer::iterator slot_iter = slots.begin();
                slot_iter != slots.end(); ++slot_iter) {
                const unsigned int w = atoi(slot_iter->c_str());
                if (w == 0) {
                    framePattern.push_back(0);
                    emptyFrames = true;
                    //cout << "Invalid input data: " << *slot_iter << endl;
                    //ifs.close();
                    ///return false;
                } else {
                    framePattern.push_back(1);
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
            map<vector<unsigned int>, unsigned int>::const_iterator it;
            it = framePatterns.find(framePattern);
            if(it == framePatterns.end()) {
                framePatterns[framePattern] = 1;
            } else {
                framePatterns[framePattern]++;
            }
            positions++;
        }
        w.push_back(unit);
    }
    U = w.size();

    ifs.close();

    cout << "Frame patterns:" << endl;
    for (map<vector<unsigned int>, unsigned int>::const_iterator it=framePatterns.begin();
            it != framePatterns.end(); ++it) {
        cout << "\t[ ";
        for (unsigned int s=1; s<=S; ++s) {
            cout << it->first[s-1] << " ";
        }
        cout << "]: " << it->second << " (" << round(100*((double)it->second / positions)) 
             << " %)" << endl;
    }

    if (F == 0) {
        cout << "F = automatic" << endl;
        infinite_F = true;
        F = 1;
    } else {
        cout << "F = " << F;
        if (reestimate_F) cout << " (will be reestimated)";
        cout << endl;
    }
    if (R == 0) {
        infinite_R = true;
        //this is untrue in the case of multiple frame patterns
        R = ceil(exp(log(F)/S)); //minimum number of roles
        cout << "R = automatic (min. " << R << ")" << endl;
    } else {
        cout << "R = " << R;
        if (reestimate_R) cout << " (will be reestimated)";
        cout << endl;
    }

    
    if (F > pow(R, S)) {
        //this is untrue in the case of multiple frame patterns
        cout << "Number of frames (F) must be lower than or equal to the number of all " <<
                "combinations of possible semantic roles (R^S)." << endl;
        return false;
    }
    if (F < framePatterns.size() && !infinite_F) {
        cout << "Number of frames (F) must be higher or equal to the number of different " <<
                "frame patterns in the input data." << endl;
        return false;
    }
    
    frameSet.setS(S);

    cout << "alpha = " << alpha0 << endl;
    cout << "beta = " << beta0 << endl;
    cout << "gamma = " << gamma0 << endl;
    cout << "delta = " << delta << endl;
    cout << "Lexical units = " << U << endl;
    cout << "Slots = " << S << endl;
    cout << "Vocabulary size = " << V << endl;


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

    for (set<unsigned int>::const_iterator it = used_frames.begin(); it != used_frames.end(); ++it) {
        for (unsigned int s=1; s<=S; ++s) {
            rfile << roles[*it-1][s-1];
            if (s != S) rfile << " ";
        }
        if (it != used_frames.end()) rfile << endl;
    }

    ffile.close();
    rfile.close();
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
    for (set<unsigned int>::const_iterator it = used_frames.begin(); it != used_frames.end(); ++it) {
        for (unsigned int s=1; s<=S; ++s) {
            cout << roles[*it-1][s-1];
            if (s != S) cout << " ";
        }
        if (it != used_frames.end()) cout << endl;
    }
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
    lfile << "Estimate number of frames:\t" << int(infinite_F) << endl;
    lfile << "Number of roles:\t" << R << endl;
    lfile << "Estimate number of roles:\t" << int(infinite_R) << endl;
    lfile << "Number of slots:\t" << S << endl;
    lfile << "Alpha:\t" << alpha0 << endl;
    lfile << "Beta:\t" << beta0 << endl;
    lfile << "Gamma:\t" << gamma0 << endl;
    lfile << "Delta:\t" << delta << endl;
    lfile << "Required number of iterations:\t" << aiter << endl;
    lfile << "Last iteration:\t" << citer << endl;
    lfile.close();

    return true;
}

bool Sampler_t::recoverParameters(string logDir) {

    string fname = logDir + "lda-frames.log";
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

    ifstream lfile(fname.c_str());
    if (!lfile.is_open()) {
        cout << "Cannot open log file '" << lfile << "\n";
        return false;
    }
    requiredIters = 0;

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
            if (lineItems.at(0) == "Estimate number of frames:") {
                unsigned int tmp = atoi(lineItems.at(1).c_str());
                if (tmp < 0 || tmp > 1) {
                    cout << "Wrong value of 'estimate number of frames' in the log file (" << lineItems.at(1) << ")." << endl;
                    return false;
                }
                infinite_F = bool(tmp);
            }
            if (lineItems.at(0) == "Number of roles:") {
                R = atoi(lineItems.at(1).c_str());
                if (R <= 0) {
                    cout << "Wrong number of roles in the log file (" << lineItems.at(1) << ")." << endl;
                    return false;
                }
            }
            if (lineItems.at(0) == "Estimate number of roles:") {
                unsigned int tmp = atoi(lineItems.at(1).c_str());
                if (tmp < 0 || tmp > 1) {
                    cout << "Wrong value of 'estimate number of roles' in the log file (" << lineItems.at(1) << ")." << endl;
                    return false;
                }
                infinite_R = bool(tmp);
            }
            if (lineItems.at(0) == "Alpha:") {
                alpha0 = atof(lineItems.at(1).c_str());
                if (alpha0 <= 0) {
                    cout << "Wrong alpha in the log file (" << lineItems.at(1) << ")." << endl;
                    return false;
                }
            }
            if (lineItems.at(0) == "Beta:") {
                beta0 = atof(lineItems.at(1).c_str());
                if (beta0 <= 0) {
                    cout << "Wrong beta in the log file(" << lineItems.at(1) << ")." << endl;
                    return false;
                }
            }
            if (lineItems.at(0) == "Gamma:") {
                gamma0 = atof(lineItems.at(1).c_str());
                if (gamma0 <= 0) {
                    cout << "Wrong gamma in the log file(" << lineItems.at(1) << ")." << endl;
                    return false;
                }
            }
            if (lineItems.at(0) == "Delta:") {
                delta = atof(lineItems.at(1).c_str());
                if (delta <= 0) {
                    cout << "Wrong delta in the log file(" << lineItems.at(1) << ")." << endl;
                    return false;
                }
            }
            if (lineItems.at(0) == "Required number of iterations:") {
                requiredIters = atoi(lineItems.at(1).c_str());
                if (requiredIters <= 0) {
                    cout << "Wrong required iteration number in the log file (" <<
                            lineItems.at(1) << ")." << endl;
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

bool Sampler_t::recoverData(string dataDir, unsigned int burn_in) {
    string ffname = dataDir + "frames.smpl";
    string rfname = dataDir + "roles.smpl";
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

    ifstream ffile(ffname.c_str());
    if (!ffile.is_open()) {
        cout << "Cannot open file '" << ffname << "'.\n";
        return false;
    }

    ifstream rfile(rfname.c_str());
    if (!rfile.is_open()) {
        cout << "Cannot open file '" << rfname << "'.\n";
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
            tau[frame] = 0;

            fc_f[frame-1]++;
            for (unsigned int s=1; s<=S; ++s) {
                if (w[u-1][t-1][s-1] != 0) {
                    fc_fsw[frame-1][s-1][w[u-1][t-1][s-1]-1]++;
                }
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

    if (!infinite_F) { //ensure that all frames are in the set
        for (unsigned int f=1; f<=F; ++f) {
            used_frames.insert(f);
        }
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
            if (role > R) {
                cout << "Wrong role number (" << *tok_iter << ")." << endl;
                ffile.close();
                rfile.close();
                return false;
            }
            //recover data
            roles[f-1][s-1] = role;
            if (role > 0) used_roles.insert(role);
        }
        if (s != S) {
            cout << "Wrong number of slots (line #" << f <<")." << endl;
            ffile.close();
            rfile.close();
            return false;
        }

    }
    
    if (!infinite_R) { //ensure that all roles are in the set
        for (unsigned int r=1; r<=R; ++r) {
            used_roles.insert(r);
        }
    }

    //recover data
    FrameKey_t fk = frameSet.makeKey(roles[f-1]);
    frameSet.insert(fk);
    if (f != F) {
        cout << "Wrong number of frames." << endl;
        ffile.close();
        rfile.close();
        return false;
    }

    cout << "...resampling phi and theta." << endl;

    initialize_beta(); 
    initialize_post_phi();
    resample_post_phi();
    initialize_post_theta();
    resample_post_theta();
    initialize_post_omega();
    if (infinite_F) {
        initialize_infinite_vars();
        cout << "...resampling tau." << endl;
        resample_tau();
    }
    if (startIter > burn_in) {
        cout << "...resampling hyperparameters." << endl;
        resample_hypers(20);
    }

    ffile.close();
    rfile.close();
    return true;
}

