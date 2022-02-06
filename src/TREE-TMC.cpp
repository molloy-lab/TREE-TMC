#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iomanip>
#include "heuristics/maxcut/burer2002.h"
#include "problem/instance.h"
#include "problem/max_cut_instance.h"


std::string input_file = "", output_file = "";
int verbose = 0, tid = 0;
std::ofstream logs[8];

const std::string help = 
"=================================== TREE-TMC ===================================\n"
"This is version 1.0.0 of TREe Embedded Triplet Max Cut (TREE-TMC).\n\n" 
"USAGE:\n"
"./TREE-QMC (-i|--input) <input file> [(-o|--output) <output file>]\n"
"           [(--maxcutseed) <integer>]\n"
"           [(-v|--verbose) <verbose mode>] [-h|--help]\n\n"
"OPTIONS:\n"
"[-h|--help]\n"
"        Prints this help message.\n"
"(-i|--input) <input file>\n"
"        Name of file containing rooted gene trees in newick format (required)\n"
"[(-o|--output) <output file>]\n"
"        Name of file for writing output rooted species tree (default: stdout)\n"
"[(--maxcutseed) <integer>]\n"
"        Seeds random number generator with <integer> prior to calling the max\n"
"        cut heuristic but after the preprocessing phase. If <integer> is set to\n"
"        -1, system time is used; otherwise, <integer> should be positive\n"
"        (default: 1).\n"
"[(-v|--verbose) <verbose mode>]\n"
"        -v 0: write no subproblem information (default)\n"
"        -v 1: write CSV with subproblem information (subproblem ID, parent\n"
"              problem ID, depth of recursion, number of taxa in subproblem,\n"
"              number of artificial taxa in the subproblem)\n"
"        -v 2: also write subproblem trees in newick format\n"
"        -v 3: also write subproblem triplet graphs in phylip matrix format\n\n"
"================================================================================\n\n";





int main(int argc, char** argv) {
    auto start = std::chrono::high_resolution_clock::now();

    std::cout << "TREE-TMC version 1.0.0\nCOMMAND: ";
    for (int i = 0; i < argc; i++) {
        std::cout << argv[i] << ' ';
    }
    std::cout << std::endl << std::endl;

    int cutseed = 1;
    int execution = 0;

    for (int i = 0; i < argc; i ++) {
        std::string opt(argv[i]);
        if (opt == "-h" || opt == "--help") { std::cout << help; return 0; }
        if (opt == "-i" || opt == "--input" && i < argc - 1) input_file = argv[++ i];
        if (opt == "-o" || opt == "--output" && i < argc - 1) output_file = argv[++ i];
        if (opt == "-v" || opt == "--verbose") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "0" && param != "1" && param != "2" && param != "3") {
                std::cout << "ERROR: invalid verbose parameter!" << std::endl;
                std::cout << help;
                return 0;
            }
            verbose = std::stoi(param);
        }
        if (opt == "--maxcutseed" && i < argc - 1) {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "-1") {
                for (int j = 0; j < param.length(); j ++) {
                    if (param[j] < '0' || param[j] > '9') {
                        std::cout << "ERROR: invalid maxcutseed parameter!" << std::endl;
                        std::cout << help;
                        return 0;
                    }
                }
            }
            cutseed = std::stoi(param);
        }
    }
    std::ifstream fin(input_file);
    if (! fin.is_open()) {
        std::cout << "ERROR: input file " << input_file << " does not exist!" << std::endl;
        std::cout << help;
        return 0;
    }

    // Read input gene trees

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Execution time: " << duration.count() << "ms" << std::endl;

    return 0;
}

