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
#include <list>
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


/*
    basic n x n x 2 matrix processing functions
    a quartet graph is represented by two matrices G and B
    storing weights of good and bad edges, respectively
*/

namespace Matrix3D {
    template <typename T> T ***new_mat(size_t nrow, size_t ncol, size_t size);
    template <typename T> void delete_mat(T ***m, size_t nrow, size_t ncol);
    template <typename T> std::string display_mat(T ***m, size_t nrow, size_t ncol, size_t size);
};


template <typename T>
T ***Matrix3D::new_mat(size_t nrow, size_t ncol, size_t size) {
    T ***m = new T**[nrow];
    for (size_t i = 0; i < nrow; i++) {
        m[i] = new T*[ncol];
        for (size_t j = 0; j < ncol; j++) {
            m[i][j] = new T[size];
            for (size_t k = 0; k < size; k++) {
                m[i][j][k] = 0;
            }
        }
    }
    return m;
}

template <typename T>
void Matrix3D::delete_mat(T ***m, size_t nrow, size_t ncol) {
    for (int i = 0; i < nrow; i ++) {
        for (int j = 0; j < ncol; j ++) {
            delete[] m[i][j];
        }
        delete[] m[i];
    }
    delete[] m;
}

template <typename T>
std::string Matrix3D::display_mat(T ***m, size_t nrow, size_t ncol, size_t size) {
    std::stringstream ss;
    for (size_t k = 0; k < size; k++) {
        ss << "Matrix2D: " << k << std::endl;
        for (int i = 0; i < size; i ++) {
            for (int j = 0; j < size; j ++) {
                ss << std::setw(12) << std::setprecision(6) << m[i][j][k];
            }
            ss << std::endl;
        }
        ss << std::endl;
    }
    return ss.str();
}


class Node {
    public:
        Node();
        Node(const std::string &name);
        ~Node();
        void add_child(Node *child);
        void remove_child(Node *child);
        void traverse_preorder(std::vector<Node*> &traversal);
        void traverse_postorder(std::vector<Node*> &traversal);
        void traverse_leaves(std::vector<Node*> &traversal);
        std::string newick(bool printindex=false);
        Node *parent;
        std::list<Node*> children;
        std::string label;
        size_t index;
        size_t size;
        bool is_leaf;
};


Node::Node() {
    parent = NULL;
    label = "";
    index = SIZE_MAX;
    size = 0;
    is_leaf = true;
}

Node::Node(const std::string &name) {
    parent = NULL;
    label = name;
    index = SIZE_MAX;
    size = 0;
    is_leaf = true;
}

Node::~Node() {
    std::list<Node*>::iterator it;
    for (it = children.begin(); it != children.end(); ++it) {
        this->remove_child(*it);
    }
    if (parent != NULL) {
        parent->remove_child(this);
    }
}

void Node::add_child(Node *child) {
    if (child == NULL) return;

    children.push_back(child);
    child->parent = this;
    is_leaf = false;
}

void Node::remove_child(Node *child) {
     if (child == NULL) return;

    children.remove(child);
    if (children.size() == 0) is_leaf = true;
}

void Node::traverse_preorder(std::vector<Node*> &traversal) {
    std::list<Node*>::iterator it;
    std::list<Node*> nlist;
    Node *node;

    traversal.clear();

    nlist.push_back(this);

    while (nlist.size() > 0) {
        node = nlist.front();
        nlist.pop_front();

        traversal.push_back(node);  // we want this to be yield...

        for (it = node->children.begin(); it != node->children.end(); ++it) {
            nlist.push_back(*it);
        }
    }
}

void Node::traverse_leaves(std::vector<Node*> &traversal) {
    std::list<Node*>::iterator it;
    std::list<Node*> nlist;
    Node *node;

    traversal.clear();

    nlist.push_back(this);

    while (nlist.size() > 0) {
        node = nlist.front();
        nlist.pop_front();

        if (node->is_leaf)
            traversal.push_back(node);  // we want this to be yield...

        for (it = node->children.begin(); it != node->children.end(); ++it) {
            nlist.push_back(*it);
        }
    }
}

void Node::traverse_postorder(std::vector<Node*> &traversal) {
    std::list<Node*>::iterator it;
    std::list<Node*> nlist1;
    std::list<Node*> nlist2;
    Node *node;

    traversal.clear();

    nlist1.push_back(this);

    while (nlist1.size() != 0) {
        node = nlist1.back();
        nlist1.pop_back();

        // Extend
        for (it = node->children.begin(); it != node->children.end(); ++it) {
             nlist1.push_back(*it);
        }

        nlist2.push_back(node);
    }

    while (nlist2.size() != 0) {
        node = nlist2.back();
        nlist2.pop_back();

        traversal.push_back(node); // we want this to be yield...
    }
}


std::string Node::newick(bool printindex) {
    // TODO: Maybe change this to pass by reference on string...

    std::list<Node*>::iterator it;
    std::string out = "";
    Node *node;
    node = this;

    if (node->is_leaf) {
        if (!node->label.empty()) {
            out += node->label;
            if (printindex) {
                out += ':';
                out += std::to_string(node->index);
            }         
        }
    } else {
        out += '(';

        for (it = node->children.begin(); it != node->children.end(); ++it) {
            out += (*it)->newick(printindex);
            out += ',';
        }
        out.pop_back();  // Drop trailing comma
        out += ')';

        if (!node->label.empty()) out += node->label;
    }

    return out;
}

class Tree {
    public:
        Tree();
        Tree(Node *myroot);
        Tree(const std::string &text);
        ~Tree();
        void traverse_preorder(std::vector<Node*> &traversal);
        void traverse_postorder(std::vector<Node*> &traversal);
        void traverse_leaves(std::vector<Node*> &traversal);
        std::string newick(bool printindex=false);
        Node *root;
};


Tree::Tree() {
    root = new Node();
}


Tree::Tree(Node *myroot) {
    root = myroot;

    //tree.suppress_unifurcations();
}


Tree::Tree(const std::string &text) {
    Node *node, *child;
    std::string label;
    size_t i;

    root = new Node();
    node = root;

    i = 0;
    while (i < text.size()) {
        if (text[i] == ';') {
            // End of Newick string
            if ((i != text.size() - 1) || (node != root)) {
                std::cerr << "Not a valid newick!" << std::endl;
                exit(1);
            }  
        } else if (text[i] == '(') {
            // Go to new child
            child = new Node();
            node->add_child(child);
            node = child;
        } else if (text[i] == ')') {
            // Go to parent
            node = node->parent;
        } else if (text[i] == ',') {
            // Go to new sibling
            node = node->parent;
            child = new Node();
            node->add_child(child);
            node = child;
        } else if (text[i] == ':') {
            // Parse edge length
            i += 1;
            while ((text[i] != ',') && (text[i] != ')') && (text[i] != ';'))
                i += 1;
            i -= 1;
        } else if ((text[i] == ' ') || (text[i] == '\t')) {
            continue;
        } else {
            // Parse node label
            label.clear();
            while ((text[i] != ':') && (text[i] != ',') && 
                    (text[i] != ')') && (text[i] != ';')) {
                label += text[i];
                i += 1;
            }
            i -= 1;
            node->label = label;
        }
        i += 1;
    }
}

Tree::~Tree() { 
    std::vector<Node*> traversal;
    root->traverse_postorder(traversal);
    for (size_t i = 0; i < traversal.size(); i++) {
        std::cout << "  Deleting " << traversal[i]->label << std::endl;
        delete traversal[i];
    }
}

std::string Tree::newick(bool printindex) {
    std::string out = this->root->newick(printindex);
    out += ";";
    return out;
}

void Tree::traverse_preorder(std::vector<Node*> &traversal) {
    traversal.clear();
    root->traverse_preorder(traversal);
}

void Tree::traverse_postorder(std::vector<Node*> &traversal) {
    traversal.clear();
    root->traverse_postorder(traversal);
}

void Tree::traverse_leaves(std::vector<Node*> &traversal) {
    traversal.clear();
    root->traverse_leaves(traversal);
}


class Forest {
    public:
        Forest();
        Forest(std::vector<Tree*> &mytrees);
        ~Forest();
        size_t num_trees();
        size_t num_labels();
        std::vector<Tree*> trees;
        std::vector<std::string> index2label;
        std::unordered_map<std::string, size_t> label2index;
};


Forest::Forest(std::vector<Tree*> &mytrees) {
    std::unordered_map<std::string, size_t>::const_iterator it;
    std::vector<Node*> traversal;
    Node *leaf;
    std::string label;
    size_t index;

    trees = mytrees;


    for (size_t i = 0; i < trees.size(); i++) {
        trees[i]->traverse_leaves(traversal);

        std::cout << "Found tree " << trees[i]->newick(true) << std::endl;

        for (size_t j = 0; j < traversal.size(); j++) {
            leaf = traversal[j];
            label = leaf->label;

            std::cout << "Found label " << label << std::endl;

            it = label2index.find(label);
            if ( it == label2index.end() ) {
                index = index2label.size();
                leaf->index = index;
                label2index.insert({label, index});
                index2label.push_back(label);
            }
            else {
                leaf->index = label2index.at(label);
            }
        }
    }
}

Forest::~Forest() {}

size_t Forest::num_trees() {
    return trees.size();
}

size_t Forest::num_labels() {
    return index2label.size();
}

namespace TripletMaxCut {
   void main(std::vector<Tree*> input);
   // void compute_good_and_bad_edges();
}


void TripletMaxCut::main(std::vector<Tree*> input) {
    double ***matrix;
    size_t n, k;
    Forest *forest = new Forest(input);

    n = forest->num_labels();
    k = forest->num_trees();

    // Make a matrix
    matrix = Matrix3D::new_mat<double>(n, n, 2);

    // Compute good and bad edges in the matrix

    // Get the cut

    // Extract the subproblems, which is a list of Tree*
    // on the subset of the problem...
    // This should be built into trees.... 
    // get_induced_subtree_copy -- 
    // given input tree and leaf set X, return a new tree
    // that is the induced subtree for the input
    
    // create_induced_subtree -- same as above but transform
    // the current tree into the induced subtree for X by deleting...

    // Clean up
    Matrix3D::delete_mat(matrix, n, n);
    for (size_t i = 0; i < input.size(); i++) {
        std::cout << "Delete tree " << i << std::endl;
        delete input[i];
    }
    input.clear();
    delete forest;

    // Now recurse on the two subproblems

    // Combine the subproblems... after the recurssion bounces back...
}


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
    std::vector<Tree*> input;
    std::string newick;

    while (std::getline(fin, newick)) {
        if (newick.find(";") == std::string::npos) break;
        Tree *tree = new Tree(newick);
        std::cout << tree->newick() << std::endl;
        input.push_back(tree);
    }
    fin.close();

    // Want to fill in your matrix based on Forest...
    TripletMaxCut::main(input);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Execution time: " << duration.count() << "ms" << std::endl;

    return 0;
}

