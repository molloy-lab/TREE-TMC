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
#include <climits>
#include <deque>
// #include <inttypes.h>
#include "heuristics/maxcut/burer2002.h"
#include "problem/instance.h"
#include "problem/max_cut_instance.h"



// Allows user to set at compile time
// Use compile flags -DUSE_SHRT or -DUSE_LONG
#ifdef USE_SHRT
    typedef uint16_t index_t;
    #define INDEX_MAX USHRT_MAX
#elif USE_LONG
    typedef uint64_t index_t;
    #define INDEX_MAX ULONG_MAX
#else
    typedef uint32_t index_t;
    #define INDEX_MAX UINT_MAX
#endif  // USE_USHRT or USE_ULONG


#define BYTES_PER_GB 1073741824


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

// Adds elements of list 2 to list 1 (does not copy!)
template <typename T> void extend(T &list1, T &list2) {
    typename T::iterator it;
    for (it = list2.begin(); it != list2.end(); ++it) {
        list1.push_back(*it);
    }
}


/*
    basic n x n x 2 matrix processing functions
    a quartet graph is represented by two matrices G and B
    storing weights of good and bad edges, respectively

    TODO: We probably want to make this 1D and then
          provide some functions or operator overloads so that it
          can be accessed like it's 3D
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
        for (int i = 0; i < nrow; i ++) {
            for (int j = 0; j < ncol; j ++) {
                ss << std::setw(2) << std::setprecision(6) << m[i][j][k];
            }
            ss << std::endl;
        }
        ss << std::endl;
    }
    return ss.str();
}


class Graph {
    public:
        Graph(std::vector<std::string> &index2label , size_t *** matrix);
        double get_cut(std::unordered_set<std::string> *A, std::unordered_set<std::string> *B);
        ~Graph();
    private:
        size_t size;
        std::unordered_map<std::string,size_t> label2index;
        std::vector<std::string> index2label;
        size_t ***graph;
        double sdp_cut(double alpha, std::unordered_set<std::string> *A, std::unordered_set<std::string> *B);
};

Graph::Graph(std::vector<std::string> &index2label, size_t *** matrix) { 
    size=index2label.size();
    for (size_t i = 0; i < size; i++) {
        label2index[index2label[i]] = i;
        this->index2label.push_back(index2label[i]);
    }
    graph = Matrix3D::new_mat<size_t>(size, size, 2);
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
            graph[i][j][0] = matrix[i][j][0];
            graph[i][j][1] = matrix[i][j][1];
        }
    }
}

double Graph::get_cut(std::unordered_set<std::string> *A, std::unordered_set<std::string> *B) {
    double positive_weight = -1.0;
    std::unordered_set<std::string> a, b;
    double lower = 0.0, upper = 6.0;
    while (lower + 0.1 < upper) {
        double alpha = (lower + upper) / 2.0;
        a.clear(); b.clear();
        double weight = sdp_cut(alpha, &a, &b);
        //if (weight < 0.001) {
        if (weight < 0.001 || a.size() == 1 || b.size() == 1) { 
            upper = alpha;
        }
        else {
            lower = alpha;
            positive_weight = weight;
            *A = a;
            *B = b;
        }
    }
    return positive_weight;
}

Graph::~Graph() {
    Matrix3D::delete_mat(graph,size,size);
}

double Graph::sdp_cut(double alpha, std::unordered_set<std::string> *A, std::unordered_set<std::string> *B) {
    std::vector<Instance::InstanceTuple> input;
    double sum = 0;
    for (int i = 0; i < size; i ++) {
        for (int j = i + 1; j < size; j ++) {
            double weight = graph[i][j][0] - alpha * graph[i][j][1];
            if (weight < 0) weight = - weight;
            sum += weight;
        }
    }
    double norm = size * (size - 1) / 2 / sum;
    for (int i = 0; i < size; i ++) {
        for (int j = i + 1; j < size; j ++) {
            double weight = (graph[i][j][0] - alpha * graph[i][j][1]) * norm;
            input.push_back(Instance::InstanceTuple(std::make_pair(i + 1, j + 1), weight));
        }
    }
    MaxCutInstance instance(input, size);
    Burer2002 heuristic(instance, -1, false, NULL);
    MaxCutSimpleSolution solution = heuristic.get_best_solution();
    std::vector<int> cut = solution.get_assignments();
    for (int i = 0; i < cut.size(); i ++) {
        if (cut[i] < 0) 
            A->insert(index2label[i]);
        else 
            B->insert(index2label[i]);
    }
    return solution.get_weight();
}

class Node {
    public:
        Node();
        Node(std::string name);
        ~Node();
        bool is_root();
        bool is_leaf();
        size_t num_children();
        Node* get_parent();
        void add_child(Node *child);
        void remove_child(Node *child);
        void contract();
        void add_children_to_list(std::list<Node*> &nodelist);
        void update_label_list(std::list<std::string> labellist);
        void update_label_list(std::string str);
        std::list<std::string> get_label_list();
        void suppress_unifurcations();
        void compute_c();
        void update_c(Node *child);
        size_t get_c();
        std::string newick(bool printindex=false);
        std::string label;  // TODO: we probably want to get rid of this at some point
                            // so we aren't storing so many copies of labels!
                            // this could be done if we read trees onto the same
                            // label set... 
        index_t index;
        index_t size;
        std::list<Node*> children;
        std::list<std::string> label_list;
    private:
        index_t c = 0;
        Node *parent;
};


namespace Traverse {
struct ToRoot
{
    // constructor that takes in a node
    ToRoot() {
        previous_node = NULL;
        current_node = NULL;
    };

    ToRoot(Node *node) {
        previous_node = NULL;
        current_node = node;
    };

    // incrementing means going to the parent node
    ToRoot &operator++() noexcept
    {
        if (current_node != NULL) {
            previous_node = current_node;
            current_node = current_node->get_parent();
        }
        return *this;
    };

    // post fixing is bad in general but it has it's usages
    ToRoot operator++(int) noexcept
    {
        ToRoot tempIter = *this;   // make a copy of the iterator
        ++*this;                   // increment
        return tempIter;           // return the copy before increment
    };

    // compare nodes
    bool operator!=(const ToRoot &other) const noexcept
    {
        return this->current_node != other.current_node;
    };

    // return the node (dereference operator)
    Node* operator*() const noexcept
    {
        return this->current_node;
    };

    // return a const pointer to the front
    ToRoot begin() const noexcept
    {
        return ToRoot(this->current_node);
    };

    // return a const pointer to the back - the back is always null
    ToRoot end() const noexcept
    {
        return ToRoot();
    };

    private:
        Node *previous_node = NULL;
        Node *current_node = NULL;
};

struct PreOrder
{
    // constructor that takes in a node
    PreOrder() {
        previous_node = NULL;
        current_node = NULL;
    };

    PreOrder(Node *node) {
        previous_node = NULL;
        current_node = node;
        current_node->add_children_to_list(nodelist);
        // TODO: Maybe be able to get rid of storage... ?
    };

    // incrementing means going to the parent node
    PreOrder &operator++() noexcept
    {
        previous_node = current_node;
        if (nodelist.size() == 0) {
            current_node = NULL;
            nodelist.clear();
        } else {
            current_node = nodelist.front();
            nodelist.pop_front();
            current_node->add_children_to_list(nodelist);
        }

        return *this;
    };

    // post fixing is bad in general but it has it's usages
    PreOrder operator++(int) noexcept
    {
        PreOrder tempIter = *this;  // make a copy of the iterator
        ++*this;                    // increment
        return tempIter;            // return the copy before increment
    };

    // compare nodes
    bool operator!=(const PreOrder &other) const noexcept
    {
        return this->current_node != other.current_node;
    };

    // return the node (dereference operator)
    Node* operator*() const noexcept
    {
        return this->current_node;
    };

    // return a const pointer to the front
    PreOrder begin() const noexcept
    {
        return PreOrder(this->current_node);
    };

    // return a const pointer to the back - the back is always null
    PreOrder end() const noexcept
    {
        return PreOrder();
    };

    private:
        Node *previous_node = NULL;
        Node *current_node = NULL;
        std::list<Node*> nodelist;
};


struct PostOrder
{
    // constructor that takes in a node
    PostOrder() {
        previous_node = NULL;
        current_node = NULL;
    };

    PostOrder(Node *node) {
        Node *tmp;

        previous_node = NULL;

        nodelist1.push_back(node);

        while (nodelist1.size() != 0) {
            tmp = nodelist1.back();
            nodelist1.pop_back();
            tmp->add_children_to_list(nodelist1);
            nodelist2.push_back(tmp);
        }

        nodelist1.clear();
        // TO DO: Maybe be able to get rid of stoarge!

        while (nodelist2.size() != 0) {
            tmp = nodelist2.back();
            nodelist2.pop_back();

            current_node = tmp;
            break;
        }
    };

    // incrementing means going to the next node in the postorder traversal
    PostOrder &operator++() noexcept
    {
        Node *node;

        previous_node = current_node;
        if (nodelist2.size() == 0) {
            current_node = NULL;
            nodelist2.clear();
        } else {
            node = nodelist2.back();
            nodelist2.pop_back();
            current_node = node;
        }

        return *this;
    };

    // post fixing is bad in general but it has it's usages
    PostOrder operator++(int) noexcept
    {
        PostOrder tempIter = *this;   // make a copy of the iterator
        ++*this;                      // increment
        return tempIter;              // return the copy before increment
    };

    // compare nodes
    bool operator!=(const PostOrder &other) const noexcept
    {
        return this->current_node != other.current_node;
    };

    // return the node (dereference operator)
    Node* operator*() const noexcept
    {
        return this->current_node;
    };

    // return a const pointer to the front
    PostOrder begin() const noexcept
    {
        return PostOrder(this->current_node);
    };

    // return a const pointer to the back - the back is always null
    PostOrder end() const noexcept
    {
        return PostOrder();
    };

    private:
        Node *previous_node = NULL;
        Node *current_node = NULL;
        std::list<Node*> nodelist1;
        std::list<Node*> nodelist2;
};


struct Leaves
{
    // constructor that takes in a node
    Leaves() {
        previous_node = NULL;
        current_node = NULL;
    };

    Leaves(Node *node) {
        Node *tmp;

        previous_node = NULL;

        nodelist.push_back(node);

        while (nodelist.size() > 0) {
            tmp = nodelist.front();
            nodelist.pop_front();
            tmp->add_children_to_list(nodelist);

            if (tmp->is_leaf()) {
                current_node = tmp;
                break;
            }
        }
    };

    // incrementing means going to the next leaf node
    Leaves &operator++() noexcept
    {
        Node *node;

        previous_node = current_node;

        if (nodelist.size() == 0) {
            current_node = NULL;
            nodelist.clear();
        } else {
            while (nodelist.size() > 0) {
                node = nodelist.front();
                nodelist.pop_front();
                node->add_children_to_list(nodelist);

                if (node->is_leaf()) {
                    current_node = node;
                    break;
                }
            }
        }

        return *this;
    };

    // post fixing is bad in general but it has it's usages
    Leaves operator++(int) noexcept
    {
        Leaves tempIter = *this;   // make a copy of the iterator
        ++*this;                   // increment
        return tempIter;           // return the copy before increment
    };

    // compare nodes
    bool operator!=(const Leaves &other) const noexcept
    {
        return this->current_node != other.current_node;
    };

    // return the node (dereference operator)
    Node* operator*() const noexcept
    {
        return this->current_node;
    };

    // return a const pointer to the front
    Leaves begin() const noexcept
    {
        return Leaves(this->current_node);
    };

    // return a const pointer to the back - the back is always null
    Leaves end() const noexcept
    {
        return Leaves();
    };

    private:
        Node *previous_node = NULL;
        Node *current_node = NULL;
        std::list<Node*> nodelist;
};

}; // namespace traversal


Node::Node() {
    parent = NULL;
    label = "";
    index = INDEX_MAX;
    size = 0;
}

Node::Node(std::string name) {
    parent = NULL;
    label = name;
    index = INDEX_MAX;
    size = 0;
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

bool Node::is_root() {
    if (parent == NULL) return true;
    return false;
}

bool Node::is_leaf() {
    if (children.size() == 0) return true;
    return false;
}

size_t Node::num_children() {
    return children.size();
}

size_t Node::get_c() {
    return c;
}

void Node::update_label_list(std::list<std::string> newnodes) {
    label_list.merge(newnodes);
}

void Node::update_label_list(std::string str) {
    label_list.push_back(str);
}

std::list<std::string> Node::get_label_list() {
    return label_list;
}

Node* Node::get_parent() {
    return parent;
}

void Node::add_child(Node *child) {
    if (child == NULL) return;

    child->parent = this;
    children.push_back(child);

}
void Node::update_c(Node *child) {
    c += child->get_c();
}

void Node::remove_child(Node *child) {
     if (child == NULL) return;

    child->parent = NULL;
    children.remove(child);
}

void Node::compute_c() {
    auto nodeItr = Traverse::PostOrder(this);
    for (; nodeItr != nodeItr.end(); ++nodeItr) {
        if ((*nodeItr)->is_leaf()) {
            (*nodeItr)->c = 1;
            //std::cout << "c["+ ((*nodeItr)->label) + "] = " + std::to_string((*nodeItr)->get_c()) << std::endl;
        continue;
        }
        std::list<Node*>::iterator it;
        for (it = (*nodeItr)->children.begin(); it != (*nodeItr)->children.end(); ++it) {
            (*nodeItr)->update_c(*it);
        }
        //std::cout << "c["+ ((*nodeItr)->label) + "] = " + std::to_string((*nodeItr)->get_c()) << std::endl;
    }
}

void Node::contract() {
    if (parent == NULL) return;

    std::list<Node*>::iterator it;
    for (it = children.begin(); it != children.end(); ++it) {
        parent->add_child(*it);
    }
    parent->remove_child(this);

    parent = NULL;
    children.clear();

    delete this;
}


void Node::add_children_to_list(std::list<Node*> &nodelist) {
    extend(nodelist, children);
}

void Node::suppress_unifurcations() {
    auto nodeItr = Traverse::PreOrder(this);
    for (; nodeItr != nodeItr.end(); ++nodeItr) {
        if ((*nodeItr)->num_children() == 1) (*nodeItr)->contract();
    }
}


std::string Node::newick(bool printindex) {
    // TODO: Maybe change to pass by reference?

    std::list<Node*>::iterator it;
    std::string out;

    if (this->is_leaf()) {
        out = "";
        if (!this->label.empty()) {
            out += this->label;
            if (printindex) {
                out += ':';
                out += std::to_string(this->index);
            }         
        }
        return out;
    }

    out = '(';

    for (it = this->children.begin(); it != this->children.end(); ++it) {
        out += (*it)->newick(printindex);
        out += ',';
    }
    out.pop_back();  // Drop trailing comma
    out += ')';

    if (!this->label.empty()) out += this->label;

    return out;
}


class Tree {
    public:
        Tree();
        Tree(Node *myroot);
        Tree(const std::string &text);
        ~Tree();
        Node* get_root();
        void suppress_unifurcations();
        void compute_c();
        std::string newick(bool printindex=false);
        Tree* get_induced_subtree_copy(std::unordered_set<std::string> taxa);
    //private:
        Node *root;
};


Tree::Tree() {
    root = new Node();
}


Tree::Tree(Node *myroot) {
    root = myroot;
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
            node = node->get_parent();
        } else if (text[i] == ',') {
            // Go to new sibling
            node = node->get_parent();
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

    root->suppress_unifurcations();
}


Tree::~Tree() { 
    auto nodeItr = Traverse::PostOrder(root);
    for (; nodeItr != nodeItr.end(); ++nodeItr) {
        delete *nodeItr;
    }
}


Node* Tree::get_root() {
    return root;
}


void Tree::suppress_unifurcations() {
    root->suppress_unifurcations();
}

void Tree::compute_c() {
    root->compute_c();
}


std::string Tree::newick(bool printindex) {
    std::string out = this->root->newick(printindex);
    out += ";";
    return out;
}

Tree* Tree::get_induced_subtree_copy(std::unordered_set<std::string> taxa) {
    std::unordered_map<std::string, Node*> label_to_leaf; 
    std::unordered_set<Node*> keep;

    auto leafItr = Traverse::Leaves(this->get_root());
    for (; leafItr != leafItr.end(); ++leafItr) {
        label_to_leaf[(*leafItr)->label] = (*leafItr);
        if (taxa.count((*leafItr)->label)) {
            keep.insert((*leafItr));
        }
    }

    for (Node* n : keep) {
        auto rootItr = Traverse::ToRoot(n);
        for (; rootItr != rootItr.end(); ++rootItr) {
            keep.insert((*rootItr));
        }
    }
    

    Tree* out = new Tree();
    out->root->label = this->root->label; //change to out->root->set_label(this->root->get_label()); 
    std::deque<Node*> q_old; 
    q_old.push_back(this->root);
    std::deque<Node*> q_new;
    q_new.push_back(out->root);
    while (q_old.size() != 0){
        Node* n_old = q_old.front(); 
        q_old.pop_front(); 
        Node* n_new = q_new.front();
        q_new.pop_front();
        std::list<Node*>::iterator c_old;
        for (c_old = (n_old)->children.begin(); c_old != (n_old)->children.end(); ++c_old) {
            if (keep.count(*c_old)) {
                Node *c_new = new Node();
                c_new->label = (*c_old)->label;
                n_new->add_child(c_new);
                q_old.push_back(*c_old);
                q_new.push_back(c_new);
            
            }   
        }   
    }
    out->suppress_unifurcations();
    return out;
}

class Forest {
    public:
        Forest();
        Forest(std::vector<Tree*> &mytrees);
        ~Forest();
        size_t num_trees();
        size_t num_labels();
        void compute_c();
        std::vector<Tree*> fetch_trees();
        Forest get_induced_subforest_copy(std::unordered_set<std::string> taxa);
    //private:
        std::vector<Tree*> trees;
        std::vector<std::string> index2label;
        std::unordered_map<std::string, index_t> label2index;
};


Forest::Forest(std::vector<Tree*> &mytrees) {
    std::unordered_map<std::string, index_t>::const_iterator mapIter;
    Node *root, *leaf;
    std::string label;
    index_t index;

    trees = mytrees;

    for (size_t i = 0; i < trees.size(); i++) {
        auto leafItr = Traverse::Leaves(trees[i]->get_root());
        for (; leafItr != leafItr.end(); ++leafItr) {
            leaf = *leafItr;
            label = leaf->label;

            mapIter = label2index.find(label);
            if ( mapIter == label2index.end() ) {
                index = index2label.size();
                leaf->index = (index_t) index;
                label2index.insert({label, index});
                index2label.push_back(label);
            }
            else {
                leaf->index = label2index.at(label);
            }
        }
    }
}


Forest::~Forest() {

}


size_t Forest::num_trees() {
    return trees.size();
}

std::vector<Tree*> Forest::fetch_trees() {
    return trees;
}


size_t Forest::num_labels() {
    return index2label.size();
}

void Forest::compute_c() {
    for (size_t i = 0; i < num_trees(); i++ ) {
        trees[i]->compute_c();
    }
}

Forest Forest::get_induced_subforest_copy(std::unordered_set<std::string> taxa) {
    std::vector<Tree*> tree_vec;
    for (Tree* t : trees) {
        auto nt = t->get_induced_subtree_copy(taxa);
        tree_vec.push_back(nt);
    }
    auto out = Forest(tree_vec);
    return out;
}


namespace TripletMaxCut {
   void main(std::vector<Tree*> input);
   // void compute_good_and_bad_edges();
}


Tree* TripletMaxCut::main(std::vector<Tree*> input) {
    size_t ***matrix;  // could probably get away with float...
    size_t n, k;
    std::vector<Tree*> trees;
    Forest *forest = new Forest(input);

    n = forest->num_labels();
    k = forest->num_trees();
    trees = forest ->fetch_trees();
    

    // Make a matrix
    matrix = Matrix3D::new_mat<size_t>(n, n, 2);

    // Get c[v] (the number of taxa beneath v) for each vertex v
    forest->compute_c();
    //Build B and G post-order
    for (size_t t = 0; t < k; t++) {
        //std::cout << std::to_string(t) << std::endl;
        auto tre = trees[t];
        auto nodeItr = Traverse::PostOrder(tre->get_root());
        for (; nodeItr != nodeItr.end(); ++nodeItr) {
            //if we are at a leaf X we should just add it's label to a singleton list [X].
            if ((*nodeItr)->is_leaf()) {
                (*nodeItr)->update_label_list((*nodeItr)->label);
                continue;
            }
            //otherwise we should iterate through the children of the node and compute B[i,j]
            size_t i = 0;
            std::list<Node*>::iterator ci;
            for (ci = (*nodeItr)->children.begin(); ci != (*nodeItr)->children.end(); ++ci) {
                size_t j = 0;
                std::list<Node*>::iterator cj;
                for (cj = (*nodeItr)->children.begin(); cj != (*nodeItr)->children.end(); ++cj) {
                    if (i < j) {
                        //for each unique pair of children we compute B[i,j] for all unique pairs.
                        //note:this is wrong as it stands. we should be looking at the label_lists
                        std::list<std::string>::iterator it1;
                        for (it1 = (*ci)->label_list.begin(); it1 != (*ci)->label_list.end(); ++it1) {
                            std::cout << *it1 << std::endl;
                            std::list<std::string>::iterator it2;
                            for (it2 = (*cj)->label_list.begin(); it2 != (*cj)->label_list.end(); ++it2) {

                                index_t labelit1 = forest->label2index[*it1];
                                index_t labelit2 = forest->label2index[*it2];
                                index_t c_lca_ij = (*nodeItr)->get_c();
                                index_t c_ci = (*ci)->get_c();
                                index_t c_cj = (*cj)->get_c();
                                //should n here be n, or should it be the number of taxa in the tree, which is <= n...?
                                //update Bad edges
                                matrix[labelit1][labelit2][1] += n - c_lca_ij;
                                matrix[labelit2][labelit1][1] += n - c_lca_ij;
                                //update Good edges
                                matrix[labelit1][labelit2][0] += (c_ci + c_cj)-2;
                                matrix[labelit2][labelit1][0] += (c_ci + c_cj)-2;
                            }

                        } 
                    }
                    j++;
                }
                //update the label list that we'll need for the next step
                (*nodeItr)->update_label_list((*ci)->get_label_list());
                //clear out information we no longer need 
                (*ci)->label_list.clear();
                i++;
            }
            

        }
    }
    //print things to make sense of it all...
    for (size_t phi = 0; phi < n; phi++){

        std:: cout << forest->index2label[phi] + ",";
    }
    std::cout << std::endl;
    std:: cout << Matrix3D::display_mat(matrix,n,n,2) << std::endl;


    // Get the cut
    //
    Graph G(forest->index2label,matrix);

    std::unordered_set<std::string> left_bipartition;
    std::unordered_set<std::string> right_bipartition;

    G.get_cut(&left_bipartition, &right_bipartition);

    Matrix3D::delete_mat(matrix, n, n);
    std::cout << "left bipartition contains: ";

    for (std::string s : left_bipartition) {
        std::cout << " " + s;
    } 
    std::cout << "\nright bipartition contains: ";

    for (std::string s : right_bipartition) {
        std::cout << " " + s;
    }
    std::cout << "\n";

    
    std::vector<Tree*> left_tree_vec;
    for (Tree* t : trees) {
        auto nt = t->get_induced_subtree_copy(left_bipartition);
        left_tree_vec.push_back(nt);
    }

    std::vector<Tree*> right_tree_vec;
    for (Tree* t : trees) {
        auto nt = t->get_induced_subtree_copy(right_bipartition);
        right_tree_vec.push_back(nt);
    }       
    //uncomment the below to make sure that the induced_subtree_code is working (it is)
    //Forest left_subproblem = Forest(left_tree_vec);
    //Forest right_subproblem = Forest(right_tree_vec);
    //std::cout << "left subtree:" + left_subproblem.trees[0]->newick() << std::endl;
    //std::cout << "right subtree:" + right_subproblem.trees[0]->newick() << std::endl;

    // Extract the subproblems, which is a list of Tree*
    // on the subset of the problem...
    // This should be built into trees.... 
    // get_induced_subtree_copy -- 
    // given input tree and leaf set X, return a new tree
    // that is the induced subtree for the input
    //
    
    // create_induced_subtree -- same as above but transform
    // the current tree into the induced subtree for X by deleting...

    // Clean up
    for (size_t i = 0; i < input.size(); i++) {
        delete input[i];
    }
    input.clear();
    //delete forest;

    // Now recurse on the two subproblems

    // Combine the subproblems... after recursion bounces back...
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

    TripletMaxCut::main(input);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Execution time: " << duration.count() << "ms" << std::endl;

    return 0;
}
