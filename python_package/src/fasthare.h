#include "graph.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <ctime>
#include <string>
#include <vector>
#include <cmath>
#include <tuple>

using namespace std;



class FastHare
{
private:
    clock_t start_time;
    ofstream f_log;
    bool logging = false;
    int n_original;

    void print_graph_test();
    void check_graph();

    bool fast_find_NS();

    void find3(int u, int v, double w1, vector< vector<int> > &links);

    bool slow_find_NS();

    bool find_NS(int k);

    void compress(int k);

public:
    // SK Ising Hamiltonian as tuple of (i, j, h_ij)
    typedef   std::vector< std::tuple<int, int, double> >    ski;
    // Output for Fasthare: SK Ising, Spin mapping, 
    // Spin reversal (-1 mean reversed), running time
    typedef   tuple< ski, vector<int>, vector<int>, double > fasthare_output;
    bool parallel_compression = true;
    bool triangle_compression = true;
    double alpha = 1.0;
    double running_time;
    graph G;

    /* Initialize the graph from an SK Ising */
    FastHare(ski sk_ising, double alpha = 1.0);

    /* Load a graph from a .net file */
    FastHare(string net_file, double alpha = 1.0);

    /* Set log file. If the file name is empty
     a string buffer will be used. */
    bool set_log(string log_file);
    double get_running_time();
    void print_graph(string out_file);
    /* Return the compressed graph and compressing information */
    fasthare_output get_output();
    void run();

    
    ~FastHare();
};
