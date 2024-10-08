#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>


using namespace std;

class edge{
public:
	edge();
    virtual ~edge();
	int neighbor;
	double w;
	int pos;
};

class node{
public:	
	node();
    virtual ~node();
    double s; // Total absolute weights of the edges incident to the node
    int l;
	vector <int> ID;
	vector <edge> edges;
};

class graph{
public:
	graph();

	// Constructor function to initialize a graph with n vertices
	// and a list of  edges  (bv[i], ev[i]) of weights w[i]
	// The three vectors u, v, w must have the same size
	void from_edge_list(int n, const vector<int> &bv, const vector<int> &ev, const vector<double> &w);

    virtual ~graph();
	int n_nodes, n_edges; // Number of nodes and edges
	vector <node> nodes;
	vector <bool> is_adj;
	vector <double> tmp;
	vector <double> tmp1;	
	vector <int> is_in_group;
	vector <bool> is_flipped;
	vector <int> pos;
	vector <int> rpos;
	bool check_slow;
	int l;
	int k;
	void load_graph(string in_file);
	void delete_node();
	void add_node(node n);
	void swap_nodes(int u, int v);	
	// void merge_nodes();
	// void merge_nodes(int u, int v);
	void merge_nodes(int u, int v, node &a, node &b);
	void restore(node &a, node &b);
	void restore_order(int u, int v);
	void flip(int u);	
	double assign_cut(vector<int> &group,int i, int s);
	double NSI(vector<int> &group);
	void merge_multiple_nodes(vector<int> &group);
};


#endif