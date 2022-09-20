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
    double s;
    int l;
	vector <int> ID;
	vector <edge> edges;
};

class graph{
public:
	graph();
    virtual ~graph();
	int n_nodes, n_edges;
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