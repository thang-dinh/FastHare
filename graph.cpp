#include "graph.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
edge::edge()
{

}

edge::~edge()
{

}

node::node()
{
    ID.clear();
    edges.clear();
    s = 0;
    l = 0;
}

node::~node()
{

}


graph::graph()
{

}

graph::~graph()
{

}

void graph::load_graph(string in_file)
{
    ifstream fi (in_file.c_str());
    int u,v;
    double k;
    fi >> n_nodes >> n_edges;    
    cout << n_nodes << "," << n_edges <<  ",";
    nodes.resize(n_nodes);        
    for (int i = 0; i < n_nodes; ++i)
        nodes[i].ID.push_back(i);
    int x,y;
    for (int i = 0; i < n_edges; ++i)
    {
        fi >> u >> v >> k;
        if (u == n_nodes) u = 0;
        if (v == n_nodes) v = 0;
        x = nodes[u].edges.size();
        y = nodes[v].edges.size();
        nodes[u].edges.push_back(edge());
        nodes[v].edges.push_back(edge());
        nodes[u].edges[x].neighbor = v;
        nodes[u].edges[x].w = k;
        nodes[u].edges[x].pos = y;
        nodes[v].edges[y].neighbor = u;
        nodes[v].edges[y].w = k;
        nodes[v].edges[y].pos = x;        
        nodes[u].s += abs(k);
        nodes[v].s += abs(k);
    }
    l = 0;
    tmp.resize(n_nodes,0);
    tmp1.resize(n_nodes,0);    
    is_adj.resize(n_nodes,false);
    is_in_group.resize(n_nodes,-1);
    is_flipped.resize(n_nodes,false);
    pos.resize(n_nodes);
    rpos.resize(n_nodes);
    check_slow = true;
}

void graph::delete_node()
{
    int u = n_nodes - 1;
    int v, pos, v_edges;
    double w;
    for (auto e: nodes[u].edges)
    {
        v = e.neighbor;
        w = e.w;
        pos = e.pos;
        v_edges = nodes[v].edges.size();
        swap(nodes[v].edges[pos],nodes[v].edges[v_edges-1]);
        nodes[v].edges.pop_back();
        nodes[v].s -= abs(w);
        if (pos < v_edges - 1)
        {            
            nodes[nodes[v].edges[pos].neighbor].edges[nodes[v].edges[pos].pos].pos = pos;
        }        
    }
    nodes[u].ID.clear();
    nodes[u].edges.clear();
    nodes[u].s = 0;
    --n_nodes;
    swap(nodes[u],nodes[n_nodes]);
    for (auto e: nodes[u].edges)
    {
        v = e.neighbor;
        w = e.w;
        pos = e.pos;
        nodes[v].edges[pos].neighbor = u;
    }    
}

void graph::add_node(node n)
{
    nodes[n_nodes] = n;
    nodes[n_nodes].l = l;
    int v;
    double w;
    int y;
    nodes[n_nodes].s = 0;
    for (int i = 0; i < int(nodes[n_nodes].edges.size()); ++i)    
    {
        v = nodes[n_nodes].edges[i].neighbor;
        w = nodes[n_nodes].edges[i].w;        
        y = nodes[v].edges.size();
        nodes[n_nodes].edges[i].pos = y;
        nodes[v].edges.push_back(edge());
        nodes[v].edges[y].neighbor = n_nodes;
        nodes[v].edges[y].w = w;
        nodes[v].edges[y].pos = i;
        nodes[v].s += abs(w);
        nodes[n_nodes].s += abs(w);
        nodes[v].l = l;
		if (nodes[n_nodes].edges.size() < sqrt(n_nodes))
			for (auto e: nodes[v].edges)
				nodes[e.neighbor].l = l;
    }
    ++n_nodes;    
}

void graph::swap_nodes(int u, int v)
{
    int k ;
    for (auto e: nodes[u].edges)
    {        
        nodes[e.neighbor].edges[e.pos].neighbor = v;
    }
    for (auto e: nodes[v].edges)
    {
        k = e.neighbor;
        if (k == v) k = u;
        nodes[k].edges[e.pos].neighbor = u;
    }
    swap(nodes[u],nodes[v]);
    swap(rpos[u],rpos[v]);
    pos[rpos[u]] = u;
    pos[rpos[v]] = v;
}


// void graph::merge_nodes(int u, int v)
// {    
//     swap_nodes(u,n_nodes-1);
//     swap_nodes(v,n_nodes-2);
//     merge_nodes();
// }

void graph::merge_nodes(int u, int v, node &a, node &b)
{    
    swap_nodes(u,n_nodes-1);
    swap_nodes(v,n_nodes-2);
    
    node c;
    u = n_nodes - 1; v = n_nodes - 2;
    // cerr << u << " " << v << endl;
    for (int i:nodes[u].ID)
        c.ID.push_back(i);
    for (int i:nodes[v].ID)
        c.ID.push_back(i);
    vector <int> adj;
    // cerr << "DONE";
    for (auto e: nodes[u].edges)
    {
        if (e.neighbor < v)
        {
            if (tmp[e.neighbor] == 0)
                adj.push_back(e.neighbor);
            tmp[e.neighbor] += e.w;
        }
    }
    // cerr << "DONE";
    for (auto e: nodes[v].edges)
    {
        if (e.neighbor < v)
        {
            if (tmp[e.neighbor] == 0)
                adj.push_back(e.neighbor);
            tmp[e.neighbor] += e.w;
        }     
    }
    // cerr << "DONE";
    int x;
    for (int k: adj)
    {
        // if (tmp[k] != 0)
        {
            x = c.edges.size();
            c.edges.push_back(edge());
            c.edges[x].neighbor = k;
            c.edges[x].w = tmp[k];
            tmp[k] = 0;
        }
    }
    // for (int j: c.ID)
    //         cerr << j << " ";
    // cerr << endl;
    // cnt = 0;
    // for (auto j: c.edges)
    // {
    //     cerr << j.neighbor << " " << j.w << " " << endl;
    //     ++cnt;
    // }
    // cerr << "--------------------------------------\n";
    // cerr << "DONE";    
    b = nodes[u];
    delete_node();
    // cerr << "DONE";
    a = nodes[v];
    delete_node();
    // cerr << "DONE";
    add_node(c);
    // cerr << "DONE";
}

void graph::restore(node &a, node &b)
{
    delete_node();
    add_node(a);
    add_node(b);    
}

void graph::restore_order(int u, int v)
{
    swap_nodes(u,n_nodes-1);
    swap_nodes(v,n_nodes-2);
}

void graph::flip(int u)
{
    for (int i: nodes[u].ID)
        is_flipped[i] = !is_flipped[i];
    // cout << "flip: ";
    // for (int i: nodes[u].ID)
    // {
    //     cout << i << " ";
    // }
    // cout << endl;
    for (int i = 0; i < int(nodes[u].edges.size()); ++i)    
    {
        nodes[u].edges[i].w = -nodes[u].edges[i].w;
        nodes[nodes[u].edges[i].neighbor].edges[nodes[u].edges[i].pos].w = nodes[u].edges[i].w;
    }
}

double graph::assign_cut(vector<int> &group,int i, int s)
{
    double res;    
    if (i == k)
    {
        if (s == 0) return 1e9;
        vector <int> adj;
        res = 0;
        int v;
        double w;
        double sum1 = 0, sum2 = 0, sum3 = 0;
        for (int u: group)
        {
            for (auto e: nodes[u].edges)
            {
                v = e.neighbor;
                w = e.w;
                if (is_in_group[v] != -1)
                {
                    if ((is_in_group[v] != is_in_group[u]) && (u < v))
                        res = res + w;
                }
                else
                {
                    if (is_adj[v] == false)
                    {
                        is_adj[v] = true;
                        adj.push_back(v);
                    }
                    if (is_in_group[u] == 0)
                    {
                        sum2 += abs(w);
                        tmp1[v] += w;
                    }
                    if (is_in_group[u] == 1)
                    {
                        sum3 += abs(w);
                        tmp1[v] -= w;
                    }
                }
            }
        }
        for (int k: adj)
        {
            sum1 = sum1 + abs(tmp1[k]);
            tmp1[k] = 0;
            is_adj[k] = false;
        }
        sum1 = sum1 / 2;
        // cout << res << " " << sum1 << " " << sum2 << " " << sum3 << endl;
        sum1 = min(sum1,sum2);
        sum1 = min(sum1,sum3);
        return res - sum1;
    }
    is_in_group[group[i]] = 0;
    res = assign_cut(group, i+1,s);
    is_in_group[group[i]] = 1;
    res = min(res, assign_cut(group, i+1,s+1));
    return res;
}

double graph::NSI(vector<int> &group)
{
    for (int i :group)
        is_in_group[i] = 0;
    k = group.size();
    if (k > 2)
    {
        int cnt = 0;
        int j = 0;
        while (cnt < k)
        {
            double sum = 0;
            for (auto e: nodes[group[j]].edges)
            {
                if (is_in_group[e.neighbor] >= 0)
                    sum += e.w;
            }
            if (sum < 0)
            {
                cnt = 0;
                flip(group[j]);
            }
            else ++cnt;
            j = (j + 1)%k;
        }
    }
    double res = assign_cut(group,1,0);
    for (int i :group)
        is_in_group[i] = -1;
    return res;   
}

void graph::merge_multiple_nodes(vector<int> &group)
{
    node a,b;
    sort(group.rbegin(), group.rend());
    int u = group[0];    
    for (int i = 1; i < int(group.size()); ++i)
    {
        merge_nodes(u,group[i],a,b);
        u = n_nodes - 1;
    }
}